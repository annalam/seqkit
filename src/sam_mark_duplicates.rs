
use crate::common::{parse_args, BamReader, BamWriter};
use std::str;
use std::collections::VecDeque;
use std::cmp::min;
use rust_htslib::bam::record::{Record, Aux};

const USAGE: &str = "
Usage:
  sam mark duplicates [options] <bam_file>

Options:
  --uncompressed    Output in uncompressed BAM format
  --ignore-umi      Ignore UMI stored in RX tag even if present

Searches BAM files for DNA fragments that were read multiple times in
sequencing. When such fragments are found, the highest quality read is
kept, and other reads are marked as duplicates.

The input BAM file must be position-sorted. Output is written to
the standard output, preserving the order and content of BAM records,
except for the duplicate flag (0x400).
";

struct Read {
	start_pos: u32,
	strand: bool,
	ready: bool,
	fraglen: u16,
	umi: Box<[u8]>,    // Consider inline string
	record: Record
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let ignore_umi = args.get_bool("--ignore-umi");

	let bam = BamReader::open(bam_path);
	let mut out = BamWriter::open("-", &bam.header(), !args.get_bool("--uncompressed"));

	// Collect statistics on how many reads were classified as duplicates.
	let mut total_reads = 0;
	let mut total_duplicates = 0;

	let mut prev_pos: u32 = 0;
	let mut prev_chr: i32 = -1;

	let mut reads: VecDeque<Read> = VecDeque::new();   // FIFO queue
	for read in bam {
		if read.is_secondary() || read.is_supplementary() {
			error!("BAM file contains secondary or supplementary reads. These are not currently supported.");
		}

		let left_pos: u32 = read.pos() as u32;
		let chr: i32 = read.tid();

		if chr != prev_chr {
			// We just entered a new chromosome, so we can finalize all
			// current read clusters, and write them out.
			find_clusters(&mut reads, u32::MAX);
			total_duplicates += flush_reads(&mut out, &mut reads);
			assert!(reads.is_empty());
			prev_chr = chr;
		} else if left_pos < prev_pos {
			error!("Input BAM file is not coordinate sorted.");
		}

		prev_pos = left_pos;

		// Calculate position of first base in the read
		let start_pos: u32 = if read.is_unmapped() { 0u32 }
			else if read.is_reverse() { read.cigar().end_pos() as u32 }
			else { left_pos };

		// Signature consists of start position, strand, and *either*
		// fragment length or UMI. Start position is stored directly
		// in the Read data structure. Unmapped reads have no signature.
		let mut umi: Vec<u8> = Vec::new();
		let mut fraglen = 0u16;
		if !read.is_unmapped() {
			if !ignore_umi {
				if let Some(Aux::String(rx)) = read.aux(b"RX") {
					umi.extend_from_slice(rx);
				}
			}

			if umi.is_empty() {
				fraglen = min(read.insert_size().abs(), u16::MAX as i64) as u16;
			}
		}

		// We guarantee that records in the BAM file are never reordered.
		// That's why even unmapped reads get added to the FIFO deque.
		reads.push_back(Read {
			start_pos, strand: !read.is_reverse(), fraglen, umi: umi.into(), 
			ready: read.is_unmapped(), record: read
		});
		total_reads += 1;

		if total_reads % 1000 == 0 {
			total_duplicates += flush_reads(&mut out, &mut reads);
			find_clusters(&mut reads, left_pos);
		}
	}

	// Make sure we flush out the last reads too.
	find_clusters(&mut reads, u32::MAX);
	total_duplicates += flush_reads(&mut out, &mut reads);
	assert!(reads.is_empty());

	eprintln!("{} / {} ({:.1}%) reads were marked as duplicates.",
		total_duplicates, total_reads,
		total_duplicates as f64 / total_reads as f64 * 100.0);
}

// Write out reads with a ready status until we find the first
// read that hasn't been processed yet. Returns the number of duplicate
// reads that were flushed, so we can keep track.
fn flush_reads(out: &mut BamWriter, reads: &mut VecDeque<Read>) -> usize {
	let mut duplicates_flushed = 0;
	while !reads.is_empty() && reads[0].ready {
		if reads[0].record.is_duplicate() { duplicates_flushed += 1; }
		out.write(&reads[0].record);
		reads.pop_front();
	}
	return duplicates_flushed;
}

// Try to identify clusters of reads sharing the same signature.
fn find_clusters(reads: &mut VecDeque<Read>, curr_pos: u32) {
	for k in 0..reads.len() {
		if reads[k].ready == true { continue; }  // Skip processed clusters

		// Skip clusters for which new reads can still be found.
		if reads[k].start_pos >= curr_pos { continue; }

		// All reads from this cluster must be in the deque, so let's
		// find them, and keep track of which duplicate has highest quality.
		let mut best = k;
		let mut best_score = reads[k].record.seq().len();
		reads[k].record.set_duplicate();
		reads[k].ready = true;
		//eprintln!("Cluster: start {}, strand {}
		for j in k+1..reads.len() {
			if reads[j].ready == true { continue; }
			if reads[j].record.pos() as u32 > reads[k].start_pos {
				break;   // Reads are sorted by left pos, so no more matches
			}
			if reads[j].start_pos != reads[k].start_pos { continue; }
			if reads[j].strand != reads[k].strand { continue; }
			if reads[j].fraglen > 0 && reads[k].fraglen > 0 &&
				reads[j].fraglen != reads[k].fraglen { continue; }
			if umi_matches(&reads[j].umi, &reads[k].umi) == false { continue; }

			reads[j].record.set_duplicate();
			reads[j].ready = true;
			let score = reads[j].record.seq().len();
			if score > best_score {
				best_score = score;
				best = j;
			}
		}

		reads[best].record.unset_duplicate();
	}
}

fn umi_matches(a: &[u8], b: &[u8]) -> bool {
	if a.is_empty() || b.is_empty() { return true; }
	if a.len() != b.len() { return false; }
	let mut mismatches = 0u16;
	for k in 0..a.len() {
		if !(a[k] == b[k] || a[k] == b'N' || b[k] == b'N') {
			mismatches += 1;
		}
	}
	return mismatches <= 1;
}

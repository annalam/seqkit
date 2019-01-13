
use crate::common::{parse_args, open_bam};
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;

// TODO: Consider using fragment signatures for paired end reads where both
// mates have aligned to the genome, but with low mapping quality. The
// reasoning is that in these situations the mapped coordinates are often
// random, and therefore duplicate identification based on mapping coordinates
// will not work properly.

const USAGE: &str = "
Usage:
  sam mark duplicates [options] <bam_file>

Options:
  --debug              Print information about duplicate marking process
  --mark-semi-aligned  Mark duplicates among mapped reads with an unmapped mate
  --mark-unaligned     Mark duplicates among unaligned reads
  --allow-erosion=N    Assume that two DNA fragments can be duplicates even if
                       one of them has eroded N bases at one end (but not
                       both ends) [default: 0]


Description:
This command identifies PCR/optical duplicates in the input BAM file, and outputs a new BAM file where duplicate reads have been flagged with the 0x400 (\"PCR or optical duplicate\") flag. Memory usage is lowest when the tool is run on position-sorted BAM files, but the tool does also work with unsorted BAM files. A BAM index is not required.

By default, the software only marks duplicates among paired end reads where both mates have aligned to the genome. For these read pairs, duplicates are identified based on fragment boundaries.

If the --mark-semi-aligned flag is provided, the software will also mark duplicates among paired end reads where only one mate has aligned to the genome. For these read pairs, a 20 bp fragment signature will be used to identify the fragment's other boundary. Two DNA fragments are then considered duplicates if the 

For unaligned reads, 20 bp from each end of a DNA fragment are used as a fragment signature.

When two or more redundant read pairs are identified, the software leaves one read pair unmarked, and marks the other redundant read pairs with the duplicate flag. If BASEQ information is available, the read pair with highest cumulative BASEQ is left unmarked. If BASEQ information is not available, the read pair with highest cumulative number of unambiguous nucleotides is left unmarked.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let debug = args.get_bool("--debug");

	let mut bam = open_bam(bam_path);
	let header = bam.header().clone();

	let mut out = bam::Writer::from_stdout(&bam::header::Header::from_template(&header)).unwrap();

	// These variables are used to collect statistics on what percentage
	// of reads were classified as duplicates.
	let mut total_aligned_reads = 0;
	let mut total_marked_duplicate = 0;

	let mut candidates: Vec<Record> = Vec::new();
	let mut cluster_chr = -1;
	let mut cluster_pos: i32 = 0;
	for r in bam.records() {
		let read = r.unwrap();

		// Unmapped reads are never marked as duplicates, so just output them
		if read.is_unmapped() {
			out.write(&read).unwrap();
			continue;
		}

		// Accumulate a cluster of reads that have all been aligned to
		// the same chromosome and the same start position.
		if read.tid() == cluster_chr && read.pos() == cluster_pos {
			candidates.push(read);
			continue;
		}

		// At this point we have a cluster of reads that have all been
		// aligned to the same chromosome and the same start position.
		// We now go through all of them and check which are real duplicates.
		total_aligned_reads += candidates.len();
		let mut duplicate = vec![false; candidates.len()];
		for c in 0..candidates.len() {
			if debug && candidates.len() > 1 {
				print_read(&candidates[c], duplicate[c]);
			}
			if duplicate[c] {
				candidates[c].set_duplicate();
				out.write(&candidates[c]).unwrap();
				total_marked_duplicate += 1;
				continue;
			}
			for d in c+1..candidates.len() {
				if duplicate[d] { continue; }  // Already flagged as duplicate
				if are_duplicates(&candidates[c], &candidates[d]) {
					duplicate[d] = true;
				}
			}
			out.write(&candidates[c]).unwrap();
		}

		if debug && candidates.len() > 1 { eprintln!("---------------"); }

		// Quick check to ensure that BAM file is position-sorted.
		assert!(read.pos() > cluster_pos || read.tid() != cluster_chr);

		// We can now start a new cluster of reads.
		cluster_chr = read.tid();
		cluster_pos = read.pos();
		candidates.clear();
		candidates.push(read);
	}

	eprintln!("{} / {} ({:.1}%) aligned reads were marked as duplicates.",
		total_marked_duplicate, total_aligned_reads,
		total_marked_duplicate as f64 / total_aligned_reads as f64 * 100.0);
}

// Two reads are considered duplicates if:
// - They are in the same chromosome
// - They have the same start position
// - They are in the same strand
// - Their mates are not in different chromosomes
// - Their mates are not in different strands
// - Their mate positions do not differ by more than 1000 bp
fn are_duplicates(a: &Record, b: &Record) -> bool {
	if a.tid() != b.tid() { return false; }
	if a.pos() != b.pos() { return false; }
	if a.is_reverse() != b.is_reverse() { return false; }
	if a.is_mate_unmapped() == false && b.is_mate_unmapped() == false {
		if a.mtid() != b.mtid() { return false; }
		if a.insert_size() != b.insert_size() { return false; }
		if a.is_mate_reverse() != b.is_mate_reverse() { return false; }
		if a.insert_size() == 0 && (a.mpos() - b.mpos()).abs() > 1000 {
			return false;
		}
		if (a.mpos() - b.mpos()).abs() > 1000 { return false; }
	}
	return true;
}

fn print_read(read: &Record, duplicate: bool) {
	let pos = if read.is_reverse() {
		read.cigar().end_pos().unwrap() - 1
	} else {
		read.pos()
	};
	eprintln!("chr?:{} ({}), insert size = {}: {}", pos,
		if read.is_reverse() { '-' } else { '+' }, read.insert_size(),
		if duplicate { "" } else { " ***" });
}

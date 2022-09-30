
use crate::common::{parse_args, BamReader, BamWriter};
use std::str;
use std::collections::{HashMap, VecDeque};
use rust_htslib::bam::Record;

const USAGE: &str = "
Usage:
  sam recalculate tlen [options] <bam_file>

Options:
  --uncompressed    Output in uncompressed BAM format
  --max-len=N       Maximum fragment length [default: 5000]

Recalculates SAM record TLEN (template length) field based on the distance
between the 5' ends of paired mates. This ensures that the TLEN represents the 
actual DNA fragment length.
";

struct Read {
	ready: bool,
	tlen: i32,
	record: Record
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let max_frag_len: i64 = args.get_str("--max-len").parse().unwrap_or(0);
	if max_frag_len <= 0 { error!("--max-len must be a positive integer."); }

	let bam = BamReader::open(&bam_path);
	let mut out = BamWriter::open("-", &bam.header(), !args.get_bool("--uncompressed"));

	// In order to preserve the order of SAM records, we need to maintain
	// a FIFO queue. We also maintain a hashmap that remembers the FIFO queue
	// position of each read's mate.
	let mut reads: VecDeque<Read> = VecDeque::new();
	let mut mates: HashMap<Box<str>, usize> = HashMap::new();
	let mut reads_written: usize = 0;

	for read in bam {
		let qname = str::from_utf8(read.qname()).unwrap();
		if read.is_paired() == false || read.is_unmapped() ||
			read.is_mate_unmapped() {
			// Template length is zero for single reads and unmapped pairs
			reads.push_back(Read { ready: true, tlen: 0, record: read });
		} else if read.is_secondary() || read.is_supplementary() {
			error!("Secondary and supplementary read alignments are not supported.");
		} else if read.tid() != read.mtid() ||
			(read.pos() - read.mpos()).abs() > max_frag_len ||
			read.is_reverse() == read.is_mate_reverse() {
			// The mates are discordantly aligned, set TLEN as zero
			reads.push_back(Read { ready: true, tlen: 0, record: read });
		} else {
			// For mapped pairs, calculate 5'-to-5' template length			
			if let Some(mate_idx) = mates.remove(qname) {
				let mate = &mut reads[mate_idx - reads_written].record;
				assert!(read.is_first_in_template() != mate.is_first_in_template());
				let start_pos = if read.is_reverse() { read.cigar().end_pos() - 1 } else { read.pos() };
				let mate_start_pos = if mate.is_reverse() { mate.cigar().end_pos() - 1 } else { mate.pos() };

				// Calculate the template length, and set the sign according to
				// the SAM specification (left-most mate has positive TLEN). 
				let mut tlen = (start_pos - mate_start_pos).abs() + 1;
				if read.pos() > mate.pos() { tlen = -tlen; }

				// If the mates are not in converging orientation, set TLEN
				// to zero for both mates.
				if (start_pos < mate_start_pos && read.is_reverse()) ||
					(start_pos > mate_start_pos && !read.is_reverse()) {
					tlen = 0;
				}

				// The TLENs of the two mates must always have opposing signs
				reads.push_back(Read { ready: true, record: read,
					tlen: tlen as i32 });
				reads[mate_idx - reads_written].tlen = -tlen as i32;
				reads[mate_idx - reads_written].ready = true;
			} else {
				mates.insert(qname.into(), reads_written + reads.len());
				reads.push_back(Read { ready: false, record: read, tlen: 0 });
			}
		}

		// Note for large BAM files the memory usage can get quite high
		// since we try to preserve the order of BAM records, and don't
		// write records into the output stream until a mate has been found.
		while !reads.is_empty() && reads[0].ready {
			let new_tlen = reads[0].tlen as i64;
			reads[0].record.set_insert_size(new_tlen);
			out.write(&reads[0].record);
			reads.pop_front();
			reads_written += 1;
		}
	}

	// If there are any remaining reads, it indicates that some mates were
	// missing from the BAM file. In that scenario, we skip the reads that
	// had a missing mate, and write out the other reads.
	for mut read in reads {
		if read.ready == false {
			eprintln!("Read {} discarded due to missing mate.", str::from_utf8(read.record.qname()).unwrap());
			continue;
		}

		let new_tlen = read.tlen as i64;
		read.record.set_insert_size(new_tlen);
		out.write(&read.record);
		//reads_written += 1;
	}
}

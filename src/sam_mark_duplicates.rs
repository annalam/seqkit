
use crate::common::{parse_args, open_bam};
use std::str;
use std::cmp::{min, max};
use rust_htslib::bam;
use rust_htslib::bam::record::{Record, Seq};
use rust_htslib::bam::{Read, ReadError};
use hashbrown::HashSet;

const USAGE: &str = "
Usage:
  sam mark duplicates [options] <bam_file>

Options:
  --uncompressed    Output in uncompressed BAM format

Description:
This command identifies PCR/optical duplicates in the input BAM file, and
outputs a new BAM file where duplicate reads have been flagged with the
0x400 (\"PCR or optical duplicate\") flag. The input BAM file must be
name-sorted.

20 bp from each end of a DNA fragment are used as a fragment signature.

When two or more redundant read pairs are identified, the software leaves one read pair unmarked, and marks the other redundant read pairs with the duplicate flag. If BASEQ information is available, the read pair with highest cumulative BASEQ is left unmarked. If BASEQ information is not available, the read pair with highest cumulative number of unambiguous nucleotides is left unmarked.
";

struct Fragment {
	signature_1: u32,	// 16 bases, 2 bits per base
	signature_2: u32    // 16 bases, 2 bits per base
}

struct FragmentUMI {
	signature_1: u16,	// 8 bases, 2 bits per base
	signature_2: u16,   // 8 bases, 2 bits per base
	umi: u32            // Max 16 bases, 2 bits per base
}

// Function for reading BAM records, with proper user-friendly messages.
// Returns false after reading the last record, or if reading fails.
fn read_record(bam: &mut bam::Reader, record: &mut bam::Record) -> bool {
	match bam.read(record) {
		Err(ReadError::NoMoreRecord) => false,
		Err(ReadError::Truncated) => error!("BAM file ended prematurely."),
		Err(ReadError::Invalid) => error!("Invalid BAM record."),
		Ok(_) => true
	}
}

// Here we pack the first 16 bases of the read into a compact
// 32-bit integer. If the read is shorter than 16 bases long, we
// mark the missing bases as 'A'.
fn mate_signature(read: &Record) -> u32 {
	let seq = read.seq();
	let mut signature = 0u32;
	if read.is_reverse() {
		let start = if seq.len() >= 16 { seq.len() - 16 } else { 0 };

		for k in (start..seq.len()).rev() {
			signature = signature * 4 + match seq.encoded_base(k) {
				1 => 3,   // A, encoded as its complement T
				2 => 2,   // C, encoded as its complement G
				4 => 1,   // G, encoded as its complement C
				_ => 0,   // T or ambiguous base, encoded as its complement A
			};
		}
	} else {
		for k in 0..min(seq.len(), 16) {
			signature = signature * 4 + match seq.encoded_base(k) {
				2 => 1,   // C
				4 => 2,   // G
				8 => 3,   // T
				_ => 0,   // A (or ambiguous base)
			};
		}
	}
	
	signature
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut bam = open_bam(bam_path);
	let header = bam.header().clone();

	let mode: &[u8] = match args.get_bool("--uncompressed") {
		true => b"wbu", false => b"wb"
	};

	let mut out = bam::Writer::new(b"-", mode,
		&bam::Header::from_template(&header)).unwrap();

	// Collect statistics on how many reads were classified as duplicates.
	let mut total_reads = 0;
	let mut total_duplicates = 0;

	let mut read_1 = Record::new();
	let mut read_2 = Record::new();

	// TODO: Use a simpler hash function (modulo would be sufficient).
	let mut seen_signatures: HashSet<u64> = HashSet::new();

	// This variable only exists so that we can make the code robust towards
	// a few missing reads in the BAM file. We can probably drop this and
	// simplify the code later.
	let mut unpaired = 0;

	loop {
		if unpaired > 10 {
			error!("Stopping early due to abnormally high number of consecutive reads with non-matching IDs. Please sort the input BAM file by read ID using 'samtools sort -n'.");
		}

		if unpaired > 0 {
			std::mem::swap(&mut read_1, &mut read_2);
			//read_1 = read_2.clone();  // Results in malformatted "read_1"
			if read_record(&mut bam, &mut read_2) == false { break; }
			total_reads += 1;
		} else {
			if read_record(&mut bam, &mut read_1) == false { break; }
			if read_record(&mut bam, &mut read_2) == false { break; }
			total_reads += 2;
		}

		if read_1.is_paired() == false || read_2.is_paired() == false {
			error!("WARNING: BAM file contains unpaired reads. Only paired end reads are currently supported.");
		}

		if read_1.is_secondary() || read_1.is_supplementary() ||
			read_2.is_secondary() || read_2.is_supplementary() {
			error!("Input BAM file contains secondary or supplementary reads. These are not currently supported.");
		}

		if read_1.qname() != read_2.qname() {
			eprintln!("WARNING: Input BAM file contains consecutive paired end reads #{} and #{} with different IDs '{}' and '{}'.",
				total_reads -1, total_reads,
				str::from_utf8(read_1.qname()).unwrap(),
				str::from_utf8(read_2.qname()).unwrap());

			unpaired += 1;
			continue;  // Jump to top and read next record in BAM file
		}

		unpaired = 0;  // Reset the 'unpaired' counter, we found a pair.

		let sig_1 = mate_signature(&read_1);
		let sig_2 = mate_signature(&read_2);
		/*if sig_1 == sig_2 {
			eprintln!("Read 1 with ID {} and signature {}:\n{} ({})",
				str::from_utf8(&read_1.qname()).unwrap(), sig_1, 
				str::from_utf8(&read_1.seq().as_bytes()).unwrap(),
				if read_1.is_reverse() { "-" } else { "+" });
			eprintln!("Read 2 with ID {} and signature {}:\n{} ({})",
				str::from_utf8(&read_2.qname()).unwrap(), sig_2, 
				str::from_utf8(&read_2.seq().as_bytes()).unwrap(),
				if read_2.is_reverse() { "-" } else { "+" });
		}*/

		// Since a double-stranded DNA fragment can be read at either
		// orientation (i.e. mate #1 can be at either one of the ends), we
		// need a DNA fragment signature that is not affected by orientation.
		// Given a paired end read, we extract 16 bp from the beginning of
		// both reads, and call these the mate-signatures A and B.
		// The total signature is then A+B if A has lower lexicographic rank,
		// or B+A if A has higher lexicographic rank.
		let signature: u64 = if sig_2 < sig_1 {
			(sig_2 as u64) | (sig_1 as u64) << 32
		} else {
			(sig_1 as u64) | (sig_2 as u64) << 32
		};

		if seen_signatures.contains(&signature) {
			read_1.set_duplicate();
			read_2.set_duplicate();
			total_duplicates += 2;
		} else {
			seen_signatures.insert(signature);
		}

		out.write(&read_1).unwrap();
		out.write(&read_2).unwrap();
	}

	eprintln!("{} / {} ({:.1}%) reads were marked as duplicates.",
		total_duplicates, total_reads,
		total_duplicates as f64 / total_reads as f64 * 100.0);
}
/*
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
*/

use crate::common::{parse_args, GzipWriter, open_bam};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::{Read, header::Header};
use rust_htslib::bam::record::Record;

const USAGE: &str = "
Usage:
  sam minimize [options] <bam_file>

Options:
  --uncompressed    Output in uncompressed BAM format

Changes read IDs into simple numeric identifiers, removes per-base qualities,
and removes all auxiliary fields.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut bam = open_bam(bam_path);
	let header = bam.header().clone();

	// These variables are used for QNAME simplification
	let mut highest_id: u32 = 0;
	let mut qname_to_id: HashMap<Vec<u8>, u32> = HashMap::new();

	let mode: &[u8] = match args.get_bool("--uncompressed") {
		true => b"wbu", false => b"wb"
	};
	//let mut out = bam::Writer::from_stdout(&Header::from_template(&header)).unwrap();
	let mut out = bam::Writer::new(b"-", mode,
		&Header::from_template(&header)).unwrap();

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));

		let mut qname = read.qname().to_vec();

		// If the read identifier is not a simple ASCII number, we simplify it
		if qname.iter().all(|c| c.is_ascii_digit()) == false {
			if let Some(slash_pos) = qname.iter().position(|x| *x == b'/') {
				qname.truncate(slash_pos);
			}

			let id = if let Some(x) = qname_to_id.remove(&qname) { x } else {
				highest_id += 1;
				qname_to_id.insert(qname, highest_id);
				highest_id
			};
			qname = format!("{}", id).into_bytes();
		}

		let cigar = read.cigar();
		let seq = read.seq().as_bytes();

		// According to BAM specification, missing per-base quality information
		// is denoted with 0xFF bytes (whose number must equal sequence length)
		let qual = vec![0xFFu8; seq.len()];

		// The call to .set() removes all AUX fields
		read.set(&qname, &cigar, &seq, &qual);
		out.write(&read).unwrap();
	}
}

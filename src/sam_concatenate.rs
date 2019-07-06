
use crate::common::{parse_args, GzipWriter, open_bam};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::{Read, header::Header};
use rust_htslib::bam::record::Record;

const USAGE: &str = "
Usage:
  sam concatenate [options] <bam_files>...

Options:
  --uncompressed    Output in uncompressed BAM format

Concatenates two or more BAM files together, while ensuring that read
identifiers do not clash. Currently we simply add e.g. a '.1' suffix to all
read identifiers found in the first BAM file.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_paths = args.get_vec("<bam_files>");
	if bam_paths.len() < 2 {
		error!("At least two BAM files must be provided for concatenation.");
	}

	let mode: &[u8] = match args.get_bool("--uncompressed") {
		true => b"wbu", false => b"wb"
	};
	let header = open_bam(bam_paths[0]).header().clone();
	let mut out = bam::Writer::new(b"-", mode,
		&Header::from_template(&header)).unwrap();

	for b in 0..bam_paths.len() {
		let bam_path = bam_paths[b];
		let mut bam = open_bam(bam_path);
		let header = bam.header().clone();
		// TODO: Sanity check that headers are consistent.

		let suffix = format!(".{}", b + 1).into_bytes();

		for r in bam.records() {
			let mut read = r.unwrap_or_else(
				|_| error!("Input BAM file ended prematurely."));

			let mut qname = read.qname().to_vec();
			qname.extend(&suffix);
			read.set_qname(&qname);
			out.write(&read).unwrap();
		}
	}
}

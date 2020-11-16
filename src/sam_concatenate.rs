
use crate::common::{parse_args, GzipWriter, BamReader};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::{Read, header::Header, record::Record, Format, CompressionLevel};

const USAGE: &str = "
Usage:
  sam concatenate [options] <bam_files>...

Options:
  --uncompressed    Output in uncompressed BAM format

Concatenates two or more BAM files together, ensuring that read
identifiers do not clash. Currently we simply add e.g. a '.1' suffix to all
read identifiers found in the first BAM file.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_paths = args.get_vec("<bam_files>");
	if bam_paths.len() < 2 {
		error!("At least two BAM files must be provided for concatenation.");
	}

	let header = BamReader::open(bam_paths[0]).header();
	let mut out = bam::Writer::from_stdout(
		&Header::from_template(&header), Format::BAM).unwrap();
	if args.get_bool("--uncompressed") {
		out.set_compression_level(CompressionLevel::Uncompressed);
	}

	for b in 0..bam_paths.len() {
		let bam_path = bam_paths[b];
		let mut bam = BamReader::open(bam_path);
		let header = bam.header().clone();
		// TODO: Sanity check that headers are consistent.

		let suffix = format!(".{}", b + 1).into_bytes();

		for mut read in bam {
			let mut qname = read.qname().to_vec();
			qname.extend(&suffix);
			read.set_qname(&qname);
			out.write(&read).unwrap();
		}
	}
}


use crate::common::{parse_args, BamReader};
use std::str;
use regex::Regex;
use rust_htslib::bam::{Header, Writer, record::Aux, Format, CompressionLevel};

const USAGE: &str = "
Usage:
  sam qname from tags [options] <bam_file>

Options:
  --uncompressed     Output in uncompressed BAM format

Finds tags (e.g. \"RX:xxxx\") in each BAM record, and appends them to the QNAME.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let tag_regex = Regex::new(r" [A-Z]+:\S*").unwrap();

	let bam = BamReader::open(&bam_path);
	let header = bam.header();

	let mut out = Writer::from_stdout(
		&Header::from_template(&header), Format::BAM).unwrap();
	if args.get_bool("--uncompressed") {
		out.set_compression_level(CompressionLevel::Uncompressed);
	}

	for mut read in bam {
		let mut qname = read.qname().to_vec();
		if let Some(Aux::String(tag)) = read.aux(b"RX") {
			qname.extend(b" RX:");
			qname.extend(tag);
			read.set_qname(&qname);
		}

		out.write(&read).unwrap();
	}
}

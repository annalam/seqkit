
use crate::common::{parse_args, BamReader};
use std::str;
use regex::Regex;
use rust_htslib::bam::{Header, Writer, record::Aux, Format, CompressionLevel};

const USAGE: &str = "
Usage:
  sam tags from qname [options] <bam_file>

Options:
  --uncompressed     Output in uncompressed BAM format

Finds tags (e.g. \"UMI:xxxx\") in the qname of each BAM record, and turns
them into actual SAM format tags.
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
		let qname = read.qname().to_vec();
		let mut parts = qname.split(|x| *x == b' ');
		let qname_without_tags = parts.next().unwrap();

		read.set_qname(&qname_without_tags);

		for tag in parts {
			if tag.starts_with(b"UMI:") {
				read.push_aux(b"RX", &Aux::String(&tag[4..]));
			} else if tag.len() >= 3 && tag[2] == b':' {
				read.push_aux(&tag[0..2], &Aux::String(&tag[3..]));
			} else {
				error!("Tag '{}' is not supported.",
					str::from_utf8(tag).unwrap());
			}
		}

		out.write(&read).unwrap();
	}
}

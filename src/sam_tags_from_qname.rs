
use crate::common::{parse_args, GzipWriter, open_bam};
use std::str;
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bam::{Read, header::Header};
use rust_htslib::bam::record::{Record, Aux};

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

	let mut bam = open_bam(bam_path);
	let header = bam.header().clone();

	let mode: &[u8] = match args.get_bool("--uncompressed") {
		true => b"wbu", false => b"wb"
	};
	let mut out = bam::Writer::new(b"-", mode,
		&Header::from_template(&header)).unwrap();

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));

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

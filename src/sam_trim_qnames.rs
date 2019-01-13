
use crate::common::parse_args;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;

const USAGE: &str = "
Usage:
  sam trim qnames [options] <bam_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut bam = if bam_path == "-" {
		bam::Reader::from_stdin().unwrap()
	} else {
		bam::Reader::from_path(&bam_path).unwrap_or_else(
			|_| error!("Cannot open BAM file '{}'", bam_path))
	};

	let header = bam.header().clone();

	let mut out = bam::Writer::from_stdout(&bam::header::Header::from_template(&header)).unwrap();

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));
		let qname = read.qname().to_vec();
		if let Some(mut trim) = qname.iter().position(|x| *x == b' ') {
			if qname[trim - 2] == b'/' && (qname[trim - 1] == b'1' || qname[trim - 1] == b'2') {
				trim -= 2;
			}
			read.set_qname(&qname[..trim]);
		}
		out.write(&read).unwrap_or_else(
			|_| error!("Output stream closed unexpectedly."));
	}
}

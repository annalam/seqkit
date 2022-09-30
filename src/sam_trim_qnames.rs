
use crate::common::{parse_args, BamReader};
use rust_htslib::bam::{Header, Writer, Format};
use std::str;

const USAGE: &str = "
Usage:
  sam trim qnames [options] <bam_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let bam = BamReader::open(&bam_path);
	let header = bam.header();

	let mut out = Writer::from_stdout(&Header::from_template(&header), Format::BAM).unwrap();

	for mut read in bam {
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

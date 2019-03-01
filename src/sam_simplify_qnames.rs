
use crate::common::{parse_args, GzipWriter, open_bam};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;

const USAGE: &str = "
Usage:
  sam repair [options] <bam_file>

Options:
  --threads=N    Number of threads to use for BAM compression [default: 1]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let threads: usize = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("Invalid value for --threads."));

	let mut bam = open_bam(bam_path);
	let header = bam.header().clone();

	let mut out = bam::Writer::from_stdout(&bam::header::Header::from_template(&header)).unwrap();
	if threads > 1 { out.set_threads(threads); }

	let mut highest_id: u32 = 0;
	let mut qname_to_id: HashMap<Vec<u8>, u32> = HashMap::new();

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));

		let mut qname = read.qname().to_owned();
		if let Some(slash_pos) = qname.iter().position(|x| *x == b'/') {
			qname.truncate(slash_pos);
		}

		let id = if let Some(x) = qname_to_id.remove(&qname) { x } else {
			highest_id += 1;
			qname_to_id.insert(qname, highest_id);
			highest_id
		};

		let new_qname = format!("{}", id).into_bytes();
		read.set_qname(&new_qname);
		out.write(&read).unwrap();
	}
}

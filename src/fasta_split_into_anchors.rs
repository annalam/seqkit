
use crate::common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta split into anchors <fastq> <anchor_len>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq>"));
	let anchor_len: usize = args.get_str("<anchor_len>").parse()
		.unwrap_or_else(|_|error!("<anchor_len> must be a positive integer."));

	let mut reads: usize = 0;
	let mut header = String::new();
	let mut seq = String::new();
	let mut qual = String::new();
	let mut line = String::new();

	while fastq.read_line(&mut header) {
		reads += 1;

		fastq.read_line(&mut seq);
		let seq_len = seq.trim_end().len();
		if seq_len < anchor_len * 2 { continue; }

		if header.starts_with('@') {
			fastq.read_line(&mut line);
			fastq.read_line(&mut qual);
			print!("@{}\n{}\n+\n{}\n", reads, &seq[0..anchor_len],
				&qual[0..anchor_len]);
			print!("@{}\n{}\n+\n{}\n", reads,
				&seq[seq_len-anchor_len..seq_len],
				&qual[seq_len-anchor_len..seq_len]);
		} else if header.starts_with('>') {
			print!(">{}\n{}\n", reads, &seq[0..anchor_len]);
			print!(">{}\n{}\n", reads, &seq[seq_len-anchor_len..seq_len]);
		} else {
			error!("Header is not valid FASTA/FASTQ:\n{}", header);
		}
	}
}

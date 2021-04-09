
use crate::common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta convert basespace <fastq_file>

Description:
FASTQ files from Illumina Basespace typically display adapter barcodes at
the end of the FASTQ header, and have read identifiers that end in /1 or /2.
This tool replaces the read identifiers by simple consecutive integers, and
places a \"BC:\" prefix in front of the barcode. An example FASTQ header in
the output could look like this: @412435 BC:TAGCTACT
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));

	let mut num_read_pairs: u64 = 0;

	let mut header = String::new();
	let mut line = String::new();
	while fastq.read_line(&mut header) == true {
		num_read_pairs += 1;
		print!("@{}", num_read_pairs);

		// TODO: Add option for keeping the original fragment ID.
		// print!("{}", header.split(' ').next().unwrap());

		let barcode = header.trim_end().split(':').last().unwrap_or("");
		if !barcode.is_empty() { print!(" BC:{}", barcode); }
		println!();

		if header.starts_with('@') {
			fastq.read_line(&mut line); print!("{}", line);
			fastq.read_line(&mut line); print!("{}", line);
			fastq.read_line(&mut line); print!("{}", line);
		} else if header.starts_with('>') {
			fastq.read_line(&mut line); print!("{}", line);
		} else {
			error!("Invalid FASTQ line:\n{}", header);
		}
	}
}

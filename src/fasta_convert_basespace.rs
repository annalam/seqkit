
use common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta convert basespace <fastq_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));

	let mut header = String::new();
	let mut line = String::new();
	while fastq.read_line(&mut header) == true {
		let barcode = header.trim_right().split(':').last().unwrap_or("");
		print!("{}", header.split(' ').next().unwrap());
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

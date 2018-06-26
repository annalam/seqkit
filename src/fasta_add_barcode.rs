
use common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta add barcode <fastq_file> <barcode_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));
	let mut barcode_file = FileReader::new(&args.get_str("<barcode_file>"));

	let mut header = String::new();
	let mut barcode = String::new();
	let mut line = String::new();
	loop {
		barcode_file.read_line(&mut header);
		if header.starts_with('@') {
			barcode_file.read_line(&mut barcode);
			barcode_file.read_line(&mut line);
			barcode_file.read_line(&mut line);
		} else if header.starts_with('>') {
			barcode_file.read_line(&mut barcode);
		}

		if fastq.read_line(&mut header) == false {
			break;    // End of file
		}

		print!("{} BC:{}\n", header.trim_right(), barcode.trim_right());

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

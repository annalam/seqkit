
use common::{parse_args, read_buffered, AsciiBufRead};
use std::str;

const USAGE: &str = "
Usage:
  fasta barcode to header <fastq_file> <barcode_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = read_buffered(&args.get_str("<fastq_file>"));
	let mut barcode_file = read_buffered(&args.get_str("<barcode_file>"));

	let mut header = String::new();
	let mut barcode = String::new();
	let mut line = String::new();
	loop {
		barcode_file.next_line(&mut barcode);
		if barcode.starts_with('@') {
			barcode_file.next_line(&mut barcode);
			barcode_file.next_line(&mut line);
			barcode_file.next_line(&mut line);
		} else if barcode.starts_with('>') {
			barcode_file.next_line(&mut barcode);
		}

		if fastq.next_line(&mut header) == false {
			break;    // End of file
		}
		if header.starts_with('@') {
			println!("{} BC:{}", header.trim_right(), barcode.trim_right());
			fastq.next_line(&mut line); print!("{}", line);
			fastq.next_line(&mut line); print!("{}", line);
			fastq.next_line(&mut line); print!("{}", line);
		} else if header.starts_with('>') {
			println!("{} BC:{}", header.trim_right(), barcode.trim_right());
			fastq.next_line(&mut line); print!("{}", line);
		} else {
			error!("Invalid FASTQ line:\n{}", header);
		}
	}
}

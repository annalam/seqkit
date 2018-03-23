
use common::{parse_args, FileReader};
use std;
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta barcode to header <fastq_file> <barcode_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));
	let mut barcode_file = FileReader::new(&args.get_str("<barcode_file>"));

	// Lock the stdout to prevent constant locking on print!().
	let stdout = std::io::stdout();
	let mut out = stdout.lock();

	let mut header = String::new();
	let mut barcode = String::new();
	let mut line = String::new();
	loop {
		barcode_file.read_line(&mut barcode);
		if barcode.starts_with('@') {
			barcode_file.read_line(&mut barcode);
			barcode_file.read_line(&mut line);
			barcode_file.read_line(&mut line);
		} else if barcode.starts_with('>') {
			barcode_file.read_line(&mut barcode);
		}

		if fastq.read_line(&mut header) == false {
			break;    // End of file
		}
		if header.starts_with('@') {
			write!(out, "{} BC:{}\n", header.trim_right(),
				barcode.trim_right());
			fastq.read_line(&mut line); print!("{}", line);
			fastq.read_line(&mut line); print!("{}", line);
			fastq.read_line(&mut line); print!("{}", line);
		} else if header.starts_with('>') {
			write!(out, "{} BC:{}\n", header.trim_right(),
				barcode.trim_right());
			fastq.read_line(&mut line); print!("{}", line);
		} else {
			error!("Invalid FASTQ line:\n{}", header);
		}
	}
}


use common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta add barcode <fastq_file> <barcode_file> <barcode_format>

Description:
Barcode format is specified as a string. Example: UUUUSSSS.
U = UMI base. S = sample barcode.
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));
	let mut barcode_file = FileReader::new(&args.get_str("<barcode_file>"));
	let barcode_format = args.get_str("<barcode_format>");
	if barcode_format.chars().any(|c| c != 'U' && c != 'S') {
		error!("Invalid barcode format: {}", barcode_format);
	}

	let mut header = String::new();
	let mut barcode = String::new();
	let mut line = String::new();
	let mut sample_barcode = String::new();
	let mut umi = String::new();
	loop {
		barcode_file.read_line(&mut header);
		if header.starts_with('@') {
			barcode_file.read_line(&mut barcode);
			barcode_file.read_line(&mut line);
			barcode_file.read_line(&mut line);
		} else if header.starts_with('>') {
			barcode_file.read_line(&mut barcode);
		}

		umi.clear(); sample_barcode.clear();
		if barcode.trim_right().len() != barcode_format.len() {
			error!("Barcode '{}' does not match with format specifier '{}'.",
				barcode, barcode_format);
		}
		for x in barcode_format.chars().zip(barcode.chars()) {
			if x.0 == 'U' { umi.push(x.1); }
			else if x.0 == 'S' { sample_barcode.push(x.1); }
		}

		if fastq.read_line(&mut header) == false {
			break;    // End of file
		}

		print!("{}", header.trim_right());
		if !sample_barcode.is_empty() { print!(" SI:{}", sample_barcode); }
		if !umi.is_empty() { print!(" UMI:{}", umi); }
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

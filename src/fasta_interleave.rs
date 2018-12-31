
use common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta interleave <fastq_1> <fastq_2>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq_1 = FileReader::new(&args.get_str("<fastq_1>"));
	let mut fastq_2 = FileReader::new(&args.get_str("<fastq_2>"));

	let mut line = String::new();
	while fastq_1.read_line(&mut line) {
		let lines = if line.starts_with('@') { 4 }
			else if line.starts_with('>') { 2 }
			else { error!("Line is not FASTA/FASTQ format: {}", line); };
		print!("{}", line);
		for k in 0..lines-1 {
			fastq_1.read_line(&mut line); print!("{}", line);
		}

		fastq_2.read_line(&mut line);
		if (lines == 4 && !line.starts_with('@')) ||
			(lines == 2 && !line.starts_with('>')) {
			error!("Input files do not share a consistent format.");
		}
		print!("{}", line);
		for k in 0..lines-1 {
			fastq_2.read_line(&mut line); print!("{}", line);
		}
	}
}

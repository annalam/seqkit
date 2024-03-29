
use crate::common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta to raw <fasta_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fasta_file>"));

	let mut line = String::new();
	while fasta_file.read_line(&mut line) {
		if line.starts_with('>') {
			fasta_file.read_line(&mut line);
			print!("{}", line);
		} else if line.starts_with('@') {
			fasta_file.read_line(&mut line);
			print!("{}", line);
			// Discard per-base qualities
			fasta_file.read_line(&mut line);
			fasta_file.read_line(&mut line);
		} else {
			error!("Invalid FASTA/FASTQ format encountered.");
		}
	}
}

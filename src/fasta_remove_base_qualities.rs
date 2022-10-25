
use crate::common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta remove base qualities <fastq_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq_file = FileReader::new(&args.get_str("<fastq_file>"));

	let mut line = String::new();
	while fastq_file.read_line(&mut line) {
		if line.starts_with('@') {
			print!(">{}", &line[1..]);   // FASTA format has '>' instead of '@'
			fastq_file.read_line(&mut line);
			print!("{}", line);
			// Discard the per-base qualities
			fastq_file.read_line(&mut line);
			fastq_file.read_line(&mut line);
		} else {
			error!("Invalid FASTQ format encountered.");
		}
	}
}

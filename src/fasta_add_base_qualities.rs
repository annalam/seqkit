
use crate::common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta add base qualities <fasta> <baseq>

Converts a FASTA file into a FASTQ file based on user-specified dummy base
quality values.
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fasta>"));
	let baseq: u8 = args.get_str("<baseq>").parse().unwrap_or_else(|_|
		error!("Base quality must be between 0 - 255."));

	let mut line = String::new();
	while fasta_file.read_line(&mut line) {
		if line.starts_with('>') {
			print!("@{}", &line[1..]);
			fasta_file.read_line(&mut line);
			let seq_len = line.len() - 1;
			print!("{}", line);
			print!("+\n{}\n", str::from_utf8(&vec!(33u8 + baseq; seq_len)).unwrap());
		} else {
			error!("Invalid FASTA format encountered.");
		}
	}
}

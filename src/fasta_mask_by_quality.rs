
use crate::common::{parse_args, FileReader};
use std::str;
use std::fmt::Write;

const USAGE: &str = "
Usage:
  fasta mask by quality <fastq_file> <min_baseq>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));
	let min_baseq: u8 = args.get_str("<min_baseq>").parse().unwrap();

	let mut seq = String::new();
	let mut base_qualities = String::new();
	let mut line = String::new();
	let mut output = String::new();
	while fasta_file.read_line(&mut line) {
		if line.starts_with('@') == false {
			error!("Invalid FASTQ format encountered.");
		}

		output.clear();
		output.push_str(&line);

		fasta_file.read_line(&mut seq);
		fasta_file.read_line(&mut line);
		fasta_file.read_line(&mut base_qualities);

		if seq.ends_with('\n') { seq.pop(); }
		if base_qualities.ends_with('\n') { base_qualities.pop(); }

		if seq.len() != base_qualities.len() {
			error!("Read sequence and base qualities are of different length.")
		}

		// TODO: Maybe printing one char at a time would be just as fast?
		for (base, qual) in seq.chars().zip(base_qualities.chars()) {
			// FIXME: We assume Sanger format base qualities here.
			output.push(if qual as u8 - 33u8 < min_baseq { 'N' } else { base });
		}
		write!(output, "\n+\n{}\n", base_qualities);
		print!("{}", output);
	}
}

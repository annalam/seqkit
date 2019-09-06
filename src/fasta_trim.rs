
use crate::common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta trim [options] <fastq_file>

Options:
  --first=N          Remove first N bases of each read [default: 0].
  --last=N           Remove last N bases of each read [default: 0].
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));
	let remove_first: usize = args.get_str("--first").parse()
		.unwrap_or_else(|_| error!("N must be a non-negative integer in --first=N."));
	let remove_last: usize = args.get_str("--last").parse()
		.unwrap_or_else(|_| error!("N must be a non-negative integer in --last=N."));

	let mut line = String::new();
	let mut seq = String::new();
	let mut qual = String::new();
	while fasta_file.read_line(&mut line) {
		if line.starts_with('>') == false && line.starts_with('@') == false {
			error!("Invalid FASTA/FASTQ format encountered.");
		}

		fasta_file.read_line(&mut seq);
		let seq_len = seq.trim_end().len();
		if remove_first + remove_last < seq_len {
			print!("{}{}\n", line, &seq[remove_first..seq_len-remove_last]);
		} else {
			print!("{}\n", line);
		}

		if line.starts_with('@') {
			fasta_file.read_line(&mut line);
			fasta_file.read_line(&mut qual);
			if remove_first + remove_last < seq_len {
				print!("+\n{}\n", &qual[remove_first..seq_len-remove_last]);
			} else {
				print!("+\n\n");
			}
		}
	}
}

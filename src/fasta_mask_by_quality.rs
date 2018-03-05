
use common::{parse_args, read_buffered, AsciiBufRead};
use std::str;
use std::process::exit;
use ascii::{AsciiString, AsciiChar};

const USAGE: &'static str = "
Usage:
  fasta mask by quality <fastq_file> <min_baseq>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = read_buffered(&args.get_str("<fastq_file>"));
	let min_baseq: i32 = args.get_str("<min_baseq>").parse().unwrap();

	let mut line = AsciiString::new();
	while fasta_file.read_ascii_line(&mut line) {
		if line[0] != '@' {
			eprintln!("Invalid FASTQ format encountered."); exit(-1);
		}
		print!("{}", line);
		fasta_file.read_ascii_line(&mut line);
		let mut seq = line.clone();
		fasta_file.read_ascii_line(&mut line);
		fasta_file.read_ascii_line(&mut line);
		let qual = &line;

		for k in 0..seq.len() {
			// FIXME: We assume Sanger format base qualities here.
			if (qual[k] as i32 - 33 as i32) < min_baseq && seq[k] != '\n' {
				seq[k] = AsciiChar::N;
			}
		}

		print!("{}+\n{}", seq, qual);
	}
}

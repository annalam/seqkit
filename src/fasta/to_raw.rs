
use parse_args;
use {read_buffered, AsciiBufRead};
use std::str;
use std::process::exit;
use ascii::AsciiString;

const USAGE: &'static str = "
Usage:
  fasta to raw <fasta_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = read_buffered(&args.get_str("<fasta_file>"));

	let mut line = AsciiString::new();
	while fasta_file.read_ascii_line(&mut line) {
		if line[0] == '>' {
			fasta_file.read_ascii_line(&mut line);
			print!("{}", line);
		} else if line[0] == '@' {
			fasta_file.read_ascii_line(&mut line);
			print!("{}", line);
			// Discard the per-base qualities, if present in the file
			fasta_file.read_ascii_line(&mut line);
			fasta_file.read_ascii_line(&mut line);
		} else {
			eprintln!("Invalid FASTA/FASTQ format encountered."); exit(-1);
		}
	}
}

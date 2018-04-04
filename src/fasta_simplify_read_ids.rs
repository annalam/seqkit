
use common::{parse_args, FileReader};
use std::str;
use std::process::exit;

const USAGE: &str = "
Usage:
  fasta simplify read ids <fastq_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));
	let mut num_reads = 0;

	let mut line = String::new();
	while fasta_file.read_line(&mut line) {
		num_reads += 1;
		if line.starts_with('@') {
			println!("@{}", num_reads);
			for _ in 0..3 {
				fasta_file.read_line(&mut line);
				print!("{}", line);
			}
		} else if line.starts_with('>') {
			println!(">{}", num_reads);
			fasta_file.read_line(&mut line);
			print!("{}", line);
		} else {
			error!("Invalid FASTA/FASTQ format encountered.");
		}
	}
}

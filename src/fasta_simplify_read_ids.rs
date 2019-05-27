
use crate::common::{parse_args, FileReader};
use std::str;
use regex::Regex;

const USAGE: &str = "
Usage:
  fasta simplify read ids [options] <fastq_file>

Options:
  --alphanumeric     Use letters a-z, A-Z and 0-9 in read identifiers
  --discard-umi      Remove \"UMI:\" tags from read identifiers, if present
";


const ALPHANUMERIC: &str =
	"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

fn alphanumeric_id(read_num: usize) -> String { "".into() }

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));
	//let alphanumeric = args.get_bool("--alphanumeric");
	let discard_umi = args.get_bool("--discard-umi");

	let umi_regex = Regex::new(r" UMI:[^\s]*").unwrap();

	let mut read_num: usize = 0;
	let mut line = String::new();

	while fasta_file.read_line(&mut line) {
		if line.is_empty() { continue; }

		let prefix = line.chars().nth(0).unwrap();
		if prefix != '@' && prefix != '>' {
			error!("Invalid FASTA/FASTQ format encountered.");
		}

		read_num += 1;
		print!("{}{}", prefix, read_num);

		// Try to extract UMI from the read identifier
		if discard_umi == false {
			if let Some(umi) = umi_regex.find(&line) {
				print!("{}", umi.as_str());
			}
		}
		println!();

		// Output the read sequence
		fasta_file.read_line(&mut line);
		print!("{}", line);

		// For FASTQ files, also output the two per-base quality lines
		if prefix == '@' {
			fasta_file.read_line(&mut line);
			println!("+");
			fasta_file.read_line(&mut line);
			print!("{}", line);
		}
	}
}

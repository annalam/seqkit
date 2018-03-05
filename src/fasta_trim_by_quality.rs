
use common::{parse_args, read_buffered, AsciiBufRead};
use std::str;
use std::process::exit;
use ascii::AsciiString;

const USAGE: &'static str = "
Usage:
  fasta trim by quality <fastq_file> <min_baseq>
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
		let seq = line.clone();
		fasta_file.read_ascii_line(&mut line);
		fasta_file.read_ascii_line(&mut line);
		let qual = &line;

		let mut total: i32 = -50;
		let mut lowest_total: i32 = total;

		let mut k: usize = qual.len() - 1;
		if qual[k] == '\n' { k -= 1; }
		let mut lowest_k = k + 1;
		loop {
			// FIXME: We assume Sanger format base qualities here.
			total += (qual[k] as i32 - 33 as i32) - min_baseq;
			if total > 0 { break; }
			if total < lowest_total {
				lowest_total = total;
				lowest_k = k;
			}
			if k == 0 { break; }
			k -= 1;
		}

		if lowest_k == 0 {
			print!("N\n+\n!\n");
		} else {
			print!("{}\n+\n{}\n", &seq[..lowest_k], &qual[..lowest_k]);
		}

	}
}

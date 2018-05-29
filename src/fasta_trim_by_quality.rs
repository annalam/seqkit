
use common::{parse_args, FileReader};
use std::str;
use std::process::exit;

const USAGE: &str = "
Usage:
  fasta trim by quality <fastq_file> <min_baseq>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));
	let min_baseq: u8 = args.get_str("<min_baseq>").parse().unwrap();
	let baseq_offset = 33u8;  // TODO: Support other base quality ASCII offsets

	let mut line = String::new();
	let mut seq = String::new();
	let mut qual = String::new();
	while fasta_file.read_line(&mut line) {
		if !line.starts_with('@') {
			error!("Invalid FASTQ format encountered.");
		}
		print!("{}", line);
		fasta_file.read_line(&mut seq);
		fasta_file.read_line(&mut line);
		fasta_file.read_line(&mut qual);

		let mut total: i32 = -50;
		let mut lowest_total: i32 = total;

		let mut k: usize = qual.trim_right().len();
		let mut lowest_k = k;
		while k > 0 {
			k -= 1;
			total += (qual.as_bytes()[k] - baseq_offset) as i32 -
				min_baseq as i32;
			if total > 0 { break; }
			if total < lowest_total {
				lowest_total = total;
				lowest_k = k;
			}
		}

		if lowest_k == 0 {
			print!("N\n+\n!\n");   // Entire read was garbage
		} else {
			print!("{}\n+\n{}\n", &seq[..lowest_k], &qual[..lowest_k]);
		}
	}
}

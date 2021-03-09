
use crate::common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  fasta extract dual umi [options] <interleaved_fastq>

Options:
  --first-bases=N   First N bases of read contain UMI bases [default: 0]
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<interleaved_fastq>"));
	let first_bases: usize = args.get_str("--first-bases").parse()
		.unwrap_or_else(|_|
		error!("N must be a non-negative integer in --first-bases=N."));

	let mut header_1 = String::new();
	let mut header_2 = String::new();
	let mut seq_1 = String::new();
	let mut seq_2 = String::new();
	let mut qual_1 = String::new();
	let mut qual_2 = String::new();
	let mut line = String::new();
	let mut umi = String::new();

	while fastq.read_line(&mut header_1) {
		umi.clear();

		let fastq_format = if header_1.starts_with('@') { true }
			else if header_1.starts_with('>') { false }
			else { error!("Header is not valid FASTA/FASTQ:\n{}", header_1); };

		if fastq_format {
			fastq.read_line(&mut seq_1);
			fastq.read_line(&mut line);
			fastq.read_line(&mut qual_1);
			fastq.read_line(&mut header_2);
			fastq.read_line(&mut seq_2);
			fastq.read_line(&mut line);
			fastq.read_line(&mut qual_2);
			if header_2.starts_with('@') == false {
				error!("Invalid FASTQ record found in input file.");
			}
		} else {
			fastq.read_line(&mut seq_1);
			fastq.read_line(&mut header_2);
			fastq.read_line(&mut seq_2);
			if header_2.starts_with('>') == false {
				error!("Invalid FASTA record found in input file.");
			}
		}

		umi += &seq_1[0..first_bases];
		umi += "+";
		umi += &seq_2[0..first_bases];

		if fastq_format {
			print!("{} RX:{}\n{}+\n{}{} RX:{}\n{}+\n{}",
				header_1.trim_end(), &umi, &seq_1[first_bases..],
				&qual_1[first_bases..], header_2.trim_end(), &umi,
				&seq_2[first_bases..], &qual_2[first_bases..])
		} else {
			print!("{} RX:{}\n{}{} RX:{}\n{}",
				header_1.trim_end(), &umi, &seq_1[first_bases..],
				header_2.trim_end(), &umi, &seq_2[first_bases..])
		}
	}
}

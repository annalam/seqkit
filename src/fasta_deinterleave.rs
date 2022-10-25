
use crate::common::{parse_args, FileReader, Compressor, GzipWriter};
use std::{str, io::Write};

const USAGE: &str = "
Usage:
  fasta deinterleave <interleaved_fastq> <out_prefix>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<interleaved_fastq>"));
	let out_prefix = args.get_str("<out_prefix>");
	let mut out_1 = GzipWriter::with_method(&format!("{}_1.fq.gz", &out_prefix),
		Compressor::GZIP);
	let mut out_2 = GzipWriter::with_method(&format!("{}_2.fq.gz", &out_prefix),
		Compressor::GZIP);

	let mut line = String::new();
	while fastq.read_line(&mut line) {
		let lines = if line.starts_with('@') { 4 }
			else if line.starts_with('>') { 2 }
			else { error!("Line is not FASTA/FASTQ format: {}", line); };
		write!(out_1, "{}", line);
		for k in 0..lines-1 {
			fastq.read_line(&mut line); write!(out_1, "{}", line);
		}

		fastq.read_line(&mut line);
		if (lines == 4 && !line.starts_with('@')) ||
			(lines == 2 && !line.starts_with('>')) {
			error!("Interleaved FASTA records are not in consistent format.");
		}
		write!(out_2, "{}", line);
		for k in 0..lines-1 {
			fastq.read_line(&mut line); write!(out_2, "{}", line);
		}
	}
}

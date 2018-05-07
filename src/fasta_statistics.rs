
use common::{parse_args, FileReader, GzipWriter};
use std::io::Write;
use std::str;
use std::collections::HashMap;
use regex::Regex;

const USAGE: &str = "
Usage:
  fasta statistics <fastq_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut fastq = FileReader::new(&args.get_str("<fastq_file>"));

	let sample_barcode_regex = Regex::new(r" SI:[ACGTNacgtn]+").unwrap();

	let mut total_records: u64 = 0;
	let mut sample_barcodes: HashMap<String, u64> = HashMap::new();

	let mut line = String::new();
	while fastq.read_line(&mut line) {
		// Find the sample barcode, formatted as SI:xxxx.
		if let Some(hit) = sample_barcode_regex.find(&line) {
			let sample_barcode = line[hit.start()+4..hit.end()].to_string();
			*sample_barcodes.entry(sample_barcode).or_insert(0) += 1;
		}

		// Read the non-header lines of this FASTQ record.
		if line.starts_with('@') {
			for _ in 0..3 { fastq.read_line(&mut line); }
		} else if line.starts_with('>') {
			fastq.read_line(&mut line);
		} else {
			error!("Invalid FASTQ header:\n{}", line);
		}

		total_records += 1;
	}

	println!("Total sequence records: {}", total_records);

	println!("Most frequent sample barcodes:");
	let mut entries: Vec<(String, u64)> =
		sample_barcodes.into_iter().collect();
	entries.sort_by_key(|x| x.1);
	entries.reverse();
	for (barcode, count) in entries {
		println!("- {}: {}", barcode, count);
	}
}

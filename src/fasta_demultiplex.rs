
use common::{parse_args, FileReader, GzipWriter, Compressor};
use std::io::Write;
use std::str;
use regex::Regex;

const USAGE: &str = "
Usage:
  fasta demultiplex [options] <sample_sheet> <fastq_1> <fastq_2>

Options:
  --parallel   Use pigz (parallel gzip) for compression

Description:
Splits a pooled FASTQ file into multiple individual FASTQ files, based on a
sample sheet. Each read in the pooled FASTQ file must carry a BC:xxxxxxxx
field in its header.
";

struct Sample {
	name: String,
	barcode: String,
	output1: GzipWriter,
	output2: GzipWriter
}

pub fn main() {
	let args = parse_args(USAGE);
	let mut sample_sheet = FileReader::new(&args.get_str("<sample_sheet>"));
	let mut fastq1 = FileReader::new(&args.get_str("<fastq_1>"));
	let mut fastq2 = FileReader::new(&args.get_str("<fastq_2>"));
	let parallel = args.get_bool("--parallel");

	let barcode_regex = Regex::new(r" BC:[ACGTNacgtn]+").unwrap();

	// Read the user-provided sample sheet into memory.
	let mut samples = Vec::new();
	let mut line = String::new();
	let mut barcode_len: usize = 0;
	while sample_sheet.read_line(&mut line) {
		if line.starts_with('#') { continue; }
		let cols: Vec<&str> = line.trim().split('\t').collect();
		if cols.len() < 2 { continue; }
		let name = cols[0];
		if cols[1].is_empty() { error!("Sample {} has no barcode.", name); }
		if barcode_len == 0 {
			barcode_len = cols[1].len();
		} else if cols[1].len() != barcode_len {
			error!("Barcodes in sample sheet must be of same length.");
		}
		let method = if parallel { Compressor::PIGZ } else { Compressor::GZIP };
		samples.push(Sample {
			name: name.into(),
			barcode: cols[1].into(),
			output1: GzipWriter::with_method(&format!("{}_1.fq.gz", name), method),
			output2: GzipWriter::with_method(&format!("{}_2.fq.gz", name), method)
		});
	}

	// Check that sample sheet did not contain samples with identical names.
	for s in 0..samples.len() {
		for k in s+1..samples.len() {
			if samples[s].name == samples[k].name {
				error!("Sample {} is listed multiple times in sample sheet.", samples[s].name);
			}
		}
	}

	// TODO: Check that UMI-including barcodes do not clash with
	// UMI-less barcodes in hybrid runs.

	// TODO: Build a lookup table for (barcode -> sample) mappings.
	/*let mut barcode_to_sample = vec![0; 4.pow(barcode_len)];
	for b in 0..4.pow(barcode_len) {
		for base in 0..barcode_len {
			
		}
	}*/

	let mut total_reads: u64 = 0;
	let mut identified_reads: u64 = 0;

	let mut line1 = String::new();
	let mut line2 = String::new();
	let mut barcode = String::new();
	let mut umi = String::new();
	while fastq1.read_line(&mut line1) && fastq2.read_line(&mut line2) {
		if line1.starts_with('@') && line2.starts_with('@') {

			// Find the sample barcode, formatted as BC:xxxx.
			// Then remove it from the header.
			// TODO: Allow user to choose if BC:xxxx should be removed.
			let (start, end) = {
				let hit = &barcode_regex.find(&line1)
					.unwrap_or_else(|| error!("No BC:xxxx field found."));
				(hit.start(), hit.end())
			};
			barcode.clear();
			barcode += &line1[(start+4)..end];
			line1.drain(start..end);
			if barcode.len() != barcode_len {
				error!("Barcode {} is of different length than barcodes in the sample sheet.", barcode);
			}

			// Remove BC:xxxx field from the second mate as well (if present).
			let (start, end) = if let Some(hit) = barcode_regex.find(&line2) {
				(hit.start(), hit.end())
			} else {
				(0, 0)
			};
			if end > 0 { line2.drain(start..end); }

			if let Some(sample) = samples.iter_mut()
				.find(|s| barcode_matches(&barcode, &s.barcode)) {

				// Extract UMI, if present
				umi.clear();
				for x in sample.barcode.chars().zip(barcode.chars()) {
					if x.0 == 'U' { umi.push(x.1); }
				}

				write!(sample.output1, "{}", line1.trim_right());
				if !umi.is_empty() { write!(sample.output1, " UMI:{}", umi); }
				write!(sample.output1, "\n");

				write!(sample.output2, "{}", line2.trim_right());
				if !umi.is_empty() { write!(sample.output2, " UMI:{}", umi); }
				write!(sample.output2, "\n");

				for _ in 0..3 {
					fastq1.read_line(&mut line1);
					fastq2.read_line(&mut line2);
					write!(sample.output1, "{}", line1);
					write!(sample.output2, "{}", line2);
				}
				identified_reads += 2;
			} else {
				// Must read all four lines even if we did not recognize the
				// barcode.
				for _ in 0..3 { 
					fastq1.read_line(&mut line1);
					fastq2.read_line(&mut line2);
				}
			}
			total_reads += 2;
		} else {
			error!("Invalid FASTQ lines:\n{}{}", line1, line2);
		}
	}

	eprintln!("{} / {} ({:.1}%) reads carried a barcode matching one of the provided samples.", identified_reads, total_reads,
		(identified_reads as f64) / (total_reads as f64) * 100.0);
}

fn barcode_matches(observed: &str, candidate: &str) -> bool {
	observed.chars().zip(candidate.chars()).all(
		|x| x.0 == x.1 || x.1 == 'N' || x.1 == 'U')
}

/*
fn barcode_to_integer(barcode: &str) -> usize {

}

fn integer_to_barcode(integer: usize) -> 
*/
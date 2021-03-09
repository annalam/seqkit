
use crate::common::{parse_args, FileReader, GzipWriter, Compressor};
use std::io::Write;
use std::str;
use regex::Regex;

const USAGE: &str = "
Usage:
  fasta demultiplex [options] <sample_sheet> <fastq_1> [<fastq_2>]

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
	output: Vec<GzipWriter>
}

pub fn main() {
	let args = parse_args(USAGE);
	let parallel = args.get_bool("--parallel");
	let barcode_regex = Regex::new(r" BC:[ACGTNacgtn+]+").unwrap();

	// Initialize the FASTQ readers
	let mut fastq: Vec<FileReader> = Vec::new();
	fastq.push(FileReader::new(&args.get_str("<fastq_1>")));
	if args.get_str("<fastq_2>").is_empty() == false {
		fastq.push(FileReader::new(&args.get_str("<fastq_2>")));
	}
	let paired_end = fastq.len() == 2;

	// Read the user-provided sample sheet into memory.
	eprintln!("Reading sample sheet...");
	let mut sample_sheet = FileReader::new(&args.get_str("<sample_sheet>"));
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

		let mut outputs: Vec<GzipWriter> = Vec::new();
		if paired_end {
			outputs.push(GzipWriter::with_method(
				&format!("{}_1.fq.gz", name), method));
			outputs.push(GzipWriter::with_method(
				&format!("{}_2.fq.gz", name), method));
		} else {
			outputs.push(GzipWriter::with_method(
				&format!("{}.fq.gz", name), method));
		}

		samples.push(Sample {
			name: name.into(),
			barcode: cols[1].into(),
			output: outputs
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
	// UMI-less barcodes in runs that contain both types of barcodes.

	// TODO: Build a lookup table for (barcode -> sample) mappings.
	/*let mut barcode_to_sample = vec![0; 4.pow(barcode_len)];
	for b in 0..4.pow(barcode_len) {
		for base in 0..barcode_len {
			
		}
	}*/

	eprintln!("Starting demultiplexing in {} end mode...",
		if paired_end { "paired" } else { "single" });
	let mut total_reads: u64 = 0;
	let mut identified_reads: u64 = 0;

	let mut line = String::new();
	let mut barcode = String::new();
	let mut umi = String::new();

	while fastq[0].read_line(&mut line) {
		if !line.starts_with('@') {
			error!("Invalid FASTQ header line:\n{}", line);
		}

		// Find the sample barcode, formatted as BC:xxxx.
		// Then remove it from the header.
		// TODO: Extract the BC:xxx search into a function.
		let (start, end) = {
			let hit = &barcode_regex.find(&line)
				.unwrap_or_else(|| error!("No BC:xxxx field found."));
			(hit.start(), hit.end())
		};
		barcode.clear();
		barcode += &line[(start+4)..end];
		line.drain(start..end);
		if barcode.len() != barcode_len {
			error!("Barcode {} is of different length than barcodes in the sample sheet.", barcode);
		}

		if let Some(sample) = samples.iter_mut()
			.find(|s| barcode_matches(&barcode, &s.barcode)) {

			// Extract UMI, if present
			umi.clear();
			for x in sample.barcode.chars().zip(barcode.chars()) {
				if x.0 == 'U' { umi.push(x.1); }
			}

			// Write the first mate into the correct output FASTQ file
			write!(sample.output[0], "{}", line.trim_right());
			if !umi.is_empty() { write!(sample.output[0], " UMI:{}", umi); }
			write!(sample.output[0], "\n");
			for _ in 0..3 {
				fastq[0].read_line(&mut line);
				write!(sample.output[0], "{}", line);
			}
			identified_reads += 1;

			// Handle the second mate (if present)
			if paired_end {
				fastq[1].read_line(&mut line);

				// Remove BC:xxx field
				let (start, end) = if let Some(hit) =
					barcode_regex.find(&line) {
					(hit.start(), hit.end())
				} else {
					(0, 0)
				};
				if end > 0 { line.drain(start..end); }

				write!(sample.output[1], "{}", line.trim_right());
				if !umi.is_empty() {
					write!(sample.output[1], " UMI:{}", umi);
				}
				write!(sample.output[1], "\n");
				for _ in 0..3 {
					fastq[1].read_line(&mut line);
					write!(sample.output[1], "{}", line);
				}
				identified_reads += 1;
			}
		} else {
			// Must read all four lines even if we did not recognize the
			// barcode.
			for _ in 0..3 { fastq[0].read_line(&mut line); }
			if paired_end {
				for _ in 0..4 { fastq[1].read_line(&mut line); }
			}
		}

		total_reads += if paired_end { 2 } else { 1 };
	}

	eprintln!("{} / {} ({:.1}%) reads carried a barcode matching one of the provided samples.", identified_reads, total_reads,
		(identified_reads as f64) / (total_reads as f64) * 100.0);
}

fn barcode_matches(observed: &str, candidate: &str) -> bool {
	observed.chars().zip(candidate.chars()).all(
		|x| x.0 == x.1 || x.1 == 'N' || x.1 == 'U')
}


use crate::common::{parse_args, FileReader, GzipWriter, Compressor};
use std::io::Write;
use std::collections::HashMap;
use std::str;
use regex::Regex;

const USAGE: &str = "
Usage:
  fasta demultiplex [options] <sample_sheet> <fastq_1> [<fastq_2>]

Options:
  --parallel      Use pigz (parallel gzip) for compression
  --index1=FASTQ  Path to FASTQ file containing the first index (optional)
  --index2=FASTQ  Path to FASTQ file containing the second index (optional)
  --dry-run=N     Analyze N reads and generate table of indexes found in the run

Splits a pooled FASTQ file into multiple individual FASTQ files, based on a
sample sheet. Each read in the pooled FASTQ file must carry a BC:xxxxxxxx
field in its header.
";

struct Sample {
	name: String,
	barcode: String,
	output: Vec<GzipWriter>,
	total_reads: u64
}

pub fn main() {
	let args = parse_args(USAGE);
	let parallel = args.get_bool("--parallel");
	let dry_run: u64 = args.get_str("--dry-run").parse().unwrap_or(0);
	if dry_run == 0 && args.get_str("--dry-run").is_empty() == false {
		error!("In --dry-run=N, N must be 64-bit positive integer.");
	}

	let barcode_regex = Regex::new(r" BC:[ACGTNacgtn+]+").unwrap();

	// Initialize the FASTQ readers
	let mut fastq: Vec<FileReader> = Vec::new();
	fastq.push(FileReader::new(&args.get_str("<fastq_1>")));
	if args.get_str("<fastq_2>").is_empty() == false {
		fastq.push(FileReader::new(&args.get_str("<fastq_2>")));
	}
	let paired_end = fastq.len() == 2;

	// Initialize index readers (if the indexes are in separate FASTQ files)
	let mut index_fastq: Vec<FileReader> = Vec::new();
	if args.get_str("--index1").is_empty() == false {
		index_fastq.push(FileReader::new(&args.get_str("--index1")));
	}
	if args.get_str("--index2").is_empty() == false {
		index_fastq.push(FileReader::new(&args.get_str("--index2")));
	}

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
			error!("Barcodes in sample sheet must all be of same length.");
		}
		let method = if parallel { Compressor::PIGZ } else { Compressor::GZIP };

		let mut outputs: Vec<GzipWriter> = Vec::new();
		if dry_run > 0 {
			// No output will be produced
		} else if paired_end {
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
			output: outputs,
			total_reads: 0
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

	eprintln!("Starting demultiplexing in {} end mode...",
		if paired_end { "paired" } else { "single" });
	let mut total_reads: u64 = 0;
	let mut identified_reads: u64 = 0;
	let mut extra_barcodes: HashMap<String, u64> = HashMap::new();

	let mut header = String::new();
	let mut line = String::new();
	let mut barcode = String::new();
	let mut umi = String::new();

	while fastq[0].read_line(&mut header) {
		if !header.starts_with('@') {
			error!("Invalid FASTQ header line:\n{}", header);
		}
		

		barcode.clear();

		// Read indexes from separate FASTQ files if requested
		if index_fastq.is_empty() == false {
			for ifq in &mut index_fastq {
				if barcode.is_empty() == false { barcode += "+"; }
				ifq.read_line(&mut line);
				assert!(line.starts_with('@'));
				ifq.read_line(&mut line);
				barcode += line.trim_end();
				ifq.read_line(&mut line);
				assert!(line.starts_with('+'));
				ifq.read_line(&mut line);
			}
		} else {
			// Search for BC:xxxx format index in the FASTQ header.
			let (start, end) = {
				let hit = &barcode_regex.find(&header)
					.unwrap_or_else(|| error!("No BC:xxxx field found."));
				(hit.start(), hit.end())
			};
			barcode += &header[(start+4)..end];
			header.drain(start..end);   // Remove BC:xxxx from header
		}

		if barcode.len() != barcode_len {
			error!("Sequenced barcode {} is of different length ({} nt) than barcodes in the sample sheet ({} nt).", barcode, barcode.len(), barcode_len);
		}

		// Find the sample sheet barcode best matching the sequenced barcode,
		// allowing for one base mismatch.
		let mut best_sample: usize = 0;
		let mut equally_fine_sample: usize = 0;
		let mut lowest_diff = usize::MAX;
		for s in 0..samples.len() {
			let diff = barcode_diff(&barcode.as_bytes(), &samples[s].barcode.as_bytes());
			if diff < lowest_diff {
				lowest_diff = diff;
				best_sample = s;
				equally_fine_sample = s;
			} else if diff == lowest_diff {
				equally_fine_sample = s;   // We have two equally good samples
			}
		}

		let MAX_BARCODE_DIFFERENCE: usize = 1;
		total_reads += 1;
		let mut write_read_out = false;

		if lowest_diff <= MAX_BARCODE_DIFFERENCE {
			if best_sample == equally_fine_sample {
				// If one sample was unambiguously the closest match to the
				// sequenced barcode, we assign the read to that sample.
				let sample = &mut samples[best_sample];
				identified_reads += 1;
				sample.total_reads += 1;
				write_read_out = !(dry_run > 0);

			} else {
				// If multiple sample sheet barcodes matched equally well,
				// report the ambiguity to the user.	
				eprintln!("WARNING: Sequenced barcode {} was an equally good match ({} mismatches) for samples {} ({}) and {} ({}), and was therefore not assigned to any sample.",
					&barcode, lowest_diff, &samples[best_sample].name,
					&samples[best_sample].barcode,
					&samples[equally_fine_sample].name,
					&samples[equally_fine_sample].barcode);
			}
		} else if dry_run > 0 {
			// If no good matches were found at all, and we are in dry run
			// mode, store the unmatched barcode.
			*extra_barcodes.entry(barcode.clone()).or_insert(0) += 1;
		}

		if write_read_out {
			let sample = &mut samples[best_sample];

			// Extract UMI, if present
			umi.clear();
			for x in sample.barcode.chars().zip(barcode.chars()) {
				if x.0 == 'U' { umi.push(x.1); }
			}

			// Write the first mate into the correct output FASTQ file
			write!(sample.output[0], "{}", header.trim_right());
			if !umi.is_empty() { write!(sample.output[0], " UMI:{}", umi); }
			write!(sample.output[0], "\n");
			for _ in 0..3 {
				fastq[0].read_line(&mut line);
				write!(sample.output[0], "{}", line);
			}

			// Handle the second mate (if present)
			if paired_end {
				fastq[1].read_line(&mut line);

				// Remove BC:xxx field, if present
				if index_fastq.is_empty() {
					let (start, end) = if let Some(hit) =
						barcode_regex.find(&line) {
						(hit.start(), hit.end())
					} else {
						(0, 0)
					};
					if end > 0 { line.drain(start..end); }
				}

				write!(sample.output[1], "{}", line.trim_right());
				if !umi.is_empty() {
					write!(sample.output[1], " UMI:{}", umi);
				}
				write!(sample.output[1], "\n");
				for _ in 0..3 {
					fastq[1].read_line(&mut line);
					write!(sample.output[1], "{}", line);
				}
			}
		} else {
			// Must read all four lines even if we did not recognize the
			// barcode.
			for _ in 0..3 { fastq[0].read_line(&mut line); }
			if paired_end {
				for _ in 0..4 { fastq[1].read_line(&mut line); }
			}
		}

		if dry_run > 0 && total_reads >= dry_run { break; };
	}

	if dry_run > 0 {
		eprintln!("Dry run completed with {} clusters. Barcodes found:",
			total_reads);
		let mut entries: Vec<(String, u64)> = samples.iter().map(|s| (s.name.clone(), s.total_reads)).collect();
		entries.extend(extra_barcodes.into_iter());
		entries.sort_by_key(|x| x.1);
		entries.reverse();
		for (barcode, count) in &entries[0..100] {
			println!("- {}: {}", barcode, count);
		}
	}

	eprintln!("{} / {} ({:.1}%) clusters carried a barcode matching one of the provided samples.", identified_reads, total_reads,
		(identified_reads as f64) / (total_reads as f64) * 100.0);
}

// Compares an observed barcode against a barcode found in the sample sheet,
// and returns the number of base positions where they differ.
fn barcode_diff(observed: &[u8], candidate: &[u8]) -> usize {
	assert!(observed.len() == candidate.len());
	let mut mismatches: usize = 0;
	for k in 0..observed.len() {
		if candidate[k] == b'N' || candidate[k] == b'U' { continue; }
		if observed[k] != candidate[k] { mismatches += 1; }
	}
	mismatches
}


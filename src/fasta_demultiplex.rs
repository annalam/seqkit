
use common::{parse_args, FileReader, GzipWriter};
use std::io::Write;
use std::str;

const USAGE: &str = "
Usage:
  fasta demultiplex <sample_sheet> <fastq_1> <fastq_2>
";

struct Sample {
	barcode: String,
	output1: GzipWriter,
	output2: GzipWriter
}

pub fn main() {
	let args = parse_args(USAGE);
	let mut sample_sheet = FileReader::new(&args.get_str("<sample_sheet>"));
	let mut fastq1 = FileReader::new(&args.get_str("<fastq_1>"));
	let mut fastq2 = FileReader::new(&args.get_str("<fastq_2>"));

	// Read the user-provided sample sheet into memory.
	let mut samples = Vec::new();
	let mut line = String::new();
	while sample_sheet.read_line(&mut line) {
		let cols: Vec<&str> = line.trim().split('\t').collect();
		if cols.len() < 2 { continue; }
		let name = cols[0];
		samples.push(Sample {
			barcode: cols[1].into(),
			output1: GzipWriter::new(&format!("{}_1.fq.gz", name)),
			output2: GzipWriter::new(&format!("{}_2.fq.gz", name))
		});
	}

	/*let barcode_len = 8;
	let mut barcode_to_sample = vec![0; 4.pow(8)];
	for b in 0..4.pow(8) {
		for base in 0..barcode_len {
			
		}
	}*/

	let mut total_reads: u64 = 0;
	let mut identified_reads: u64 = 0;

	let mut line1 = String::new();
	let mut line2 = String::new();
	while fastq1.read_line(&mut line1) && fastq2.read_line(&mut line2) {
		if line1.starts_with('@') && line2.starts_with('@') {
			let barcode = line1.split(':').last().unwrap().trim().to_string();
			if let Some(sample) = samples.iter_mut()
				.find(|s| barcode_matches(&barcode, &s.barcode)) {

				write!(sample.output1, "{}", line1);
				write!(sample.output2, "{}", line2);

				for _ in 0..3 {
					fastq1.read_line(&mut line1);
					fastq2.read_line(&mut line2);
					write!(sample.output1, "{}", line1);
					write!(sample.output2, "{}", line2);
				}
				identified_reads += 2;
			} else {
				// Must read all four lines even if we do not recognize the
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

fn barcode_matches(observed: &str, provided: &str) -> bool {
	assert!(observed.len() == provided.len());
	observed.chars().zip(provided.chars()).all(|x| x.0 == x.1 || x.1 == 'N')
}

/*
fn barcode_to_integer(barcode: &str) -> usize {

}

fn integer_to_barcode(integer: usize) -> 
*/

use crate::common::{parse_args, FileReader};
use std::str;
use std::collections::HashMap;
use bio::io::fasta;

const USAGE: &str = "
Usage:
  fasta gc content <genome.fa> <regions.bed>

Description:
Calculates the GC content percentage of FASTA file regions listed in the input
BED file. Ambiguous N nucleotides are omitted from both the numerator and the
denominator.
";

pub fn main() {
	let args = parse_args(USAGE);
	let fasta_path = args.get_str("<genome.fa>");
	let bed_path = args.get_str("<regions.bed>");

	eprintln!("Reading reference genome into memory...");
	let mut fasta = fasta::Reader::from_file(&fasta_path)
		.unwrap_or_else(|_| error!("Input FASTA file {} could not be read.", &fasta_path));
	let mut genome = HashMap::new();
	for entry in fasta.records() {
		let chr = entry.unwrap();
		genome.insert(chr.id().to_owned(), chr.seq().to_owned());
	}

	let mut regions_bed = FileReader::new(&bed_path);
	let mut line = String::new();
	while regions_bed.read_line(&mut line) {
		let cols: Vec<&str> = line.trim().split('\t').collect();
		if cols.len() < 3 {
			eprintln!("WARNING: Input BED file contains line with less than 3 columns:\n{}\n", &line);
		}

		if let Some(chr_seq) = genome.get(&cols[0] as &str) {
			let start: usize = cols[1].parse().unwrap_or_else(|_| error!("Invalid region:\n{}\n", &line));
			let stop: usize = cols[2].parse().unwrap_or_else(|_| error!("Invalid region:\n{}\n", &line));
			let seq = chr_seq.get(start..stop).unwrap_or_else(|| error!("Invalid region:\n{}\n", &line));

			// TODO: Check for ambiguous nucleotides other than N?
			let gc_nucs = seq.iter().filter(|&&b| b == b'C' || b == b'G' || b == b'c' || b == b'g').count();
			let total_nucs = seq.iter().filter(|&&b| b != b'N' && b != b'n').count();
			println!("{}\t{}\t{:.3}", gc_nucs, total_nucs, gc_nucs as f32 / total_nucs as f32);
		}
	}
}

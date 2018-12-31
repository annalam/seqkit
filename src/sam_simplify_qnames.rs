
use common::{parse_args, GzipWriter};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;

const USAGE: &str = "
Usage:
  sam simplify qnames <bam_file>

These commands convert BAM files into FASTQ, FASTA, or raw sequence-per-line
format. Both name-sorted and position-sorted BAM files are supported,
but memory usage can reach several GB for position-sorted BAM files.

Output is written into files whose name is derived based on output prefix
and format. For example, with output format FASTQ and prefix \"sample\",
paired end reads are written into files sample_1.fq.gz and sample_2.fq.gz,
and orphan reads are written into sample.fq.gz.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut bam = if bam_path == "-" {
		bam::Reader::from_stdin().unwrap()
	} else {
		bam::Reader::from_path(&bam_path).unwrap_or_else(
			|_| error!("Cannot open BAM file '{}'", bam_path))
	};

	let header = bam.header().clone();

	let mut out = bam::Writer::from_stdout(&bam::header::Header::from_template(&header)).unwrap();

	let mut reads_1: HashMap<Box<str>, u32> = HashMap::new();
	let mut reads_2: HashMap<Box<str>, u32> = HashMap::new();

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));
		if read.is_secondary() || read.is_supplementary() { continue; }

		let qname = str::from_utf8(read.qname()).unwrap();
		let mut read_seq = sequence(&read, 10);

		if read.is_paired() == false {
			write_read(&mut out_single, output_format, qname, &read_seq);
		} else if read.is_first_in_template() {
			if let Some(mate_seq) = reads_2.remove(qname) {
				write_read(&mut out_1, output_format, qname, &read_seq);
				write_read(&mut out_2, output_format, qname, &mate_seq);
			} else {
				reads_1.insert(qname.into(), read_seq.into());
			}
		} else if read.is_last_in_template() {
			if let Some(mate_seq) = reads_1.remove(qname) {
				write_read(&mut out_1, output_format, qname, &mate_seq);
				write_read(&mut out_2, output_format, qname, &read_seq);
			} else {
				reads_2.insert(qname.into(), read_seq.into());
			}
		}
	}

	// If we are left with any orphan reads for which a pair was not found,
	// we add those to the prefix.fq.gz output file.
	for (qname, seq) in reads_1.iter().chain(reads_2.iter()) {
		write_read(&mut out_single, output_format, qname, &seq);
	}
}


use common::{parse_args, GzipWriter};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;

const USAGE: &str = "
Usage:
  sam to fastq <bam_file> <out_prefix>
";

fn sequence(read: &Record) -> String {
	let seq = read.seq();
	let qual = read.qual();
	let mut ret = String::with_capacity(seq.len());
	for k in 0..seq.len() {
		if qual[k] < 10 {
			ret.push('N');
		} else {
			ret.push(match seq.encoded_base(k) {
				1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', _ => 'N'
			});
		}
	}
	ret
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let out_prefix = args.get_str("<out_prefix>");

	let mut bam = if bam_path == "-" {
		bam::Reader::from_stdin().unwrap()
	} else {
		bam::Reader::from_path(&bam_path).unwrap_or_else(
			|_| error!("Cannot open BAM file '{}'", bam_path))
	};

	let header = bam.header().clone();

	let mut fastq_1 = GzipWriter::new(&format!("{}_1.fa.gz", out_prefix));
	let mut fastq_2 = GzipWriter::new(&format!("{}_2.fa.gz", out_prefix));
	let mut fastq_single = GzipWriter::new(&format!("{}.fa.gz", out_prefix));

	let mut reads_1: HashMap<String, Box<str>> = HashMap::new();
	let mut reads_2: HashMap<String, Box<str>> = HashMap::new();

	let mut total_reads = 0;

	for r in bam.records() {
		let mut read = r.unwrap_or_else(
			|_| error!("Input BAM file ended prematurely."));
		if read.is_secondary() || read.is_supplementary() { continue; }

		let qname = str::from_utf8(read.qname()).unwrap();

		if read.is_paired() == false {
			total_reads += 1;
			write!(fastq_single, ">{}\n{}\n", total_reads, sequence(&read));
		} else if read.is_first_in_template() {
			if let Some(mate_seq) = reads_2.remove(qname) {
				total_reads += 1;
				write!(fastq_1, ">{}\n{}\n", total_reads, sequence(&read));
				write!(fastq_2, ">{}\n{}\n", total_reads, mate_seq);
			} else {
				reads_1.insert(qname.into(), sequence(&read).into());
			}
		} else if read.is_last_in_template() {
			if let Some(mate_seq) = reads_1.remove(qname) {
				total_reads += 1;
				write!(fastq_1, ">{}\n{}\n", total_reads, mate_seq);
				write!(fastq_2, ">{}\n{}\n", total_reads, sequence(&read));
			} else {
				reads_2.insert(qname.into(), sequence(&read).into());
			}
		}
	}
}

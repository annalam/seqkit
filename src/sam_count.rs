
use crate::common::parse_args;
use std::str;
use std::thread;
use std::fs::{File, remove_file};
use std::process::{Command, Stdio};
use std::io::{BufReader, BufWriter, BufRead, Write};
use rust_htslib::bam;
use rust_htslib::bam::Read;

const USAGE: &str = "
Usage:
  sam count [options] <bam_file> <regions.bed>

Options:
  --frac-inside=F   Minimum overlap between read and region [default: 0.0]
  --min-mapq=N      Only count reads with MAPQ ≥ threshold [default: 0]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>").to_string();
	let bed_path = args.get_str("<regions.bed>");
	let min_mapq: u8 = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be an integer between 0 - 255."));

	let bam = bam::Reader::from_path(&bam_path).unwrap();
	let mut chr_names: Vec<String> = Vec::new();
	for name in bam.header().target_names() {
		chr_names.push(str::from_utf8(name).unwrap().to_string());
	}
	let index_path = format!(".{}.idx", bam_path);
	let mut index_file = File::create(&index_path).unwrap();
	for c in 0..chr_names.len() {
		write!(&mut index_file, "{}\t{}\n", chr_names[c], bam.header().target_len(c as u32).unwrap());
	}

    let mut cmd = Command::new("bedtools");
    cmd.arg("coverage"); cmd.arg("-sorted"); cmd.arg("-counts");
    cmd.arg("-g"); cmd.arg(&index_path);
    cmd.arg("-a"); cmd.arg(&bed_path);
    cmd.arg("-b"); cmd.arg("stdin");
    let bedtools = cmd.stdin(Stdio::piped()).stdout(Stdio::piped())
        .spawn().expect("Could not start bedtools process.");

	let mut bedtools_in = BufWriter::new(bedtools.stdin.unwrap());
	let bedtools_out = BufReader::new(bedtools.stdout.unwrap());
	thread::spawn(move || {
		let mut bam = bam::Reader::from_path(&bam_path).unwrap();
		for r in bam.records() {
			let read = r.unwrap();
			if read.is_paired() == false { continue; }
			if read.is_unmapped() || read.is_mate_unmapped() { continue; }
			if read.is_duplicate() || read.is_secondary() { continue; }
			if read.is_supplementary() { continue; }
			if read.tid() != read.mtid() { continue; }
			if read.mapq() < min_mapq { continue; }

			// Only count each DNA fragment once.
			if read.pos() > read.mpos() || (read.pos() == read.mpos() && read.is_first_in_template() == false) { continue; }

			let frag_size = read.insert_size().abs();
			if frag_size > 5000 { continue; }

			write!(bedtools_in, "{}\t{}\t{}\n", chr_names[read.tid() as usize], read.pos(), read.pos() + frag_size);
	    }
    });

	for l in bedtools_out.lines() {
		let line = l.unwrap();
		println!("{}", line.split('\t').last().unwrap());
    }

    remove_file(index_path);
}

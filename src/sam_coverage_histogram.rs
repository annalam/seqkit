
use crate::common::parse_args;
use std::str;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead};

const USAGE: &str = "
Usage:
  sam coverage histogram [options] <bam_file>

Options:
  --regions=BED   Regions to calculate coverage in [default: everywhere]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let regions_bed = args.get_str("--regions");
	let max_coverage = 10_000;

	let mut cmd = Command::new("samtools");
	cmd.arg("depth").arg("-a");
	if regions_bed != "everywhere" {
		cmd.arg("-b").arg(regions_bed);
	}
	if bam_path == "-" {
		cmd.arg("-"); cmd.stdin(Stdio::inherit());
	} else {
		cmd.arg(bam_path);
	}
	cmd.stdout(Stdio::piped());
	let samtools = cmd.spawn().unwrap_or_else(
		|_| error!("Could not start 'samtools depth'."));

	let samtools_out = BufReader::new(samtools.stdout.unwrap());

	let mut hist = vec![0; max_coverage + 1];
	for l in samtools_out.lines() {
		let line = l.unwrap();
		let mut cols = line.split('\t');
		cols.next(); cols.next();
		let count: usize = cols.next().unwrap().parse().unwrap();
		if count >= hist.len() { continue; }
		hist[count] += 1;
	}

	for k in 0..hist.len() {
		println!("{}\t{}", k, hist[k]);
	}
}

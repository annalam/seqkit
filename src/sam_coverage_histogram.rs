
use crate::common::parse_args;
use std::str;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead};

const USAGE: &str = "
Usage:
  sam coverage histogram [options] <bam_file>

Options:
  --region=REGION   Region to calculate coverage in [default: everywhere]
  --regions=BED     BED file of regions to calculate coverage in
                    [default: everywhere]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let region = args.get_str("--region");
	let regions_bed = args.get_str("--regions");
	let max_coverage = 10_000;

	if region != "everywhere" && regions_bed != "everywhere" {
		error!("Only one of --region or --regions can be provided.");
	}

	let mut cmd = Command::new("samtools");
	cmd.arg("depth").arg("-a");
	if region != "everywhere" {
		cmd.arg("-r").arg(region);
	} else if regions_bed != "everywhere" {
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

	let mut hist = vec![0u64; max_coverage + 1];
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

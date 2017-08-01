
use parse_args;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;

const USAGE: &'static str = "
Usage:
  sam statistics [options] <bam_file>

Options:
  --target-regions=BED   Count on-target% for regions in BED file [optional]
";

struct Region {
	chrom: String,
	start: u64,
	end: u64
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut target_regions: Vec<Region> = Vec::new();
	if args.get_str("--target-regions").is_empty() == false {
		//let bed = buffered_reader(args.get_str("--target-regions"));
		
	}

	let mut total_reads: u64 = 0;
	let mut aligned_reads: u64 = 0;
	let mut duplicate_reads: u64 = 0;

	let bam = bam::Reader::from_path(&bam_path).unwrap();
	for r in bam.records() {
		let read = r.unwrap();
		if read.is_secondary() || read.is_supplementary() { continue; }
		total_reads += 1;
		if read.is_unmapped() { continue; }

		aligned_reads += 1;
		if read.is_duplicate() { duplicate_reads += 1; }

		// TODO: Check if read is on-target (i.e. inside the targeted regions)

    }

    println!("Total reads: {}", total_reads);
    println!("Aligned reads: {} ({:.1}% of all reads)", aligned_reads, aligned_reads as f64 / total_reads as f64 * 100.0);
    println!("Duplicate reads: {} ({:.1}% of aligned reads)", duplicate_reads, duplicate_reads as f64 / aligned_reads as f64 * 100.0);
}

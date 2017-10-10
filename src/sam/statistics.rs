
use parse_args;
use ErrorHelper;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::bed;

const USAGE: &'static str = "
Usage:
  sam statistics [options] <bam_file>

Options:
  --on-target=BED   Count on-target% for regions in BED file [optional]
";

struct Region {
	chr: i32,    // Numeric chromosome ID, as used in the BAM file
	start: i32,
	end: i32
}

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let bam = bam::Reader::from_path(&bam_path)
		.on_error("Could not open BAM file.");
	let bam_header = bam.header();

	let mut target_regions: Vec<Region> = Vec::new();
	if args.get_str("--target-regions").is_empty() == false {
		let mut bed = bed::Reader::from_file(args.get_str("--on-target"))
			.on_error("Cannot open BED file.");
		for r in bed.records() {
			let record = r.unwrap();
			target_regions.push(Region {
				chr: bam_header.tid(record.chrom().as_bytes()).unwrap() as i32,
				start: (record.start() + 1) as i32, end: record.end() as i32
			});
		}
	}

	let mut total_reads: u64 = 0;
	let mut aligned_reads: u64 = 0;
	let mut duplicate_reads: u64 = 0;
	let mut on_target_reads: u64 = 0;

	for r in bam.records() {
		let read = r.unwrap();
		if read.is_secondary() || read.is_supplementary() { continue; }
		total_reads += 1;
		if read.is_unmapped() { continue; }

		aligned_reads += 1;
		if read.is_duplicate() { duplicate_reads += 1; }

		if target_regions.is_empty() == false {
			// TODO: Check if read overlaps one of the target regions
			let start = read.pos() + 1;
			let end = read.cigar().end_pos().unwrap() + 1;
			for region in &target_regions {
				on_target_reads += (read.tid() == region.chr && start <= region.end && end >= region.start) as u64;
			}
		}
    }

    println!("Total reads: {}", total_reads);
    println!("Aligned reads: {} ({:.1}% of all reads)", aligned_reads, aligned_reads as f64 / total_reads as f64 * 100.0);
    println!("Duplicate reads: {} ({:.1}% of aligned reads)", duplicate_reads, duplicate_reads as f64 / aligned_reads as f64 * 100.0);

    if target_regions.is_empty() == false {
    	println!("On-target reads: {} ({:.1}% of aligned reads)", on_target_reads, on_target_reads as f64 / aligned_reads as f64 * 100.0);
    }
}

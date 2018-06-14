
use common::{parse_args, PathArgs};
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::bed;

const USAGE: &str = "
Usage:
  sam statistics [options] <bam_file>

Options:
  --on-target=BED   Count on-target% for regions in BED file [optional]
";

struct Region { start: i32, end: i32 }

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_path("<bam_file>");
	let targets_path = args.get_path("--on-target");

	let mut bam = bam::Reader::from_path(&bam_path)
		.unwrap_or_else(|_| error!("Could not open BAM file {}.", &bam_path));
	let bam_header = bam.header().clone();

	// Each chromosome has its own vector of target regions. Chromosomes
	// are identified by their numeric IDs in the BAM file.
	let mut target_regions: Vec<Vec<Region>> = Vec::new();
	if targets_path.is_empty() == false {
		eprintln!("Reading target regions into memory...");
		for _ in 0..bam_header.target_count() {
			target_regions.push(Vec::new());
		}

		let mut bed = bed::Reader::from_file(&targets_path).unwrap_or_else(
			|_| error!("Could not open BED file {}.", targets_path));
		for r in bed.records() {
			let record = r.unwrap();
			let chr = record.chrom();
			let tid = bam_header.tid(chr.as_bytes()).unwrap_or_else(||
				error!("Chromosome {} is listed in target region BED file, but is not found in BAM file.", chr));
			target_regions[tid as usize].push(Region { 
				start: (record.start() + 1) as i32, end: record.end() as i32
			});
		}

		// Sort the regions by start coordinate so we can search faster
		for tid in 0..target_regions.len() {
			target_regions[tid].sort_unstable_by_key(|r| r.start);
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
			// TODO: Handle spliced regions in reads
			let start = read.pos() + 1;
			let end = read.cigar().end_pos().unwrap() + 1;
			for region in &target_regions[read.tid() as usize] {
				if start <= region.end && end >= region.start {
					on_target_reads += 1;
					break;
				}

				// Can stop searching since later regions cannot overlap
				// with the alignment.
				if region.start > end { break; }
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

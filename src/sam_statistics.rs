
use crate::common::{parse_args, PathArgs, FileReader, BamReader};
use std::str;

const USAGE: &str = "
Usage:
  sam statistics [options] <bam_file>

Options:
  --on-target=BED   Count on-target% for regions in BED file [optional]
";

struct Region { start: i64, end: i64 }

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_path("<bam_file>");
	let targets_path = args.get_path("--on-target");
	let max_frag_len = 5000;

	let mut bam = BamReader::open(&bam_path);
	let bam_header = bam.header();

	// Each chromosome has its own vector of target regions. Chromosomes
	// are identified by their numeric IDs in the BAM file.
	let mut target_regions: Vec<Vec<Region>> = Vec::new();
	if targets_path.is_empty() == false {
		eprintln!("Reading target regions into memory...");
		for _ in 0..bam_header.target_count() {
			target_regions.push(Vec::new());
		}


		let mut bed_file = FileReader::new(&targets_path);
		let mut line = String::new();
		while bed_file.read_line(&mut line) {
			if line.trim().is_empty() || line.starts_with('#') { continue; }
			let cols: Vec<&str> = line.trim().split('\t').collect();
			if cols.len() < 3 {
				error!("Invalid line in BED file {}:\n{}", targets_path, line);
			}
			let tid = bam_header.tid(cols[0].as_bytes()).unwrap_or_else(||
				error!("Chromosome {} is listed in target region BED file, but is not found in BAM file.", cols[0]));
			target_regions[tid as usize].push(Region {
				start: cols[1].parse::<i64>().unwrap() + 1,
				end: cols[2].parse::<i64>().unwrap()
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

	let mut total_fragments: u64 = 0;
	let mut on_target_fragments: u64 = 0;

	for read in bam {
		if read.is_secondary() || read.is_supplementary() { continue; }
		total_reads += 1;
		if read.is_unmapped() { continue; }

		aligned_reads += 1;
		if read.is_duplicate() { duplicate_reads += 1; }

		// Check whether the fragment is on-target
		if target_regions.is_empty() { continue; }

		let mut start = 0; let mut end = 0;
		if read.is_paired() {
			if read.is_mate_unmapped() { continue; }
			if read.tid() != read.mtid() { continue; }

			// Make sure that this read represents the left-most one. This also
			// ensures that we only count each DNA fragment once.
			if read.pos() > read.mpos() || (read.pos() == read.mpos() && read.is_first_in_template() == false) { continue; }

			let tlen = read.insert_size().abs();
			if tlen > max_frag_len { continue; }

			start = read.pos() + 1;
			end = start + tlen;
		} else {
			// TODO: Handle spliced reads
			start = read.pos() + 1;
			end = read.cigar().end_pos() + 1;
		}

		total_fragments += 1;

		// TODO: Use an interval tree (faster)
		for region in &target_regions[read.tid() as usize] {
			if start <= region.end && end >= region.start {
				on_target_fragments += 1;
				break;
			}

			// Can stop searching since later regions cannot overlap
			// with the alignment.
			if region.start > end { break; }
		}
    }

    println!("Total reads: {}", total_reads);
    println!("Aligned reads: {} ({:.1}% of all reads)", aligned_reads, aligned_reads as f64 / total_reads as f64 * 100.0);
    println!("Duplicate reads: {} ({:.1}% of aligned reads)", duplicate_reads, duplicate_reads as f64 / aligned_reads as f64 * 100.0);

    if target_regions.is_empty() == false {
    	println!("On-target: {:.1}%", on_target_fragments as f64 / total_fragments as f64 * 100.0);
    }
}

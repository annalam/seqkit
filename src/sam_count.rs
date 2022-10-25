
use std::collections::VecDeque;
use std::str;
use crate::common::{parse_args, BamReader, read_regions};

const USAGE: &str = "
Usage:
  sam count [options] <bam_file> <regions.bed>

Options:
  --min-mapq=N      Only count reads with MAPQ ≥ threshold [default: 0]
  --max-frag-len=N  Maximum allowed DNA fragment length [default: 5000]
  --single-end      Count individual reads, rather than DNA fragments
  --center          Only count fragments whose center is within a region

Counts the number of DNA fragments (or single reads) in the input BAM file
that overlap each region described in the input BED file. The BAM file must
be position-sorted.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>").to_string();
	let min_mapq: u8 = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be an integer between 0 - 255."));
	let max_frag_len: u32 = args.get_str("--max-frag-len").parse().unwrap_or_else(|_| error!("--max-frag-len must be an integer."));
	let single_end = args.get_bool("--single-end");
	let count_centers = args.get_bool("--center");

	// Read target regions into memory
	eprintln!("Reading target regions from BED file...");
	let mut regions = read_regions(args.get_str("<regions.bed>"));
	let mut region_frags = vec![0u32; regions.len()];

	eprintln!("Counting {}...",
		if single_end { "reads" } else { "DNA fragments" });
	let bam = BamReader::open(&bam_path);
	let chr_names: Vec<String> = bam.header().target_names().iter()
		.map(|name| str::from_utf8(name).unwrap().into()).collect();

	let mut prev_chr: i32 = -1;
	let mut prev_pos: i64 = 0;

	let mut chr_regions: VecDeque<usize> = VecDeque::new();

	for read in bam {
		if read.is_unmapped() { continue; }
		if read.is_duplicate() || read.is_secondary() { continue; }
		if read.is_supplementary() { continue; }
		if read.mapq() < min_mapq { continue; }

		if read.tid() != prev_chr {
			// We entered a new chromosome in the position-sorted BAM file
			prev_chr = read.tid();
			let chr = &chr_names[read.tid() as usize];

			// Generate a vector containing the indices of all target regions
			// within this new chromosome. Sort the indices according to the
			// left edge positions of the target regions, to enable efficient
			// interval searches.
			let mut regions_vec: Vec<usize> = Vec::new();
			for r in 0..regions.len() {
				if regions[r].chr == *chr { regions_vec.push(r); }
			}
			regions_vec.sort_by_key(|r| regions[*r].start);
			chr_regions = regions_vec.into();   // Convert into deque

		} else if read.pos() < prev_pos {
			error!("Input BAM file is not coordinate sorted.");
		}
		prev_pos = read.pos();

		let mut start = read.pos() as u32;      // 0-based inclusive start
		let mut end = if single_end {
			read.cigar().end_pos() as u32       // 0-based exclusive end
		} else {
			if read.is_paired() == false { continue; }
			if read.is_mate_unmapped() { continue; }
			if read.tid() != read.mtid() { continue; }

			// In paired end sequencing, each DNA fragment produces two reads.
			// Fragment boundaries can be inferred based on the BAM record of
			// either read, based on the template length (TLEN) field.
			// To prevent double counting and to ensure that we process
			// fragments in positional order, we ignore the BAM record that
			// represents the right-most read from the fragment. If both reads
			// have the same left-edge position, we ignore the second mate.
			if read.pos() > read.mpos() || (read.pos() == read.mpos() && read.is_first_in_template() == false) { continue; }

			let insert_size = read.insert_size().abs() as u32;
			if insert_size < 20 { continue; }

			start + insert_size                // 0-based exclusive end
		};

		if end - start > max_frag_len { continue; }

		// If the user only wants to count fragments whose center position lies
		// within a target region, we collapse each fragment to a center point.
		if count_centers {
			let len = end - start;
			start += len / 2;   
			end = start + 1;
		}

		// Now that we know the boundaries of the DNA fragment, we need to
		// figure out which target regions it overlaps with. If it overlaps
		// with multiple target regions, we count it towards each of them.

		// Remove regions from the beginning of the deque if their right edge
		// is already behind us.
		while chr_regions.is_empty() == false &&
			regions[chr_regions[0]].end < prev_pos as u32 {
			chr_regions.pop_front();
		} 

		// Calculate which region(s) this fragment overlaps with. Note that
		// the fragments and target regions are both described using semi-open
		// 0-based coordinates (inclusive start, exclusive end).
		for r in &chr_regions {
			if regions[*r].start >= end { break; }
			if regions[*r].end <= start { continue; }
			region_frags[*r] += 1;
		}
    }

    for r in 0..regions.len() {
    	println!("{}", region_frags[r]);
    }
}


use crate::common::{parse_args, BamReader};
use std::str;

const USAGE: &str = "
Usage:
  sam fragment lengths [options] <bam_file>

Options:
  --max-frag-size=F     Maximum fragment size [default: 5000]
  --reads=N             Finish after analyzing this many reads [default: Inf]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let max_frag_size: usize =
		args.get_str("--max-frag-size").parse().unwrap();
	let sufficient_sample_size: usize = if args.get_str("--reads") == "Inf" {
		std::usize::MAX
	} else {
		args.get_str("--reads").parse().unwrap()
	};

	let mut histogram: Vec<usize> = vec![0; max_frag_size + 1];

	let mut total_reads: usize = 0;
	let mut bam = BamReader::open(&bam_path);
	for read in bam {
		if read.is_paired() == false { continue; }
		if read.is_first_in_template() == false { continue; }
		if read.is_unmapped() || read.is_mate_unmapped() { continue; }
		if read.is_duplicate() || read.is_secondary() { continue; }
		if read.is_supplementary() { continue; }
		if read.tid() != read.mtid() { continue; }

		let frag_size = read.insert_size().abs() as usize;
		if frag_size > max_frag_size { continue; }

		total_reads += 1;
		histogram[frag_size] += 1;
		if total_reads >= sufficient_sample_size { break; }
    }

    for size in 1..max_frag_size+1 {
    	println!("{}\t{}", size, histogram[size]);
    }
}

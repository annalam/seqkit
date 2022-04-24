
use crate::common::{parse_args, BamReader};
use std::str;

const USAGE: &str = "
Usage:
  sam fragments [options] <bam_file>

Options:
  --min-size=N     Minimum fragment size [default: 0]
  --max-size=N     Maximum fragment size [default: 5000]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let min_frag_size: i64 = args.get_str("--min-size").parse().unwrap();
	let max_frag_size: i64 = args.get_str("--max-size").parse().unwrap();

	let mut bam = BamReader::open(&bam_path);

	let mut chr_names: Vec<String> = Vec::new();
	for name in bam.header().target_names() {
		chr_names.push(str::from_utf8(name).unwrap().to_string());
	}

	for read in bam {
		if read.is_paired() == false { continue; }
		if read.is_unmapped() || read.is_mate_unmapped() { continue; }
		if read.is_duplicate() || read.is_secondary() { continue; }
		if read.is_supplementary() { continue; }
		if read.tid() != read.mtid() { continue; }
		if read.is_reverse() || !read.is_mate_reverse() { continue; }
		if read.is_quality_check_failed() { continue; }
		// TODO: add support for non-converging read pair orientation

		let frag_size = read.insert_size().abs();
		if frag_size > max_frag_size || frag_size < min_frag_size { continue; }

		// Output in BED format (0-based half-open segments)
		println!("{}\t{}\t{}", chr_names[read.tid() as usize], read.pos(), read.pos() + frag_size);
    }
}


use common::parse_args;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;

const USAGE: &str = "
Usage:
  sam filter by sequence [options] <bam_file> <sequence>...
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let sequences = args.get_vec("<sequence>");

	let mut bam = bam::Reader::from_path(&bam_path).unwrap();

	for r in bam.records() {
		let read = r.unwrap();
		if read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		/*
		if read.tid() != read.mtid() { continue; }

		if read.pos() > read.mpos() || (read.pos() == read.mpos() && !read.is_first_in_template()) {
			continue;
		}

		let frag_size = read.insert_size().abs();
		if frag_size > max_frag_size || frag_size < min_frag_size { continue; }

		// Output in BED format (0-based half-open segments)
		println!("{}\t{}\t{}", chr_names[read.tid() as usize], read.pos(), read.pos() + frag_size);*/
    }
}

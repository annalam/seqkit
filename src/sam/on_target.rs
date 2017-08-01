
use parse_args;
use std::str;
use rust_htslib::bam;
use rust_htslib::bam::Read;

const USAGE: &'static str = "
Usage:
  sam on target [options] <bam_file> <targets.bed>

Options:
  --max-frag-size=F     Maximum fragment size [default: 5000]
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let max_frag_size: usize =
		args.get_str("--max-frag-size").parse().unwrap();

	let mut curr_chr_tid: i32 = -1;
	let bam = bam::Reader::from_path(&bam_path).unwrap();
	for r in bam.records() {
		let read = r.unwrap();
		if read.is_unmapped() { continue; }
		if read.is_duplicate() || read.is_secondary() { continue; }
		if read.is_supplementary() { continue; }
		if read.tid() != read.mtid() { continue; }

		let frag_size = read.insert_size().abs() as usize;
		if frag_size > max_frag_size { continue; }
		histogram[frag_size] += 1;
    }

    for size in 1..max_frag_size+1 {
    	println!("{}\t{}", size, histogram[size]);
    }
}


extern crate docopt;
extern crate bio;
extern crate rust_htslib;
extern crate ascii;
extern crate rand;

use std::env;

#[macro_use] mod common;
mod sam_count; mod sam_fragments; mod sam_fragment_lengths; mod sam_statistics;
mod sam_trim_qnames;
mod sam_mark_duplicates;
mod sam_to_fastq; mod sam_subsample;

const USAGE: &str = "
Usage:
  sam count <bam_file> <regions.bed>
  sam fragments <bam_file>
  sam fragment lengths <bam_file>
  sam coverage histogram <bam_file>
  sam statistics <bam_file>
  sam mark duplicates <bam_file>
  sam trim qnames <bam_file>
  sam subsample <bam_file> <fraction>

Extract reads from BAM files:
  sam to raw <bam_file> <out_prefix>
  sam to fasta <bam_file> <out_prefix>
  sam to fastq <bam_file> <out_prefix>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "count" {
		sam_count::main();
	} else if args.len() >= 2 && args[1] == "fragments" {
		sam_fragments::main();
	} else if args.len() >= 2 && args[1] == "statistics" {
		sam_statistics::main();
	} else if args.len() >= 3 && args[1..3] == ["fragment", "lengths"] {
		sam_fragment_lengths::main();
	} else if args.len() >= 3 && args[1..3] == ["mark", "duplicates"] {
		sam_mark_duplicates::main();
	//} else if args.len() >= 3 && args[1..3] == ["coverage", "histogram"] {
	//	coverage_histogram::main();
	} else if args.len() >= 3 && args[1..3] == ["trim", "qnames"] {
		sam_trim_qnames::main();
	/*} else if args.len() >= 4 && args[1..4] == ["discard", "tail", "artifacts"] {
		sam_discard_tail_artifacts::main();*/
	} else if args.len() >= 3 && args[1] == "to" &&
		(args[2] == "raw" || args[2] == "fasta" || args[2] == "fastq") {
		sam_to_fastq::main();
	} else if args.len() >= 2 && args[1] == "subsample" {
		sam_subsample::main();
	} else {
		eprintln!("{}", USAGE);
	}
}

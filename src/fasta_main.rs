
extern crate docopt;
extern crate flate2;
extern crate ascii;
extern crate bio;
extern crate num_traits;

use std::env;

#[macro_use] mod common;
mod fasta_to_raw; mod fasta_trim_by_quality; mod fasta_mask_by_quality;
mod fasta_mappability_track;
mod fasta_demultiplex;

const USAGE: &str = "
Usage:
  fasta to raw <fasta/fastq>
  fasta trim by quality <fastq_file> <min_baseq>
  fasta mask by quality <fastq_file> <min_baseq>
  fasta mappability track <genome>
  fasta demultiplex <sample_sheet> <fastq_1> <fastq_2>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 3 && args[1..3] == ["to", "raw"] {
		fasta_to_raw::main();
	} else if args.len() >= 4 && args[1..4] == ["trim", "by", "quality"] {
		fasta_trim_by_quality::main();
	} else if args.len() >= 4 && args[1..4] == ["mask", "by", "quality"] {
		fasta_mask_by_quality::main();
    } else if args.len() >= 3 && args[1..3] == ["mappability", "track"] {
		fasta_mappability_track::main();
	} else if args.len() >= 2 && args[1] == "demultiplex" {
		fasta_demultiplex::main();
	} else {
		eprintln!("{}", USAGE);
	}
}

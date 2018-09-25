
extern crate docopt;
extern crate bio;
extern crate regex;

use std::env;

#[macro_use] mod common;
mod fasta_to_raw; mod fasta_simplify_read_ids;
mod fasta_trim_by_quality; mod fasta_mask_by_quality;
mod fasta_mappability_track; mod fasta_add_barcode;
mod fasta_demultiplex; mod fasta_convert_basespace;
mod fasta_statistics; mod fasta_remove_base_qualities;

const USAGE: &str = "
Usage:
  fasta to raw <fasta/fastq>
  fasta remove base qualities <fastq>
  fasta simplify read ids <fastq_file>
  fasta trim by quality <fastq_file> <min_baseq>
  fasta mask by quality <fastq_file> <min_baseq>
  fasta mappability track <genome>
  fasta add barcode <fastq_file> <barcode_file> <barcode_format>
  fasta convert basespace <fastq_file>
  fasta demultiplex <sample_sheet> <fastq_1> <fastq_2>
  fasta statistics <fastq_file>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 3 && args[1..3] == ["to", "raw"] {
		fasta_to_raw::main();
	} else if args.len() >= 4 && args[1..4] == ["remove", "base", "qualities"] {
		fasta_remove_base_qualities::main();
	} else if args.len() >= 4 && args[1..4] == ["simplify", "read", "ids"] {
		fasta_simplify_read_ids::main();
	} else if args.len() >= 4 && args[1..4] == ["trim", "by", "quality"] {
		fasta_trim_by_quality::main();
	} else if args.len() >= 4 && args[1..4] == ["mask", "by", "quality"] {
		fasta_mask_by_quality::main();
    } else if args.len() >= 3 && args[1..3] == ["mappability", "track"] {
		fasta_mappability_track::main();
	} else if args.len() >= 3 && args[1..3] == ["add", "barcode"] {
		fasta_add_barcode::main();
	} else if args.len() >= 3 && args[1..3] == ["convert", "basespace"] {
		fasta_convert_basespace::main();
	} else if args.len() >= 2 && args[1] == "demultiplex" {
		fasta_demultiplex::main();
	} else if args.len() >= 2 && args[1] == "statistics" {
		fasta_statistics::main();
	} else {
		eprintln!("{}", USAGE);
	}
}

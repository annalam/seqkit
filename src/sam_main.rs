
#![allow(unused_must_use)]

use std::env;

#[macro_use] mod common;
mod sam_count; mod sam_fragments; mod sam_fragment_lengths;
mod sam_statistics;
mod sam_to_fastq; mod sam_subsample;
mod sam_coverage_histogram; mod sam_merge;
mod sam_minimize; mod sam_tags_from_qname; mod sam_qname_from_tags;
mod sam_trim_qnames;
mod sam_mark_duplicates;
mod sam_consensus;

const USAGE: &str = "
Usage:
  sam merge <bam_files>...
  sam consensus <bam_file>
  sam count <bam_file> <regions.bed>
  sam coverage histogram <bam_file>
  sam fragments <bam_file>
  sam fragment lengths <bam_file>
  sam mark duplicates <bam_file>
  sam minimize <bam_file>
  sam statistics <bam_file>
  sam subsample <bam_file> <fraction>  
  sam tags from qname <bam_file>
  sam qname from tags <bam_file>
  sam trim qnames <bam_file>  

Extract reads from BAM files:  
  sam to fasta <bam_file> <out_prefix>
  sam to fastq <bam_file> <out_prefix>  
  sam to interleaved fasta <bam_file>
  sam to interleaved fastq <bam_file>
  sam to interleaved raw <bam_file>
  sam to raw <bam_file> <out_prefix>
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
	} else if args.len() >= 3 && args[1..3] == ["coverage", "histogram"] {
		sam_coverage_histogram::main();
	} else if args.len() >= 3 && args[1] == "to" && 
		(args[2] == "raw" || args[2] == "fasta" || args[2] == "fastq") {
		sam_to_fastq::main();
	} else if args.len() >= 4 && args[1..3] == ["to", "interleaved"] &&
		(args[3] == "raw" || args[3] == "fasta" || args[3] == "fastq") {
		sam_to_fastq::main();
	} else if args.len() >= 2 && args[1] == "subsample" {
		sam_subsample::main();
	} else if args.len() >= 2 && args[1] == "merge" {
		sam_merge::main();
	} else if args.len() >= 2 && args[1] == "minimize" {
		sam_minimize::main();
	} else if args.len() >= 4 && args[1..4] == ["tags", "from", "qname"] {
		sam_tags_from_qname::main();
	} else if args.len() >= 4 && args[1..4] == ["qname", "from", "tags"] {
		sam_qname_from_tags::main();
	} else if args.len() >= 3 && args[1..3] == ["trim", "qnames"] {
		sam_trim_qnames::main();
	} else if args.len() >= 3 && args[1..3] == ["mark", "duplicates"] {
		sam_mark_duplicates::main();
	} else if args.len() >= 2 && args[1] == "consensus" {
		sam_consensus::main();
	} else {
		eprintln!("{}", USAGE);
	}
}

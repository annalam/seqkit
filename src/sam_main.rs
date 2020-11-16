
use std::env;

#[macro_use] mod common;
mod sam_count; mod sam_fragments; mod sam_fragment_lengths;
mod sam_statistics;
mod sam_to_fastq; mod sam_subsample;
mod sam_coverage_histogram; mod sam_concatenate;
mod sam_minimize; mod sam_tags_from_qname;
mod sam_trim_qnames;
mod sam_determine_sex;

const USAGE: &str = "
Usage:
  sam count <bam_file> <regions.bed>
  sam fragments <bam_file>
  sam fragment lengths <bam_file>
  sam coverage histogram <bam_file>
  sam statistics <bam_file>
  sam subsample <bam_file> <fraction>
  sam concatenate <bam_files>...
  sam minimize <bam_file>
  sam tags from qname <bam_file>
  sam trim qnames <bam_file>
  sam determine sex <bam_file>

Extract reads from BAM files:
  sam to raw <bam_file> <out_prefix>
  sam to fasta <bam_file> <out_prefix>
  sam to fastq <bam_file> <out_prefix>
  sam to interleaved raw <bam_file>
  sam to interleaved fasta <bam_file>
  sam to interleaved fastq <bam_file>
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
	} else if args.len() >= 2 && args[1] == "concatenate" {
		sam_concatenate::main();
	} else if args.len() >= 2 && args[1] == "minimize" {
		sam_minimize::main();
	} else if args.len() >= 4 && args[1..4] == ["tags", "from", "qname"] {
		sam_tags_from_qname::main();
	} else if args.len() >= 3 && args[1..3] == ["trim", "qnames"] {
		sam_trim_qnames::main();
	} else if args.len() >= 3 && args[1..3] == ["determine", "sex"] {
		sam_determine_sex::main();
	} else {
		eprintln!("{}", USAGE);
	}
}

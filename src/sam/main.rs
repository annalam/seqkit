
extern crate docopt;
extern crate bio;
extern crate rust_htslib;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};

mod count; mod fragment_lengths; mod statistics;

const USAGE: &'static str = "
Usage:
  sam count <bam_file> <regions.bed>
  sam fragment lengths <bam_file>
  sam statistics <bam_file>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "count" {
		count::main();
	} else if args.len() >= 2 && args[1] == "statistics" {
		statistics::main();
	} else if args.len() >= 3 && args[1..3] == ["fragment", "lengths"] {
		fragment_lengths::main();
	}
	//else if args.len() >= 2 && args[1] == "somatic" { somatic::main(); }
	else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		eprintln!("Invalid arguments.\n{}", usage); exit(-1);
	})
}
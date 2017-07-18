
extern crate docopt;
//extern crate bio;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};

mod to_raw;

const USAGE: &'static str = "
Usage:
  fasta <subcommand>

Available subcommands:
  to raw             Convert to raw (headerless) format
  trim by quality    Trim by quality
  mask by quality    Mask by quality
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 3 && args[1..3] == ["to", "raw"] {
		to_raw::main();
	//} else if args.len() >= 4 && args[1..4] == ["trim", "by", "quality"] {
	//	trim_by_quality::main();
	//} else if args.len() >= 4 && args[1..4] == ["mask", "by", "quality"] {
	//	mask_by_quality::main();
	} else {
		println!("{}", USAGE); exit(-1);
	}
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		println!("Invalid arguments.\n{}", usage); exit(-1);
	})
}
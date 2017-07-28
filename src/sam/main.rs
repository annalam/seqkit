
extern crate docopt;
extern crate bio;
extern crate rust_htslib;

use std::env;
use std::process::exit;
use std::io::{Write, stderr};
use docopt::{Docopt, ArgvMap};

mod sam;

const USAGE: &'static str = "
Usage:
  sam <subcommand>

Available subcommands:
  count      Count reads that overlap BED file regions
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "count" { sam::main(); }
	//else if args.len() >= 2 && args[1] == "somatic" { somatic::main(); }
	else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		writeln!(stderr(), "Invalid arguments.\n{}", usage); exit(-1);
	})
}
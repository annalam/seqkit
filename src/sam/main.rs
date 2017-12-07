
extern crate docopt;
extern crate bio;
extern crate rust_htslib;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};

mod count; mod fragments; mod fragment_lengths; mod statistics;
mod coverage_histogram; mod mark_duplicates;

const USAGE: &'static str = "
Usage:
  sam count <bam_file> <regions.bed>
  sam fragments <bam_file>
  sam fragment lengths <bam_file>
  sam coverage histogram <bam_file>
  sam statistics <bam_file>
  sam mark duplicates <bam_file>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "count" {
		count::main();
	} else if args.len() >= 2 && args[1] == "fragments" {
		fragments::main();
	} else if args.len() >= 2 && args[1] == "statistics" {
		statistics::main();
	} else if args.len() >= 3 && args[1..3] == ["fragment", "lengths"] {
		fragment_lengths::main();
	} else if args.len() >= 3 && args[1..3] == ["mark", "duplicates"] {
		mark_duplicates::main();
	} else if args.len() >= 3 && args[1..3] == ["coverage", "histogram"] {
		coverage_histogram::main();
	} else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		eprintln!("Invalid arguments.\n{}", usage); exit(-1);
	})
}


// Helper methods for error reporting
trait ErrorHelper<T> {
	fn on_error(self, msg: &str) -> T;
}

impl<T> ErrorHelper<T> for Option<T> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Some(x) => x,
			None => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}

impl<T, E> ErrorHelper<T> for Result<T, E> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Ok(x) => x,
			Err(_) => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}

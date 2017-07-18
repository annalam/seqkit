
use parse_args;
use std::str;
use std::io::{BufReader, BufWriter, BufRead, Write, stderr};

const USAGE: &'static str = "
Usage:
  fasta to raw <fasta_file>
";

pub fn main() {
	let args = parse_args(USAGE);
}


extern crate docopt;
extern crate flate2;
extern crate ascii;
extern crate bio;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};
use std::fs::File;
use std::io::{stdin, BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use ascii::AsciiString;
use std::mem;

mod to_raw; mod trim_by_quality; mod mask_by_quality;
mod mapq_track;

const USAGE: &'static str = "
Usage:
  fasta to raw <fasta/fastq>
  fasta trim by quality <fastq_file> <min_baseq>
  fasta mask by quality <fastq_file> <min_baseq>
  fasta mapq track <genome>
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 3 && args[1..3] == ["to", "raw"] {
		to_raw::main();
	} else if args.len() >= 4 && args[1..4] == ["trim", "by", "quality"] {
		trim_by_quality::main();
	} else if args.len() >= 4 && args[1..4] == ["mask", "by", "quality"] {
		mask_by_quality::main();
    } else if args.len() >= 3 && args[1..3] == ["mapq", "track"] {
		mapq_track::main();
	} else {
		eprintln!("{}", USAGE); exit(-1);
	}
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		eprintln!("Invalid arguments.\n{}", usage); exit(-1);
	})
}

pub fn read_buffered(path: &str) -> Box<BufRead> {
	if path == "-" {
		Box::new(BufReader::new(stdin()))
	} else {
		let file = File::open(path).on_error(&format!(
			"Cannot open file {}.", path));
		if path.ends_with(".gz") {
			Box::new(BufReader::new(MultiGzDecoder::new(file).unwrap()))
		} else {
			Box::new(BufReader::new(file))
		}
	}
}

// Helper method for reading ASCII format files
trait AsciiBufRead {
	fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool;
}

impl<T: BufRead> AsciiBufRead for T {
	fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool {
		line.clear();
		let v = unsafe { mem::transmute::<&mut AsciiString, &mut Vec<u8>>(line) };
		self.read_until('\n' as u8, v).unwrap() > 0
	}
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

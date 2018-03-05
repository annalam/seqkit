
use docopt::{Docopt, ArgvMap};
use std::io::{stdin, BufRead, BufReader};
use std::fs::File;
use flate2::read::MultiGzDecoder;
use std::mem;
use ascii::AsciiString;

macro_rules! error {
	($($arg:tt)+) => ({
		use std::process::exit;
		eprint!("ERROR: "); eprintln!($($arg)+); exit(-1);
	})
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		error!("Invalid arguments.\n{}", usage);
	})
}

pub fn read_buffered(path: &str) -> Box<BufRead> {
	if path == "-" {
		Box::new(BufReader::new(stdin()))
	} else {
		let file = File::open(path).unwrap_or_else(
			|_| error!("Cannot open file {}.", path));
		if path.ends_with(".gz") {
			Box::new(BufReader::new(MultiGzDecoder::new(file)))
		} else {
			Box::new(BufReader::new(file))
		}
	}
}

// Helper method for reading ASCII format files
pub trait AsciiBufRead {
	fn next_line(&mut self, line: &mut String) -> bool;
	fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool;
}

impl<T: BufRead> AsciiBufRead for T {
	fn next_line(&mut self, line: &mut String) -> bool {
		line.clear();
		match self.read_line(line) {
			Ok(len) => len > 0,
			_ => { error!("I/O error while reading from file."); }
		}
	}

	fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool {
		line.clear();
		let v = unsafe { mem::transmute::<&mut AsciiString, &mut Vec<u8>>(line) };
		self.read_until(b'\n', v).unwrap() > 0
	}
}

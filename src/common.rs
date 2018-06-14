
use docopt::{Docopt, ArgvMap};
use std::process::{Command, Stdio, ChildStdin};
use std::fmt::Arguments;
use std;
use std::io::{stdin, BufRead, BufReader, Write};
use std::fs::File;
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

pub trait PathArgs {
	fn get_path(&self, arg: &str) -> String;
}

impl PathArgs for ArgvMap {
	fn get_path(&self, arg: &str) -> String {
		let path = self.get_str(arg);
		// TODO: Check path validity?
		if path.starts_with('~') {
			if let Some(home) = std::env::home_dir() {
				return format!("{}{}", home.display(), &path[1..]);
			}
		}
		path.into()
	}
}

pub struct GzipWriter {
	//gzip: Child
	gzip: ChildStdin
}

impl GzipWriter {
	pub fn new(path: &str) -> GzipWriter {
		let file = File::create(path).unwrap_or_else(
			|_| error!("Cannot open file {} for writing.", path));
		GzipWriter {
			gzip: Command::new("gzip").arg("-c")
				.stdin(Stdio::piped()).stdout(file).spawn()
				.unwrap_or_else(|_| error!("Cannot start gzip process."))
				.stdin.unwrap()
		}
	}

	pub fn write_fmt(&mut self, fmt: Arguments) -> std::io::Result<()> {
		self.gzip.write_fmt(fmt)
	}
}

pub struct FileReader {
	bufread: Box<BufRead>
}

impl FileReader {
	pub fn new(path: &str) -> FileReader {
		let bufread: Box<BufRead> = if path == "-" {
			Box::new(BufReader::new(stdin()))
		} else {
			let file = File::open(path).unwrap_or_else(
				|_| error!("Cannot open file {} for reading.", path));
			if path.ends_with(".gz") {
				Box::new(BufReader::new(Command::new("gunzip").arg("-c")
					.stdout(Stdio::piped()).stdin(file).spawn()
					.unwrap_or_else(|_| error!("Cannot start gunzip process."))
					.stdout.unwrap()))
			} else {
				Box::new(BufReader::new(file))
			}
		};
		FileReader { bufread }
	}

	pub fn read_line(&mut self, line: &mut String) -> bool {
		line.clear();
		match self.bufread.read_line(line) {
			Ok(len) => len > 0,
			_ => { error!("I/O error while reading from file."); }
		}
	}

	pub fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool {
		line.clear();
		let v = unsafe { std::mem::transmute::<&mut AsciiString, &mut Vec<u8>>(line) };
		self.bufread.read_until(b'\n', v).unwrap() > 0
	}
}

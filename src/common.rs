
use docopt::{Docopt, ArgvMap};
use std::process::{Command, Stdio, ChildStdin};
//use std::fmt::Arguments;
use std;
use std::fs::File;
use std::io::{stdin, BufRead, BufReader, Write};
use std::os::unix::io::{FromRawFd, AsRawFd};
use rust_htslib::bam::{self, Read, Header, HeaderView, Format};
use rust_htslib::bam::CompressionLevel;

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

pub fn dirname(path: &str) -> String {
	path.rsplitn(2, '/').last().unwrap().into()
}

pub fn filename(path: &str) -> String {
	path.split('/').last().unwrap().into()
}

pub struct GzipWriter {
	//gzip: Child
	gzip: ChildStdin      // TODO: Add BufWriter for improved performance
}

#[derive(Copy, Clone)]
pub enum Compressor { GZIP, PIGZ }

impl GzipWriter {
	pub fn new(path: &str) -> GzipWriter {
		GzipWriter::with_method(path, Compressor::GZIP)
	}

	pub fn with_method(path: &str, method: Compressor) -> GzipWriter {
		let file = File::create(path).unwrap_or_else(
			|_| error!("Cannot open file {} for writing.", path));
		let compressor = match method {
			Compressor::GZIP => Command::new("gzip").arg("-c")
				.stdin(Stdio::piped()).stdout(file).spawn()
				.unwrap_or_else(|_| error!("Cannot start gzip process.")),
			Compressor::PIGZ => Command::new("pigz").arg("-c")
				.stdin(Stdio::piped()).stdout(file).spawn()
				.unwrap_or_else(|_| error!("Cannot start pigz process."))
		};
		GzipWriter { gzip: compressor.stdin.unwrap() }
	}
}

impl Write for GzipWriter {
	fn write(&mut self, buf: &[u8]) -> std::io::Result<usize>
		{ self.gzip.write(buf) }
	fn flush(&mut self) -> std::io::Result<()> { self.gzip.flush() }
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

	/*pub fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool {
		line.clear();
		let v = unsafe { std::mem::transmute::<&mut AsciiString, &mut Vec<u8>>(line) };
		self.bufread.read_until(b'\n', v).unwrap() > 0
	}*/
}

pub struct BamReader {
	reader: bam::Reader,
	record: bam::Record
}

impl BamReader {
	pub fn open(path: &str) -> BamReader {
		let reader = if path == "-" {
			bam::Reader::from_stdin().unwrap_or_else(
				|_| error!("Failed to read BAM file from standard input."))
		} else {
			bam::Reader::from_path(&path).unwrap_or_else(
				|_| error!("Cannot open BAM file '{}'", path))
		};
		BamReader { reader, record: bam::Record::new() }
	}

	pub fn header(&self) -> bam::HeaderView {
		self.reader.header().clone()
	}
}

impl Iterator for BamReader {
	type Item = bam::Record;

	fn next(&mut self) -> Option<bam::Record> {
		match self.reader.read(&mut self.record) {
			Ok(true) => Some(self.record.clone()),
			Ok(false) => None,
			Err(bam::Error::TruncatedRecord) =>
				error!("BAM file ended prematurely."),
			Err(bam::Error::InvalidRecord) =>
				error!("Invalid BAM record."),
			Err(e) => error!("Failed reading BAM record: {}", e)
		}
	}
}

pub struct BamWriter {
	writer: bam::Writer
}

impl BamWriter {
	pub fn open(path: &str, header: &HeaderView, compressed: bool) -> BamWriter {

		let mut writer = if path == "-" {
			bam::Writer::from_stdout(&Header::from_template(&header),
				Format::BAM).unwrap_or_else(|_|
				error!("Could not write BAM records into standard output."))
		} else {
			bam::Writer::from_path(&path, &Header::from_template(&header),
				Format::BAM).unwrap_or_else(
				|_| error!("Cannot open BAM file '{}' for writing.", path))
		};
		if compressed == false {
			writer.set_compression_level(CompressionLevel::Uncompressed);
		}

		BamWriter { writer }
	}

	pub fn write(&mut self, record: &bam::Record) {
		self.writer.write(&record).unwrap_or_else(|_|
			error!("Failed at writing BAM record."));
	}
}



// CODE FOR READING TARGET REGIONS INTO MEMORY
pub struct Region {
	pub chr: Box<str>,
	pub start: usize,   // 0-based position (inclusive)
	pub end: usize,     // 0-based position (exclusive)
	pub name: Box<str>
}

// Read genomic regions into memory from a BED file
pub fn read_regions(bed_path: &str) -> Vec<Region> {
	let mut regions: Vec<Region> = Vec::new();
	let mut regions_file = FileReader::new(bed_path);
	let mut line = String::new();
	while regions_file.read_line(&mut line) {
		if line.starts_with('#') { continue; }
		let cols: Vec<&str> = line.trim().split('\t').collect();
		if cols.len() < 3 {
			error!("Invalid region in BED file:\n{}", &line);
		}

		// Allow the fourth column to be missing.
		let name = if cols.len() >= 4 { cols[3].into() } else { "".into() };
		
		regions.push(Region {
			chr: cols[0].into(), name: name,
			start: cols[1].parse().unwrap(), end: cols[2].parse().unwrap()
		});
	}
	// Sort the regions by chromosome and then start position.
	//regions.sort_unstable_by(|a, b| {
	//	let ord = a.chr.cmp(&b.chr);
	//	if ord == Ordering::Equal { a.start.cmp(&b.start) } else { ord }
	//});
	regions
}


use crate::common::{parse_args, FileReader};
use std::str;
use std::collections::VecDeque;

const USAGE: &str = "
Usage:
  fasta check <fasta/fastq>

Description:
Checks that the input FASTA or FASTQ file is correctly formatted, and reports
the line number if any malformatted lines are found.
";

struct ReaderWithMemory {
	file: FileReader,
	lines_read: usize,
	prev_lines: VecDeque<String>
}

impl ReaderWithMemory {
	fn new(path: &str) -> ReaderWithMemory {
		ReaderWithMemory {
			file: FileReader::new(path),
			prev_lines: VecDeque::new(),
			lines_read: 0
		}
	}

	fn read_line(&mut self, line: &mut String) -> bool {
		if self.file.read_line(line) == false { return false; }
		self.prev_lines.push_back(line.clone());
		if self.prev_lines.len() > 10 {
			self.prev_lines.pop_front();
		}
		self.lines_read += 1;
		return true;
	}

	fn history(&self) -> String {
		let mut history = String::new();
		for line in &self.prev_lines {
			history += &line; history += "\n";
		}
		history
	}
}

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta = ReaderWithMemory::new(&args.get_str("<fasta/fastq>"));

	let mut line = String::new();
	while fasta.read_line(&mut line) {
		if line.starts_with('>') {
			fasta.read_line(&mut line);
		} else if line.starts_with('@') {
			fasta.read_line(&mut line);
			fasta.read_line(&mut line);
			if !line.starts_with('+') {
				error!("Missing quality header prefix '+' on line {}:\n{}\n",
					fasta.lines_read, fasta.history());
			}
			fasta.read_line(&mut line);
		} else {
			error!("Missing header prefix '>' or '@' on line {}:\n{}\n",
				fasta.lines_read, fasta.history());
		}
	}
}

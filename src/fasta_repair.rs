
use common::{parse_args, FileReader};
use std::str;
use std::io::Write;

const USAGE: &str = "
Usage:
  sam repair <sam_file>
";

pub fn main() {
	let args = parse_args(USAGE);
	let mut sam = FileReader::new(&args.get_str("<sam_file>"));

	let mut line = String::new();
	while sam.read_line(&mut line) {
		if let Some(pos) = line.find("BC:") {
			print!("{}UMI:{}{}", &line[0..pos], &line[pos+3..pos+7], &line[pos+11..]);
		} else {
			print!("{}", line);
		}
	}
}

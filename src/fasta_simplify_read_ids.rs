
use common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta simplify read ids <fastq_file>
";

static IDCHARS: [u8;62] = [  b'0',b'1',b'2',b'3',b'4',b'5',b'6',b'7',b'8',b'9',
						     b'A',b'B',b'C',b'D',b'E',b'F',b'G',b'H',b'I',b'J',b'K',b'L',b'M',b'N',b'O',b'P',b'Q',b'R',b'S',b'T',b'U',b'V',b'W',b'X',b'Y',b'Z', 
						     b'a',b'b',b'c',b'd',b'e',b'f',b'g',b'h',b'i',b'j',b'k',b'l',b'm',b'n',b'o',b'p',b'q',b'r',b's',b't',b'u',b'v',b'w',b'x',b'y',b'z' ];
static NCHARS: usize = 62;

pub fn next_id( id : &mut Vec<usize>, str_id : &mut Vec<u8>) {

	let mut i = id.len(); // index of running identifier char
	
    if i != 0 {
    
        i -= 1;
    
    	loop  {
    
    		let running_num = id[ i];
    
    		if running_num == NCHARS-1 {
    		    // Reset running char and examine previous
    			id[ i] = 0;
    			str_id[ i] = IDCHARS[ 0];
    			if i == 0 { break; }
    			i -= 1;
    			continue;
    		}
    		
    		// Replace running char with next
    		id[ i] = running_num+1;
    		str_id[ i] = IDCHARS[ running_num+1];
    		return;
    	}
	}

    id.insert(0, 0);
    str_id.insert(0, IDCHARS[ 0]);
}

pub fn main() {
	let args = parse_args(USAGE);
	let mut fasta_file = FileReader::new(&args.get_str("<fastq_file>"));	

	let mut line = String::new();

	let mut str_id : Vec<u8> = Vec::new();
	let mut running_id : Vec<usize> = Vec::new();

	while fasta_file.read_line(&mut line) {
		
		next_id( &mut running_id, &mut str_id);

		if line.starts_with('@') {
			println!("@{}", str::from_utf8(&str_id).unwrap());
			for _ in 0..3 {
				fasta_file.read_line(&mut line);
				print!("{}", line);
			}
		} else if line.starts_with('>') {
			println!(">{}", str::from_utf8(&str_id).unwrap());
			fasta_file.read_line(&mut line);
			print!("{}", line);
		} else {
			error!("Invalid FASTA/FASTQ format encountered.");
		}
	}
}

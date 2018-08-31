
use common::{parse_args, FileReader};
use std::str;

const USAGE: &str = "
Usage:
  fasta simplify read ids <fastq_file>
";


static IDCHARS: [char;62] = ['0','1','2','3','4','5','6','7','8','9',
						     'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z', 
						     'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z' ];
static NCHARS: usize = 62;

pub fn next_id( id : &mut Vec<usize>, str_id : &mut Vec<char>) {

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

	let mut str_id : Vec<char> = Vec::new();
	let mut running_id : Vec<usize> = Vec::new();

	while fasta_file.read_line(&mut line) {
		
		next_id( &mut running_id, &mut str_id);
        let s: String = str_id.iter().collect();		

		if line.starts_with('@') {
			println!("@{}", s);
			for _ in 0..3 {
				fasta_file.read_line(&mut line);
				print!("{}", line);
			}
		} else if line.starts_with('>') {
			println!(">{}", s);
			fasta_file.read_line(&mut line);
			print!("{}", line);
		} else {
			error!("Invalid FASTA/FASTQ format encountered.");
		}
	}
}

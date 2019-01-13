use crate::common::parse_args;
use std::str;
use std::process::exit;
use std::ascii::AsciiExt;
use std::time::{Duration, Instant};
 
use rust_htslib::bam;
use rust_htslib::bam::{Records, Read};
use rust_htslib::bam::record::{Cigar, CigarStringView, Record, Seq};

use genome_reader::RefGenomeReader;

// Ref genome:
// /home/annalam/homo_sapiens/hg38.fa

const USAGE: &str = "
Usage:
  sam discard tail artifacts [options] <ref_genome.fa> <bam_file>
  sam discard tail artifacts --test

Options:  
  --tail-length=N    The tail length of the read that is examined [default: 20].
                     The tail length specifies the number of bases that are 
                     examined from each end of the read. Indels encountered do
                     not affect this number.
  --threshold=FLOAT  Threshold ratio (0.0-1.0, inclusive) of variations that
                     will cause the read to be discarded [default: 0.15].
  --mismatches=N     Number of mismatches (inclusive) in a tail that will cause 
                     the read to be discarded.
                     This number will be converted to a threshold ratio.
  --debug            Print information and save discarded reads to a separate 
                     file.
  --test             Run tests and exit.

Description:

This script processes all the reads in a bam-file and removes reads that have 
a number of mismatching bases in the tails (ends) of the reads that is above 
the specified threshold ratio. Output (new bam file) is written to stdout.

";


pub fn main() {

	let version = "0.5";	
	let args = parse_args(USAGE);
	if args.get_bool("--test") { exit( test()); }

	let genome_path = args.get_str("<ref_genome.fa>");
	let bam_path = args.get_str("<bam_file>");
	let debug = args.get_bool("--debug");	
	let tail_length = args.get_str("--tail-length").parse::<usize>().unwrap_or_else(|_| error!("--tail-length requires a positive integer argument."));
	let mut mm_treshold_ratio = args.get_str("--threshold").parse::<f32>().unwrap_or_else(|_| error!( "--threshold requires a floating point number argument (default: 0.10)."));
	let mismatches = args.get_str("--mismatches");


	//Sanity checks
	if tail_length < 1 {
		error!("--tail-length requires a positive integer argument.");
	}

	if !mismatches.is_empty() { 				
		let mm	= mismatches.parse::<u32>().unwrap_or_else(|_| error!( "--mismatches requires a positive integer argument."));
		mm_treshold_ratio = (mm as f32) / (tail_length as f32);
		eprintln!( "INFO: Mismatch threshold ratio set to {:.2}", mm_treshold_ratio);
	}

	if mm_treshold_ratio < 0.0 || mm_treshold_ratio > 1.0 {
		error!("--threshold requires a floating point number between 0.0 and 1.0.");
	}

	eprintln!( "INFO: {} mismatches in {} read tail bases will cause read to be discarded.", check_threshold( tail_length, mm_treshold_ratio), tail_length);


	//Copy bam header and add comment
	let mut bam_reader = bam::Reader::from_path(&bam_path)
		.expect(&format!("Cannot open BAM file {}", &bam_path));
	let header_view = bam_reader.header().clone();
	
	let mut header = bam::header::Header::from_template( &header_view);
	header.push_comment( &format!("Processed with discard tail artifacts (ver {}) TAIL_LEN:{} THRESHOLD:{:.1}", version, tail_length, mm_treshold_ratio).as_bytes()); //@CO
	let mut out = bam::Writer::from_stdout( &header).unwrap_or_else(|_| error!("Failed to create output stream."));


	let debug_filename = &bam_path.replace(".bam", "_tail_discards_debug.bam");

	let mut discarded_out = if debug { 								
		eprintln!("DEBUG: Writing discarded reads to file '{}'. (debug)", &debug_filename);
		bam::Writer::from_path( &debug_filename, &bam::header::Header::from_template(&header_view)).unwrap_or_else(
			|_| error!("Cannot create file '{}'", &debug_filename))
	} else { 
		bam::Writer::from_path( "/dev/null", &bam::header::Header::from_template(&header_view)).unwrap_or_else(
			|_| error!("Failed to create dummy bam::Writer.")) 
	};


	// Load the chromosome indices from FASTA file into memory. 
	// RefGenomeReader is used later to read chromosome sequenes into memory.
	let mut ref_genome_reader = RefGenomeReader::new( genome_path);	

	let mut chr_names: Vec<String> = Vec::new();
	let mut chr_lens: Vec<usize> = Vec::new();

	//Read the names and lengths of chromosomes that reads
	//have been aligned to
	for tid in 0..header_view.target_names().len() {

		let name = String::from_utf8( header_view.target_names()[ tid].to_vec()).unwrap_or_else(|_| error!("Failed to read target names."));
		let chr_len = header_view.target_len( tid as u32)
			.unwrap_or_else(|| error!("Failed to read target length.")) as usize;

		chr_names.push( name);
		chr_lens.push( chr_len);
	}

	let mut chr_index = 0;
	//Load first chromosome sequence to memory
	let mut chr_seq = ref_genome_reader.load_chromosome_seq( &chr_names[ chr_index]);

	let mut records_total = 0i32;
	let mut records_passed = 0i32;
	let mut records_failed = 0i32;
	let mut report_progress = 0;

	eprintln!( "INFO: Running DISCARD TAIL ARTIFACTS...");
	let start_time = Instant::now(); //save processing start time

	//Loop through all records
	loop {

		let mut read = Record::new();

		if read_record(&mut bam_reader, &bam_path, &mut read) == false {
			break; //No more records
		}

		records_total += 1;
		
		//Do not analyze unmapped reads (just output), analyze secondary or supplementary reads
		if read.is_unmapped() { //read.is_secondary() || read.is_supplementary

			//Unmapped reads have tid() == -1
			out.write( &read).unwrap_or_else(|_| error!("Failed to write output stream."));
			continue;
		}

		//Need to load new chromosome?
		if read.tid() != chr_index as i32 {					
			chr_index = read.tid() as usize;

			//DEBUG printing
			if debug {
				//eprintln!( "CHR index: {}", read.tid());
				//eprintln!( "CHR name:'{}'", &chr_names[ chr_index]);
			}

			if chr_index >= chr_names.len() {
				error!("Chromosome index error!");
			}
			chr_seq = ref_genome_reader.load_chromosome_seq( &chr_names[ chr_index]);
		}

		//DEBUG
		//let mapq = read.mapq() as i32; //Mapping quality not used
		//eprintln!( "Read:  suppl: {}, secondary: {}", read.is_secondary(), read.is_supplementary());
		//eprintln!( "QNAME: {}", str::from_utf8( &read.qname()).unwrap());
		//eprintln!( "CIGAR: {}", read.cigar());
				
		let mut tail_left_mismatches = 0;
		let mut tail_right_mismatches = 0;

		//Cannot call mutable and immutable functions of Record in same scope
		//Extra brackets limit the scope of read.seq() return value
		{ 
			let seq = read.seq();
			let read_len = seq.len();	

			let mut ref_offset = read.pos() as usize; // 0-based position within chromosome
			//Iterate cigar from left to right
			tail_left_mismatches = count_mismatches( tail_length, read.cigar().iter(), &seq, 0usize, ref_offset, &chr_seq, 1 );

			//To iterate through the seq from right to left, we need to know 
			//the offset between the indices of seq and refseq due to Indels
			//at the right end of the read 
			//NOTE: it would be possible to skip reading the right tail if
			//      left tail missmatches are above the threshold
			let offset = count_right_end_offset( read.cigar());			
			ref_offset = ((ref_offset as i32) + offset) as usize;
			//Iterate cigar from right to left
			tail_right_mismatches = count_mismatches( tail_length, read.cigar().iter().rev(), &seq, read_len-1, ref_offset, &chr_seq, -1 );

		} //end scope for read.seq()

		let left_tail_mm_ratio = tail_left_mismatches as f32  / tail_length as f32;
		let right_tail_mm_ratio = tail_right_mismatches as f32 / tail_length as f32;

		if left_tail_mm_ratio >= mm_treshold_ratio {			
			//No need to set quality check failed
			//if only passed reads are outputted
			//read.set_quality_check_failed(); 
			records_failed += 1;
			if debug { discarded_out.write( &read).unwrap_or_else(|_| error!( "Failed to write debug output")); }
			//eprintln!( "DISCARDING: {}, read: {}-{}", &chr_names[ chr_index], read.pos(), read.pos() + read.seq().len() as i32);
		}
		else if right_tail_mm_ratio >= mm_treshold_ratio {
			records_failed += 1;
			//read.set_quality_check_failed();
			if debug { discarded_out.write( &read).unwrap_or_else(|_| error!("Failed to write debug output")); }
			//eprintln!( "DISCARDING: {}, read: {}-{}", &chr_names[ chr_index], read.pos(), read.pos() + read.seq().len() as i32);
		}
		else {
			records_passed += 1;
			out.write( &read).unwrap_or_else(|_| error!("Failed to write output stream."));			
		}

		if debug {
			//Print out progress
			report_progress += 1;
			if report_progress >= 100000 {
				report_progress = 0;
				let fail_pct = ( records_failed as f32  / records_total as f32) * 100.0;
				eprintln!( "Records processed: {}/{} ({:.1}%) [FAILED/TOTAL]", records_failed, records_total, fail_pct );			
			}
		}
	} //loop records

	//Print statistics
	eprintln!( "\nDISCARD TAIL ARTIFACT RESULTS");	
	eprintln!( "Total number of reads processed: {}", records_total);
	eprintln!( "Number of reads with {:.1}% or more tail artifacts: {} [FAILED]", mm_treshold_ratio * 100.0, records_failed);
	eprintln!( "Number of reads with less than {:.1}% tail artifacts: {} [PASSED]", mm_treshold_ratio * 100.0, records_passed);
	eprintln!( "");
	eprintln!( "Percentage of processed records discarded {:.1}%", (records_failed as f32  / records_total as f32) * 100.0);	
	eprintln!( "");
	eprintln!("Processing took {} seconds.", start_time.elapsed().as_secs());
	
	if debug { eprintln!("Reads with tail artifacts written to file: '{}'", debug_filename ); }
	eprintln!( "ALL DONE.");
}

//Count combined length of indels
fn count_right_end_offset( csv : CigarStringView ) -> i32 {

	let mut offset = 0i32;

	for c in csv.iter() {	
		match *c {
			Cigar::Ins(len) => { offset -= len as i32 },
			Cigar::Del(len) => { offset += len as i32 },
			_ => {},
		}
	}
	return offset;
}

fn count_mismatches<'a, I>( tail_length : usize, cigar_iter: I, seq : &Seq, mut seq_index : usize, mut ref_offset : usize, chr_seq : &Vec<u8>, step : i32 ) -> i32
where I : ::std::iter::Iterator<Item = &'a Cigar>,
{

	let mut n_bases_examined  = 0usize;	
	let mut tail_mismatches   = 0usize;
	let mut tail_insertions   = 0usize;
	let mut tail_deletions    = 0usize;

	let mut print_after_indel = 0usize; //Debug

	'cigar_loop: for c in cigar_iter {		

		if n_bases_examined >= tail_length { break; }

		match *c {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {

				let l = len as usize;

				for k in 0..l {

					n_bases_examined += 1;
					//eprintln!( "seqindex: {}, ref_offset: {}", seq_index, ref_offset);

					//Get record nucleotide
					let read_acgt = encoded_base_to_char( seq.encoded_base( seq_index as usize));
					// Get reference genome nucleotide
					let ref_acgt = (chr_seq[ ref_offset+seq_index] as char).to_ascii_uppercase();

					//'N' counts as a match always
					let equal = read_acgt == ref_acgt || ref_acgt == 'N' || read_acgt == 'N';

					if !equal {
						tail_mismatches += 1;
						//eprintln!( "SI: {} '{}' not equal to RI: {} '{}' ({})", seq_index, read_acgt, ref_offset+seq_index, ref_acgt, n_bases_examined);
					}
					/*
					else {
						eprintln!( "SI: {} equal to RI: {} ({})", seq_index, ref_offset+seq_index, n_bases_examined);						
					}
					*/

					seq_index = ((seq_index as i32) + step) as usize;

					//DEBUG
					if false && print_after_indel > 0 {
						eprintln!( "{}: ref {} vs {} = {}", seq_index, ref_acgt, read_acgt, equal );
						print_after_indel -= 1;
					}			

					//When enough (tail_length) bases have been examined, we can exit 
					if n_bases_examined >= tail_length { break 'cigar_loop; }
				}
			},
			Cigar::Ins(len) => {
				
				//Count indels in tails as a single mismatch
				tail_mismatches += 1;
				tail_insertions += 1;					
				//n_bases_examined += 1; //Indels count as one
				
				let l = len as i32;
				seq_index = ((seq_index as i32) + (l * step)) as usize;				
				ref_offset = ((ref_offset as i32) - (l * step)) as usize; //add offset

				print_after_indel = 20; //DEBUG
			},
			Cigar::Del(len) => {

				//eprintln!( "DEL at SI: {} RI: {}", seq_index, ref_offset+seq_index);

				//Count indels in tails as a single mismatch				
				tail_mismatches += 1;
				tail_deletions += 1;
				//n_bases_examined += 1; //Indels count as one				

				let l = len as i32;
				ref_offset = ((ref_offset as i32) + (l * step)) as usize; //add offset

				print_after_indel = 20; //DEBUG
			},
			Cigar::RefSkip(len)  => panic!("Unexpected CIGAR type: N"),
			Cigar::SoftClip(len) => panic!("Unexpected CIGAR type: S"),
			Cigar::HardClip(len) => panic!("Unexpected CIGAR type: H"),
			Cigar::Pad(len)      => panic!("Unexpected CIGAR type: P"),
			Cigar::Back(len)     => panic!("Unexpected CIGAR type: 9")
		};
	} //cigar

	return tail_mismatches as i32;
}

#[derive(PartialEq)]
enum TailPos {
	None,
    Left,
    Right,
}

fn in_tail( index : usize, read_len : usize, tail_size : usize) -> TailPos {
	if index < tail_size { return TailPos::Left; }
	if index >= read_len - tail_size { return TailPos::Right; }
	return TailPos::None;
}

fn first_index_in_tail( index : usize, cigar_len: usize, tail_size : usize, total : usize) -> i32 {
	if index < tail_size { return index as i32; }
	if index + cigar_len > total - tail_size { return (index + ((index + cigar_len) - (total - tail_size))) as i32; }
	-1
}

fn read_record(bam: &mut bam::Reader, bam_path: &str, read: &mut Record) -> bool {
	if let Err(e) = bam.read(read) {
		match e {
		bam::ReadError::Truncated => {
			eprintln!("ERROR: Input file '{}' is truncated.", bam_path);
			exit(-1);
		},
		bam::ReadError::Invalid => {
			eprintln!("ERROR: Input file '{}' is corrupted.", bam_path);
			exit(-1);
		},
		bam::ReadError::NoMoreRecord => { return false; }
		}
	}
	return true;
}

fn check_threshold( tail_length : usize, threshold_ratio : f32) -> usize {

	let mut n_mismatches = 0usize;

	while n_mismatches <= tail_length {

		let cur_ratio = (n_mismatches as f32)  / (tail_length as f32);
		if cur_ratio >= threshold_ratio { break; }
		n_mismatches += 1;
	}

	return n_mismatches;
}

fn encoded_base_to_acgt_index(encoded: u8) -> usize {
	match encoded {
		1 => 0,   // A
		2 => 1,   // C
		4 => 2,   // G
		8 => 3,   // T
		_ => 4    // N
	}
}

fn encoded_base_to_char(encoded: u8) -> char {
	match encoded { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', _ => 'N' }
}

fn from_slice(bytes: &[u8]) -> [u8; 32] {
    let mut a = [0; 32];
    for i in 0..a.len() {
        // Panics if not enough input
        a[i] = bytes[i];
    }
    a
}

fn print_read(read: &Record, duplicate: bool) {
	let pos = if read.is_reverse() {
		read.cigar().end_pos().unwrap() - 1
	} else {
		read.pos()
	};
	eprintln!("chr?:{} ({}), insert size = {}: {}", pos,
		if read.is_reverse() { '-' } else { '+' }, read.insert_size(),
		if duplicate { "" } else { " ***" });
}

//Number of failed tests for flag "--test"
static mut FAILED_TESTS: i32 = 0;

//Run tests to see if discard tail artifacts is working as it
//is supposed to. For these tests, insertions and deletions count
//as a single mismatch. The tail length specifies the number of
//bases that are examined from each end of the read. Indels
//do not affect this number.
fn test() -> i32 {

	use rust_htslib::bam::record::{Cigar, CigarString};	

	//Create test data
	let qname = [b't',b'e',b's',b't']; //dummy
	let quality_arr = [b'!'; 50]; //dummy

	//Create identical ref and test sequences
	let mut seq_arr = [b'A'; 50];	
	let mut ref_seq = vec![b'A'; 50]; 

	//Create REcord for testing
	let mut test_rec = Record::new();
	test_rec.set( &qname, &CigarString( vec![Cigar::Match(50)]), &seq_arr, &quality_arr);

	eprintln!( "[{}] IDENTICAL LEFT: ", ftr( 0i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));	
	eprintln!( "[{}] IDENTICAL RIGHT", ftr( 0i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	eprintln!( "[{}] LONG TAIL LEFT", ftr( 0i32, count_mismatches( 100, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));	
	eprintln!( "[{}] LONG TAIL RIGHT", ftr( 0i32, count_mismatches( 100, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));
	
	ref_seq[ 0] = b'N';
	eprintln!( "[{}] REF N", ftr( 0i32, count_mismatches( 1, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));	
	//Insert N to test seq
	seq_arr[ 1] = b'N'; test_rec.set( &qname, &CigarString( vec![Cigar::Match(50)]), &seq_arr, &quality_arr);
	eprintln!( "[{}] SEQ N", ftr( 0i32, count_mismatches( 1, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	//Restore test seq back to original state
	seq_arr[ 1] = b'A'; test_rec.set( &qname, &CigarString( vec![Cigar::Match(50)]), &seq_arr, &quality_arr);

	ref_seq[ 0] = b'a'; ref_seq[ 49] = b'a';
	eprintln!( "[{}] LOWER CASE LEFT", ftr( 0i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));	
	eprintln!( "[{}] LOWER CASE RIGHT", ftr( 0i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	ref_seq[ 0] = b'C'; ref_seq[ 49] = b'C';
	eprintln!( "[{}] MISMATCH LEFT", ftr( 1i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] MISMATCH RIGHT", ftr( 1i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	ref_seq[ 19] = b'G'; ref_seq[ 30] = b'T';
	eprintln!( "[{}] 2 MISMATCH LEFT", ftr( 2i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] 2 MISMATCH RIGHT", ftr( 2i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	ref_seq[ 20] = b'G'; ref_seq[ 29] = b'T';
	eprintln!( "[{}] NON-TAIL MISMATCH LEFT", ftr( 2i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] NON-TAIL MISMATCH RIGHT", ftr( 2i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	ref_seq = vec![b'F'; 50];
	eprintln!( "[{}] NONSENSE LEFT", ftr( 20i32, count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] NONSENSE RIGHT", ftr( 20i32, count_mismatches( 20, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1)));

	eprintln!( "[{}] OVERLAP LEFT", ftr( 50i32, count_mismatches( 50, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] OVERLAP RIGHT",ftr( 50i32, count_mismatches( 50, test_rec.cigar().iter().rev(), &test_rec.seq(),49, 0, &ref_seq, -1))); 

	//New ref_seq with "mismatching matches" in the middle
	let mut ref_seq = vec![b'A'; 20];
	let mut piece2 = vec![ b'C';20];
	let mut piece3 = vec![ b'A';20];
	ref_seq.append( &mut piece2);
	ref_seq.append( &mut piece3);
	 
	let mut offset = count_right_end_offset( test_rec.cigar()) as usize;
	eprintln!( "[{}] REF OFFSET MATCH", ftr( 0i32, offset as i32));

	//Deletion of 10
	test_rec.set( &qname, &CigarString( vec![Cigar::Match(20), Cigar::Del(10), Cigar::Match(20) ]), &seq_arr, &quality_arr);
	offset = count_right_end_offset( test_rec.cigar()) as usize;
	eprintln!( "[{}] REF OFFSET DEL", ftr( 10i32, offset as i32));
	//eprintln!( "{}", test_rec.cigar() );

	eprintln!( "[{}] DELETION LEFT", ftr( 0i32,  count_mismatches( 20, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	eprintln!( "[{}] DELETION LEFT (2)", ftr( 2i32,  count_mismatches( 21, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));

	//Ten mismatches + 1 from del
	eprintln!( "[{}] INDEL RIGHT", ftr( 11i32, count_mismatches( 40, test_rec.cigar().iter().rev(), &test_rec.seq(), 49, offset, &ref_seq, -1)));

	test_rec.set( &qname, &CigarString( vec![Cigar::Match(20), Cigar::Ins(10), Cigar::Match(20) ]), &seq_arr, &quality_arr);
	offset = count_right_end_offset( test_rec.cigar()) as usize;
	eprintln!( "[{}] REF OFFSET INS", ftr( -10i32, offset as i32));
	eprintln!( "[{}] INSERTION LEFT", ftr( 2i32,  count_mismatches( 21, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	//20 from mismatches (at right end of tail) + 1 from insert
	eprintln!( "[{}] INSERTION RIGHT", ftr( 21i32, count_mismatches( 50, test_rec.cigar().iter().rev(), &test_rec.seq(), 49, offset, &ref_seq, -1)));

	test_rec.set( &qname, &CigarString( vec![Cigar::Match(20), Cigar::Ins(5), Cigar::Del(5), Cigar::Match(20) ]), &seq_arr, &quality_arr);
	offset = count_right_end_offset( test_rec.cigar()) as usize;
	eprintln!( "[{}] REF OFFSET INS-DEL", ftr( 0i32, offset as i32));
	//eprintln!( "{}", test_rec.cigar() );

	eprintln!( "[{}] INS-DEL LEFT", ftr( 5i32,  count_mismatches( 23, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1)));
	//15 mismatches + 1 Ins + 1 from del
	eprintln!( "[{}] INS-DEL RIGHT", ftr( 17i32, count_mismatches( 50, test_rec.cigar().iter().rev(), &test_rec.seq(), 49, offset, &ref_seq, -1)));

	//eprintln!( "{:?}", ref_seq );

	//Very simple performance test
	let start_time = Instant::now();
	let mut dummy = 0;
	let mut timeout = false;

	for i in 0..2000000 {
		dummy += count_mismatches( 50, test_rec.cigar().iter(), &test_rec.seq(), 0,  0, &ref_seq,  1);		
		dummy += count_mismatches( 50, test_rec.cigar().iter().rev(), &test_rec.seq(), 49, offset, &ref_seq, -1);
		if dummy >= 1000000 {
			dummy = 0;
			if start_time.elapsed().as_secs() >= 5 {
				timeout = true;
				break; 
			}
		}
	}

	let elapsed = start_time.elapsed().as_secs();
	eprintln!( "[{}] PERFORMANCE ({} seconds)", ftr( 1, (elapsed < 5 && !timeout) as i32), elapsed);
	let mut retval = 0;

	unsafe {
		if FAILED_TESTS == 0 { eprintln!( "INFO: All tests passed. [OK]"); }
		else { eprintln!( "ERROR: Number of failed tests: {}", FAILED_TESTS); }
		retval = -1*FAILED_TESTS
	};

	return retval; //TESTS OK == 0
}

//format_test_result
fn ftr( expect:i32, value:i32 ) -> String {

	use ansi_term::Colour::Red;
	use ansi_term::Colour::Green;

	if value == expect { return Green.paint("PASS").to_string();}
	unsafe { FAILED_TESTS += 1 };
	return Red.paint(format!("FAIL - {} vs. {}", value, expect)).to_string();
}

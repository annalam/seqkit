use parse_args;
use std::str;
use std::fs::File;
use std::process::exit;
use std::ascii::AsciiExt;

use bio::io::fasta;

use rust_htslib::bam;
use rust_htslib::bam::{Records, Read};
use rust_htslib::bam::record::{Cigar, Record};

use ErrorHelper;

pub struct RefGenomeReader {

	filepath: String,
	genome_reader: fasta::IndexedReader<File>,

}

impl RefGenomeReader {

    pub fn new( genome_fasta_path: &str) -> RefGenomeReader {

        RefGenomeReader {
            filepath: genome_fasta_path.to_string(),
            genome_reader: fasta::IndexedReader::from_file(&genome_fasta_path).on_error( &format!("Could not open genome FASTA file '{}'.", &genome_fasta_path))
        }
    }

	pub fn load_chromosome_seq( &mut self, chromosome_name : &str ) -> Vec<u8> {

    	let mut chr_seq: Vec<u8> = Vec::new();        
    	self.genome_reader.read_all( &chromosome_name, &mut chr_seq);
        eprintln!("INFO: Loaded chromosome '{}', len: {}", &chromosome_name, chr_seq.len());
    	return chr_seq; //transfer ownership
    }

}

/*
trait LoadChromosomeSeq {

	//fn new( &self, genome_fasta_path: &str) -> GenomeReader;
	fn load_chromosome_seq( &mut self, chromosome_name : &str) -> Vec<u8>;
}
*/

/*
impl LoadChromosomeSeq for RefGenomeReader {

    fn load_chromosome_seq( &mut self, chromosome_name : &str ) -> Vec<u8> {

    	let mut chr_seq: Vec<u8> = Vec::new();
    	self.genome_reader.read_all( &chromosome_name, &mut chr_seq);
    	return chr_seq; //transfer ownership
    }
}
*/
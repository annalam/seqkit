
use parse_args;
use {read_buffered, AsciiBufRead};
use std::str;
use std::process::exit;
use ascii::{AsciiString, AsciiChar};
use ErrorHelper;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufWriter};
use std::collections::HashMap;
use bio::io::fasta;


const USAGE: &'static str = "
Usage:
  fasta mapq track [options] <genome> <min_baseq>

 Options:
    --win-size=N     window size for read aligment [default: 48]
    --sliding        sliding window mode
";

pub fn main() {
	let args = parse_args(USAGE);
    let genome_path = args.get_str("<genome>");
    let win_size: usize = args.get_str("<win_size>").parse().unwrap();
    let slide_window = args.get_bool("--sliding");


    let fasta = fasta::Reader::from_file(format!("{}.fa", genome_path))
		.on_error(&format!("Genome FASTA file {}.fa could not be read.", genome_path));
	eprintln!("Reading reference genome into memory...");

	let mut genome = HashMap::new();
	for entry in fasta.records() {
		let chr = entry.unwrap();
		genome.insert(chr.id().to_owned(), chr.seq().to_owned());
	}


    // to make chromosome number order sequence
    let chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr16",
                        "chr17", "chr18","chr19","chr20","chr21","chr22","chrX","chrY"];



    // starting from chromosome 1, towards 22 then X and Y
    for ch  in chromosomes.iter() {
        for (chr, seq) in &genome {
            let chrm = chr.trim();
            if &ch.to_string() == chrm {
                let read = String::from_utf8(seq.to_owned()).unwrap();

                let mut qual = Vec::new();
                if slide_window {
                    qual = sliding(&read, win_size);
                } else {
                    qual = moving(&read, win_size, &genome_path);
                }
                println!("{}\tafter parse \t{}", chrm, read.len());
            }
        }
    }


    /*
	let mut line = AsciiString::new();
	while fasta_file.read_ascii_line(&mut line) {
		if line[0] != '@' {
			eprintln!("Invalid FASTQ format encountered."); exit(-1);
		}
		print!("{}", line);
		fasta_file.read_ascii_line(&mut line);
		let mut seq = line.clone();
		fasta_file.read_ascii_line(&mut line);
		fasta_file.read_ascii_line(&mut line);
		let qual = &line;

		for k in 0..seq.len() {
			// FIXME: We assume Sanger format base qualities here.
			if (qual[k] as i32 - 33 as i32) < min_baseq && seq[k] != '\n' {
				seq[k] = AsciiChar::N;
			}
		}

		print!("{}+\n{}", seq, qual);
	}
    */
}

/// Compute mapping quality for reference genome
/// using sliding window of given size
fn sliding(seq: &str, window_size: usize) -> Vec<usize> {
    let mut qual: Vec<usize> = Vec::new();
    let mut strt: usize = 0;
    let mut endn  = strt + window_size;
    let mut read = seq.get(strt..endn).unwrap();

    println!("{:?}", read.len());
    unimplemented!()
}


/// Compute mapping quality for reference genome
/// using moving window of given size
fn moving(seq: &str, window_size: usize, genome_path: &str) -> Vec<usize> {
    let mut qual: Vec<usize> = Vec::new();
    let mut strt: usize = 0;
    let chrsize = seq.len();

    loop {
        let mut endn  = strt + window_size;

        if strt >= 0 {
            let mut read = seq.get(strt..endn).unwrap();
            // println!("{:?}\t{}", read.len(), &read);
            // println!("{:?}\t{}\t{}", strt,  endn, chrsize);
            let q = align(read, &genome_path);
            qual.push(q);
            strt += window_size;
        } if strt >= chrsize {
            break;
        }
    }
    qual
}

/// Compute mapping quality for given read
/// using bowtie1
fn align(read: &str, genome_path: &str) -> usize {
    let mut q: usize = 255;
    eprintln!("Aligning {} bp windows against the genome...", read.len());
    let bowtie = Command::new("bowtie")
        .args(&["-f", "-p1", "-v0", "-m1", "-B1", "--suppress", "5,6,7,8", &genome_path, "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start Bowtie process.");

    // let mut bowtie_in = BufWriter::new(bowtie.stdin.unwrap());
    let bowtie_out = BufReader::new(bowtie.stdout.unwrap());

    println!("{:?}", bowtie_out);
    q
}

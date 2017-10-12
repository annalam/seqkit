
use parse_args;
use std::str;
use ErrorHelper;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Write};
use std::collections::HashMap;
use bio::io::fasta;


const USAGE: &'static str = "
Usage:
  fasta mapq track [options] <genome>

 Options:
    --win-size=N     window size for read aligment [default: 48]
    --sliding        sliding window mode
";


pub fn main() {
    let args = parse_args(USAGE);
    let genome_path = args.get_str("<genome>");
    let win_size: i32 = args.get_str("--win-size").parse().unwrap();
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
    let mut qual = HashMap::new();

    for ch  in chromosomes.iter() {
        for (chr, seq) in &genome {
            let chrm = chr.trim();

            if &ch.to_string() == chrm {
                let read   = String::from_utf8(seq.to_owned()).unwrap();

                if slide_window {
                    qual = sliding(&chr, &seq, win_size);
                } else {
                    let quals = moving(&chr, &seq, win_size, &genome_path);
                    for (key, val) in quals {
                        qual.insert(key, val);
                    }
                }
            }

        } // end-for-inner
    }
    // end-for-outer
}

/// Compute mapping quality for given read
/// using bowtie1
fn align(read: String, genome_path: &str) -> usize {
    let mut q: usize = 0;

    let bowtie = Command::new("bowtie")
        .args(&["-p1", "-a", &genome_path, "-f", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start Bowtie process.");


    match bowtie.stdin.unwrap().write_all(read.as_bytes()) {
        Err(why) => panic!("input not sent to bowtie"),
        Ok(_)    => println!("Bowtie recived input"),
    }

    let bowtie_out = BufReader::new(bowtie.stdout.unwrap());

    for l in bowtie_out.lines() {
        let line = l.unwrap();
        q += 1;
    }

    q
}


/// Compute mapping quality for reference genome
/// using moving window of given size
fn moving(chr: &str, seq: &[u8], window_size: i32, genome_path: &str) -> HashMap<String, usize> {
    let mut qual = HashMap::new();
    let mut strt: usize = 0;
    let mut endn: usize = 0;
    let chrsize = seq.len();

    while strt <= chrsize {
        endn = strt + window_size as usize;

        let tmp = seq.get(strt..endn).unwrap();
        let header  = ">".to_string() + chr + ":" + &strt.to_string() + "-" + &endn.to_string();
        let read = header.to_owned() + "\n"+ &String::from_utf8(tmp.to_vec()).unwrap();
        let window  = chr.to_owned() + ":" + &strt.to_string() + "-" + &endn.to_string();
        let q = align(read, &genome_path);

        println!("{:?}\t{:?}", &window, &q);
        qual.insert(window, q);
        strt += window_size as usize;
    }
    qual // return quality HashMap
}


/// Compute mapping quality for reference genome
/// using sliding window of given size
fn sliding(chr: &str,  seq: &[u8], window_size: i32) -> HashMap<String, usize> {
    /*
    let mut qual: Vec<usize> = Vec::new();
    let mut strt: usize = 0;
    let mut endn  = strt + window_size as usize;
    let read = seq.get(strt..endn).unwrap();

    println!("{:?}", read.len());
    */
    unimplemented!()
}


use parse_args;
use std::str;
use std::thread;
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
    //let mut qual = HashMap::new();
    for ch in &chromosomes {
        let s   = genome.get(&ch.to_string());
        let seq = s.unwrap();
        let ref_genome = genome_path.to_owned();
        if slide_window {
            sliding(&ch, &seq, win_size);
        } else {
            moving(&ch, &seq, win_size, ref_genome);
        }
    }
}


/// Compute mapping quality for given read
/// using bowtie1
fn align(read: String, genome_path: &str) -> usize {
    let read = read.to_owned();
    let genome_path  = genome_path.to_owned();
    let bowtie = Command::new("bowtie")
        .args(&["-p1", "-a", &genome_path, "-f", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start Bowtie process.");


    match bowtie.stdin.unwrap().write_all(read.as_bytes()) {
        Err(_why) => panic!("input not sent to bowtie"),
        Ok(_)    => print!(""),
    }

    let bowtie_out = BufReader::new(bowtie.stdout.unwrap());
    let q = bowtie_out.lines().count();
    q
}


/// Compute mapping quality for reference genome
/// using moving window of given size
fn moving(chr: &str, seq: &[u8], window_size: i32, genome_path: String) {
    let mut strt: usize = 0;
    let mut endn: usize = 0;
    let ref_genome = genome_path.to_owned();
    while strt <= seq.len() {
        endn = strt + window_size as usize;

        let tmp = seq.get(strt..endn).unwrap();
        let header  = ">".to_string() + chr + ":" + &strt.to_string() + "-" + &endn.to_string();
        let read = header.to_owned() + "\n"+ &String::from_utf8(tmp.to_vec()).unwrap();
        let window  = chr.to_owned() + ":" + &strt.to_string() + "-" + &endn.to_string();

        let q = align(read, &ref_genome);
        println!("{}\t{:?}", &window, &q);

        strt += window_size as usize;
    }
}


/// Compute mapping quality for reference genome
/// using sliding window of given size
fn sliding(chr: &str,  seq: &[u8], window_size: i32){
    /*
    let mut qual: Vec<usize> = Vec::new();
    let mut strt: usize = 0;
    let mut endn  = strt + window_size as usize;
    let read = seq.get(strt..endn).unwrap();

    println!("{:?}", read.len());
    */
    unimplemented!()
}

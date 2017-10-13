
use parse_args;
use std::str;
use std::thread;
use std::sync::{Arc, Mutex};
use ErrorHelper;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::collections::HashMap;
use std::fs::File;
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
    let win_size: usize = args.get_str("--win-size").parse().unwrap();
    let sliding = args.get_bool("--sliding");

    let genome = genome_path.to_owned();

    let bowtie = Command::new("bowtie")
        .args(&["-p1", "-a", &genome_path, "-f", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start Bowtie process.");

    let mut bowtie_in = BufWriter::new(bowtie.stdin.unwrap());
    let bowtie_out = BufReader::new(bowtie.stdout.unwrap());

    thread::spawn(move || {
        eprintln!("spawning a new thread for bowtie!");
        eprintln!("Reading reference genome into memory...");
        let fasta = fasta::Reader::from_file(format!("{}.fa", &genome))
    		.on_error(&format!("Genome FASTA file {}.fa could not be read.", &genome));

        if sliding {
                println!("{}", "Sliding window mode selected!");
    			send_sliding_windows(fasta, &mut bowtie_in, win_size);
    		} else {
                println!("{}", "Moving window mode selected!");
    			send_moving_windows(fasta, &mut bowtie_in, win_size);
    		}
    });
}


fn send_moving_windows(fasta: fasta::Reader<File>, aligner_in: &mut Write, win_size: usize) {
    println!("running send_moving_windows function");
    let mut num_reads_sent = 0;
    let mut strt: usize = 0;
    let mut endn: usize = 0;
    println!("Not going beyond this point!");
    println!("{}\t////////////////////////", fasta.records().count());
    /*
    for entry in fasta.records() {
        let chr = entry.unwrap();
        println!("{}\t{}", chr.id(), chr.seq().len());

        let ch  = chr.id().to_owned();
        let seq = chr.seq().to_owned();

        while strt <= seq.len() {
            endn = strt + win_size as usize;
            let read = seq.get(strt..endn).unwrap();
            let window  = ch.to_owned() + ":" + &strt.to_string() + "-" + &endn.to_string();

            println!(">{:?}\n{:?}", window, read);
            write!(aligner_in, ">{}\n", window);
    		aligner_in.write_all(read);
            writeln!(aligner_in);
            num_reads_sent += 1;
            strt += win_size as usize;
        }
    }
    */
    println!("{:?}\t", num_reads_sent);
}

/// Compute mapping quality for reference genome
/// using sliding window of given size
fn send_sliding_windows(chr: fasta::Reader<File>, aligner_in: &mut Write, win_size: usize){
    unimplemented!()
}

/*
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
*/

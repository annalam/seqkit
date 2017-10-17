
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
    --win-size=N     window size for bowtie aligment (3-1023) [default: 48]
    --sliding        enable sliding window mode
";


pub fn main() {
    let args = parse_args(USAGE);
    let genome_path = args.get_str("<genome>");
    let win_size: usize = args.get_str("--win-size").parse().unwrap();
    let sliding = args.get_bool("--sliding");

    let genome = genome_path.to_owned();
    let bowtie = Command::new("bowtie")
        .args(&["-p1", "-k10", &genome_path, "-f", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start Bowtie process.");

    let mut bowtie_in = BufWriter::new(bowtie.stdin.unwrap());
    let bowtie_out = BufReader::new(bowtie.stdout.unwrap());

    let child = thread::spawn(move || {
        eprintln!("spawning a new thread for bowtie!");
        eprintln!("Reading reference genome into memory...");
        let fasta = fasta::Reader::from_file(format!("{}.fa", &genome))
    		.on_error(&format!("Genome FASTA file {}.fa could not be read.", &genome));

        send_seq_slices(fasta, &mut bowtie_in, win_size, sliding);
    });

    let mut prev = String::new();
    let mut prev_read_count = 0;

    for l in bowtie_out.lines() {
        let line = l.unwrap();
        let mut cols = line.split('\t');
        let window = cols.nth(0).unwrap().to_string();
        let mut reads_count = 0;

        if prev.is_empty() && prev_read_count == 0 {
            reads_count += 1;
            prev = window;
            prev_read_count = reads_count;
        } else if window == prev && prev_read_count > 0 {
            prev_read_count += 1;
            prev = window;
        } else if window != prev && prev_read_count > 0 {
            println!("{}\t{}", &window, &prev_read_count);
            reads_count += 1;
            prev_read_count = reads_count;
            prev = window;
        } else {
            println!("{}\tSomthing else happend", line);
        }
    }

    let res = child.join();
}


/// Compute mapping quality for reference genome
/// using moving window of given size
fn send_seq_slices(fasta: fasta::Reader<File>, aligner_in: &mut Write, win_size: usize, sliding: bool) {
    if sliding {
        eprintln!("running sliding-window mode");
    } else {
        eprintln!("running moving-window mode");
    }

    for entry in fasta.records() {
        let mut strt: usize = 0;
        let mut endn: usize = 0;
        let chr = entry.unwrap();
        eprintln!("{}\t{}", chr.id(), chr.seq().len());
        let ch  = chr.id().to_owned();
        let seq = chr.seq().to_owned();

        while strt + 48 <= seq.len() + 1 {
            endn = strt + win_size as usize + 1;
            let read = seq.get(strt..endn).unwrap();
            let window  = ch.to_owned() + ":" + &strt.to_string() + "-" + &endn.to_string();

            //println!("{:?}\n{:?}", &window, &read);
        	write!(aligner_in, ">{}:\n", &window);
    		aligner_in.write_all(&read);
            if sliding {
                strt += 1;
            } else {
                strt += win_size as usize;
            }
        }
        eprintln!("Processing {}\tcompleted!", &ch);
    }
}

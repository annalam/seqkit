
use parse_args;
use std::str;
use std::thread;
use ErrorHelper;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::collections::HashSet;
use std::fs::File;
use bio::io::fasta;


const USAGE: &'static str = "
Usage:
  fasta pairwise similarity <sv_file> <genome> <bed>
";


pub fn main() {
    let args = parse_args(USAGE);
    let sv_path         = args.get_str("<sv_file>");
    let genome_path     = args.get_str("<genome>");
    let bed_path        = args.get_str("<bed>");

    let genome = genome_path.to_owned();

    let blast = Command::new("blast")
        .args(&["-query", &genome_path, "-subject", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start BLAST process.");

    let mut blast_in = BufWriter::new(blast.stdin.unwrap());
    let blast_out = BufReader::new(blast.stdout.unwrap());

    let child = thread::spawn(move || {
        eprintln!("spawning a new thread for blast!");
        eprintln!("Reading reference genome into memory...");
        let fasta = fasta::Reader::from_file(format!("{}.fa", &genome))
    		.on_error(&format!("Genome FASTA file {}.fa could not be read.", &genome));
        if !sv_path.is_empty() {
            send_list_slices(fasta, &mut blast_in);
        } else {
            send_seq_slices(fasta,  &mut blast_in);
        }
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
            println!("{}\t{}", &window.trim_right_matches(':'), &prev_read_count);
            reads_count += 1;
            prev_read_count = reads_count;
            prev = window;
        } else {
            println!("{}\tSomthing else happend", line);
        }
    }

    let _res = child.join();
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
        let chr = entry.unwrap();
        eprintln!("{}\t{}", chr.id(), chr.seq().len());
        let ch  = chr.id().to_owned();
        let seq = chr.seq().to_owned();

        while strt + win_size as usize <= seq.len() + 1 {
            let endn = strt + win_size as usize;
            let read = seq.get(strt..endn).unwrap();
            // printing genome slice into 1-based co-ordinates
            let mut window = String::new();
            if sliding {
                window  = ch.to_owned() + ":" + &strt.to_string();
            } else {
                window  = ch.to_owned() + ":" + &strt.to_string() + "-" + &endn.to_string();
            }

            //println!(">{}\n{}", &window, String::from_utf8(read.to_owned()).unwrap());
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


/// compute mappbility score for list of input chromosome positions
/// each line in input line should be in chrX:NNNNNNN format
fn send_list_slices(fasta: fasta::Reader<File>, aligner_in: &mut Write, list_pos: Vec<String>, win_size: usize) {

    for entry in fasta.records() {
        let chr = entry.unwrap();
        eprintln!("{}\t{}", chr.id(), chr.seq().len());
        let ch  = chr.id().to_owned();
        let seq = chr.seq().to_owned();

        for line in list_pos.iter() {
            let chr         = line.split(':').nth(0).unwrap().to_string();
            let pos: usize  = line.split(':').nth(1).unwrap().parse().unwrap();
            let mut strt: usize = pos - (win_size / 2) ;

            if ch == chr && strt + win_size <= seq.len() {
                let endn = strt + win_size as usize;
                let read = seq.get(strt..endn).unwrap();
                // printing genome slice into 1-based co-ordinates
                // output id will have input position
                let window  = ch.to_owned() + ":" + &pos.to_string();
                //println!(">{}\n{}", &window, String::from_utf8(read.to_owned()).unwrap());
                write!(aligner_in, ">{}:\n", &window);
                aligner_in.write_all(&read);
            }
        }
        eprintln!("Processing {}\tcompleted!", &ch);
    }
}


use parse_args;
use std::io;
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
  fasta pairwise blast <genome> <bed> <fusions_table>
";


pub fn main() {
    let args = parse_args(USAGE);
    let genome_path     = args.get_str("<genome>");
    let bed_path        = args.get_str("<bed>");
    let fusion_path     = args.get_str("<fusions_table>");

    let genome = genome_path.to_owned();

    let mut list_bed = Vec::<String>::new();
    if !bed_path.is_empty() {
        let ls = BufReader::new(File::open(&bed_path).on_error(
            &format!("Could not open BED file '{}'.", bed_path)));

        for l in ls.lines() {
            let line = l.unwrap().to_string();
            list_bed.push(line);
        }
    }

    let mut list_fusions = Vec::<String>::new();
    if !fusion_path.is_empty() {
        let ls = BufReader::new(File::open(&fusion_path).on_error(
            &format!("Could not open fusions table file '{}'.", fusion_path)));

        for l in ls.lines() {
            let line = l.unwrap().to_string();
            list_fusions.push(line);
        }
    }

    println!("{}\t{}", &list_fusions.len(), &list_bed.len());

    let blast = Command::new("blastn")
        .args(&["-query", &genome_path, "-subject", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped()).spawn()
        .on_error("Could not start BLAST process.");

    let mut blast_in = BufWriter::new(blast.stdin.unwrap());
    let blast_out = BufReader::new(blast.stdout.unwrap());

    let child = thread::spawn(move || {
        eprintln!("spawning a new thread for blast!");
        eprintln!("Reading reference genome into memory...");

        let mut fasta_reader = fasta::IndexedReader::from_file(&format!("{}.fa", &genome))
    		.on_error(&format!("Genome FASTA file {}.fa could not be read.", &genome));

        send_gene_seqs(&mut fasta_reader, /*&mut blast_in,*/ list_fusions, list_bed);

    });

    let mut prev = String::new();
    let mut prev_read_count = 0;

    for l in blast_out.lines() {
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


fn gene_boundaries(gene: String, bedlines: &Vec<String>) -> (u64, u64){
    let mut start: u64 = 0;
    let mut end:   u64 = 0;

    for line in bedlines {
        let cols: Vec<&str> = line.split_whitespace().collect();
        let name = cols[3].to_string();

        if gene == name {
            start = cols[1].to_string().parse::<u64>().unwrap();
            end   = cols[2].to_string().parse::<u64>().unwrap();
        }
    }
    (start, end)
}


/// compute mappbility score for list of input chromosome positions
/// each line in input line should be in chrX:NNNNNNN format
fn send_gene_seqs(fasta_reader: &mut fasta::IndexedReader<File>, /*aligner_in: &mut Write,*/ list_fusions: Vec<String>, list_bed: Vec<String>) {

    for line in list_fusions.iter() {
        let cols: Vec<&str> = line.split('\t').collect();
        let pair: String = cols[0].to_string();
        let pos1: String = cols[1].to_string();
        let pos2: String = cols[2].to_string();

        let genes: Vec<&str> = pair.split(':').collect();
        let gene_a: String = genes[0].to_string();
        let gene_b: String = genes[1].to_string();

        let mut tmp1  =  pos1.split(':');
        let chr1      = tmp1.next().unwrap();

        let mut tmp2  =  pos2.split(':');
        let chr2      = tmp2.next().unwrap();

        let header_a  = &gene_a.to_owned();
        let header_b  = &gene_b.to_owned();

        let (gene_a_start, gene_a_end) = gene_boundaries(gene_a, &list_bed);
        let (gene_b_start, gene_b_end) = gene_boundaries(gene_b, &list_bed);

        let mut seq_a: Vec<u8> = Vec::new();
        let mut seq_b: Vec<u8> = Vec::new();
        fasta::IndexedReader::read(fasta_reader, chr1, gene_a_start, gene_a_end, &mut seq_a);

        fasta::IndexedReader::read(fasta_reader, chr2, gene_b_start, gene_b_end, &mut seq_b);

        print!(">{}\t{}\t", header_a, seq_a.len());
        println!(">{}\t{}", header_b, seq_b.len());
    }
}

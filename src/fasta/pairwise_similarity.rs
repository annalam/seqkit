
use parse_args;
use std::str;
use std::thread;
use ErrorHelper;
use std::io::{BufReader, BufRead};
use std::fs::File;
use bio::io::fasta;
use parasailors::*;

const USAGE: &'static str = "
Compute pairwise sequence similarity between fusion partners
using vectorized version of Smith-Waterman algorithm

Usage:
  fasta pairwise similarity <genome> <fusions_table>
";


pub fn main() {
    let args = parse_args(USAGE);
    let genome_path     = args.get_str("<genome>");
    let fusion_path     = args.get_str("<fusions_table>");

    let genome = genome_path.to_owned();

    let mut list_fusions = Vec::<String>::new();
    if !fusion_path.is_empty() {
        let ls = BufReader::new(File::open(&fusion_path).on_error(
            &format!("Could not open fusions table file '{}'.", fusion_path)));

        for l in ls.lines() {
            let line = l.unwrap().to_string();
            list_fusions.push(line);
        }
    }

    let child = thread::spawn(move || {
        eprintln!("spawning a new thread for local alignment!");
        eprintln!("Reading reference genome into memory...");

        let mut fasta_reader = fasta::IndexedReader::from_file(&format!("{}.fa", &genome))
    		.on_error(&format!("Genome FASTA file {}.fa could not be read.", &genome));

        send_gene_seqs(&mut fasta_reader, list_fusions);

    });
    let _res = child.join();
}


/// compute sequence similarity score
/// each line in input file should be in below format
/// A2ML1:KLRG1	chr12:8950032 (+)	chr12:9017889 (+)	TCGA-HNSC	TCGA-CR-7364-01A-11R-2016-07
fn send_gene_seqs(fasta_reader: &mut fasta::IndexedReader<File>, list_fusions: Vec<String>) {

    for line in list_fusions.iter() {
        let cols: Vec<&str> = line.split('\t').collect();
        let pair: String = cols[0].to_string();
        let pos1: String = cols[1].to_string();
        let pos2: String = cols[2].to_string();

        let genes: Vec<&str> = pair.split(':').collect();
        let gene_a: String = genes[0].to_string();
        let gene_b: String = genes[1].to_string();

        let mut tmp1_chr  =  pos1.split(':');
        let chr1          =  tmp1_chr.next().unwrap();
        if  chr1.contains("chrM") { continue; }
        let pos_a         =  tmp1_chr.next().unwrap().to_string();

        let mut bp1       = pos_a.split_whitespace();
        let breakpoint1   =    bp1.next().unwrap().parse::<u64>().unwrap();

        let gene_a_start: u64       =   breakpoint1 - 200;
        let gene_a_end: u64         =   breakpoint1 + 200;

        let mut tmp2_chr  =  pos2.split(':');
        let chr2          =  tmp2_chr.next().unwrap();
        if  chr2.contains("chrM") { continue; }
        let pos_b         =  tmp2_chr.next().unwrap().to_string();

        let mut bp2       = pos_b.split_whitespace();
        let breakpoint2   =    bp2.next().unwrap().parse::<u64>().unwrap();

        let gene_b_start: u64       =   breakpoint2 - 200;
        let gene_b_end: u64         =   breakpoint2 + 200;

        // let header_a  = &gene_a.to_owned();
        // let header_b  = &gene_b.to_owned();

        let mut seq_a = Vec::new();
        let mut seq_b = Vec::new();
        fasta::IndexedReader::read(fasta_reader, chr1, gene_a_start, gene_a_end, &mut seq_a);
        fasta::IndexedReader::read(fasta_reader, chr2, gene_b_start, gene_b_end, &mut seq_b);

        let seq_size: f64 = seq_a.len() as f64;
        let identity_matrix = Matrix::new(MatrixType::Identity);
        let profile = Profile::new(&seq_a[..], &identity_matrix);

        let score: f64 = local_alignment_score(&profile, &seq_b[..], 1, 1).to_string().parse::<f64>().unwrap();

        let similarity = format!("{:.*}", 2, 100.0 * (score / seq_size));
        println!("{}\t{}", &line, &similarity);
    }
}

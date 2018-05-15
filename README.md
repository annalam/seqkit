# Seqkit

Seqkit is a suite of software utilities for manipulating and analyzing common genome sequencing data types (FASTA, SAM). Seqkit is written in Rust, and uses rust-htslib for reading and writing BAM files. Seqkit is divided into two utilities: `fasta` and `sam`. Each utility provides various useful subcommands. For a complete listing type the command name without arguments into your shell.

Features
--------

For FASTA/FASTQ files, Seqkit can:
- Convert between FASTA, FASTQ and raw sequence-per-line formats
- Extract sample barcodes and UMIs from multiplexed FASTQ sequencing data
- Demultiplex FASTQ sequencing data
- Trim FASTQ sequencing data by per-base quality (BASEQ) values
- Mask low quality bases in FASTQ sequencing data
- Replace FASTQ read identifiers with compact numeric IDs

For BAM files, Seqkit can:
- Calculate a histogram of fragment lengths
- Calculate statistics about unaligned, aligned and duplicate-flagged reads
- Extract reads from position-sorted, name-sorted or unsorted BAM files

Installation
------------

Install Rust (version 1.24 or later). Download the Seqkit source code from Github, then run the following command to install the binaries:
```
cargo install --path <seqkit_source_directory>
```

Examples
--------

Extract UMIs and demultiplex Illumina sequencing data where both the sample barcode and UMI are stored in the adapter:
```
fasta demultiplex sample_sheet.tsv
  <(fasta simplify read ids multiplexed_R1.fq.gz | fasta add barcode - multiplexed_I1.fq.gz UUUUSSSS)
  <(fasta simplify read ids multiplexed_R2.fq.gz | fasta add barcode - multiplexed_I1.fq.gz UUUUSSSS)
```

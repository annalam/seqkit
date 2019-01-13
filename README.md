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
- Extract reads from name-sorted or position-sorted BAM files
- Calculate a histogram of fragment lengths
- Calculate statistics about unaligned, aligned and duplicate-flagged reads

Installation
------------

Install Rust (version 1.31 or later). Then run the following command:
```
cargo install --force --git https://github.com/annalam/seqkit
```

Examples
--------

#### Extracting reads
Extract reads from a name-sorted or position-sorted BAM file called `tumor.bam`. Paired end reads are written in gzip-compressed FASTQ format into output files `tumor_1.fq.gz` and `tumor_2.fq.gz`, which are automatically created. Orphan reads are written into output file `tumor.fq.gz`. The second parameter specifies the prefix used for the output file names.
```
sam to fastq tumor.bam tumor
```

#### Demultiplex a pooled sequencing run
Extract UMIs and demultiplex Illumina sequencing data where both the sample barcode and UMI are stored in the adapter:
```
fasta demultiplex sample_sheet.tsv
  <(fasta add barcode multiplexed_R1.fq.gz multiplexed_I1.fq.gz)
  <(fasta add barcode multiplexed_R2.fq.gz multiplexed_I1.fq.gz)
```

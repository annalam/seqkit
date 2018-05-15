# Seqkit

Seqkit is a suite of software utilities for manipulating and analyzing common genome sequencing data types (FASTA, SAM). Seqkit is written in Rust, and uses rust-htslib for reading and writing BAM files. Seqkit is divided into two utilities: `fasta` and `sam`. Each utility provides various useful subcommands - for a complete listing type the command name without arguments into your shell:
```
$ fasta

Usage:
  fasta to raw <fasta/fastq>
  fasta simplify read ids <fastq_file>
  fasta trim by quality <fastq_file> <min_baseq>
  fasta mask by quality <fastq_file> <min_baseq>
  fasta mappability track <genome>
  fasta add barcode <fastq_file> <barcode_file> <barcode_format>
  fasta convert basespace <fastq_file>
  fasta demultiplex <sample_sheet> <fastq_1> <fastq_2>
  fasta statistics <fastq_file>
```


Examples
--------

Extract UMIs and demultiplex Illumina sequencing data where both the sample barcode and UMI are stored in the adapter:
```
fasta demultiplex sample_sheet.tsv
  <(fasta simplify read ids multiplexed_R1.fq.gz |
    fasta add barcode multiplexed_R1.fq.gz multiplexed_I1.fq.gz UUUUSSSS)
  <(fasta simplify read ids multiplexed_R2.fq.gz |
    fasta add barcode - multiplexed_I1.fq.gz UUUUSSSS)
```

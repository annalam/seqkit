# Seqkit

Seqkit is a suite of software utilities for manipulating and analyzing common genome sequencing data types (FASTA, SAM). Seqkit is written in Rust, and uses rust-htslib for reading and writing BAM files. Seqkit is divided into two utilities: `fasta` and `sam`. Each utility provides various useful subcommands. For a complete listing type the command name without arguments into your shell.


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

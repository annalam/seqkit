# seqkit

**Bioinformatics sequence manipulation utilities written in Rust**

**Dependencies**

Install [Rust](https://www.rust-lang.org/en-US/)

**Installation**

  ```sh
    git clone https://github.com/annalam/seqkit.git
    cd seqkit
    cargo install --root=/your/tools/path/seqkit
  ```


sam
===

**utils for working with SAM, BAM and CRAM files**

Usage:
 ```sh
  sam count <bam_file> <regions.bed>
  sam fragments <bam_file>
  sam fragment lengths <bam_file>
  sam statistics <bam_file>
  sam mark duplicates <bam_file>

 ```

fasta
=====
**utils for working with FASTA and FASTQ files**

Usage:
```sh
  fasta to raw <fasta/fastq>
  fasta trim by quality <fastq_file> <min_baseq>
  fasta mask by quality <fastq_file> <min_baseq>
  fasta mappability track <genome>
```

**fasta mappability track**

Usage:
```sh
  fasta mappability track [options] <genome>

 Options:
    --win-size=N     window size for bowtie aligment (4-1024) [default: 48]
    --sliding        enable sliding window mode
```

- Building mappability track for a given reference genome in FASTA format
- Moving windows of sequence slices aligned against the genome `default`
- `--sliding` will enable sliding window mode (`will take longer time`)
- Bowtie1 is used for alignment
- window size can be adjusted at running time with `--win-size`
- window size can range from `4` to `1024` (limitation by `bowtie1`)

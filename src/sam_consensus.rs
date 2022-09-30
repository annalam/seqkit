
use crate::common::{parse_args, BamReader, BamWriter};
use std::str;
use std::io::{self, Write};
use std::collections::VecDeque;
use std::collections::HashMap;
use std::cmp::{min, max};
use rust_htslib::bam::record::{Record, Aux, Cigar, CigarString};

const USAGE: &str = "
Usage:
  sam consensus [options] <bam_file>

Options:
  --uncompressed      Output in uncompressed BAM format
  --ignore-umi        Ignore UMI stored in RX tag even if present
  --min-evidence=N    Minimum evidence level [default: 1]
  --max-len=N         Maximum fragment length [default: 5000]
  --min-mapq=N        MapQ threshold for mates for being eligible
                      to consensus processing [default: 5]
  --keep-discordant   Include discordant reads that could not be merged
                      into consensus fragments in the output BAM file.
                      They will be marked with a QC_FAIL (0x200) flag.
  --human-readable    Output consensus alignments as human readable text
                      instead of the default BAM format
                      
Generates consensus DNA fragments based on redundant paired end reads.
Currently only Illumina-style sequencing data (i.e. converging orientation)
is supported. 

Any fragment bases with a below 80% consensus across the redundant reads
are replaced with an ambiguous N nucleotide.

The input BAM file must be position-sorted. Consensus fragments are written
to the standard output in BAM format. An auxiliary field DP is added to each
consensus fragment, describing the number of duplicate read pairs that were
used in generating the consensus.

Each called consensus base is assigned an evidence level:
  1: read in one direction, no duplicates
  2: read in one direction, with duplicates
  3: read in both directions, no duplicates
  4: read in both directions, has duplicates
  5: read in both directions, has duplicates, duplex strand consensus

When viewing consensus BAM files in IGV, the consensus read lengths may
cause IGV to automatically enable a visualization mode meant for nanopore
sequencing reads. You can disable this by choosing 'Experiment type -> Other'
and 'Quick consensus mode -> Off' under the right-click menu.
";

// The input BAM file is position-sorted, and we want the output BAM to
// preserve that property. To achieve that, we serially read through all reads
// in the input BAM file, adding new reads to the end of FIFO queue, and
// releasing and outputting reads from the beginning of the FIFO queue once
// they have been processed.
//
// The ReadPair objects can be in one of several states:
// 1) A single read from a pair that is not amenable for consensus generation.
//    Stored in FIFO to preserve ordering, and written out as-is.
// 2) An incomplete read pair containing only the leftmost mate. The right mate
//    will be added to this pair when it is encountered in the input BAM file.
// 3) A complete read pair amenable for consensus generation. Merged into a
//    consensus fragment once we can prove that no further duplicates can be
//    encountered in the BAM file. Only the resulting consensus fragment
//    gets written into the output BAM file, not the original pairs themselves. 
#[derive(Debug)]
struct ReadPair {
	r1: Record,
	r2: Record,
	left_pos: u32,  // Leftmost 1-based inclusive coordinate of fragment
	right_pos: u32, // Rightmost 1-based inclusive coordinate of fragment
	umi: Box<[u8]>,
	strand: u8,   // Strand of the sequenced ssDNA fragment (true = plus)
}

impl ReadPair {
	fn invalid(r: Record) -> ReadPair {
		ReadPair { r1: r, r2: Record::new(),
			left_pos: u32::MAX, right_pos: 0, strand: b'+', umi: Box::new(*b"") }
	}
	fn is_ready(&self) -> bool { self.left_pos > 0 }
	fn is_invalid(&self) -> bool { self.left_pos == u32::MAX }
	fn is_merged(&self) -> bool { self.left_pos == u32::MAX - 1 }
	fn mark_invalid(&mut self) { self.left_pos = u32::MAX; }
	fn mark_merged(&mut self) { self.left_pos = u32::MAX - 1; }
}

// Altref represent a single genomic coordinate
#[derive(Default, Clone)]
struct Altref  {
	n_total: u32,  // includes N's
	alts: HashMap<String, u32>,
	umi_fwd: u32,  // Number of reads with non-flipped umi
	umi_rev: u32,  // Number of reads with flipped umi (in comparison to first dup found)
	fwd : u32,     // Evidence from fwd strand reads (excludes N)
	rev : u32     // Evidence from rev strand reads (excludes N)
}

#[derive(Default)]
struct Settings  {
	uncompressed: bool,
	print_alignment: bool,
	max_frag_len: u32,
	ignore_umi: bool,
	min_evidence: usize,
	min_mapq: u8,
	keep_discordant: bool,
	chr_names: Vec<String>
}

#[derive(Default)]
struct Stats {
	total_reads: usize,
	concordant: usize,

	// Different reasons for omitting a read from consensus generation
	unpaired: usize,
	low_mapq: usize,
	not_converging: usize,
	unmapped: usize,
	too_long: usize,
	diff_chr: usize
}


pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");

	let mut settings = Settings::default();
	let mut stats = Stats::default();

	settings.max_frag_len = args.get_str("--max-len").parse().unwrap_or_else(
		|_| error!("--max-len must be a positive integer."));

	settings.min_mapq = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be an integer 0-255."));
	eprintln!("Only reads with mapping quality {} or higher are processed for consensus.", settings.min_mapq);

	settings.ignore_umi = args.get_bool("--ignore-umi");
	settings.print_alignment = args.get_bool("--human-readable");
	settings.min_evidence = args.get_str("--min-evidence").parse().unwrap_or_else(|_| error!("--min-evidence must be a non-negative integer."));
	settings.uncompressed = args.get_bool("--uncompressed");
	settings.keep_discordant = args.get_bool("--keep-discordant");

	eprintln!("Writing {}compressed output.",
		if settings.uncompressed { "un" } else { "" });
 
	let bam = BamReader::open(&bam_path);
	let header = bam.header();
	settings.chr_names = header.target_names().iter()
		.map(|x| str::from_utf8(x).unwrap().into()).collect();

	let mut out = BamWriter::open(
		if settings.print_alignment { "/dev/null" } else { "-" },
		&header, !settings.uncompressed);

	// In order to preserve the order of BAM records, we need to maintain
	// a FIFO queue. We also maintain a hashmap that maps read identifiers into
	// FIFO queue positions. We use a variable "n_pairs_completed" to
	// keep track of how many items have been released from the beginning
	// of the FIFO queue, so we can always correctly index into the FIFO queue.
	let mut fifo: VecDeque<ReadPair> = VecDeque::new();
	let mut mates: HashMap<Box<str>, usize> = HashMap::new();
	let mut n_pairs_completed: usize = 0;

	let mut prev_chr: i32 = -1;
	let mut prev_pos: u32 = 0;

	for read in bam {
		if read.is_secondary() || read.is_supplementary() { continue; }

		stats.total_reads += 1;

		let chr: i32 = read.tid();
		let pos: u32 = read.pos() as u32 + 1;  // 1-based position

		if chr == -1 {
			// Chromosome ID can be -1 if read is unmapped. If so, ignore.
		} else if chr != prev_chr {
			// We just entered a new chromosome, so all reads in the FIFO queue
			// can be processed. If there were any incomplete read pairs in
			// this chromosome, we know they will never be amenable for
			// consensus generation, so we mark them as invalid.
			for pair in &mut fifo {
				if pair.is_ready() == false { pair.mark_invalid(); }
			}
			n_pairs_completed += write_consensus(&mut out, &mut fifo,
				u32::MAX, &settings);
			prev_chr = chr;

			eprintln!("Processing {}...", settings.chr_names[chr as usize]);
		} else if pos < prev_pos {
			error!("Input BAM file is not sorted by coordinate.");
		}

		prev_pos = pos;

		// Check if the read is valid for consensus generation
		let mut valid_for_consensus = if !read.is_paired() {
			stats.unpaired += 1; false
		} else if read.is_unmapped() || read.is_mate_unmapped() {
			stats.unmapped += 1; false
		} else if read.tid() != read.mtid() {
			stats.diff_chr += 1; false
		} else if read.is_reverse() == read.is_mate_reverse() {
			stats.not_converging += 1; false
		} else {
			true
		};

		// If the read is not valid for consensus generation, we store it
		// in the FIFO queue as an "invalid" ReadPair, which will simply be
		// written into the output BAM file when it reaches the start of the
		// FIFO queue.
		if valid_for_consensus == false {
			fifo.push_back(ReadPair::invalid(read));
			continue;
		} 

		// At this point we have confirmed that the read looks potentially
		// amenable for consensus generation. We now check if we have already
		// stored its mate in the FIFO queue.
		let qname = str::from_utf8(read.qname()).unwrap();
		if let Some(mate_idx) = mates.remove(qname) {
			// Mate was previously stored in the FIFO queue, so we add this
			// read to the pair, and mark as ready for processing.
			let pair = &mut fifo[mate_idx - n_pairs_completed];

			// Calculate fragment boundaries. We assume converging paired
			// end read orientation (i.e. Illumina paired end sequencing).
			// Filling in the "left_pos" and "right_pos" fields implicitly
			// marks the read pair as ready for processing.
			if read.is_reverse() {
				pair.left_pos = pair.r1.pos() as u32 + 1;
				pair.right_pos = read.cigar().end_pos() as u32;
			} else {
				pair.left_pos = read.pos() as u32 + 1;
				pair.right_pos = pair.r1.cigar().end_pos() as u32;
			}

			// Perform final confirmation that the read pair is amenable
			// for consensus generation.
			let frag_len = pair.left_pos.abs_diff(pair.right_pos) + 1;
			if read.pos() as u32 + 1 < pair.left_pos ||
				pair.r1.pos() as u32 + 1 < pair.left_pos ||
				read.cigar().end_pos() as u32 > pair.right_pos ||
				pair.r1.cigar().end_pos() as u32 > pair.right_pos {
				valid_for_consensus = false;
				stats.not_converging += 2;
			} else if frag_len > settings.max_frag_len {
				valid_for_consensus = false;
				stats.too_long += 2;
			} else if min(pair.r1.mapq() as u8, read.mapq() as u8) < settings.min_mapq {
				valid_for_consensus = false;
				stats.low_mapq += 2;
			}

			// If the read pair was not amenable for consensus generation,
			// we store it as an invalid ReadPair in the FIFO queue, and mark
			// the mate also as invalid in the FIFO queue.
			if valid_for_consensus == false {
				pair.mark_invalid();
				fifo.push_back(ReadPair::invalid(read));
	            continue;
			}

			// Calculate the physical strand of the sequenced ssDNA fragment
			pair.strand = if read.is_first_in_template() == read.is_reverse()
				{ b'-' } else { b'+' };

			assert!(pair.umi == umi_for_read(&read, settings.ignore_umi));
			pair.r2 = read;
			stats.concordant += 2;
			
		} else {
			// First half of a mate pair found. Save it and wait for its
			// paired mate to show up.
			mates.insert(qname.into(), n_pairs_completed + fifo.len());
			fifo.push_back(ReadPair {
				umi: umi_for_read(&read, settings.ignore_umi),
				r1: read, r2: Record::new(),
				left_pos: 0, right_pos: 0, strand: b'+'
			});
		}		

		n_pairs_completed += write_consensus(&mut out, &mut fifo, pos, &settings);
	}

	n_pairs_completed += write_consensus(&mut out, &mut fifo, u32::MAX, &settings);

	// Print statistics
	eprintln!("\nReads used for consensus generation: {} / {} ({:.1}%)",
		stats.concordant, stats.total_reads,
		(stats.concordant as f32 / stats.total_reads as f32 * 100.0));

	let bad = (stats.total_reads - stats.concordant) as f32;

	eprintln!("\nBreakdown of discordant reads:");
	if stats.unpaired > 0 { eprintln!("- Unpaired: {} ({:.1}%)", stats.unpaired,
		stats.unpaired as f32 / bad * 100.0); }
	eprintln!("- Low MAPQ: {} ({:.1}%)", stats.low_mapq,
		stats.low_mapq as f32 / bad * 100.0);
	eprintln!("- Non-converging: {} ({:.1}%)", stats.not_converging,
		stats.not_converging as f32 / bad * 100.0);
	eprintln!("- Unaligned: {} ({:.1}%)", stats.unmapped,
		stats.unmapped as f32 / bad * 100.0);
	eprintln!("- Too long: {} ({:.1}%)", stats.too_long,
		stats.too_long as f32 / bad * 100.0);
	eprintln!("- Interchromosomal: {} ({:.1}%)", stats.diff_chr,
		stats.diff_chr as f32 / bad * 100.0);
	eprintln!();
}


// This function identifies clusters of duplicate fragments by studying the
// FIFO queue of read pairs. Duplicate identification is only carried out for
// complete read pairs (i.e. both mates have been seen) for which we can prove
// that no further duplicates can be found in the BAM file.
fn write_consensus(out: &mut BamWriter, fifo: &mut VecDeque<ReadPair>, cur_pos: u32, settings: &Settings) -> usize {

	let mut n_pairs_completed: usize = 0;

	while !fifo.is_empty() && fifo[0].is_ready() && cur_pos > fifo[0].right_pos {

		let mut pair = fifo.pop_front().unwrap();
		n_pairs_completed += 1;

		// If the fragment was merged into consensus, do not write into BAM
		if pair.is_merged() { continue; }

		// Reads that were not amenable for consensus generation get marked
		// with a QC_FAIL flag and written out.
		if pair.is_invalid() {
			if settings.keep_discordant {
				pair.r1.set_quality_check_failed();
				out.write(&pair.r1);
			}
			continue;
		}

		// This read pair is valid for consensus generation. Find all of its
		// duplicates.
		let mut duplicates: Vec<usize> = Vec::new();
		for j in 0..fifo.len() {
			// TODO: Is this ok for earlier stop?
			//if fifo[j].r1.pos() + 1 > pair.left_pos { break; }

			if fifo[j].is_invalid() || fifo[j].is_merged() { continue; }
			if fifo[j].is_ready() == false { continue; }

			// We can stop looking if we see a ready fragment with a higher
			// left edge coordinate, since we know we won't see any
			// fragments with compatible boundaries after that.
			if fifo[j].left_pos > pair.left_pos { break; }

			// Not a duplicate if fragment boundaries don't match
			if pair.left_pos != fifo[j].left_pos { continue; }
			if pair.right_pos != fifo[j].right_pos { continue; }
			
			let umi_match = umi_diff(&pair.umi, &fifo[j].umi) <= 1;
			/*eprintln!("{} and {} distance: {}",
				str::from_utf8(&pairs[0].umi).unwrap(),
				str::from_utf8(&pairs[j].umi).unwrap(),
				umi_dst(&pairs[0].umi, &pairs[j].umi));*/

			// Reject the duplicate if the UMI does not match
			if umi_match == false && !settings.ignore_umi { continue; }

			// Looks like a valid duplicate
			duplicates.push(j);
		}

		let mut duplicate_refs: Vec<&ReadPair> = vec![&pair];
		for d in &duplicates { duplicate_refs.push(&fifo[*d]); }

		let consensus_record = build_consensus_for_duplicates(
			duplicate_refs, pair.left_pos, pair.right_pos, &settings);

		for d in &duplicates { fifo[*d].mark_merged(); }

		// WRITE OUTPUT
		out.write(&consensus_record);
	} 

	return n_pairs_completed;
}


fn build_consensus_for_duplicates(dups: Vec<&ReadPair>, from: u32, to: u32, settings: &Settings) -> Record {

	// Size of the reference genome region overlapped by the consensus fragment
	let ref_len: usize = (to - from + 1) as usize;

	
	// Add both mates for all duplicates into the pileup vector
	let mut pileups = vec![Altref::default(); ref_len];
	for dup in &dups {
		add_read_to_consensus(&dup.r1, dup.strand, from, &mut pileups);
		add_read_to_consensus(&dup.r2, dup.strand, from, &mut pileups);
	}

	// Calculate average mapping quality
	let mut mapq_sum = 0u32;
	for dup in &dups {
		mapq_sum += dup.r1.mapq() as u32 + dup.r2.mapq() as u32;
	}
	let avg_mapq = (mapq_sum as f32 / dups.len() as f32 / 2.0).round() as u8;

	let mut consensus_str: Vec<String> = Vec::new();
	for c in 0..ref_len {
		let cs = do_calc_consensus(&pileups[c]);
		consensus_str.push(cs);
	}

	// Create synthetic BAM record with consensus seq
	let new_name = dups[0].r1.qname().to_owned(); // Use name of first duplicate

	let mut cons_rec = Record::new();	
	cons_rec.set_flags(0u16); // Not paired, aligned, forward strand...
	cons_rec.set_tid(dups[0].r1.tid());
	cons_rec.set_pos(from as i64 - 1);
	cons_rec.set_mapq(avg_mapq);

	let new_cigar = do_compose_cigar(&consensus_str);
	let new_seq = do_compose_seq(&consensus_str);
	
	// BASEQ	
	let quals = do_compose_base_qualities(&pileups, &settings);
	// Without BASEQ: vec![ 255u8; new_seq.len()]	
	assert!(quals.len() == new_seq.len());
	
	// Note: Using set() invalidates AUX data
	cons_rec.set(new_name.as_ref(), Some(&new_cigar), new_seq.as_ref(), quals.as_ref());

	// Set mate info to non-existent
	cons_rec.set_mpos(-1i64);
	cons_rec.set_mtid(-1i32); 
	
	// TLEN, always positive for new (unpaired) segments	
	cons_rec.set_insert_size((dups[0].right_pos - dups[0].left_pos + 1) as i64);
	cons_rec.set_bin(reg2bin(from - 1, cons_rec.cigar().end_pos() as u32) as u16);

	// Add the number of read pairs included in the consensus
	cons_rec.push_aux(b"DP", &Aux::Integer(dups.len() as i64)); 
	if dups[0].umi.is_empty() == false {
		cons_rec.push_aux(b"RX", &Aux::String(&dups[0].umi))
	}

	if settings.print_alignment {
		print_aligned_seqs(&cons_rec, &dups, from, to, &new_cigar, &settings);
	}
	return cons_rec;
}



fn add_read_to_consensus(read: &Record, strand: u8, start_coord: u32, consensus: &mut Vec<Altref>) {

	let seq = read.seq().as_bytes();
	let cigar = read.cigar();
	let seqpos = read.pos() as u32 + 1;  // 1-based start position (inclusive)

	// Index bounds sanity check
	assert!(seqpos >= start_coord);
	assert!(read.cigar().end_pos() <= start_coord as i64 + consensus.len() as i64);

	let mut seq_idx: usize = 0;
	let mut ref_idx: usize = (seqpos - start_coord) as usize;		

	for s in cigar.iter() {
		match *s {
			// Cigar::Equal and Cigar::Diff are newer additions and not
			// supported by all aligners
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
				for _ in 0..len {
					let pileup = &mut consensus[ref_idx];

					// Add to count or insert new key
					*(pileup.alts).entry((seq[seq_idx] as char).to_string()).or_insert(0) += 1;
					pileup.n_total += 1;

					// Fill consensus base quality information
					if seq[seq_idx] != b'N' {
						if read.is_reverse() { pileup.rev += 1; }
						else { pileup.fwd += 1; }
	
						if strand == b'-' { pileup.umi_rev += 1; }
						else { pileup.umi_fwd += 1; }					
					}

					seq_idx += 1; 						
					ref_idx += 1;
				}
			},
			Cigar::Ins(len) => {
				if seq_idx == 0 { 
					error!("Insertion as first CIGAR element." ); 
				}

				// Create insertion string attached to previous ref-coord base
				let mut ins = String::with_capacity(len as usize + 3);
				ins.push(seq[seq_idx - 1] as char);
				for _ in 0..len {
					ins.push(seq[seq_idx] as char);
					seq_idx += 1;
				}
				// Only accept insertions without N bases
				if ins[1..].contains('N') == false {
					*(consensus[ref_idx-1].alts).entry(ins).or_insert(0) += 1;
				}
				// Insertions do not count towards total because they are
				// attached to previous base in seq.
			},
			Cigar::Del(len) => {
				for _ in 0..len {
					let pileup = &mut consensus[ref_idx];
					*(pileup.alts).entry('-'.to_string()).or_insert(0) += 1;
					pileup.n_total += 1;

					// Fill consensus base quality information
					if read.is_reverse() { consensus[ref_idx].rev += 1; }
					else { consensus[ref_idx].fwd += 1; }		

					if strand == b'-' {
						consensus[ref_idx].umi_rev += 1;
					} else {
						consensus[ref_idx].umi_fwd += 1;
					}	

					ref_idx += 1;
				}
			},
			Cigar::SoftClip(_len) | Cigar::HardClip(_len) => 
				error!("Unexpected hard/soft clip in CIGAR."),
			Cigar::RefSkip(_len) =>
				error!("Unexpected CIGAR type: N"),
			Cigar::Pad(_len) =>
				error!("Unexpected CIGAR type: P")
		}
	}
}



fn do_complete_cigar(cig_vec: &mut Vec<Cigar>, cig_type: char, count: u32) {	

	if count == 0 { return; }
	match cig_type {
		'M' => { cig_vec.push(Cigar::Match(count)); }
		'D' => { cig_vec.push(Cigar::Del(count)); }
		'I' => { cig_vec.push(Cigar::Ins(count)); }
		 _  => { error!("{}", format!("Unexpected CIGAR value '{}' encountered", cig_type)); }
	};
}



fn do_compose_cigar(cons_str : &Vec<String>) -> CigarString {

	let mut cig_vec : Vec<Cigar> = Vec::new();
	
	let mut cur_cig : char = '?';
	let mut prev_cig : char = '?';
	let mut prev_count : u32 = 0;

	for i in 0..cons_str.len() {

		let cur_str = &cons_str[ i];
		let len = cur_str.len();
		
		// Insertions:
		// Separate first base from insertion and complete previous cigar element
		if len > 1 {

			let attached_base_type = if cur_str.chars().nth(0).unwrap() != '-' { 'M' } else { 'D' };			

			if attached_base_type == prev_cig || prev_count == 0 {
				// Add to previous cigar elem
				do_complete_cigar( &mut cig_vec, attached_base_type, prev_count+1);
			} else {
				// Complete previous cigar elem separately
				do_complete_cigar( &mut cig_vec, prev_cig, prev_count);
				do_complete_cigar( &mut cig_vec, attached_base_type, 1);	
			}

			do_complete_cigar( &mut cig_vec, 'I', (len-1) as u32);
			prev_cig = '?';
			prev_count = 0;
			continue;
		}		

		assert!(len > 0);			

		if cur_str == " " { continue; }

		// Update current cigar
		cur_cig = if cur_str.chars().nth(0).unwrap() != '-' { 'M' } else { 'D' };
		let flush = prev_cig != cur_cig && prev_count > 0;

		if flush {
			do_complete_cigar( &mut cig_vec, prev_cig, prev_count);
			prev_count = 0;
		}

		prev_cig = cur_cig;
		prev_count += 1;
	}

	if prev_count > 0 {
		do_complete_cigar( &mut cig_vec, prev_cig, prev_count);
	}
	return CigarString(cig_vec);
}

fn do_compose_seq(cons_str: &Vec<String>) -> Vec<u8> {
	let mut retval: Vec<u8> = Vec::new();
	for s in cons_str {
		for c in s.chars() {
			if c == '-' { continue; } // Deletions do not go to record seq
			retval.push(c as u8);
		}
	}
	return retval;
}


fn do_compose_base_qualities(altrefs: &Vec<Altref>, settings: &Settings) -> Vec<u8> {

	// Five levels of confidence:
	// 1: one direction, no duplicates
	// 2: has duplicates, one direction
	// 3: both directions, no duplicates
	// 4: both directions, has duplicates
	// 5: both directions, has duplicates, flipped UMIs

	let mut retvals: Vec<u8> = Vec::new();
	for a in 0..altrefs.len() {
		let altref = &altrefs[a];
		let cons_str = do_calc_consensus(altref);
		let total = altref.fwd + altref.rev; // Excludes N's

		let both_dirs = altref.fwd > 0 && altref.rev > 0;

		let evidence_level = if total == 1 { 1 }
		else if altref.umi_fwd > 1 && altref.umi_rev > 1 && both_dirs { 5 }	
		else if (altref.fwd > 1 || altref.rev > 1) && !both_dirs { 2 }
		else if total == 2 && both_dirs { 3 }
		else if total > 2 && both_dirs { 4 }			
		else { 0 };

		for c in cons_str.bytes() {
			if c == b'-' { continue; }  // Deletions not part of SEQ or BASEQ
			retvals.push(if evidence_level >= settings.min_evidence { 42 } else { 0 });
		}
	}
	
	retvals
}



fn do_print_seq(seq: &Vec<u8>, cigar: &CigarString, slot_vec: &Vec<usize>, seqpos: u32, min_coord: u32) {

	let mut seq_idx: usize = 0;
	let mut ref_idx: usize = (seqpos - min_coord) as usize;

	// Left align new sequence accounting for ref slot sizes
	print!("{:1$}", "", slot_vec[0..ref_idx].iter().sum());

	let mut filler: usize = 0;

	for s in cigar.iter() { 	 
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
				if filler > 0 { print!("{:.1$}", "", filler); }
			    filler = 0;

				for k in 0..len {
					print!("{}", seq[seq_idx] as char);
					if k != len - 1 {
						print!("{:.1$}", "", slot_vec[ref_idx] - 1);
					} else {
						filler = slot_vec[ref_idx] - 1;
					}

					seq_idx += 1;
					ref_idx += 1;
				}
			},
			Cigar::Ins(len) => {
				for k in 0..len {
					print!("{}", (seq[seq_idx] as char).to_ascii_lowercase());
					seq_idx += 1;
				}
	
				if ref_idx > 0 { filler -= len as usize; }
				else { continue; }

				print!("{:.1$}", "", filler);
				filler = 0;
			},
			Cigar::Del(len) => {
				print!("{:.1$}", "", filler);
				filler = 0;
				
				for k in 0..len {	 
					print!("-");
					filler = slot_vec[ref_idx] - 1;
					if k != len - 1 {
						print!("{:.1$}", "", filler);
						filler = 0;
					}
					ref_idx += 1;
				}
			}
			_ => { continue; }
		};
	}	

	print!("{:.1$}\n", "", filler);
}

fn print_aligned_seqs(cons: &Record, dups: &Vec<&ReadPair>, min_coord: u32, max_coord: u32, ref_cigar: &CigarString, settings: &Settings) {

	let ref_len = (max_coord - min_coord + 1) as usize;

	// Calculate how much horizontal space is needed for the read IDs and UMIs
	let mut max_name_len: usize = "Reference".len();
	let mut max_umi_len: usize = 0;
	for dup in dups {
		let qname = str::from_utf8(dup.r1.qname()).unwrap();
		max_name_len = max(max_name_len, qname.len());
		max_umi_len = max(max_umi_len, dup.umi.len());
	}

	// Leave a one space gap between the UMI and sequence if the UMI is shown
	if max_umi_len > 0 { max_umi_len += 1; }

	// For each reference genome position, figure out what is the longest
	// insertion found there.
	// TODO: Surely we don't need the consensus sequence here?
	let mut slot_vec = vec![1usize; ref_len];
	for d in -1i32..(dups.len() as i32 * 2) {
		let record = if d < 0 { &cons }
			else if d % 2 == 0 { &dups[(d / 2) as usize].r1 }
			else { &dups[(d / 2) as usize].r2 };
		
		let seqpos = record.pos() as u32 + 1;
		assert!(seqpos >= min_coord);
		let mut ref_idx: usize = (seqpos - min_coord) as usize;

		for s in record.cigar().iter() { 
			match *s {
				Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
					ref_idx += len as usize;
				},
				Cigar::Ins(len) => {
					assert!(ref_idx != 0);   // Insertion at the beginning
					slot_vec[ref_idx - 1] = max((len + 1) as usize, slot_vec[ ref_idx - 1]);
				},
				Cigar::Del(len) => { ref_idx += len as usize; },
				_ => { error!("Unsupported CIGAR element."); }
			}
		}
	}

	// Print the fragment information and consensus sequence
	let qname = str::from_utf8(dups[0].r1.qname()).unwrap();
	println!("\n Fragment {} ({}:{}-{})", &qname,
		settings.chr_names[dups[0].r1.tid() as usize], min_coord, max_coord);
	print!(" {:<len$}", "Consensus", len=max_name_len);
	print!(" {:<len$}", "", len=max_umi_len);
	do_print_seq(&cons.seq().as_bytes(), &cons.cigar(), &slot_vec, cons.pos() as u32 + 1, min_coord);

	// Print the reads that were used to form the consensus fragment
	for d in 0..(dups.len() * 2) {
		let record = if d % 2 == 0 { &dups[(d / 2) as usize].r1 }
			else { &dups[(d / 2) as usize].r2 };

		let seq = record.seq().as_bytes();
		let cigar = record.cigar();
		let seqpos = record.pos() as u32 + 1;
		let qname = str::from_utf8(record.qname()).unwrap();

		let dir = if record.is_reverse() { '-' } else { '+' };
		let umi = str::from_utf8(&dups[(d / 2) as usize].umi).unwrap();
		print!("{}{:<len$}", dir, &qname, len=max_name_len);
		print!(" {:<len$}", umi, len=max_umi_len);
		
		do_print_seq(&seq, &cigar, &slot_vec, seqpos, min_coord);
	}
	print!("\n");
}

// Calculates consensus for a single reference genome coordinate
fn do_calc_consensus(altref: &Altref) -> String {

	let mut retval: String = 'N'.to_string();
	if altref.n_total == 0 { return retval; }

	let threshold: u32 = ((altref.n_total as f32) * 0.8).ceil() as u32; 
	let mut longest: usize = 0;
	
	// Longest threshold reaching alt is preferred to give
	// insertions a chance to win over single base substitions
    for (key, value) in &altref.alts {
        if *value >= threshold && key.len() > longest {
        	retval = key.to_string();
        	longest = key.len();
        }
    }

    // If no base amongst the duplicates reaches threshold, return 'N'
    return retval;
}

fn umi_for_read(read: &Record, ignore_umi: bool) -> Box<[u8]> {
	let mut umi: Vec<u8> = Vec::new();
	if !ignore_umi {
		if let Some(Aux::String(rx)) = read.aux(b"RX") {
			umi.extend_from_slice(rx);
		}
	}
	return umi.into();
}

fn umi_diff(a: &[u8], b: &[u8]) -> u8 {

	if a.is_empty() || b.is_empty() { return 0; }
	if a.len() != b.len() { return 255; }
	
	let a_flippable = a.iter().position(|&x| x == b'+');
	let b_flippable = b.iter().position(|&x| x == b'+');
	assert!(a_flippable == b_flippable);  // TODO: Add support for asymmetry?
	//if a_flippable != b_flippable { return 255; }

	// Single-stranded UMI
	if a_flippable == None {
		let mut mismatches = 0u8;
		for k in 0..a.len() {
			if !(a[k] == b[k] || a[k] == b'N' || b[k] == b'N') {
				mismatches += 1;
			}
		}
		return mismatches;
	}

	// Duplex UMI
	let sep_pos = a_flippable.unwrap();
	let mut fwd_miss = 0u8;
	let mut flip_miss = 0u8;

	for f in 0..sep_pos {
		let r = sep_pos+f+1;

		if !(a[f] == b[f] || a[f] == b'N' || b[f] == b'N') { fwd_miss += 1; }	
		if !(a[r] == b[r] || a[r] == b'N' || b[r] == b'N') { fwd_miss += 1; }

		if !(a[f] == b[r] || a[f] == b'N' || b[r] == b'N') { flip_miss += 1; }
		if !(a[r] == b[f] || a[r] == b'N' || b[f] == b'N') { flip_miss += 1; }
	}		

	return min(fwd_miss, flip_miss);
}


// calculate BAM bin-field given an fragment [beg,end)
// This function calculates the bin number for a BAM record (used in
// BAM indexing). The calculation is defined in SAM specification section 5.3.
// The beginning and end coordinates must be given as a [beg, end) zero-based
// semi-open interval. If you need to provide a bin number for a BAM record
// with POSITION = -1 (i.e. an unmapped read), the bin should be set to 4680
// according to the specification. 
pub fn reg2bin(beg: u32, end: u32) -> u32 {
	let end = end - 1;
	if beg >> 14 == end >> 14 { return ((1 << 15) - 1) / 7 + (beg >> 14); }
	if beg >> 17 == end >> 17 { return ((1 << 12) - 1) / 7 + (beg >> 17); }
	if beg >> 20 == end >> 20 { return ((1 << 9) - 1) / 7 + (beg >> 20); }
	if beg >> 23 == end >> 23 { return ((1 << 6) - 1) / 7 + (beg >> 23); }
	if beg >> 26 == end >> 26 { return ((1 << 3) - 1) / 7 + (beg >> 26); }
	return 0;
}

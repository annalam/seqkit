
use crate::common::{parse_args, BamReader, BamWriter};
use std::{str, cmp::Ordering};
use std::collections::BinaryHeap;
use rust_htslib::bam;
use rust_htslib::bam::{header::Header, Format, CompressionLevel};

const USAGE: &str = "
Usage:
  sam merge [options] <bam_files>...

Options:
  --suffix          Add a suffix to read identifiers to avoid clashes
  --uncompressed    Output in uncompressed BAM format

Merges two or more position-sorted BAM files together, ensuring that the
resulting output BAM file is also position-sorted.
";

// A wrapper type for position sorted BAM records, implementing Ord.
// We want to use this type inside a BinaryHeap to implement a min-heap,
// so we implement the Ord in reverse (BinaryHeap is a max-heap by default).
struct SortedRecord { record: bam::Record, bam_idx: usize }
impl PartialEq for SortedRecord {
	fn eq(&self, other: &Self) -> bool {
		self.record.tid() == other.record.tid() &&
			self.record.pos() == other.record.pos()
	}
}
impl Eq for SortedRecord {}
impl Ord for SortedRecord {
	fn cmp(&self, other: &Self) -> Ordering {
		// We convert the TIDs to u32 to make -1 overflow to u32::MAX, since
		// in sorted BAM files unmapped reads have TID = -1 and come last.
		if (self.record.tid() as u32) < (other.record.tid() as u32) {
			Ordering::Greater
		} else if self.record.tid() as u32 > other.record.tid() as u32 {
			Ordering::Less
		} else if self.record.pos() < other.record.pos() {
			Ordering::Greater
		} else if self.record.pos() > other.record.pos() {
			Ordering::Less
		} else if self.record.tid() == other.record.tid() &&
			self.record.pos() == other.record.pos() {
			Ordering::Equal
		} else {
			unreachable!()
		}
	}
}
impl PartialOrd for SortedRecord {
	fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		Some(self.cmp(other))
	}
}


pub fn main() {
	let args = parse_args(USAGE);
	let bam_paths = args.get_vec("<bam_files>");
	let add_suffix = args.get_bool("--suffix");
	if bam_paths.len() < 2 {
		error!("At least two BAM files must be provided for concatenation.");
	}

	// Open all of the BAM files for reading and check that they contain
	// identical SQ fields.
	let mut bams: Vec<_> = bam_paths.iter().map(|bam_path| BamReader::open(bam_path)).collect();
	let header = bams[0].header();
	let chr_names = header.target_names();
	for b in 1..bams.len() {
		if bams[b].header().target_names() != chr_names {
			error!("Input BAM files {} and {} have different SQ fields.",
				&bam_paths[0], &bam_paths[b]);
		}
	}

	let mut out = BamWriter::open("-", &header, !args.get_bool("--uncompressed"));

	// Initialize the position-sorted binary heap that we will use for merging
	// BAM records while maintaining the position-sorted property.
	let mut heap: BinaryHeap<SortedRecord> = BinaryHeap::new();
	for b in 0..bams.len() {
		if let Some(record) = bams[b].next() {
			heap.push(SortedRecord { record, bam_idx: b });
		}
	}

	loop {
		if let Some(mut first) = heap.pop() {
			if let Some(record) = bams[first.bam_idx].next() {
				heap.push(SortedRecord { record, bam_idx: first.bam_idx });
			}
			if add_suffix {
				let suffix = format!(".{}", first.bam_idx + 1).into_bytes();
				let mut qname = first.record.qname().to_vec();
				qname.extend(&suffix);
				first.record.set_qname(&qname);
			}
			out.write(&first.record);
		} else { break; }
	}
}

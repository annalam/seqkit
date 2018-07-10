
use common::{parse_args, PathArgs};
use std::str;
use std::collections::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rand::random;

const USAGE: &str = "
Usage:
  sam subsample <bam_file> <fraction>

If your BAM file has been duplicate-flagged, remember to re-run duplicate
flagging after subsampling, otherwise random subsampling can delete the only
non-duplicate-flagged DNA fragment in a duplicate cluster.
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_path("<bam_file>");
	let keep_frac: f32 = args.get_str("<fraction>").parse().unwrap_or(-1.0);
	if !(keep_frac >= 0.0 && keep_frac <= 1.0) {
		error!("Subsampling fraction must be between 0 - 1.");
	}

	let mut bam = bam::Reader::from_path(&bam_path)
		.unwrap_or_else(|_| error!("Could not open BAM file {}.", &bam_path));
	let bam_header = bam.header().clone();

	let mut out = bam::Writer::from_stdout(
		&bam::header::Header::from_template(&bam_header)).unwrap();

	let mut total_reads: u64 = 0;
	let mut kept_reads: u64 = 0;

	let mut keep_mate: HashMap<Vec<u8>, bool> = HashMap::new();

	for r in bam.records() {
		let read = r.unwrap();
		if read.is_supplementary() { continue; }  // TODO: Warn user about this
		
		if read.is_paired() {
			let qname = read.qname();
			let mut keep = false;
			if let Some(x) = keep_mate.remove(qname) {
				keep = x;
			} else {
				// We haven't seen this DNA fragment yet, so we roll the dice
				// to decide if we keep the fragment or not.
				keep = random::<f32>() <= keep_frac;
				keep_mate.insert(qname.into(), keep);
			}

			if keep {
				out.write(&read).unwrap_or_else(
					|_| error!("Output stream closed unexpectedly."));
				kept_reads += 1;
			}
			total_reads += 1;
		} else {
			error!("Only paired end sequencing data supported for now.");
		}
    }

    eprintln!("Total reads: {}", total_reads);
    eprintln!("Kept reads: {} ({:.1}% of all reads)", kept_reads, kept_reads as f64 / total_reads as f64 * 100.0);
}

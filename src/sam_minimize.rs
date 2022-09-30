
use crate::common::{parse_args, BamReader, BamWriter};
use std::str;
use std::collections::HashMap;

const USAGE: &str = "
Usage:
  sam minimize [options] <bam_file>

Options:
  --uncompressed    Output in uncompressed BAM format
  --read-ids        Minimize read identifiers (i.e. QNAME fields)
  --base-qualities  Remove per-base qualities
  --tags            Remove all aux fields (tags)
  --baseq-fill=N    Base quality value to fill in as placeholder [default: 255]

Changes read IDs into simple numeric identifiers, removes per-base qualities,
and removes all auxiliary fields (tags).
";

pub fn main() {
	let args = parse_args(USAGE);
	let bam_path = args.get_str("<bam_file>");
	let minimize_qnames = args.get_bool("--read-ids");
	let remove_baseq = args.get_bool("--base-qualities");
	let remove_tags = args.get_bool("--tags");
	let baseq_fill: u8 = args.get_str("--baseq-fill").parse().unwrap_or_else(
		|_| error!("--baseq-fill must be an integer between 0 and 255."));

	if !minimize_qnames && !remove_baseq && !remove_tags {
		error!("One of --read-ids, --base-qualities, or --tags must be given.");
	}

	if remove_baseq && !remove_tags {
		error!("Running 'sam minimize' with --base-qualities but without the --tags flag is not yet supported.");
	}

	// These variables are used for QNAME simplification
	let mut highest_id: u32 = 0;
	let mut qname_to_id: HashMap<Vec<u8>, u32> = HashMap::new();

	let bam = BamReader::open(&bam_path);
	let mut out = BamWriter::open("-", &bam.header(), !args.get_bool("--uncompressed"));

	for mut read in bam {
		let mut qname = read.qname().to_vec();

		if minimize_qnames {
			if let Some(slash_pos) = qname.iter().position(|x| *x == b'/') {
				qname.truncate(slash_pos);
			}

			let id = if let Some(x) = qname_to_id.remove(&qname) { x } else {
				highest_id += 1;
				qname_to_id.insert(qname.clone(), highest_id);
				highest_id
			};
			qname = format!("{}", id).into_bytes();
		}

		let cigar = read.cigar();
		let seq = read.seq().as_bytes();

		let qual = if remove_baseq {
			// If the user requested to remove base quality information, we
			// replace the BASEQ field with enough placeholder bytes to match
			// the sequence length. By default, the placeholder byte 255 (0xFF)
			// is used, which according to the BAM specification indicates a
			// missing quality value. The user can override the default
			// placeholder byte using the --baseq-fill argument.
			vec![baseq_fill; seq.len()]
		} else {
			read.qual().to_vec()
		};

		if minimize_qnames && !remove_baseq && !remove_tags {
			read.set_qname(&qname);
		} else {
			// The call to .set() removes all AUX fields
			read.set(&qname, Some(&cigar), &seq, &qual);
		}
		out.write(&read);
	}
}

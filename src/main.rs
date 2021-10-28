use std::{cell::RefCell, fs::File, io::BufWriter};

use debruijn::dna_string::DnaString;
use fastq_dedup_ec::DeduplicationLookup;
use seq_io::fastq::Reader;
use clap::{Parser};

#[derive(Parser)]
#[clap(version = "0.1", author = "Daniel Cameron")]
struct Opts {
    #[clap(long)]
    in1: String,
    #[clap(long)]
    in2: String,
    #[clap(long)]
    out1: String,
    #[clap(long)]
    out2: String,
    #[clap(short, long)]
    mapping: String,
    #[clap(long, default_value="8")]
    max_hamming_distance: u32,
    #[clap(long, default_value="150")]
    max_phred_distance: u32,
    #[clap(long, default_value="2.0")]
    max_per_base_phred_distance: f32,
    #[clap(long, default_value="1024")]
    uninformative_kmer_threshold: usize,
    // Maximum edit distance between two different consensuses for them to be merged
    // Consensuses will only be merged into consensuses that have greater read support
    //#[clap(long, default_value="8")]
    //consensus_merge_max_hamming_distance: usize,
    /// File to output all distance results to.
    /// Warning: this file is extremely verbose.
    #[clap(long)]
    debug_distances: Option<String>,
}

fn go(args : Opts) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut lookup = DeduplicationLookup::new();
    lookup.max_hamming_distance = args.max_hamming_distance;
    lookup.max_phred_distance = args.max_phred_distance;
    lookup.max_per_base_phred_distance = args.max_per_base_phred_distance;
    lookup.uninformative_kmer_threshold = args.uninformative_kmer_threshold;
    if let Some(filename) = args.debug_distances {
        lookup.debug_distances = Some(RefCell::new(Box::new(BufWriter::new(File::create(&filename)?))));
    }
    let mut mapping_writer = BufWriter::new(File::create(&args.mapping)?);
    let mut reader1 = Reader::from_path(args.in1)?;
    let mut reader2 = Reader::from_path(args.in2)?;
    for (r1, r2) in reader1.records().zip(reader2.records()) {
        let r1 = r1?;
        let r2 = r2?;
        let seq1 = DnaString::from_acgt_bytes(&r1.seq);
        let seq2 = DnaString::from_acgt_bytes(&r2.seq);
        let qual1 = &r1.qual;
        let qual2 = &r2.qual;
        let id = lookup.add_read_pair(seq1, qual1, seq2, qual2);
        write_mapping(&mut mapping_writer, id, &r1.head)?;
    }
    //lookup.merge_consenses(args.consensus_merge_max_hamming_distance);
    // Output consensus
    let mut writer1 = BufWriter::new(File::create(&args.out1)?);
    let mut writer2 = BufWriter::new(File::create(&args.out2)?);
    for (i, cons) in lookup.seq.iter().enumerate() {
        if cons.read_count() > 0 {
            let name = format!("cons_{}_r_{}_", i, cons.read_count());
            let seq = cons.sequence();
            let qual = cons.paired_fastq_qual_score();
            seq_io::fastq::write_to(&mut writer1, name.as_bytes(), &seq.0.to_ascii_vec(), &qual.0)?;
            seq_io::fastq::write_to(&mut writer2, name.as_bytes(), &seq.1.to_ascii_vec(), &qual.1)?;
        }
    }
    Ok(())
}
fn write_mapping<W: std::io::Write>(
    writer: &mut W,
    id: u32, // ConsensusSequenceId
    read_name: &[u8]
) -> std::io::Result<()> {
    writer.write_all(id.to_string().as_bytes())?;
    writer.write_all(b",")?;
    writer.write_all(read_name)?;
    writer.write_all(b"\n")?;
    Ok(())
}
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
   go(Opts::parse())
}
#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn test_command_line_small_consensus_distances() {
        let args = Opts{
            in1: "SRR13320978_1.fastq".to_string(),
            in2: "SRR13320978_2.fastq".to_string(),
            out1: "out.SRR13320978_1.fq".to_string(),
            out2: "out.SRR13320978_2.fq".to_string(),
            mapping: "out.SRR13320978.mapping.csv".to_string(),
            max_hamming_distance: 8,
            max_phred_distance: 150,
            max_per_base_phred_distance: 2.0,
            uninformative_kmer_threshold: 1024,
            debug_distances: Some("out.debug_distances.tsv".to_string()),
        };
        go(args).unwrap();
        let args = Opts{
            in1: "out.SRR13320978_1.fq".to_string(),
            in2: "out.SRR13320978_2.fq".to_string(),
            out1: "out.SRR13320978_1.cons.fq".to_string(),
            out2: "out.SRR13320978_2.cons.fq".to_string(),
            mapping: "out.SRR13320978.cons.mapping.csv".to_string(),
            max_hamming_distance: 1000,
            max_phred_distance: 100000,
            max_per_base_phred_distance: 0.0,
            uninformative_kmer_threshold: 50000,
            debug_distances: Some("out.debug_distances.cons.tsv".to_string()),
        };
        go(args).unwrap();
    }
    #[test]
    fn test_command_line_large() {
        let args = Opts{
            in1: "SRR13320955_1.fastq".to_string(),
            in2: "SRR13320955_2.fastq".to_string(),
            out1: "out.SRR13320955_1.fq".to_string(),
            out2: "out.SRR13320955_2.fq".to_string(),
            mapping: "out.SRR13320955.mapping.csv".to_string(),
            max_hamming_distance: 8,
            max_phred_distance: 150,
            max_per_base_phred_distance: 2.0,
            uninformative_kmer_threshold: 1024,
            debug_distances: None,
        };
        go(args).unwrap();
    }
}
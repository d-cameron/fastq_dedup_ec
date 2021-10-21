use fastq_dedup_ec::DeduplicationLookup;

fn main() {
    let lookup = DeduplicationLookup::new();
    // iterator over read pairs
    {
        let id = lookup.add_read_pair(r1, qual1, r2, qual2);
        // write read -> consensus mapping to file
    }
    // Output consensus sequences:
    for cons in &lookup.seq {
        let (r1, r2) = cons.sequence();
        let (q1, q2) = cons.paired_fastq_qual_score();
        let reads = cons.read_count();
        // write r1/r2
    }
}

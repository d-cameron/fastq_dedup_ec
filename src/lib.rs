use std::collections::HashMap;

use debruijn::{Mer, Vmer, dna_string::DnaString, kmer::Kmer16};

// Enhancements:
// - Ignore lookups with many kmer
// - Collapse similar consensus sequences

type BaseQualAccumulation = u32;
pub struct ConsensusSequence {
    seq: [DnaString; 2],
    basequals: [Vec<[BaseQualAccumulation;4]>; 2],
    reads: u32,
}
// Maximum number of unique sequences
type ConsensusSequenceId = u32;
type Kmer = Kmer16; // UPDATE NEXT LINE IF YOU CHANGE THIS
const KMER_SIZE: usize = 16; // TODO: how do we access Kmer::_k()? DnaString::from_acgt_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").get_kmer::<Kmer>(0).len();
pub struct DeduplicationLookup {
    pub seq: Vec<ConsensusSequence>,
    lookup: [Vec<HashMap<Kmer, ConsensusSequenceId>>; 2],
}
const MAX_EDIT_DISTANCE: usize = 8;
const MAX_PHRED_DISTANCE: usize = MAX_EDIT_DISTANCE * 42;
impl DeduplicationLookup {
    pub fn new() -> DeduplicationLookup {
        DeduplicationLookup {
            seq: vec![],
            lookup: [vec![], vec![]],
        }
    }
    pub fn add_read_pair(&mut self, r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> ConsensusSequenceId {
        // find best candidate
        match self.find_closest_consensus(r1, qual1, r2, qual2) {
            Some((bestid, best_edit_distance, best_phred_distance)) => {
                if self.seq[bestid as usize].add_read_pair(r1, qual1, r2, qual2) {
                    self.update_lookup(bestid);
                }
                bestid
            }
            None => {
                self.seq.push(ConsensusSequence::new(r1, qual1, r2, qual2));
            let newid = (self.seq.len() - 1) as ConsensusSequenceId;
            self.update_lookup(newid);
            newid
            }
        }
    }
    pub fn find_closest_consensus(&self, r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> Option<(ConsensusSequenceId, usize, usize)> {
        let mut kmer_matches = HashMap::new();
        // matching kmers
        // ===================================
        // ===================================
        //   TODO: count matching kmers per consens
        // ===================================
        // ===================================

        // ===================================
        // Find the one with the lowest edit/phred distance
        // ===================================
    }
    /// Update the lookup for this consensus
    /// Implementation note: we don't know what's already in the lookup
    /// I've chosen not to remove the 'old' consensus kmers as it's simpler and doesn't really change the outcome
    /// The consensus shouldn't change often and more error kmers actually gives a slightly better chance
    /// of the kmer lookup matching
    fn update_lookup(&mut self, id: ConsensusSequenceId) {
        for rindex in [0, 1] {
            DeduplicationLookup::update_lookup_read(&mut self.lookup[rindex], &self.seq[id as usize].seq[rindex], id);
        }
    }
    fn update_lookup_read(lookup: &mut Vec<HashMap<Kmer16, ConsensusSequenceId>>, seq: &DnaString, id: ConsensusSequenceId) {
        while lookup.len() * KMER_SIZE < seq.len() {
            lookup.push(HashMap::new());
        }
        let mut i = 0;
        while (i * KMER_SIZE) <= seq.len() {
            let read_offset = i * KMER_SIZE;
            let kmer = seq.get_kmer(read_offset);
            lookup[i].insert(kmer, id);
            i += 1;
        }
    }
}
const FASTQ_PHRED_ENCODING: u8 = 33;
impl ConsensusSequence {
    fn fastq_qual_score(seq: &DnaString, basequals: &Vec<[BaseQualAccumulation;4]>) -> Vec<u8> {
        let mut out = Vec::with_capacity(seq.len());
        for i in 0..out.len() {
            let encoded_base = seq.get(i);
            let score= 2 * basequals[i][encoded_base as usize] - basequals[i].iter().sum::<BaseQualAccumulation>();
            // 0-93 allowed: https://en.wikipedia.org/wiki/FASTQ_format
            out.push(FASTQ_PHRED_ENCODING + std::cmp::max(0, std::cmp::min(93, score)) as u8);
        }
        out
    }
    pub fn sequence(&self) -> (&DnaString, &DnaString) {
        (&self.seq[0], &self.seq[1])
    }
    pub fn read_count(&self) -> u32 { 
        self.reads
    }
    pub fn paired_fastq_qual_score(&self) -> (Vec<u8>, Vec<u8>) {
        (
            ConsensusSequence::fastq_qual_score(&self.seq[0], &self.basequals[0]),
            ConsensusSequence::fastq_qual_score(&self.seq[1], &self.basequals[1])
        )
    }
    pub fn new(r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> ConsensusSequence {
        let quals = [ConsensusSequence::basequals(&r1, qual1), ConsensusSequence::basequals(&r2, qual2)];
        ConsensusSequence {
            seq: [r1, r2],
            basequals: quals,
            reads: 1,
        }
    }
    fn basequals(seq: &DnaString, qual : &[u8]) -> Vec<[BaseQualAccumulation;4]> {
        let mut out = vec![[0 as BaseQualAccumulation; 4]; qual.len()];
        for i in 0..out.len() {
            out[i][seq.get(i) as usize] += (qual[i] - FASTQ_PHRED_ENCODING) as BaseQualAccumulation;
        }
        out
    }
    pub fn add_read_pair(&mut self, r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> bool {
        self.reads += 1;
        // not using || shortcut since we always want to add both reads
        self.add_read(0, r1, qual1) | self.add_read(1, r2, qual2)
    }
    pub fn add_read(&mut self, rindex: usize, r: DnaString, qual: &[u8]) -> bool {
        let seq = &mut self.seq[rindex];
        let q = &mut self.basequals[rindex];
        let mut has_changed = false;
        if r.len() > seq.len() {
            seq.extend(r.slice(seq.len(), r.len()).iter());
            for _ in seq.len()..r.len() {
                q.push([0; 4]);
            }
            assert_eq!(seq.len(), q.len());
            has_changed = true;
        }
        for i in 0..r.len() {
            let existing_base = seq.get(i);
            let existing_base_score = q[i][existing_base as usize];
            let read_base = r.get(i);
            q[i][read_base as usize] += (qual[i] - FASTQ_PHRED_ENCODING) as BaseQualAccumulation;
            let read_base_score = q[i][read_base as usize];
            if read_base != existing_base && existing_base_score < read_base_score {
                seq.set_mut(i, read_base);
                has_changed = true;
            }
        }
        has_changed
    }
}
/*
// Decided I didn't like the needletail API (e.g. newline termination) and went with seq_io
pub struct NeedleTailFastqIterator {
    reader1: Box<dyn FastxReader>,
    reader2: Box<dyn FastxReader>,
}
impl NeedleTailFastqIterator {
    pub fn new_paired(file1 : &str, file2 : &str) -> Result<NeedleTailFastqIterator, ParseError> {
        let reader1 = parse_fastx_file(file1)?;
        let reader2 = parse_fastx_file(file2)?;
        Ok(NeedleTailFastqIterator { reader1, reader2 })
    }
    pub fn next(&mut self) -> Option<Result<(SequenceRecord, SequenceRecord), ParseError>> {
        let rec1 = self.reader1.next();
        let rec2 = self.reader2.next();
        if rec1.is_none() && rec2.is_none() {
            return None;
        }
        let seq1 = match rec1 {
            None => {
                let err_pos = ErrorPosition { id : None, line : self.reader1.position().line() };
                let err = needletail::errors::ParseError::new_unexpected_end(err_pos, needletail::parser::Format::Fastq);
                return Some(Err(err));
            },
            Some(Ok(x)) => x,
            Some(Err(x)) => return Some(Err(x)),
        };
        let seq2 = match rec2 {
            None => {
                let err_pos = ErrorPosition { id : None, line : self.reader2.position().line() };
                let err = needletail::errors::ParseError::new_unexpected_end(err_pos, needletail::parser::Format::Fastq);
                return Some(Err(err));
            },
            Some(Ok(x)) => x,
            Some(Err(x)) => return Some(Err(x)),
        };
        Some(Ok((seq1, seq2)))
    }
}
*/
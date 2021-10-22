use std::{cmp::{self, Ordering}, collections::HashMap};

use multimap::MultiMap;

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
    lookup: [Vec<MultiMap<Kmer, ConsensusSequenceId>>; 2],
    max_hamming_distance: HammingDistance,
    max_phred_distance: PhredDistance,
    max_per_base_phred_distance: PerBasePhredDistance,
}
type HammingDistance = u32;
type PhredDistance = u32;
type PerBasePhredDistance = f32;
impl DeduplicationLookup {
    pub fn new() -> DeduplicationLookup {
        DeduplicationLookup {
            seq: vec![],
            lookup: [vec![], vec![]],
            max_hamming_distance: 8,
            max_phred_distance: 6 * 42,
            max_per_base_phred_distance: 2.0,
        }
    }
    pub fn add_read_pair(&mut self, r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> ConsensusSequenceId {
        // find best candidate
        match self.find_closest_consensus(&r1, qual1, &r2, qual2) {
            Some((bestid, _, _, _)) => {
                if self.seq[bestid as usize].add_read_pair(&r1, qual1, &r2, qual2) {
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
    pub fn find_closest_consensus(&self, r1: &DnaString, qual1: &[u8], r2: &DnaString, qual2: &[u8]) -> Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)> {
        let mut kmer_matches = HashMap::new();
        DeduplicationLookup::add_matching_kmers(&self.lookup[0], r1, &mut kmer_matches);
        DeduplicationLookup::add_matching_kmers(&self.lookup[1], r2, &mut kmer_matches);
        let mut best_match = None;
        for id in kmer_matches.keys() {
            // TODO: add fast exits
            let edit_distance = DeduplicationLookup::calc_edit_distances(*id, &self.seq[*id as usize], &r1, qual1, &r2, qual2);
            if edit_distance.1 < self.max_hamming_distance
                    && edit_distance.2 < self.max_phred_distance
                    && edit_distance.3 < self.max_per_base_phred_distance
                    && DeduplicationLookup::best_id_cmp(best_match, Some(edit_distance)) == Ordering::Less {
                best_match = Some(edit_distance);
            }
        }
        best_match
    }
    fn calc_edit_distances(id: ConsensusSequenceId, consensus: &ConsensusSequence, r1: &DnaString, qual1: &[u8], r2: &DnaString, qual2: &[u8]) -> (ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance) {
        let (edit_distance1, phred_distance1, shared_bases1) = DeduplicationLookup::calc_edit_distance(&consensus.seq[0], r1, qual1);
        let (edit_distance2, phred_distance2, shared_bases2) = DeduplicationLookup::calc_edit_distance(&consensus.seq[1], r2, qual2);
        (id, edit_distance1 + edit_distance2, phred_distance1 + phred_distance2, ((phred_distance1 + phred_distance2) as PerBasePhredDistance) / ((shared_bases1 + shared_bases2) as PerBasePhredDistance))
    }
    fn calc_edit_distance(seq: &DnaString, r: &DnaString, qual: &[u8]) -> (HammingDistance, PhredDistance, usize) {
        let shared_bases = cmp::min(seq.len(), r.len());
        let hamming = seq.slice(0, shared_bases).hamming_dist(&r.slice(0, shared_bases));
        let mut phred_distance = 0;
        for i in 0..shared_bases {
            if seq.get(i) != r.get(i) {
                phred_distance += (qual[i] - FASTQ_PHRED_ENCODING) as PhredDistance;
            }
        }
        (hamming, phred_distance, shared_bases)
    }
    fn best_id_cmp(a: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>, b: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>) -> Ordering {
        match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(y)) => Ordering::Less,
            (Some(y), None) => Ordering::Greater,
            (Some((_, _, _, a_per_base_phred_errors)), Some((_, _, _, b_per_base_phred_errors))) =>
            b_per_base_phred_errors.partial_cmp(&a_per_base_phred_errors).unwrap()
        }
    }
    fn add_matching_kmers(lookup: &Vec<MultiMap<Kmer, ConsensusSequenceId>>, seq: &DnaString, counts: &mut HashMap<ConsensusSequenceId, u16>) {
        let mut i = 0;
        static EMPTY_LOOKUP : Vec<ConsensusSequenceId> = Vec::new();
        while (i * KMER_SIZE) <= seq.len() {
            let read_offset = i * KMER_SIZE;
            let kmer: Kmer = seq.get_kmer(read_offset);
            let consensuses_containing_kmer_at_position = lookup[i].get_vec(&kmer).unwrap_or(&EMPTY_LOOKUP);
            // TODO: exclude consensus with too many matching kmers
            for id in consensuses_containing_kmer_at_position {
                if let Some(x) = counts.get_mut(id) {
                    *x += 1;
                } else {
                    counts.insert(*id, 1);
                }
            }
            i += 1;
        }
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
    fn update_lookup_read(lookup: &mut Vec<MultiMap<Kmer, ConsensusSequenceId>>, seq: &DnaString, id: ConsensusSequenceId) {
        while lookup.len() * KMER_SIZE < seq.len() {
            lookup.push(MultiMap::new());
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
    pub fn add_read_pair(&mut self, r1: &DnaString, qual1: &[u8], r2: &DnaString, qual2: &[u8]) -> bool {
        self.reads += 1;
        // not using || shortcut since we always want to add both reads
        self.add_read(0, r1, qual1) | self.add_read(1, r2, qual2)
    }
    pub fn add_read(&mut self, rindex: usize, r: &DnaString, qual: &[u8]) -> bool {
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
#[cfg(test)]
mod tests {
    use crate::*;
    fn s(seq: &str) -> DnaString {
        DnaString::from_acgt_bytes(seq.as_bytes())
    }
    #[test]
    fn test_edit_distances() {
        assert_eq!((1, 1, 4), DeduplicationLookup::calc_edit_distance(&s("ACGT"), &s("ACGA"), &[33,33,33,34]));
        assert_eq!((1, 2, 4), DeduplicationLookup::calc_edit_distance(&s("ACGT"), &s("ACGA"), &[33,33,33,35]));
        assert_eq!((2, 5, 4), DeduplicationLookup::calc_edit_distance(&s("ACGT"), &s("ACTA"), &[33,33,34,37]));
        assert_eq!((1, 7, 3), DeduplicationLookup::calc_edit_distance(&s("ACG"), &s("ACTA"), &[33,33,40,37]));
    }
}
use std::{cmp::{self, Ordering}, collections::{hash_map::DefaultHasher}, hash::{Hash, Hasher}};
use multimap::MultiMap;
use rustc_hash::FxHashMap;
use debruijn::{Mer, Vmer, dna_string::DnaString, kmer::{Kmer12}};

// Enhancements:
// - [done] Ignore lookups with many kmer
// - [done] fast exact matching
// - Collapse similar consensus sequences

type BaseQualAccumulation = u32;
pub struct ConsensusSequence {
    seq: [DnaString; 2],
    basequals: [Vec<[BaseQualAccumulation;4]>; 2],
    reads: u32,
}
// Maximum number of unique sequences
type ConsensusSequenceId = u32;
type Kmer = Kmer12; // UPDATE NEXT LINE IF YOU CHANGE THIS
const KMER_SIZE: usize = 12; // TODO: how do we access Kmer::_k()? DnaString::from_acgt_bytes(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").get_kmer::<Kmer>(0).len();
const KMER_STRIDE: usize = KMER_SIZE;
pub struct DeduplicationLookup {
    pub seq: Vec<ConsensusSequence>,
    exact_lookup: FxHashMap<u64, ConsensusSequenceId>,
    lookup: [Vec<MultiMap<Kmer, ConsensusSequenceId>>; 2],
    pub max_hamming_distance: HammingDistance,
    pub max_phred_distance: PhredDistance,
    pub max_per_base_phred_distance: PerBasePhredDistance,
    pub uninformative_kmer_threshold: usize,
}
type HammingDistance = u32;
type PhredDistance = u32;
type PerBasePhredDistance = f32;
impl DeduplicationLookup {
    pub fn new() -> DeduplicationLookup {
        DeduplicationLookup {
            seq: vec![],
            exact_lookup: FxHashMap::default(),
            lookup: [vec![], vec![]],
            max_hamming_distance: 8,
            max_phred_distance: 6 * 42,
            max_per_base_phred_distance: 2.0,
            uninformative_kmer_threshold: 512,
        }
    }
    pub fn add_read_pair(&mut self, r1: DnaString, qual1: &[u8], r2: DnaString, qual2: &[u8]) -> ConsensusSequenceId {
        assert_eq!(r1.len(), qual1.len());
        assert_eq!(r2.len(), qual2.len());
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
    fn exact_hash_key(r1: &DnaString, r2: &DnaString) -> u64 {
        let mut hasher = DefaultHasher::new();
        r1.hash(&mut hasher);
        r2.hash(&mut hasher);
        hasher.finish()
    }
    fn min_hits_required(&self, r1len: usize, r2len: usize) -> u16 {
        let r1kmers = 1 + (r1len - KMER_SIZE) / KMER_STRIDE;
        let r2kmers = 1 + (r2len - KMER_SIZE) / KMER_STRIDE;
        let kmers = r1kmers + r2kmers;
        let allowable_errors = self.max_hamming_distance as usize;
        let mismatched_kmers_per_error: usize = KMER_SIZE / std::cmp::min(KMER_STRIDE, KMER_SIZE);
        let min_kmers = kmers - allowable_errors * mismatched_kmers_per_error;
        min_kmers as u16
    }
    pub fn find_closest_consensus(&self, r1: &DnaString, qual1: &[u8], r2: &DnaString, qual2: &[u8]) -> Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)> {
        // Exact match
        match self.exact_lookup.get(&DeduplicationLookup::exact_hash_key(r1, r2)) {
            Some(exact_id) => {
                let edit_distance = self.calc_edit_distances(*exact_id, &self.seq[*exact_id as usize], &r1, qual1, &r2, qual2);
                if edit_distance.1 == 0 {
                    return Some(edit_distance);
                }
            },
            None => {},
        };
        // Find a match with errors allowed
        let mut kmer_matches = FxHashMap::default();
        let mut uninformative_hits = 0;
        uninformative_hits += self.add_matching_kmers(&self.lookup[0], r1, &mut kmer_matches);
        uninformative_hits += self.add_matching_kmers(&self.lookup[1], r2, &mut kmer_matches);
        let mut best_match = None;
        let min_hits_required = self.min_hits_required(r1.len(), r2.len());
        //println!("ConsensusesWithHits={}", kmer_matches.len());
        for (id, hits) in kmer_matches {
            let hits = hits + uninformative_hits; // just assume everything matches for positions with many matches
            if hits >= min_hits_required {
                let edit_distance = self.calc_edit_distances(id, &self.seq[id as usize], &r1, qual1, &r2, qual2);
                // println!("edit_distance={},{},{}", edit_distance.1, edit_distance.2, edit_distance.3);
                if edit_distance.1 <= self.max_hamming_distance
                        && edit_distance.2 <= self.max_phred_distance
                        && edit_distance.3 <= self.max_per_base_phred_distance
                        && self.cmp_distance_abundance(best_match, Some(edit_distance)) == Ordering::Less {
                    best_match = Some(edit_distance);
                }
            }
        }
        best_match
    }
    fn calc_edit_distances(&self, id: ConsensusSequenceId, consensus: &ConsensusSequence, r1: &DnaString, qual1: &[u8], r2: &DnaString, qual2: &[u8]) -> (ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance) {
        let hamming = hamming_distance_common_bases(&consensus.seq[0], r1) + hamming_distance_common_bases(&consensus.seq[1], r2);
        if hamming > self.max_hamming_distance {
            return (id, hamming, PhredDistance::MAX, MAX_PHRED_ENCODABLE as PerBasePhredDistance);
        }
        let (phred_distance1, shared_bases1) = DeduplicationLookup::phred_distance(&consensus.seq[0], r1, qual1);
        let (phred_distance2, shared_bases2) = DeduplicationLookup::phred_distance(&consensus.seq[1], r2, qual2);
        (id, hamming, phred_distance1 + phred_distance2, ((phred_distance1 + phred_distance2) as PerBasePhredDistance) / ((shared_bases1 + shared_bases2) as PerBasePhredDistance))
    }
    fn phred_distance(seq: &DnaString, r: &DnaString, qual: &[u8]) -> (PhredDistance, usize) {
        let shared_bases = cmp::min(seq.len(), r.len());
        let mut phred_distance = 0;
        for i in 0..shared_bases {
            if seq.get(i) != r.get(i) {
                phred_distance += (qual[i] - FASTQ_PHRED_ENCODING) as PhredDistance;
            }
        }
        (phred_distance, shared_bases)
    }
    /*
    fn cmp_distance_per_base_phred(&self, a: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>, b: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>) -> Ordering {
        match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Less,
            (Some(_), None) => Ordering::Greater,
            (Some((_a_id, _a_hamming, _a_phred, a_per_base_phred_errors)),
             Some((_b_id, _b_hamming, _b_phred, b_per_base_phred_errors))) =>
            b_per_base_phred_errors.partial_cmp(&a_per_base_phred_errors).unwrap()
        }
    }
    */
    fn cmp_distance_abundance(&self, a: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>, b: Option<(ConsensusSequenceId, HammingDistance, PhredDistance, PerBasePhredDistance)>) -> Ordering {
        match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Less,
            (Some(_), None) => Ordering::Greater,
            (Some((a_id, _a_hamming, _a_phred, _a_per_base_phred_errors)),
             Some((b_id, _b_hamming, _b_phred, _b_per_base_phred_errors))) =>
                self.seq[a_id as usize].read_count().partial_cmp(&self.seq[b_id as usize].read_count()).unwrap()
        }
    }
    // returns the number of uninformative positions that were skipped because there were too many matching consensuses
    fn add_matching_kmers(&self, lookup: &Vec<MultiMap<Kmer, ConsensusSequenceId>>, seq: &DnaString, counts: &mut FxHashMap<ConsensusSequenceId, u16>) -> u16 {
        let mut uninformative = 0;
        let mut i = 0;
        static EMPTY_LOOKUP : Vec<ConsensusSequenceId> = Vec::new();
        while i * KMER_STRIDE + KMER_SIZE <= seq.len() && i < lookup.len() {
            let read_offset = i * KMER_STRIDE;
            let kmer: Kmer = seq.get_kmer(read_offset);
            let consensuses_containing_kmer_at_position = lookup[i].get_vec(&kmer).unwrap_or(&EMPTY_LOOKUP);
            // TODO: exclude consensus with too many matching kmers
            if consensuses_containing_kmer_at_position.len() > self.uninformative_kmer_threshold {
                uninformative += 1;
            } else {
                for id in consensuses_containing_kmer_at_position {
                    match counts.entry(*id) {
                        std::collections::hash_map::Entry::Occupied(mut entry) => {
                            *entry.get_mut() += 1;
                        },
                        std::collections::hash_map::Entry::Vacant(entry) => {
                            entry.insert(1);
                        },
                    }
                }
            }
            i += 1;
        }
        uninformative
    }
    /// Update the lookup for this consensus
    /// Implementation note: we don't know what's already in the lookup
    /// I've chosen not to remove the 'old' consensus kmers as it's simpler and doesn't really change the outcome
    /// The consensus shouldn't change often and more error kmers actually gives a slightly better chance
    /// of the kmer lookup matching
    fn update_lookup(&mut self, id: ConsensusSequenceId) {
        for rindex in [0, 1] {
            DeduplicationLookup::update_lookup_read(&mut self.lookup[rindex], &self.seq[id as usize].seq[rindex], id, false);
        }
        self.exact_lookup.insert(DeduplicationLookup::exact_hash_key(&self.seq[id as usize].seq[0], &self.seq[id as usize].seq[0]), id);
    }
    fn update_lookup_read(lookup: &mut Vec<MultiMap<Kmer, ConsensusSequenceId>>, seq: &DnaString, id: ConsensusSequenceId, ensure_no_double_counting: bool) {
        while lookup.len() * KMER_STRIDE + KMER_SIZE <= seq.len() {
            lookup.push(MultiMap::new());
        }
        let mut i = 0;
        while i * KMER_STRIDE + KMER_SIZE <= seq.len() {
            let read_offset = i * KMER_STRIDE;
            let kmer: Kmer = seq.get_kmer(read_offset);
            match lookup[i].get_vec_mut(&kmer) {
                Some(vector) if !ensure_no_double_counting || !vector.contains(&id) => {
                    vector.push(id);
                },
                None => lookup[i].insert(kmer, id),
                // already in lookup and don't want to double-count it
                _ => {},
            }
            i += 1;
        }
    }
}
struct DirtyHackDnaString {
    storage: Vec<u64>,
    len: usize,
}
fn hamming_distance_common_bases(a: &DnaString, b: &DnaString) -> HammingDistance {
    if a.len() > b.len() {
        // ensure a is the shorter of the two sequences
        return hamming_distance_common_bases(b, a);
    }
    let ahack: &DirtyHackDnaString = unsafe {
        std::mem::transmute(a)
    };
    let bhack: &DirtyHackDnaString = unsafe {
        std::mem::transmute(b)
    };
    if ahack.len == 0 {
        return 0;
    }
    let mut edit_distance = 0;
    let final_bin = ahack.storage.len() - 1;
    if final_bin > 0 {
        // all bins except the last are full
        for i in 0..final_bin {
            edit_distance += count_diff_2_bit_packed(ahack.storage[i], bhack.storage[i]);
        }
    }
    // ignore all bases past the length of a:
    let final_bin_base_count = ((a.len() - 1) & 0x1F) + 1;
    let final_bin_a_bases = ahack.storage[final_bin] >> (64 - 2 * final_bin_base_count);
    let final_bin_b_bases = bhack.storage[final_bin] >> (64 - 2 * final_bin_base_count);
    edit_distance += count_diff_2_bit_packed(final_bin_a_bases, final_bin_b_bases);
    edit_distance
}
/// count Hamming distance between 2 2-bit DNA packed u64s
#[inline]
fn count_diff_2_bit_packed(a: u64, b: u64) -> u32 {
    let bit_diffs = a ^ b;
    let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & 0x5555555555555555;
    two_bit_diffs.count_ones()
}
const FASTQ_PHRED_ENCODING: u8 = 33;
const MAX_PHRED_ENCODABLE: u8 = 93;
impl ConsensusSequence {
    fn fastq_qual_score(seq: &DnaString, basequals: &Vec<[BaseQualAccumulation;4]>) -> Vec<u8> {
        let mut out = Vec::with_capacity(seq.len());
        for i in 0..seq.len() {
            let encoded_base = seq.get(i);
            let score: i32 = 2 * basequals[i][encoded_base as usize] as i32 - basequals[i].iter().sum::<BaseQualAccumulation>() as i32;
            // 0-93 allowed: https://en.wikipedia.org/wiki/FASTQ_format
            out.push(FASTQ_PHRED_ENCODING + std::cmp::max(0, std::cmp::min(MAX_PHRED_ENCODABLE as i32, score)) as u8);
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
            for _ in seq.len()..r.len() {
                q.push([0; 4]);
            }
            seq.extend(r.slice(seq.len(), r.len()).iter());
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
    fn to_p(seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|x| x + 33).collect()
    }
    fn from_p(seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|x| x - 33).collect()
    }
    #[test]
    fn test_hamming_distance_common_bases_simple() {
        assert_eq!(1, hamming_distance_common_bases(&s("ACGT"), &s("ACGA")));
        assert_eq!(2, hamming_distance_common_bases(&s("ACGT"), &s("ACTA")));
        assert_eq!(1, hamming_distance_common_bases(&s("ACG"), &s("ACTA")));
        assert_eq!(33, hamming_distance_common_bases(&s("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), &s("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")));
        assert_eq!(32, hamming_distance_common_bases(&s("AAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAA"), &s("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")));
    }
    #[test]
    fn test_edit_distances() {
        assert_eq!((1, 4), DeduplicationLookup::phred_distance(&s("ACGT"), &s("ACGA"), &[33,33,33,34]));
        assert_eq!((2, 4), DeduplicationLookup::phred_distance(&s("ACGT"), &s("ACGA"), &[33,33,33,35]));
        assert_eq!((5, 4), DeduplicationLookup::phred_distance(&s("ACGT"), &s("ACTA"), &[33,33,34,37]));
        assert_eq!((7, 3), DeduplicationLookup::phred_distance(&s("ACG"), &s("ACTA"), &[33,33,40,37]));
    }
    #[test]
    fn test_consensus_update() {
        let mut con = ConsensusSequence::new(s("ACGT"), &[33,34,35,36], s("AACCGGTT"), &[37,38,39,40,41,42,43,44]);
        assert_eq!("ACGT", con.seq[0].to_string());
        assert_eq!("AACCGGTT", con.seq[1].to_string());
        assert_eq!(vec![0, 1, 2, 3], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![4, 5, 6, 7, 8, 9, 10, 11], from_p(&con.paired_fastq_qual_score().1));
        assert_eq!(false, con.add_read_pair(&s("ACGT"), &[33,34,35,36], &s("AACCGGTT"), &[37,38,39,40,41,42,43,44]));
        assert_eq!("ACGT", con.seq[0].to_string());
        assert_eq!("AACCGGTT", con.seq[1].to_string());
        assert_eq!(vec![0, 2, 4, 6], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![8, 10, 12, 14, 16, 18, 20, 22], from_p(&con.paired_fastq_qual_score().1));
        assert_eq!(true, con.add_read_pair(&s("ACGC"), &[34,34,34,63], &s("TTGGAAGG"), &[34,34,34,34,34,34,34,34]));
        assert_eq!("ACGC", con.seq[0].to_string());
        assert_eq!("AACCGGTT", con.seq[1].to_string());
        assert_eq!(vec![1, 3, 5, 30-6], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![7, 9, 11, 13, 15, 17, 19, 21], from_p(&con.paired_fastq_qual_score().1));
    }
    #[test]
    fn max_phred_is_93() {
        let mut con = ConsensusSequence::new(s("A"), to_p(&[90]).as_slice(), s("T"), &[33+90]);
        con.add_read_pair(&s("A"), &[50], &s("T"), &[50]);
        assert_eq!(vec![93], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![93], from_p(&con.paired_fastq_qual_score().1));
    }
    #[test]
    fn min_phred_is_0() {
        let mut con = ConsensusSequence::new(s("A"), &[34], s("T"), &[33]);
        con.add_read_pair(&s("C"), &[34], &s("C"), &[33]);
        con.add_read_pair(&s("G"), &[35], &s("G"), &[33]);
        con.add_read_pair(&s("T"), &[34], &s("A"), &[33]);
        assert_eq!("G", con.seq[0].to_string());
        assert_eq!("T", con.seq[1].to_string());
        assert_eq!(vec![0], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![0], from_p(&con.paired_fastq_qual_score().1));
    }
    #[test]
    fn fastq_qual_score() {
        assert_eq!(vec![0,3,93,4], from_p(&ConsensusSequence::fastq_qual_score(&s("ACGT"), &vec![
            [0,1,2,3],
            [1,6,1,1],
            [0,100,200,0],
            [1,4,1,10],
            ])));
    }
    #[test]
    fn consensus_add_append() {
        let mut con = ConsensusSequence::new(s("A"), &[63], s("CG"), &[63, 34]);
        con.add_read_pair(&s("AT"), &[63, 63], &s("CGG"), &[64, 64, 64]);
        assert_eq!("AT", con.seq[0].to_string());
        assert_eq!("CGG", con.seq[1].to_string());
        assert_eq!(vec![60, 30], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![30+31, 1+31, 31], from_p(&con.paired_fastq_qual_score().1));
    }
    #[test]
    fn consensus_add_shorter() {
        let mut con = ConsensusSequence::new(s("AT"), &[63, 63], s("CG"), &[63, 63]);
        con.add_read_pair(&s("G"), &[34], &s("T"), &[36]);
        con.add_read_pair(&s("C"), &[35], &s("G"), &[37]);
        assert_eq!("AT", con.seq[0].to_string());
        assert_eq!("CG", con.seq[1].to_string());
        assert_eq!(vec![30-1-2, 30], from_p(&con.paired_fastq_qual_score().0));
        assert_eq!(vec![30-3-4, 30], from_p(&con.paired_fastq_qual_score().1));
    }
    #[test]
    fn test_consensus_lookup() {
        let mut lookup = DeduplicationLookup::new();
        // SRR13320978.1
        lookup.add_read_pair(
            s("GCATGAGAACAAGAAGAACAGTACGTCCTGCAGTGAGTCGAAGTCAAGATAGGAAGAACACACGTCTGAACAACAGTCACA"),
            "/E//6///EA//6A/EE<E/E/EE</A6//A</6//////E///////E6//////A/6EE/6///////<//6A//6///".as_bytes(),
            s("CTTCAAGTTCTCAAGCTGCTAGAGATTTTTCCACACTGACTAAATGTTCTGAGGGATCTCTAGTTACCAGAGTCAGTCGTCTTGTAAATCAAAATAGAAACAAAAAAAAAAAAAAAAAAA"),
            "AAA6/////A<///</AAAA///<////66/A///////<A666///</////E6/////////<//////////////<////6/////A6///////////////6A///////A///".as_bytes());
        lookup.add_read_pair(
            s("TTCCATTACATTCAATTCTATTCCATTCCATTCAAATCCAATCAGTTCAATTCCATACCATTACAATCTATTCAGTATAAT"),
            "A/E/EAE//6A///AE//A<A/A<A6//<6E//AE////6/<<///A<//A//<<</</</6/6/666//6////6////<".as_bytes(),
            s("ATTAACCACAATAGAATGAAATAGAATAAAATAGAACAAAATAAAATAAAATAAAATAAAATAAAATAAAATCAAATAAAAAACTAAAAATAATTCAACAAAAAATAAATAAAATAAAAA"),
            "//////////////////////////////////////////A/////////A////A/////A////6/////////A///6///</<//////////////////////////<////".as_bytes());
        assert_eq!(2, lookup.seq.len());
        // exact match
        lookup.add_read_pair(
            s("TTCCATTACATTCAATTCTATTCCATTCCATTCAAATCCAATCAGTTCAATTCCATACCATTACAATCTATTCAGTATAAT"),
            "A/E/EAE//6A///AE//A<A/A<A6//<6E//AE////6/<<///A<//A//<<</</</6/6/666//6////6////<".as_bytes(),
            s("ATTAACCACAATAGAATGAAATAGAATAAAATAGAACAAAATAAAATAAAATAAAATAAAATAAAATAAAATCAAATAAAAAACTAAAAATAATTCAACAAAAAATAAATAAAATAAAAA"),
            "//////////////////////////////////////////A/////////A////A/////A////6/////////A///6///</<//////////////////////////<////".as_bytes());
        assert_eq!(2, lookup.seq.len());
        // 1 mismatch
        lookup.add_read_pair(
            s("ATCCATTACATTCAATTCTATTCCATTCCATTCAAATCCAATCAGTTCAATTCCATACCATTACAATCTATTCAGTATAAT"),
            "A/E/EAE//6A///AE//A<A/A<A6//<6E//AE////6/<<///A<//A//<<</</</6/6/666//6////6////<".as_bytes(),
            s("ATTAACCACAATAGAATGAAATAGAATAAAATAGAACAAAATAAAATAAAATAAAATAAAATAAAATAAAATCAAATAAAAAACTAAAAATAATTCAACAAAAAATAAATAAAATAAAAA"),
            "//////////////////////////////////////////A/////////A////A/////A////6/////////A///6///</<//////////////////////////<////".as_bytes());
        assert_eq!(2, lookup.seq.len());
        // 2 mismatches (r1:2, r2:2)
        lookup.add_read_pair(
            s("TGCCATTACATTCAATTCTATTCCATTCCATTCAAATCCAATCAGTTCAATTCCATACCATTACAATCTATTCAGTATAAT"),
            "A/E/EAE//6A///AE//A<A/A<A6//<6E//AE////6/<<///A<//A//<<</</</6/6/666//6////6////<".as_bytes(),
            s("AGTAACCACAATAGAATGAAATAGAATAAAATAGAACAAAATAAAATAAAATAAAATAAAATAAAATAAAATCAAATAAAAAACTAAAAATAATTCAACAAAAAATAAATAAAATAAAAA"),
            "//////////////////////////////////////////A/////////A////A/////A////6/////////A///6///</<//////////////////////////<////".as_bytes());
        assert_eq!(2, lookup.seq.len());
        // 10 mismatches
        lookup.add_read_pair(
            s("TGCCAATACATTCAATGCTATTCCAGTCCATTCAAAACCAATCAGTACAATTCCAAACCATAACAATCTATTCAGAATAAT"),
            "A/E/EAE//6A///AE//A<A/A<A6//<6E//AE////6/<<///A<//A//<<</</</6/6/666//6////6////<".as_bytes(),
            s("AGTAACCACAATAGAATGAAATAGAATAAAATAGAACAAAATAAAATAAAATAAAATAAAATAAAATAAAATCAAAAAAAAAACTAAAAATAATTCAACAAAAAAAAAATAAAATAAAAA"),
            "//////////////////////////////////////////A/////////A////A/////A////6/////////A///6///</<//////////////////////////<////".as_bytes());
        assert_eq!(3, lookup.seq.len());
    }
    #[test]
    fn test_hamming_distance_common_bases() {
        let mut a  = String::new();
        for _ in 0..=130 {
            a.push('A');
            let dna = DnaString::from_acgt_bytes(a.as_bytes());
            assert_eq!(0, hamming_distance_common_bases(&dna, &dna));
            let mut b = String::new();
            for _ in 0..=130 {
                b.push('T');
                let dnb = DnaString::from_acgt_bytes(b.as_bytes());
                let len = std::cmp::min(a.len(), b.len());
                let _stra = &a[0..len];
                let _strb = &b[0..len];
                assert_eq!(len as u32, hamming_distance_common_bases(&dna, &dnb));
                assert_eq!(len as u32, hamming_distance_common_bases(&dnb, &dna));
                assert_eq!(dna.slice(0, len).hamming_dist(&dnb.slice(0, len)), hamming_distance_common_bases(&dna, &dnb));
            }
        }
    }
    #[test]
    fn cmp_distance_abundance() {
        let mut lookup = DeduplicationLookup::new();
        lookup.max_hamming_distance = 1;
        lookup.add_read_pair(
            s("AAAAAAAAAAAAAAAAT"),
            "AAAAAAAAAAAAAAAAT".as_bytes(),
            s("AAAAAAAAAAAAAAAAT"),
            "AAAAAAAAAAAAAAAAT".as_bytes());
        lookup.add_read_pair(
            s("AAAAAAAAAAAAAAAAG"),
            "AAAAAAAAAAAAAAAAA".as_bytes(),
            s("AAAAAAAAAAAAAAAAG"),
            "AAAAAAAAAAAAAAAAA".as_bytes());
        lookup.add_read_pair(
            s("AAAAAAAAAAAAAAAAG"),
            "AAAAAAAAAAAAAAAAA".as_bytes(),
            s("AAAAAAAAAAAAAAAAG"),
            "AAAAAAAAAAAAAAAAA".as_bytes());
        assert_eq!(2, lookup.seq.len());
        lookup.add_read_pair(
            s("AAAAAAAAAAAAAAAAT"),
            "AAAAAAAAAAAAAAAAA".as_bytes(),
            s("AAAAAAAAAAAAAAAAG"),
            "AAAAAAAAAAAAAAAAA".as_bytes());
        assert_eq!(Ordering::Less, lookup.cmp_distance_abundance(Some((0, 0, 0, 0.0)), Some((1, 1, 1, 1.0))));
        assert_eq!(3, lookup.seq[1].read_count());
    }
}
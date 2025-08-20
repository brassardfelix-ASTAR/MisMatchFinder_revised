#[derive(Clone)]
pub struct Fragment {
    reads: Vec<GappedRead>,
    chrom: String,
    insertions: Vec<Mismatch>,
    // we store this so we know if we only need to look at read1 instead of both
    only_overlap: bool,
}

use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};

use lazy_static::lazy_static;
use regex::Regex;

use crate::bamreader::mismatch::MismatchType;

use super::cigar::GappedRead;
use super::mismatch::Mismatch;

fn is_real_sbs(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    // Must be a tri-nucleotide context with no gaps on neighbors
    // and the center base must differ.
    ref_seq.len() == 3
        && alt_seq.len() == 3
        && ref_seq[0] != b'-' && ref_seq[2] != b'-'
        && alt_seq[0] != b'-' && alt_seq[2] != b'-'
        && ref_seq[1] != alt_seq[1]
}

impl Fragment {

    pub fn make_fragment(
        mut read1: GappedRead,
        mut read2: GappedRead,
        only_overlap: bool,
        strict: bool,
        seqname: &String,
    ) -> Option<Fragment> {
        //calculate the offset and the overlap for both reads
        let read1_start = read1.start();
        let read1_end = read1.end();

        let read2_start = read2.start();
        let read2_end = read2.end();

        let overlap = max(read1_start, read2_start)..min(read1_end, read2_end);

        //if there is no overlap, but the only_overlap flag is set, we can just return None here
        let full_fragment;
        if only_overlap && overlap.is_empty() {
            return None;
        } else {
            //otherwise we continue building the full fragment
            full_fragment = min(read1_start, read2_start)..max(read1_end, read2_end);
        }

        for ref_pos in full_fragment {
            let read1_pos = read1.get_read_pos_from_ref(ref_pos);
            let read2_pos = read2.get_read_pos_from_ref(ref_pos);

            let read1_pos = match read1_pos {
                Some(v) => *v as i64,
                None => -1,
            };
            let read2_pos = match read2_pos {
                Some(v) => *v as i64,
                None => -1,
            };

            // if one of them is <0 means we have no overlap
            if read1_pos < 0 || read2_pos < 0 {
                //if we only want to analyse the overlap, we just reduce the quality of these positions to 0
                if only_overlap {
                    if read1_pos >= 0 {
                        read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                    } else if read2_pos >= 0 {
                        read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                    } else {
                        // gap in both reads
                    }
                }
                continue;
            }

            //now both reads are overlapping and we can make a consensus

            let r1_seq = read1.get_read_seq()[read1_pos as usize];
            let r2_seq = read2.get_read_seq()[read2_pos as usize];

            if r1_seq == r2_seq {
                // we update the quality to reflect the additional support in read1
                read1.get_read_qual_mut()[read1_pos as usize] =
                    read1.get_read_qual()[read1_pos as usize] + read2.get_read_qual()[read2_pos as usize];
                // we set the quality to 0 in the second read, to not get the mismatch twice
                read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
            } else if strict {
                // if strict and there is an overlap, we set both quals to 0 to drop this site
                read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
            } else {
                // pick the better-supported base; downweight the other one
                let r1_qual = read1.get_read_qual()[read1_pos as usize];
                let r2_qual = read2.get_read_qual()[read2_pos as usize];

                if r1_qual > r2_qual {
                    read1.get_read_qual_mut()[read1_pos as usize] = r1_qual - (r2_qual / 2);
                    read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                } else if r2_qual > r1_qual || read1.is_reverse() {
                    read1.get_read_qual_mut()[read1_pos as usize] = 0u8;
                    read2.get_read_qual_mut()[read2_pos as usize] = r2_qual - (r1_qual / 2);
                } else {
                    // equal qual, prefer forward read (read1 is forward here)
                    read1.get_read_qual_mut()[read1_pos as usize] = r1_qual - (r2_qual / 2);
                    read2.get_read_qual_mut()[read2_pos as usize] = 0u8;
                }
            }
        }

        // and now we check the inserts for overlaps
        // create final storage
        let mut fragment_ins: HashMap<u32, Mismatch> = HashMap::new();
        let mut insert_pos: HashSet<u32> = HashSet::new();

        for mm in read1.get_insertions().iter().cloned() {
            if only_overlap && !overlap.contains(&(mm.position as i64)) {
                continue;
            } else {
                insert_pos.insert(mm.position);
                fragment_ins.insert(mm.position, mm);
            }
        }

        //going through read2 and checking if we already had it
        for mm in read2.get_insertions().iter().cloned() {
            if only_overlap && !overlap.contains(&(mm.position as i64)) {
                continue;
            } else {
                if fragment_ins.contains_key(&mm.position) {
                    let entry = fragment_ins.get_mut(&mm.position).unwrap();
                    entry.quality = entry.quality + mm.quality;
                    insert_pos.remove(&mm.position);
                } else {
                    if strict {
                        continue;
                    } else {
                        fragment_ins.insert(mm.position, mm);
                    }
                }
            }
        }

        for pos in insert_pos {
            fragment_ins.remove(&pos);
        }

        Some(Fragment {
            reads: Vec::from([read1, read2]),
            chrom: seqname.to_string(),
            insertions: fragment_ins.into_values().collect(),
            only_overlap,
        })
    }

    pub fn strand(&self) -> Option<char> {
        self.reads.first().map(|r| if r.is_reverse() { '-' } else { '+' })
    }

    pub fn is_overlapping(&self) -> bool {
        if self.reads.len() < 2 {
            return false;
        }
        let r1 = &self.reads[0];
        let r2 = &self.reads[1];
        let start1 = r1.start();
        let end1 = r1.end();
        let start2 = r2.start();
        let end2 = r2.end();

        start1 <= end2 && start2 <= end1
    }

    pub fn reads_agree(&self, mm: &Mismatch) -> Option<bool> {
        if self.reads.len() < 2 {
            return None;
        }
        let pos = mm.pos_in_read?;
        let r1 = &self.reads[0];
        let r2 = &self.reads[1];

        let base1 = r1.get_read_seq().get(pos)?;
        let base2 = r2.get_read_seq().get(pos)?;

        Some(base1 == base2)
    }

    pub fn make_se_fragment(read: GappedRead, chrom: &String) -> Self {
        Fragment {
            insertions: read.get_insertions().to_vec(),
            reads: Vec::from([read]),
            chrom: chrom.to_string(),
            only_overlap: false,
        }
    }

    pub fn get_mismatches(&self, min_base_qual: u8) -> Vec<Mismatch> {
        let mut ret: Vec<Mismatch> = Vec::new();

        for (i, read) in self.reads.iter().enumerate() {
            //if we only want the overlap, we can just skip read 2
            if self.only_overlap && i > 0 {
                break;
            }

            let md = read.get_md_str();
            //if we have a full integer as MD string we can skip, because we have no changes
            let md_parse: Result<i64, _> = md.parse();
            match md_parse {
                Ok(_) => return ret,
                Err(v) => v, //we have to work
            };

            lazy_static! {
                static ref RE: Regex =
                    Regex::new(r"(?P<digit>[0-9]+)(\^?)(?P<change>[A-Z]+)?").unwrap();
            }

            let ref_pos = read.get_read_pos();
            let mut ref_index = 0;

            //because we cant make the iterator peekable, we store the intermediate result, and only insert at the end;
            let mut prev_mismatch: Option<Mismatch> = None;
            let mut prev_pos: usize = 0;

            let upper_bound = read.len() - 1;

            for elem in RE.captures_iter(md) {
                //skip all positions which are just matched
                ref_index += &elem["digit"].parse().unwrap();

                if ref_index >= upper_bound || ref_index == 0 {
                    //there is no way we can make a tri nucleotide here
                    ref_index += 1;
                    continue;
                }

                //if we have a center group, we have a deletion
                let deletion = match elem.get(2) {
                    None => false,
                    Some(v) => !v.as_str().is_empty(),
                };

                //if the quality of the base is low, we stop here
                if !deletion && read.get_read_qual()[ref_index] < min_base_qual {
                    ref_index += 1;
                    continue;
                }

                //store the change
                let change = match elem.get(3) {
                    None => {
                        break; //nothing more to do
                    }
                    Some(v) => v.as_str().as_bytes(),
                };

                if deletion {
                    // if we have something stored, we push that as well, as there is no DBS with deletions
                    if let Some(v) = prev_mismatch {
                        ret.push(v);
                    }

                    //we use the quality of the neighbouring bases as an indication of the qual of the deletion
                    let qual = (read.get_read_qual()[ref_index - 1]
                        + read.get_read_qual()[ref_index + change.len()])
                        / 2;
                    if qual > min_base_qual {
                        // we take the base just before the deletion as the alt
                        let alt_seq = read.get_read_seq()[(ref_index - 1)..ref_index].to_vec();
                        // and for the reference, we have to append what the read is missing
                        let mut ref_seq = alt_seq.clone();
                        ref_seq.append(change.to_vec().as_mut());

                        // deletions remove bases from the read; there is not a single base in read at the deleted run
                        // safest is to mark pos_in_read as unknown here
                        let pos_in_read = None;

                        // DEL: position is leftmost base (1-based)
                        ret.push(Mismatch {
                            chromosome: self.chrom.to_string(),
                            rid: read.tid(),
                            position: (read.start() + ref_index as i64) as u32,
                            reference: ref_seq,
                            alternative: alt_seq,
                            quality: qual,
                            typ: MismatchType::DEL,

                            // New fields
                            fragment_length: Some(self.fragment_length()),
                            read_orientation: Some(self.read_orientation()),
                            region_tag: None,
                            pos_in_read,
                        });
                    }
                    prev_mismatch = None;
                    prev_pos = 0;
                    //step over the deletion to be at the base after
                    ref_index += change.len();
                } else {
                    // otherwise we have single base substitutions
                    let change = change[0];

                    // require continuous reference stretch
                    if ref_pos[ref_index - 1] + 1 == ref_pos[ref_index]
                        && ref_pos[ref_index] + 1 == ref_pos[ref_index + 1]
                    {
                        //get the read sequence (which is the alternative)
                        let mut alt_seq =
                            read.get_read_seq()[(ref_index - 1)..(ref_index + 2)].to_vec();

                        // drop cases where neighbors are gaps
                        if alt_seq[0] == b'-' || alt_seq[2] == b'-' {
                            ref_index += 1;
                            continue;
                        }

                        let mut ref_seq = alt_seq.clone();
                        //and change according to MD to get the ref string
                        ref_seq[1] = change;

                        // ROOT FIX: if the center base matches the reference, this is NOT a mismatch.
                        if !is_real_sbs(&ref_seq, &alt_seq) {
                            // If a previous mismatch was pending, flush it (it was valid),
                            // and skip this position.
                            if let Some(v) = prev_mismatch.take() {
                                ret.push(v);
                            }
                            // advance past this MD position
                            ref_index += 1;
                            continue;
                        }

                        //if we have that, we need to check if there is another base change just behind
                        if let Some(mm) = prev_mismatch {
                            if mm.position == (read.start() + ref_index as i64 - 1) as u32 {
                                if mm.typ == MismatchType::DBS {
                                    prev_mismatch = None;
                                    prev_pos = ref_index;
                                } else {
                                    //alt_seq will be shortened to only 2 bases
                                    alt_seq.pop();
                                    // but we need to update the ref_seq with the center from the previous
                                    ref_seq[0] = mm.reference[1];
                                    // and shorten as well
                                    ref_seq.pop();

                                    // GUARD: the second (current) base must be a real mismatch
                                    if ref_seq[1] == alt_seq[1] {
                                        // new site isn't a true mismatch — keep the previous SBS
                                        ret.push(mm);
                                        prev_mismatch = None;
                                        ref_index += 1;
                                        continue;
                                    }

                                    // DBS harmonization:
                                    // - position: LEFTMOST base (1-based) => start + ref_index
                                    // - pos_in_read: LEFTMOST base
                                    let left_ref_pos = ref_pos[ref_index - 1];
                                    let pos_in_read = read.ref_to_read_pos(left_ref_pos);

                                    prev_mismatch = Some(Mismatch {
                                        quality: (mm.quality + read.get_read_qual()[ref_index]) / 2,
                                        chromosome: mm.chromosome,
                                        rid: read.tid(),
                                        position: (read.start() + ref_index as i64) as u32,
                                        reference: ref_seq,
                                        alternative: alt_seq,
                                        typ: MismatchType::DBS,

                                        // New fields
                                        fragment_length: Some(self.fragment_length()),
                                        read_orientation: Some(self.read_orientation()),
                                        region_tag: None,
                                        pos_in_read,
                                    });

                                    prev_pos = ref_index;
                                }
                            } else {
                                // if it isnt, then we push the old and build a new mismatch
                                ret.push(mm);

                                let pos_in_read = read.ref_to_read_pos(ref_pos[ref_index]);

                                // SBS: position is CENTER base (1-based) => start + ref_index + 1
                                prev_mismatch = Some(Mismatch {
                                    quality: read.get_read_qual()[ref_index],
                                    chromosome: self.chrom.to_string(),
                                    rid: read.tid(),
                                    position: (read.start() + ref_index as i64 + 1) as u32,
                                    reference: ref_seq,
                                    alternative: alt_seq,
                                    typ: MismatchType::SBS,

                                    // New fields
                                    fragment_length: Some(self.fragment_length()),
                                    read_orientation: Some(self.read_orientation()),
                                    region_tag: None,
                                    pos_in_read,
                                })
                            }
                        } else {
                            // no previous mismatch; create one if not immediately contiguous with a previous pos
                            if ref_index - 1 == prev_pos {
                                prev_pos = ref_index;
                            } else {
                                //get the read sequence
                                let mut ref_seq =
                                    read.get_read_seq()[ref_index - 1..ref_index + 2].to_vec();
                                // change according to MD to get the ref string
                                ref_seq[1] = change;

                                let pos_in_read = read.ref_to_read_pos(ref_pos[ref_index]);

                                // SBS: position is CENTER base (1-based) => start + ref_index + 1
                                prev_mismatch = Some(Mismatch {
                                    quality: read.get_read_qual()[ref_index],
                                    chromosome: self.chrom.to_string(),
                                    rid: read.tid(),
                                    position: (read.start() + ref_index as i64 + 1) as u32,
                                    reference: ref_seq,
                                    alternative: alt_seq,
                                    typ: MismatchType::SBS,

                                    // New fields
                                    fragment_length: Some(self.fragment_length()),
                                    read_orientation: Some(self.read_orientation()),
                                    region_tag: None,
                                    pos_in_read,
                                })
                            }
                        }
                    } else {
                        // reference not continuous; flush any previous mismatch
                        if let Some(v) = prev_mismatch {
                            ret.push(v);
                        }
                        prev_mismatch = None;
                    }

                    // step over the mismatch
                    ref_index += 1;
                }
            }
            // and if we have something left in the end, we push that as well
            if let Some(v) = prev_mismatch {
                ret.push(v);
            }
        }

        // go through all of the insertions we already collected before and add them to the final result if they are of high enough qual
        for mut mm in self.insertions.iter().cloned() {
            if mm.quality >= min_base_qual {
                mm.fragment_length = Some(self.fragment_length());
                mm.read_orientation = Some(self.read_orientation());
                mm.region_tag = None;
                // Map the mismatch position to *any* read in this fragment
                mm.pos_in_read = self.ref_to_read_pos_any(mm.position as i64);
                ret.push(mm);
            }
        }

        // Final belt & suspenders: never emit no-op variants
        ret.retain(|m| m.reference != m.alternative);

        ret
    }

    pub fn average(numbers: &[u8]) -> f32 {
        numbers.iter().map(|x| *x as u64).sum::<u64>() as f32 / numbers.len() as f32
    }

    pub fn fragment_length(&self) -> u32 {
        let start = self.reads.first().unwrap().start();
        let end = self.reads.last().unwrap().end();
        (end - start) as u32
    }

    pub fn read_orientation(&self) -> String {
        self.reads
            .iter()
            .map(|r| if r.is_reverse() { "R" } else { "F" })
            .collect::<Vec<&str>>()
            .join("")
    }

    #[inline]
    fn ref_to_read_pos_any(&self, pos: i64) -> Option<usize> {
        for r in &self.reads {
            if let Some(p) = r.ref_to_read_pos(pos) {
                return Some(p);
            }
        }
        None
    }
}


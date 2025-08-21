// We generate a gapped read how it is aligned to the reference, without insertions
// but the deletions will be shown as '-' in the sequence
#[derive(Clone)]
pub struct GappedRead {
    gapped_seq: Vec<u8>,
    gapped_qual: Vec<u8>,
    gapped_pos: Vec<i64>,
    tid: i32,
    start: i64,
    end: i64,
    ref_to_read_map: HashMap<i64, usize>,
    is_reverse: bool,
    md_str: String,
    inserts: Vec<Mismatch>,
}

impl GappedRead {
    pub fn start(&self) -> i64 {
        self.start
    }

    pub fn end(&self) -> i64 {
        self.end
    }

    pub fn len(&self) -> usize {
        self.gapped_pos.len()
    }

    pub fn get_read_pos_from_ref(&self, ref_pos: i64) -> Option<&usize> {
        self.ref_to_read_map.get(&ref_pos)
    }

    pub fn get_read_pos(&self) -> &Vec<i64> {
        &self.gapped_pos
    }

    pub fn get_read_qual(&self) -> &Vec<u8> {
        &self.gapped_qual
    }

    pub fn get_read_qual_mut(&mut self) -> &mut Vec<u8> {
        &mut self.gapped_qual
    }

    pub fn get_read_seq(&self) -> &Vec<u8> {
        &self.gapped_seq
    }

    // we dont really need this one, as we dont change the sequence (so far)
    pub fn _get_read_seq_mut(&mut self) -> &mut Vec<u8> {
        &mut self.gapped_seq
    }

    pub fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    pub fn get_md_str(&self) -> &String {
        &self.md_str
    }

    pub fn tid(&self) -> i32 {
        self.tid
    }

    pub fn get_insertions(&self) -> &Vec<Mismatch> {
        &self.inserts
    }

    pub fn get_insertions_mut(&mut self) -> &mut Vec<Mismatch> {
        &mut self.inserts
    }
}

use std::{collections::HashMap, panic};

use rust_htslib::bam::{self, record::Cigar};

use super::mismatch::Mismatch;

pub fn parse_cigar_str(read: &bam::Record, chrom: &String) -> GappedRead {
    let mut gapped_seq: Vec<u8> = Vec::new();
    let mut gapped_qual: Vec<u8> = Vec::new();
    let mut gapped_pos: Vec<i64> = Vec::new();
    let mut ref_to_read_map: HashMap<i64, usize> = HashMap::new();

    let mut i: u32 = 0;
    let cigars = read.cigar();

    let mut ref_pos = read.pos();
    let mut internal_pos: usize = 0;

    let mut mismatches: Vec<Mismatch> = Vec::new();

    for op in cigars.into_iter() {
        match op {
            Cigar::Match(n) | Cigar::Diff(n) | Cigar::Equal(n) => {
                // if there is a (mis)match, we use the sequence as is
                for x in 0..*n {
                    let idx = (i + x) as usize;
                    gapped_seq.push(read.seq()[idx]);
                    gapped_qual.push(read.qual()[idx]);
                    gapped_pos.push(idx as i64);

                    //store internal mapping
                    ref_to_read_map.insert(ref_pos as i64, internal_pos);
                    internal_pos += 1;
                    ref_pos += 1;
                }
                i += n;
            }
            Cigar::SoftClip(n) => {
                // if there is a softclip we remove that part of the read
                i += n;
            }
            Cigar::Ins(n) => {
                // if there is an insertion we store the insertion as a mismatch, which we can report on later

                //however we can only do that, if its not the VERY first thing, because we wont know the
                // base before, so here we just drop the mismatch
                if i > 0 {
                    //start from the position before
                    let idx = i as usize - 1;
                    let ref_seq: Vec<u8> = Vec::from([read.seq()[idx]]);
                    let mut alt_seq = Vec::new();
                    let mut qual: u32 = 0;
                    // this is one longer because we need to ref start
                    for x in 0..=*n {
                        let x = idx + x as usize;
                        alt_seq.push(read.seq()[x]);
                        qual += read.qual()[x] as u32;
                    }
                    mismatches.push(Mismatch {
                        chromosome: chrom.to_string(),
                        position: ref_pos as u32,
                        reference: ref_seq,
                        alternative: alt_seq,
                        quality: (qual / (n + 1)) as u8,
                        typ: super::mismatch::MismatchType::INS,
                        rid: read.tid(),
                    });
                }
                // and we remove that part of the read
                i += n;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                //in a deletion, we add dashes to show the deletion
                for _ in 0..*n {
                    gapped_seq.push(b'-');
                    gapped_qual.push(0u8);
                    gapped_pos.push(-1);

                    //store internal mapping
                    ref_to_read_map.insert(ref_pos as i64, internal_pos);
                    internal_pos += 1;
                    ref_pos += 1;
                }
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                //do absolutely nothing
            }
        }
    }

    let md = match read.aux(b"MD") {
        Ok(bam::record::Aux::String(md)) => md.to_string(),
        _ => {
            panic!(
                "Read ({})did not contain MD string",
                std::str::from_utf8(read.qname()).unwrap()
            );
        }
    };

    return GappedRead {
        gapped_seq,
        gapped_qual,
        gapped_pos,
        tid: read.tid(),
        start: read.pos(),
        end: ref_pos,
        ref_to_read_map,
        is_reverse: read.is_reverse(),
        md_str: md,
        inserts: mismatches,
    };
}

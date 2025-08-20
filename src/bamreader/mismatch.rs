use std::{cmp::Ordering, fmt};

#[derive(Debug, Hash, Clone, Eq, Ord, PartialOrd)]
pub enum MismatchType {
    SBS,
    DBS,
    INS,
    DEL,
}

impl PartialEq for MismatchType {
    fn eq(&self, other: &Self) -> bool {
        core::mem::discriminant(self) == core::mem::discriminant(other)
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

#[derive(Debug, Hash, Clone)]
pub struct Mismatch {
    pub chromosome: String,
    pub position: u32,
    pub reference: Vec<u8>,
    pub alternative: Vec<u8>,
    pub quality: u8,
    pub typ: MismatchType,
    pub rid: i32,

    // New fields
    pub fragment_length: Option<u32>,
    pub read_orientation: Option<String>,
    pub region_tag: Option<String>,
    pub pos_in_read: Option<usize>,

}

impl fmt::Display for Mismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{:?}",
            self.chromosome,
            self.position,
            String::from_utf8(self.reference.to_vec()).unwrap(),
            String::from_utf8(self.alternative.to_vec()).unwrap(),
            self.quality,
            self.typ
        )
    }
}

impl PartialOrd for Mismatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Mismatch {
    fn cmp(&self, other: &Self) -> Ordering {
        // check first if the chromsosomes are the same
        if self.chromosome == other.chromosome {
            if self.position < other.position {
                return Ordering::Less;
            } else if self.position > other.position {
                return Ordering::Greater;
            } else {
                if self.reference == other.reference && self.alternative == other.alternative {
                    return Ordering::Equal;
                } else {
                    // if they are the same position, we just compare the length of the variant
                    let len_cmp = (self.reference.len() + self.alternative.len())
                        .cmp(&(other.reference.len() + other.alternative.len()));

                    match len_cmp {
                        Ordering::Equal => self.reference.cmp(&other.reference),
                        Ordering::Greater => Ordering::Greater,
                        Ordering::Less => Ordering::Less,
                    }
                }
            }
        } else {
            self.chromosome.cmp(&other.chromosome)
        }
    }

    fn max(self, other: Self) -> Self {
        // HACK(fee1-dead): go back to using `self.max_by(other, Ord::cmp)`
        // when trait methods are allowed to be used when a const closure is
        // expected.
        match self.cmp(&other) {
            Ordering::Less | Ordering::Equal => other,
            Ordering::Greater => self,
        }
    }

    fn min(self, other: Self) -> Self {
        // HACK(fee1-dead): go back to using `self.min_by(other, Ord::cmp)`
        // when trait methods are allowed to be used when a const closure is
        // expected.
        match self.cmp(&other) {
            Ordering::Less | Ordering::Equal => self,
            Ordering::Greater => other,
        }
    }

    fn clamp(self, min: Self, max: Self) -> Self {
        assert!(min <= max);
        if self < min {
            min
        } else if self > max {
            max
        } else {
            self
        }
    }
}

impl PartialEq for Mismatch {
    fn eq(&self, other: &Mismatch) -> bool {
        //we start with the easy things, this will make the method a bit bloaty, but more performant
        if self.typ != other.typ {
            return false;
        }
        // we do the position before the chromosome, as we wont compare mismatches between chromosomes
        if self.position != other.position {
            return false;
        }
        if self.chromosome != other.chromosome {
            return false;
        }
        // each of the things needs to fit
        for (n_self, n_other) in self.alternative.iter().zip(other.alternative.iter()) {
            if n_self != n_other {
                return false;
            }
        }
        for (n_self, n_other) in self.reference.iter().zip(other.reference.iter()) {
            if n_self != n_other {
                return false;
            }
        }

        //we dont really care about the quality being the same, so we ignore that and return true if we got to here
        return true;
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

// our partial eq already satisfies this
impl Eq for Mismatch {}

impl Mismatch {
    /// Build the INFO string for VCF.
    pub fn to_vcf_info(
        &self,
        rp_count: Option<usize>,
        overlap: Option<bool>,
        agree: Option<bool>,
        mapq_minmax: Option<(u8, u8)>,
        frag_count: Option<usize>,
        nm_minmax: Option<(u8, u8)>,
        bq: Option<u8>,
        pos_in_read: Option<(usize, usize)>,
        rp_names_preview: Option<String>,
        strand_counts: Option<(usize, usize)>, // (fwd, rev)
        fraglen_minmax: Option<(u32, u32)>,
    ) -> String {
        let mut info = String::new();

        // OVERLAP as a flag
        if overlap.unwrap_or(false) {
            info.push_str(";OVERLAP");
        }

        // AGREE
        match agree {
            Some(val) => info.push_str(&format!(";AGREE={}", val)),
            None => info.push_str(";AGREE=NA"),
        }

        // RP_COUNT
        if let Some(n) = rp_count {
            info.push_str(&format!(";RP_COUNT={}", n));
        }

        // FRAG_COUNT
        if let Some(n) = frag_count {
            info.push_str(&format!(";FRAG_COUNT={}", n));
        }

        // MAPQ min..max (collapse if equal)
        if let Some((lo, hi)) = mapq_minmax {
            if lo == hi {
                info.push_str(&format!(";MAPQ={}", lo));
            } else {
                info.push_str(&format!(";MAPQ={}-{}", lo, hi));
            }
        }

        // NM min..max (collapse if equal)
        if let Some((lo, hi)) = nm_minmax {
            if lo == hi {
                info.push_str(&format!(";NM={}", lo));
            } else {
                info.push_str(&format!(";NM={}-{}", lo, hi));
            }
        }

        // BQ (aggregate)
        if let Some(b) = bq {
            info.push_str(&format!(";BQ={}", b));
        }

        // FRAG_LEN min..max (collapse if equal)
        if let Some((lo, hi)) = fraglen_minmax {
            if lo == hi {
                info.push_str(&format!(";FRAG_LEN={}", lo));
            } else {
                info.push_str(&format!(";FRAG_LEN={}-{}", lo, hi));
            }
        }

        // STRAND counts
        if let Some((fwd, rev)) = strand_counts {
            info.push_str(&format!(";STRAND=FWD={},REV={}", fwd, rev));
        }

        // POS_IN_READ min..max (collapse if equal)
        if let Some((lo, hi)) = pos_in_read {
            if lo == hi {
                info.push_str(&format!(";POS_IN_READ={}", lo));
            } else {
                info.push_str(&format!(";POS_IN_READ={}-{}", lo, hi));
            }
        }

        // RP names (preview)
        if let Some(names) = rp_names_preview {
            if !names.is_empty() {
                info.push_str(&format!(";RP={}", names));
            }
        }

        if info.starts_with(';') {
            info.remove(0);
        }
        info
    }
}

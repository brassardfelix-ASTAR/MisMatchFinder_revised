use std::{
    cmp::{max, min},
    collections::{BTreeMap, HashMap},
    ops::RangeInclusive,
};

use log::{debug, info, warn};
use rust_htslib::bam::{self, record::Aux, Read, Record};

use crate::bamreader::{
    cigar::parse_cigar_str,
    fragments::Fragment,
};

use self::filter::region::BedObject;

use crate::bamreader::mismatch::{Mismatch, MismatchType};

pub mod cigar;
pub mod filter;
pub mod fragments; // src/bamreader/fragments.rs
pub mod mismatch;

/// Normalize a Mismatch for use as a map key so fragments at the same locus coalesce
/// without altering allele strings expected by downstream filters (e.g., germline).
/// - Keep REF/ALT exactly as constructed upstream (SBS 3-mer, DBS 2-mer, indels as-is)
/// - Strip per-fragment/read fields that should NOT shard grouping
fn canonical_key(mut m: Mismatch) -> Mismatch {
    // DO NOT collapse REF/ALT here; germline filter expects SBS as 3-mer.
    // Just neutralize fields that are read/fragment-specific.
    m.rid = 0;
    m.quality = 0;
    m.fragment_length = None;
    m.read_orientation = None;
    m.region_tag = None;
    m.pos_in_read = None;
    m
}

#[derive(Debug, Default)]
pub struct MismatchMeta {
    pub rp_names: Vec<String>,
    pub mapq: Vec<u8>,
    pub nm: Vec<u8>,
    pub frag_len: Vec<u32>,
    pub strand: Vec<String>,
    pub pos_in_read: Vec<usize>,
    pub overlap: Option<bool>,
    pub agree: Option<bool>,
}

// Make the thresholds u8 to match NM getter and avoid i8/u8 clashes
pub fn find_mismatches(
    bam: &mut bam::Reader,
    white_list: &Option<BedObject>,
    fragment_length_intervals: &Vec<RangeInclusive<i64>>,
    min_edit_distance_per_read: u8,
    max_edit_distance_per_read: u8,
    min_mapping_quality: u8,
    min_avg_base_quality: f32,
    min_base_quality: u8,
    only_overlap: bool,
    strict_overlap: bool,
) -> BTreeMap<Mismatch, MismatchMeta> {
    // store all mismatches found
    let mut mismatch_store: BTreeMap<Mismatch, MismatchMeta> = BTreeMap::new();

    // map tid -> chromosome
    let mut tid_map: HashMap<i32, String> = HashMap::new();
    {
        let header = bam.header();
        for (id, chr) in header.target_names().into_iter().enumerate() {
            tid_map.insert(
                id as i32,
                String::from_utf8(chr.to_vec()).expect("Could not convert tid to String"),
            );
        }
    }

    let mut counter = 0;
    let mut last_chr = &String::from("*");
    let mut last_pos = -1;

    // simple counters
    let mut fragments_total = 0;
    let mut fragments_analysed = 0;
    let mut fragments_wrong_length = 0;
    let mut fragments_low_base_quality = 0;

    // cache mates by qname
    let mut read_cache: HashMap<String, Record> = HashMap::with_capacity(10000);

    for r in bam.records() {
        let record = match r {
            Ok(rec) => rec,
            Err(e) => {
                warn!("Failed to parse BAM record: {}", e);
                continue;
            }
        };

        let qname = std::str::from_utf8(record.qname()).unwrap().to_owned();

        counter += 1;
        if counter % 500000 == 0 {
            info!(
                "Read through {counter} reads - last position: {}:{}",
                last_chr, last_pos
            );
        }

        if !record.is_paired() {
            // skip low quality / non-primary / unmapped
            if record.is_secondary()
                || record.is_supplementary()
                || record.is_unmapped()
                || record.is_duplicate()
                || record.is_quality_check_failed()
                || record.mapq() < min_mapping_quality
            {
                continue;
            }

            debug!("Working on single end read: {}", qname);

            // XA tag indicates alt mappings
            if let Ok(_) = record.aux(b"XA") {
                continue;
            }

            fragments_total += 1;

            let chrom = tid_map.get(&record.tid()).unwrap();
            last_chr = chrom;
            last_pos = record.pos();

            let edit_dist: u8 = get_edit_distance(&record);
            if edit_dist >= min_edit_distance_per_read && edit_dist <= max_edit_distance_per_read {
                // approximate fragment size by read length
                let frag_size = record.seq_len() as i64;
                let mut skip = true;
                for ivl in fragment_length_intervals {
                    if ivl.contains(&frag_size) { skip = false; }
                }
                if skip {
                    debug!("Discarded read due to wrong fragment size");
                    fragments_wrong_length += 1;
                    continue;
                }

                // average base qual
                if Fragment::average(record.qual()) < min_avg_base_quality {
                    fragments_low_base_quality += 1;
                    continue;
                }

                // cigar -> gapped read (no insertions)
                let read = parse_cigar_str(&record, chrom);

                // whitelist check
                let analyse = match white_list {
                    Some(wl) => wl.has_overlap(chrom, read.start() as usize, read.end() as usize),
                    None => true,
                };

                if analyse {
                    fragments_analysed += 1;

                    let frag = Fragment::make_se_fragment(read, chrom);
                    let mismatches = frag.get_mismatches(min_base_quality);

                    debug!("Found {} mismatches in fragment", mismatches.len());

                    for mm in mismatches {
                        // Coalesce by canonical key (alleles unchanged to keep SBS 3-mer)
                        let key = canonical_key(mm.clone());
                        let entry = mismatch_store.entry(key).or_default();

                        // avoid accidental dup names
                        if !entry.rp_names.contains(&qname) {
                            entry.rp_names.push(qname.clone());
                        }
                        entry.mapq.push(record.mapq());
                        entry.nm.push(get_edit_distance(&record));
                        entry.frag_len.push(record.seq_len() as u32);
                        if let Some(strand) = frag.strand() {
                            entry.strand.push(strand.to_string());
                        }
                        if let Some(pos) = mm.pos_in_read {
                            entry.pos_in_read.push(pos);
                        }
                    }
                } else {
                    debug!("Discarded read due to whitelist check");
                }
            }
        } else if read_cache.contains_key(&qname) {
            // mate available: process PE fragment

            // skip non-primary before pulling mate
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }

            let mate = read_cache.remove(&qname).unwrap();

            debug!("Working on paired end read: {}", qname);

            // fragment-level quality checks
            if !(record.is_unmapped()
                || mate.is_unmapped()
                || record.is_duplicate()
                || mate.is_duplicate()
                || record.is_quality_check_failed()
                || mate.is_quality_check_failed()
                || record.tid() != record.mtid()
                || (record.mapq() + mate.mapq()) / 2 < min_mapping_quality)
            {
                let chrom = tid_map.get(&record.tid()).unwrap();
                last_chr = chrom;
                last_pos = record.pos();

                if let Ok(_) = record.aux(b"XA") { continue; }
                if let Ok(_) = mate.aux(b"XA")   { continue; }

                fragments_total += 1;

                let read1_edit: u8 = get_edit_distance(&record);
                let read2_edit: u8 = get_edit_distance(&mate);

                if (read1_edit >= min_edit_distance_per_read && read1_edit <= max_edit_distance_per_read)
                    || (read2_edit >= min_edit_distance_per_read && read2_edit <= max_edit_distance_per_read)
                {
                    let frag_size = record.insert_size().abs();
                    let mut skip = true;
                    for ivl in fragment_length_intervals {
                        if ivl.contains(&frag_size) { skip = false; }
                    }
                    if skip {
                        debug!("Discarded read-pair due to wrong fragment size");
                        fragments_wrong_length += 1;
                        continue;
                    }

                    if (Fragment::average(record.qual()) + Fragment::average(mate.qual())) / 2.0
                        < min_avg_base_quality
                    {
                        debug!("Discarded read-pair due to low average base quality");
                        fragments_low_base_quality += 1;
                        continue;
                    }

                    let read1 = parse_cigar_str(&record, chrom);
                    let read2 = parse_cigar_str(&mate, chrom);

                    let analyse = match white_list {
                        Some(wl) => wl.has_overlap(
                            chrom,
                            min(read1.start() as usize, read2.start() as usize),
                            max(read1.end() as usize, read2.end() as usize),
                        ),
                        None => true,
                    };

                    if analyse {
                        debug!("Analysing read(pair)");
                        fragments_analysed += 1;

                        let frag = if read1.get_read_pos() <= read2.get_read_pos() {
                            Fragment::make_fragment(read1, read2, only_overlap, strict_overlap, chrom)
                        } else {
                            Fragment::make_fragment(read2, read1, only_overlap, strict_overlap, chrom)
                        };

                        let mismatches = match &frag {
                            Some(v) => v.get_mismatches(min_base_quality),
                            None => Vec::new(),
                        };

                        debug!("Found {} mismatches in fragment", mismatches.len());

                        for mm in mismatches {
                            // Coalesce by canonical key (alleles unchanged to keep SBS 3-mer)
                            let key = canonical_key(mm.clone());
                            let entry = mismatch_store.entry(key).or_default();

                            if !entry.rp_names.contains(&qname) {
                                entry.rp_names.push(qname.clone());
                            }
                            entry.mapq.push(record.mapq());
                            entry.nm.push(get_edit_distance(&record));
                            entry.frag_len.push(record.seq_len() as u32);

                            if let Some(v) = &frag {
                                if let Some(strand) = v.strand() {
                                    entry.strand.push(strand.to_string());
                                }
                                // Fragment-level overlap / agree
                                entry.overlap = Some(v.is_overlapping());
                                entry.agree = v.reads_agree(&mm);
                            }
                            if let Some(pos) = mm.pos_in_read {
                                entry.pos_in_read.push(pos);
                            }
                        }
                    } else {
                        debug!("Discarded read pair due to whitelist check");
                    }
                } else {
                    debug!("Discarded read pair due to wrong edit distance");
                }
            } else {
                debug!("Discarded read pair due to mapping quality filters");
            }

            debug!("Done with read {}", &qname);
        } else if !(record.is_supplementary() || record.is_secondary()) {
            // cache this read until its mate arrives
            read_cache.insert(qname, record);
        } else {
            // discard non-primary
        }
    }

    // should have no reads left in cache
    let mut break_out = 10;
    let un_paired_reads = read_cache.len();
    if un_paired_reads != 0 {
        warn!("Read cache contained unpaired read at the end of the analysis, this shouldnt happen with a well formed bam");
        for (qname, _) in read_cache {
            warn!("{qname}");
            if break_out == 0 {
                warn!("... and {} more", un_paired_reads - break_out);
                break;
            }
            break_out -= 1;
        }
    }

    info! {"Analysed {fragments_total} fragments and {fragments_wrong_length} were excluded due to wrong length, leaving {fragments_analysed} after whitelist and base quality (excluded: {fragments_low_base_quality}) check"};
    mismatch_store
}

fn get_edit_distance(read: &Record) -> u8 {
    match read.aux(b"NM") {
        Ok(Aux::U8(v)) => v as u8,
        Ok(Aux::I8(v)) => if v < 0 { 0 } else { v as u8 },
        Ok(Aux::U16(v)) => v.min(u8::MAX as u16) as u8,
        Ok(Aux::I16(v)) => if v < 0 { 0 } else { (v as u16).min(u8::MAX as u16) as u8 },
        Ok(Aux::U32(v)) => v.min(u8::MAX as u32) as u8,
        Ok(Aux::I32(v)) => if v < 0 { 0 } else { (v as u32).min(u8::MAX as u32) as u8 },
        Ok(Aux::String(s)) => s.parse::<u32>().map(|x| x.min(u8::MAX as u32) as u8).unwrap_or(0),
        _ => 0,
    }
}

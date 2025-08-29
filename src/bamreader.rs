use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, HashMap},
    ops::RangeInclusive,
};

use log::{debug, info, warn};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar},
    Read, Record,
};

#[derive(Clone, Copy, Debug, Default)]
pub struct SupportStats {
    pub count: usize,
    // keep means precise
    pub sum_edit_dist: f64,
    pub sum_frag_size: i64,
    pub sum_avg_base_qual: f64,
    pub sum_mapq: i64,
    pub overlap_count: usize,
    pub rp_agree_count: usize,
}

use crate::bamreader::{cigar::parse_cigar_str, fragment::Fragment};

use self::{filter::region::BedObject, mismatch::Mismatch};

pub mod cigar;
pub mod filter;
pub mod fragment;
pub mod mismatch;

pub fn find_mismatches(
    bam: &mut bam::Reader,
    white_list: &Option<BedObject>,
    fragment_length_intervals: &Vec<RangeInclusive<i64>>,
    min_edit_distance_per_read: i8,
    max_edit_distance_per_read: i8,
    min_mapping_quality: u8,
    min_avg_base_quality: f32,
    min_base_quality: u8,
    only_overlap: bool,
    strict_overlap: bool,
) -> BTreeMap<Mismatch, SupportStats> {
    let mut mismatch_store: BTreeMap<Mismatch, SupportStats> = BTreeMap::new();

    // build tid -> chrom name map
    let mut tid_map: HashMap<i32, String> = HashMap::new();
    {
        let header = bam.header();
        for (id, chr) in header.target_names().into_iter().enumerate() {
            tid_map.insert(
                id as i32,
                String::from_utf8(chr.to_vec()).expect("We couldnt fonvert the tid to String"),
            );
        }
    }

    let mut counter = 0;
    let mut last_chr = &String::from("*");
    let mut last_pos = -1;

    // counters
    let mut fragments_total = 0;
    let mut fragments_analysed = 0;
    let mut fragments_wrong_length = 0;
    let mut fragments_low_base_quality = 0;

    // read cache for mates
    let mut read_cache: HashMap<String, Record> = HashMap::with_capacity(10000);

    for r in bam.records() {
        let record = r.unwrap();
        let qname = std::str::from_utf8(record.qname()).unwrap().to_owned();

        counter += 1;
        if counter % 500000 == 0 {
            info!(
                "Read through {counter} reads - last position: {}:{}",
                last_chr, last_pos
            );
        }

        if !record.is_paired() {
            // filter low-quality SE
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

            // skip XA (alt mappings)
            if let Ok(_) = record.aux(b"XA") {
                continue;
            }

            fragments_total += 1;

            let chrom = tid_map.get(&record.tid()).unwrap();
            last_chr = chrom;
            last_pos = record.pos();

            let edit_dist = get_edit_distance(&record);
            if edit_dist >= min_edit_distance_per_read && edit_dist <= max_edit_distance_per_read {
                let frag_size = record.seq_len() as i64;
                let mut skip = true;
                for ivl in fragment_length_intervals {
                    if ivl.contains(&frag_size) {
                        skip = false;
                    }
                }
                if skip {
                    debug!("Discarded read due to wrong fragment size");
                    fragments_wrong_length += 1;
                    continue;
                }

                if Fragment::average(record.qual()) < min_avg_base_quality {
                    fragments_low_base_quality += 1;
                    continue;
                }

                let read = parse_cigar_str(&record, chrom);

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
                        // fragment-level stats (SE)
                        let edit = get_edit_distance(&record) as f64; // whole-read NM
                        let frag_size = record.seq_len() as i64;
                        let avg_bq = Fragment::average(record.qual()) as f64;
                        // (4) MAPQ_AVG: ignore 255 (unknown)
                        let mapq_avg = if record.mapq() != 255 {
                            record.mapq() as i64
                        } else {
                            0
                        };

                        let stats = mismatch_store.entry(mm).or_insert(SupportStats::default());
                        stats.count += 1;
                        stats.sum_edit_dist += edit;
                        stats.sum_frag_size += frag_size;
                        stats.sum_avg_base_qual += avg_bq;
                        stats.sum_mapq += mapq_avg;
                        // overlap / agreement not applicable to SE
                    }
                } else {
                    debug!("Discarded read due to whitelist check");
                }
            }
        } else if read_cache.contains_key(&qname) {
            // skip secondary/supplementary before fetching mate
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }

            let mate = read_cache.remove(&qname).unwrap();

            debug!("Working on paired end read: {}", qname);

            // basic pair quality / co-location
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

                if let Ok(_) = record.aux(b"XA") {
                    continue;
                }
                if let Ok(_) = mate.aux(b"XA") {
                    continue;
                }

                fragments_total += 1;

                let read1_edit = get_edit_distance(&record);
                let read2_edit = get_edit_distance(&mate);

                if (read1_edit >= min_edit_distance_per_read && read1_edit <= max_edit_distance_per_read)
                    || (read2_edit >= min_edit_distance_per_read && read2_edit <= max_edit_distance_per_read)
                {
                    let frag_size = record.insert_size().abs();
                    let mut skip = true;
                    for ivl in fragment_length_intervals {
                        if ivl.contains(&frag_size) {
                            skip = false;
                        }
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

                    // aligned segments (no I)
                    let r1_ns = parse_cigar_str(&record, chrom);
                    let r2_ns = parse_cigar_str(&mate, chrom);

                    // whitelist window uses the non-strict build’s bounds
                    let analyse = match white_list {
                        Some(wl) => wl.has_overlap(
                            chrom,
                            min(r1_ns.start() as usize, r2_ns.start() as usize),
                            max(r1_ns.end() as usize, r2_ns.end() as usize),
                        ),
                        None => true,
                    };

                    if analyse {
                        debug!("Analysing read(pair)");
                        fragments_analysed += 1;

                        // overlap window from the non-strict build
                        let r1s = r1_ns.start();
                        let r1e = r1_ns.end();
                        let r2s = r2_ns.start();
                        let r2e = r2_ns.end();
                        let ovl_start = std::cmp::max(r1s, r2s);
                        let ovl_end = std::cmp::min(r1e, r2e);

                        // Build fragment with the user’s flags (may be strict=false)
                        let res_ns = if r1_ns.get_read_pos() <= r2_ns.get_read_pos() {
                            Fragment::make_fragment(r1_ns, r2_ns, only_overlap, strict_overlap, chrom)
                        } else {
                            Fragment::make_fragment(r2_ns, r1_ns, only_overlap, strict_overlap, chrom)
                        };

                        let mismatches = match res_ns {
                            Some(v) => v.get_mismatches(min_base_quality),
                            None => Vec::new(),
                        };

                        // Build a second fragment with strict_overlap=true just to know which sites have mate agreement.
                        // Parse again to avoid ownership moves.
                        let r1_st = parse_cigar_str(&record, chrom);
                        let r2_st = parse_cigar_str(&mate, chrom);
                        let res_strict = if r1_st.get_read_pos() <= r2_st.get_read_pos() {
                            Fragment::make_fragment(r1_st, r2_st, only_overlap, true, chrom)
                        } else {
                            Fragment::make_fragment(r2_st, r1_st, only_overlap, true, chrom)
                        };
                        let strict_set: BTreeSet<Mismatch> = match res_strict {
                            Some(v) => v.get_mismatches(min_base_quality).into_iter().collect(),
                            None => BTreeSet::new(),
                        };

                        debug!("Found {} mismatches in fragment", mismatches.len());

                        for mm in mismatches {
                            // ----- site-in-overlap (indel/DBS-aware) -----
                            // mm.position is 1-based VCF POS (anchor base).
                            // Convert to 0-based anchor coordinate.
                            let anchor0 = mm.position as i64 - 1;

                            let ref_len = mm.reference.len() as i64;
                            let alt_len = mm.alternative.len() as i64;
                            let del_len = ref_len - alt_len;

                            // (2) DBS overlap gating:
                            //   - DEL (ref>alt): deleted span AFTER anchor -> [anchor0+1, anchor0+del_len+1)
                            //   - SUBSTITUTION with len>1 (DBS): require full ref span -> [anchor0, anchor0+ref_len)
                            //   - SNV or INS: gate by anchor base -> [anchor0, anchor0+1)
                            let (site_lo, site_hi) = if del_len > 0 {
                                (anchor0 + 1, anchor0 + del_len + 1)
                            } else if ref_len == alt_len && ref_len > 1 {
                                (anchor0, anchor0 + ref_len)
                            } else {
                                (anchor0, anchor0 + 1)
                            };

                            // Is the (site/span) contained within the mates' overlap window?
                            let site_in_overlap = site_lo >= ovl_start && site_hi <= ovl_end;

                            // Enforce --only_overlaps at the SITE level
                            if only_overlap && !site_in_overlap {
                                continue;
                            }

                            // ---------- fragment-level stats ----------
                            // EDITDIST_AVG = mean of whole-read NM tags
                            let edit = (read1_edit as f64 + read2_edit as f64) / 2.0;

                            let frag_size_i = frag_size as i64;
                            let avg_bq = ((Fragment::average(record.qual())
                                + Fragment::average(mate.qual())) / 2.0) as f64;

                            // (4) MAPQ_AVG: ignore 255s
                            let mut mapq_sum: i64 = 0;
                            let mut mapq_n: i64 = 0;
                            let m1 = record.mapq();
                            let m2 = mate.mapq();
                            if m1 != 255 {
                                mapq_sum += m1 as i64;
                                mapq_n += 1;
                            }
                            if m2 != 255 {
                                mapq_sum += m2 as i64;
                                mapq_n += 1;
                            }
                            let mapq_avg = if mapq_n > 0 { mapq_sum / mapq_n } else { 0 };

                            let stats = mismatch_store
                                .entry(mm.clone())
                                .or_insert(SupportStats::default());
                            stats.count += 1;
                            stats.sum_edit_dist += edit;
                            stats.sum_frag_size += frag_size_i;
                            stats.sum_avg_base_qual += avg_bq;
                            stats.sum_mapq += mapq_avg;

                            if site_in_overlap {
                                stats.overlap_count += 1;

                                // (8) Count pair-agreement even when strict_overlap was OFF for the main call:
                                // A mismatch "agrees" if it also appears in the strict build (same fragment),
                                // AND the site is in the pair-overlap (by definition of RP_AGREE_COUNT field).
                                if strict_set.contains(&mm) {
                                    stats.rp_agree_count += 1;
                                }
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
            read_cache.insert(qname, record);
        } else {
            // discard extra alignments
        }
    }

    // final cache check
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

// ---- NM helper with robust fallback ----
// Prefer NM tag if present; otherwise compute a fallback.
// (7) Avoid double-counting: if CIGAR::Diff was seen, do NOT also count MD mismatches.
fn get_edit_distance(read: &Record) -> i8 {
    // Prefer NM tag when available
    if let Ok(aux) = read.aux(b"NM") {
        let v = match aux {
            Aux::I8(n) => n as i32,
            Aux::I16(n) => n as i32,
            Aux::I32(n) => n as i32,
            Aux::U8(n) => n as i32,
            Aux::U16(n) => n as i32,
            Aux::U32(n) => n as i32,
            _ => 0,
        };
        return if v > i8::MAX as i32 { i8::MAX } else { v as i8 };
    }

    // Fallback from CIGAR (+MD only if Diff not present)
    let mut nm_fallback: i32 = 0;
    let mut saw_diff = false;
    for op in read.cigar().iter() {
        match *op {
            Cigar::Ins(len) | Cigar::Del(len) => nm_fallback += len as i32,
            Cigar::Diff(len) => {
                nm_fallback += len as i32;
                saw_diff = true;
            }
            _ => {}
        }
    }
    if !saw_diff {
        if let Ok(Aux::String(md_str)) = read.aux(b"MD") {
            // Count mismatching letters outside deletion blocks
            let mut in_del = false;
            for ch in md_str.chars() {
                match ch {
                    '^' => in_del = true,         // deletion (already counted via CIGAR)
                    '0'..='9' => in_del = false,  // run of matches
                    'A' | 'C' | 'G' | 'T' | 'N' => {
                        if !in_del {
                            nm_fallback += 1;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    if nm_fallback > i8::MAX as i32 { i8::MAX } else { nm_fallback as i8 }
}


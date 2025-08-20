use std::{
    cmp::{max, min},
    collections::{BTreeMap, HashMap},
    ops::RangeInclusive,
};

use log::{debug, info, warn};

use rust_htslib::bam::{self, record::Aux, Read, Record};

use crate::bamreader::{
    cigar::parse_cigar_str,
    fragments::Fragment, // <-- use `fragments` module
};

use self::{
    filter::region::BedObject,
    mismatch::Mismatch,
};

pub mod cigar;
pub mod filter;
pub mod fragments; // <-- ensure your file is `src/bamreader/fragments.rs`
pub mod mismatch;

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
    //store all mismatches found and how often they were found (we also know its going to be a big hash)
    let mut mismatch_store: BTreeMap<Mismatch, MismatchMeta> = BTreeMap::new();

    //get the header for the mapping of tid to chromosome name
    let mut tid_map: HashMap<i32, String> = HashMap::new();
    {
        let header = bam.header();
        let seqnames = header.target_names().into_iter();

        for (id, chr) in seqnames.enumerate() {
            tid_map.insert(
                id as i32,
                String::from_utf8(chr.to_vec()).expect("Could not convert tid to String"),
            );
        }
    }

    let mut counter = 0;
    let mut last_chr = &String::from("*");
    let mut last_pos = -1;

    // count stats
    let mut fragments_total = 0;
    let mut fragments_analysed = 0;
    let mut fragments_wrong_length = 0;
    let mut fragments_low_base_quality = 0;

    //we create a read cache and set the initial capacity to 10k to reduce the reallocation operations
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
            // we skipp anything that isnt high quality
            if record.is_secondary()
                || record.is_supplementary()
                || record.is_unmapped()
                || record.is_duplicate()
                || record.is_quality_check_failed()
                || record.mapq() < min_mapping_quality
            {
                continue;
            }

            debug!(
                "Working on single end read: {} ",
                std::str::from_utf8(record.qname()).unwrap()
            );

            // if the read has any additional mapping locations, we cant really trust the alignment as much
            if let Ok(_) = record.aux(b"XA") {
                continue;
            }

            // from here on we think this is a proper fragment that could be in the analysis
            fragments_total += 1;

            // get the chromosome the record is on
            let chrom = tid_map.get(&record.tid()).unwrap();
            last_chr = chrom;
            last_pos = record.pos();

            let edit_dist: u8 = get_edit_distance(&record);
            if edit_dist >= min_edit_distance_per_read && edit_dist <= max_edit_distance_per_read {
                // we cant really do a fragment size check here, so we have to check the read length instead
                let frag_size = record.seq_len() as i64;
                // but we can only really estimate it
                let mut skip = true;
                for ivl in fragment_length_intervals {
                    if ivl.contains(&frag_size) {
                        //we are good here
                        skip = false;
                    }
                }
                if skip {
                    debug!("Discarded read due to wrong fragment size");
                    fragments_wrong_length += 1;
                    continue;
                }

                //then we check for the average base quality of the read
                if Fragment::average(record.qual()) < min_avg_base_quality {
                    fragments_low_base_quality += 1;
                    continue;
                }

                // get only the aligned part of the read, without insertions
                let read = parse_cigar_str(&record, chrom);

                //check if the read is in the whitelist or if no white list was supplied
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
                        let entry = mismatch_store.entry(mm.clone()).or_default();

                        entry.rp_names.push(qname.clone());
                        entry.mapq.push(record.mapq());
                        entry.nm.push(get_edit_distance(&record));
                        entry.frag_len.push(record.seq_len() as u32); // or insert_size().abs() for PE
                        if let Some(strand) = frag.strand() {
                            entry.strand.push(strand.to_string()); // char -> String
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
            // we have to skip out of this if the read is a secondary (before we get the read from the cache)
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }

            // get and delete the record from the cache, because we dont want to bloat the storage
            let mate = read_cache.remove(&qname).unwrap();

            debug!(
                "Working on paired end read: {} ",
                std::str::from_utf8(record.qname()).unwrap()
            );

            // we check all the quality of the read AND a few for the mate, because we need to be sure we have a proper pair
            if !(record.is_unmapped()
                || mate.is_unmapped()
                || record.is_duplicate()
                || mate.is_duplicate()
                || record.is_quality_check_failed()
                || mate.is_quality_check_failed()
                || record.tid() != record.mtid()
                // we use the average mapping quality of the read as the mapping quality of the fragment
                || (record.mapq() + mate.mapq()) / 2 < min_mapping_quality)
            {
                // get the chromosome the record is on (because we know they are both on the same)
                let chrom = tid_map.get(&record.tid()).unwrap();
                last_chr = chrom;
                last_pos = record.pos();

                //if the read has any additional mapping locations, we cant really trust the alignment as much
                if let Ok(_) = record.aux(b"XA") {
                    continue;
                }
                if let Ok(_) = mate.aux(b"XA") {
                    continue;
                }

                // from here on we think this is a proper fragment that could be in the analysis
                fragments_total += 1;

                // we check if the reads actually have any changes (edit distance), this contains both mismatches and indels
                let read1_edit: u8 = get_edit_distance(&record);
                let read2_edit: u8 = get_edit_distance(&mate);

                //both reads need to be within the edit distance requirements
                if (read1_edit >= min_edit_distance_per_read
                    && read1_edit <= max_edit_distance_per_read)
                    || (read2_edit >= min_edit_distance_per_read
                        && read2_edit <= max_edit_distance_per_read)
                {
                    //now we check if the fragment has the right size
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

                    //then we check for the average base quality of the reads
                    if (Fragment::average(record.qual()) + Fragment::average(mate.qual())) / 2.
                        < min_avg_base_quality
                    {
                        debug!("Discarded read-pair due to low average base quality");
                        fragments_low_base_quality += 1;
                        continue;
                    }

                    // get only the aligned part of the read, without insertions
                    let read1 = parse_cigar_str(&record, chrom);
                    let read2 = parse_cigar_str(&mate, chrom);

                    //check if the read is in the whitelist or if no white list was supplied
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
                            let entry = mismatch_store.entry(mm.clone()).or_default();

                            entry.rp_names.push(qname.clone());
                            entry.mapq.push(record.mapq());
                            entry.nm.push(get_edit_distance(&record));
                            entry.frag_len.push(record.seq_len() as u32); // or insert_size().abs() for PE
                            if let Some(v) = &frag {
                                if let Some(strand) = v.strand() {
                                    entry.strand.push(strand.to_string()); // char -> String
                                }
                                entry.overlap = Some(v.is_overlapping());
                                entry.agree = v.reads_agree(&mm); // Option<bool>
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
        // we dont check if the read is mapped here, so we get mapped and unmapped pairs as well
        } else if !(record.is_supplementary() || record.is_secondary()) {
            read_cache.insert(qname, record);
        } else {
            // we do discard any additional reads, which are not primary alignment or not mapped to the same chromosome
        }
    }

    // if we have done everything correctly, we should have no reads in the cache anymore
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


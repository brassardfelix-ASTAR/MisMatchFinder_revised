use std::io::Write;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::path::PathBuf;

use bgzip::BGZFWriter;
use crate::bamreader::MismatchMeta;
use crate::bamreader::mismatch::{Mismatch, MismatchType};

/// Return (0-based VCF POS, REF, ALT) with SBS context collapsed to the middle base for printing.
/// Falls back to original alleles for DBS/indels/edge cases.
fn vcf_pos_and_alleles(m: &Mismatch) -> (i64, Vec<u8>, Vec<u8>) {
    match m.typ {
        MismatchType::SBS if m.reference.len() == 3 && m.alternative.len() == 3 => {
            // In fragments.rs we set SBS position to the CENTER base (1-based).
            // VCF POS in print is 1-based; the writer here expects to print m.position as given.
            // But many downstream tools use 0-based internally; keep conversion local as needed.
            // We return 0-based POS for internal use below.
            let pos0 = (m.position as i64) - 1;
            let r = vec![m.reference[1]];
            let a = vec![m.alternative[1]];
            (pos0, r, a)
        }
        _ => {
            // Keep original representation; convert POS to 0-based for internal indexing if needed.
            let pos0 = (m.position as i64) - 1;
            (pos0, m.reference.clone(), m.alternative.clone())
        }
    }
}

// normal string literal so \t and \n are interpreted as tabs/newlines
pub const VCF_HEADER: &str = "\
##fileformat=VCFv4.2
##source=MisMatchFinder
##FILTER=<ID=PASS,Description=\"Mismatch of high quality\">
##INFO=<ID=RP,Number=.,Type=String,Description=\"Read-pair IDs supporting the mismatch (unique per pair)\">
##INFO=<ID=RP_COUNT,Number=1,Type=Integer,Description=\"Number of read pairs supporting the mismatch (unique)\">
##INFO=<ID=OVERLAP,Number=1,Type=Flag,Description=\"Mismatch in overlapping region of read pair\">
##INFO=<ID=AGREE,Number=1,Type=String,Description=\"Do both reads agree on mismatch (true/false/NA)\">
##INFO=<ID=MAPQ,Number=.,Type=Integer,Description=\"Mapping qualities of supporting reads (or range)\">
##INFO=<ID=FRAG_COUNT,Number=1,Type=Integer,Description=\"Number of fragments supporting mismatch (keep MULTI)\">
##INFO=<ID=NM,Number=.,Type=Integer,Description=\"Edit distances (NM tag) of supporting reads (or range)\">
##INFO=<ID=BQ,Number=.,Type=Integer,Description=\"Base qualities at mismatch position (or aggregate)\">
##INFO=<ID=FRAG_LEN,Number=.,Type=Integer,Description=\"Fragment lengths (or range)\">
##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand orientation of supporting reads\">
##INFO=<ID=POS_IN_READ,Number=.,Type=Integer,Description=\"Position(s) of mismatch within read(s)\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

pub fn write_vcf(
    mismatches: &BTreeMap<Mismatch, MismatchMeta>,
    file: PathBuf,
    _include_small_vars: bool,
    somatic: bool,
) -> std::io::Result<()> {
    let mut vcf_fh = File::create(file)?;
    let mut writer = BGZFWriter::new(&mut vcf_fh, flate2::Compression::default());

    // write header
    writer.write_all(VCF_HEADER.as_bytes())?;

    for (mm, meta) in mismatches {
        // RP preview (already capped during aggregation; add ellipsis if total > preview)
        let mut rp_preview = if meta.rp_names.is_empty() {
            String::new()
        } else {
            meta.rp_names.join(",")
        };
        if meta.rp_set.len() > meta.rp_names.len() {
            if !rp_preview.is_empty() {
                rp_preview.push_str(",...");
            } else {
                rp_preview.push_str("...");
            }
        }

        // MAPQ, NM min..max
        let (min_mapq, max_mapq) = (
            meta.mapq.iter().min().cloned().unwrap_or(0),
            meta.mapq.iter().max().cloned().unwrap_or(0),
        );
        let (min_nm, max_nm) = (
            meta.nm.iter().min().cloned().unwrap_or(0),
            meta.nm.iter().max().cloned().unwrap_or(0),
        );

        // FRAG_LEN min..max
        let (min_fraglen, max_fraglen) = (
            meta.frag_len.iter().min().cloned().unwrap_or(0),
            meta.frag_len.iter().max().cloned().unwrap_or(0),
        );

        // STRAND counts (FWD='+', REV='-')
        let (mut fwd, mut rev) = (0usize, 0usize);
        for s in &meta.strand {
            if s == "+" { fwd += 1; }
            else if s == "-" { rev += 1; }
        }

        // POS_IN_READ min..max
        let pos_in_read_minmax = if meta.pos_in_read.is_empty() {
            None
        } else {
            let lo = *meta.pos_in_read.iter().min().unwrap();
            let hi = *meta.pos_in_read.iter().max().unwrap();
            Some((lo, hi))
        };

        // Build INFO via your existing helper on Mismatch
        let rp_count_unique = meta.rp_set.len();

        let info = mm.to_vcf_info(
            Some(rp_count_unique),                 // RP_COUNT (unique)
            meta.overlap,                          // OVERLAP flag
            meta.agree,                            // AGREE
            Some((min_mapq, max_mapq)),            // MAPQ range
            Some(meta.frag_len.len()),             // FRAG_COUNT (kept as number of fragments seen)
            Some((min_nm, max_nm)),                // NM range
            Some(mm.quality),                      // BQ (aggregate carried on the Mismatch)
            pos_in_read_minmax,                    // POS_IN_READ min..max
            if rp_preview.is_empty() { None } else { Some(rp_preview) }, // RP preview
            Some((fwd, rev)),                      // STRAND counts
            Some((min_fraglen, max_fraglen)),      // FRAG_LEN range
        );

        // Collapse SBS context for REF/ALT when printing, but keep POS as 1-based in VCF
        let (_pos0, ref_allele_bytes, alt_allele_bytes) = vcf_pos_and_alleles(mm);
        let ref_str = std::str::from_utf8(&ref_allele_bytes).unwrap();
        let alt_str = std::str::from_utf8(&alt_allele_bytes).unwrap();

        // If after collapsing REF == ALT, skip (sanity)
        if ref_str == alt_str {
            continue;
        }

        // Write line (POS printed as m.position which is 1-based already)
        let line = format!(
            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}{}",
            mm.chromosome,
            mm.position,
            ref_str,
            alt_str,
            mm.quality,
            info,
            if somatic { ";SOMATIC\n" } else { "\n" }
        );

        writer.write_all(line.as_bytes())?;
    }

    writer.close()?;
    Ok(())
}

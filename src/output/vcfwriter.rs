use std::io::Write;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::PathBuf;

use bgzip::BGZFWriter;
use crate::bamreader::MismatchMeta;
use std::collections::HashSet;
use crate::bamreader::mismatch::{Mismatch, MismatchType};

fn vcf_pos_and_alleles_pos1(m: &Mismatch) -> (u32, Vec<u8>, Vec<u8>) {
    match m.typ {
        MismatchType::SBS if m.reference.len() == 3 && m.alternative.len() == 3 => {
            let pos1 = m.position;           // CENTER already
            let r = vec![m.reference[1]];
            let a = vec![m.alternative[1]];
            (pos1, r, a)
        }
        _ => (m.position, m.reference.clone(), m.alternative.clone()),
    }
}

// normal string literal so \t and \n are interpreted as tabs/newlines
pub const VCF_HEADER: &str = "\
##fileformat=VCFv4.2
##source=MisMatchFinder
##FILTER=<ID=PASS,Description=\"Mismatch of high quality\">
##INFO=<ID=RP,Number=.,Type=String,Description=\"Names of read pairs supporting the mismatch (e.g., read1/read2)\">
##INFO=<ID=RP_COUNT,Number=1,Type=Integer,Description=\"Number of read pairs supporting the mismatch\">
##INFO=<ID=OVERLAP,Number=1,Type=Flag,Description=\"Mismatch in overlapping region of read pair\">
##INFO=<ID=AGREE,Number=1,Type=String,Description=\"Do both reads agree on mismatch (true/false/NA)\">
##INFO=<ID=MAPQ,Number=.,Type=Integer,Description=\"Mapping qualities of supporting reads (or average)\">
##INFO=<ID=FRAG_COUNT,Number=1,Type=Integer,Description=\"Number of fragments supporting mismatch (keep MULTI)\">
##INFO=<ID=NM,Number=.,Type=Integer,Description=\"Edit distances (NM tag) of supporting reads (or average)\">
##INFO=<ID=BQ,Number=.,Type=Integer,Description=\"Base qualities at mismatch position (or average)\">
##INFO=<ID=FRAG_LEN,Number=.,Type=Integer,Description=\"Fragment lengths (or average)\">
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

    // write header (tabs/newlines are real because VCF_HEADER is a normal string)
    writer.write_all(VCF_HEADER.as_bytes())?;

    for (mm, meta) in mismatches {
        // RP preview (up to 10)
        let rp_preview = if meta.rp_names.is_empty() {
            String::new()
        } else {
            let mut s = meta
                .rp_names
                .iter()
                .take(10)
                .cloned()
                .collect::<Vec<String>>()
                .join(",");
            if meta.rp_names.len() > 10 {
                s.push_str(",...");
            }
            s
        };

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


        // Unique RP count (dedupe by read-pair name)
        let rp_count_unique = {
            let mut set: HashSet<&str> = HashSet::new();
            for n in &meta.rp_names {
                set.insert(n.as_str());
            }
            set.len()
        };

        let info = mm.to_vcf_info(
            Some(rp_count_unique),               // RP_COUNT
            meta.overlap,                            // OVERLAP flag
            meta.agree,                              // AGREE
            Some((min_mapq, max_mapq)),              // MAPQ min..max
            Some(meta.frag_len.len()),               // FRAG_COUNT
            Some((min_nm, max_nm)),                  // NM min..max
            Some(mm.quality),                        // BQ (aggregate)
            pos_in_read_minmax,                      // POS_IN_READ min..max
            if rp_preview.is_empty() { None } else { Some(rp_preview) }, // RP preview
            Some((fwd, rev)),                        // STRAND counts
            Some((min_fraglen, max_fraglen)),        // FRAG_LEN min..max
        );

        // Collapse SBS REF/ALT and adjust POS when appropriate
        let (pos1, ref_allele, alt_allele) = vcf_pos_and_alleles_pos1(mm);
        let ref_str = std::str::from_utf8(&ref_allele).unwrap();
        let alt_str = std::str::from_utf8(&alt_allele).unwrap();

        let line = format!(
            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}{}",
            mm.chromosome,
            pos1,
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

use crate::bamreader::mismatch::Mismatch;
use crate::bamreader::SupportStats;
use bgzip::BGZFWriter;
use std::io::Write;
use std::{collections::BTreeMap, fs::File, path::PathBuf};

const VCF_HEADER: &[u8] = "##fileformat=VCFv4.2
##source==MisMatchFinder
##FILTER=<ID=PASS,Description=\"Mismatch of high quality\">
##INFO=<ID=MULTI,Number=1,Type=Integer,Description=\"(Deprecated) Number of reads supporting the mismatch\">
##INFO=<ID=READSUPPORT,Number=1,Type=Integer,Description=\"Number of fragments supporting the mismatch\">
##INFO=<ID=EDITDIST_AVG,Number=1,Type=Float,Description=\"Average edit distance (NM) across supporting fragments\">
##INFO=<ID=FRAGSIZE_MEAN,Number=1,Type=Float,Description=\"Mean fragment size across supporting fragments\">
##INFO=<ID=AVG_BQ,Number=1,Type=Float,Description=\"Mean average base quality across supporting fragments\">
##INFO=<ID=MAPQ_AVG,Number=1,Type=Float,Description=\"Average mapping quality across supporting fragments\">
##INFO=<ID=OVERLAP_COUNT,Number=1,Type=Integer,Description=\"Number of supporting fragments where the site is in read overlap\">
##INFO=<ID=RP_AGREE_COUNT,Number=1,Type=Integer,Description=\"Number of supporting fragments where overlapping read pairs agree\">
##INFO=<ID=SOMATIC,Number=1,Type=String,Description=\"This variant is somatic after germline check\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    .as_bytes();

pub fn write_vcf(
    mismatches: &BTreeMap<Mismatch, SupportStats>,
    file: PathBuf,
    include_small_vars: bool,
    somatic: bool,
) -> std::io::Result<()> {
    let mut vcf_fh = File::create(file).expect("Could not open file to write vcf");
    let mut writer = BGZFWriter::new(&mut vcf_fh, flate2::Compression::default());

    // write header
    writer.write_all(VCF_HEADER)?;

    for (mm, stats) in mismatches {
        let read_support = stats.count as i64;
        let multi = read_support; // backward compatible
        let edit_avg = if stats.count > 0 { (stats.sum_edit_dist as f64) / (stats.count as f64) } else { 0.0 };
        let frag_mean = if stats.count > 0 { (stats.sum_frag_size as f64) / (stats.count as f64) } else { 0.0 };
        let avg_bq = if stats.count > 0 { stats.sum_avg_base_qual / (stats.count as f64) } else { 0.0 };
        let mapq_avg = if stats.count > 0 { (stats.sum_mapq as f64) / (stats.count as f64) } else { 0.0 };

        let info_tail = format!(
            "MULTI={};READSUPPORT={};EDITDIST_AVG={:.2};FRAGSIZE_MEAN={:.2};AVG_BQ={:.2};MAPQ_AVG={:.2};OVERLAP_COUNT={};RP_AGREE_COUNT={}",
            multi, read_support, edit_avg, frag_mean, avg_bq, mapq_avg, stats.overlap_count, stats.rp_agree_count
        );

        let line = match mm.typ {
            crate::bamreader::mismatch::MismatchType::SBS => {
                if include_small_vars {
                    let mut line = format!(
                        "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}",
                        mm.chromosome,
                        // original code appears to use 1-based or as stored; keep as-is
                        mm.position,
                        std::str::from_utf8(&mm.reference[1..2]).unwrap(),
                        std::str::from_utf8(&mm.alternative[1..2]).unwrap(),
                        mm.quality,
                        info_tail
                    );
                    if somatic {
                        line.push_str(";SOMATIC\n");
                    } else {
                        line.push('\n');
                    }
                    Some(line)
                } else {
                    None
                }
            }
            crate::bamreader::mismatch::MismatchType::DBS => {
                if include_small_vars {
                    let mut line = format!(
                        "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}",
                        mm.chromosome,
                        mm.position,
                        std::str::from_utf8(&mm.reference[..]).unwrap(),
                        std::str::from_utf8(&mm.alternative[..]).unwrap(),
                        mm.quality,
                        info_tail
                    );
                    if somatic {
                        line.push_str(";SOMATIC\n");
                    } else {
                        line.push('\n');
                    }
                    Some(line)
                } else {
                    None
                }
            }
            crate::bamreader::mismatch::MismatchType::INS
            | crate::bamreader::mismatch::MismatchType::DEL => {
                let mut line = format!(
                    "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}",
                    mm.chromosome,
                    mm.position,
                    std::str::from_utf8(&mm.reference[..]).unwrap(),
                    std::str::from_utf8(&mm.alternative[..]).unwrap(),
                    mm.quality,
                    info_tail
                );
                if somatic {
                    line.push_str(";SOMATIC\n");
                } else {
                    line.push('\n');
                }
                Some(line)
            }
        };
        if let Some(l) = line {
            writer.write_all(l.as_bytes())?;
        }
    }

    writer.close()?;
    Ok(())
}

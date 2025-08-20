use std::{
    fs::{self, File},
    ops::RangeInclusive,
    path::PathBuf,
    collections::BTreeMap,
};

use clap::{Parser, ValueHint};
use log::{debug, error, info};
use mismatchfinder::{bamreader, output};
use mismatchfinder::bamreader::mismatch::Mismatch;
use rust_htslib::bam;
use simple_logger::SimpleLogger;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Options {
    /// Output folder to write files to
    #[clap(short='o', long="output", value_hint = ValueHint::FilePath)]
    output_folder: PathBuf,

    /// Bed file for genomic regions to ignore in the analysis. White list bed file regions overwrite black list regions
    #[clap(long="blacklist_bed", value_hint = ValueHint::FilePath)]
    blacklist_file: Option<PathBuf>,

    /// Bed file for genomic regions to include in the analysis. White list bed file regions overwrite black list regions [default: all regions]
    #[clap(long="whitelist_bed", value_hint = ValueHint::FilePath)]
    whitelist_file: Option<PathBuf>,

    /// File to read germline information from default is an echtvar file
    #[clap(long="germline_file", value_hint = ValueHint::FilePath)]
    germline_file: Option<PathBuf>,

    /// Mimimum mapping quality for a read to be considered
    #[clap(
        short = 'q',
        long = "minimum_mapping_quality",
        value_parser,
        default_value_t = 20
    )]
    min_mapping_quality: u8,

    /// Mimimum base quality of the mismatch (BQ is summed for a readpair if reads agree)
    #[clap(
        short = 'Q',
        long = "minimum_base_quality",
        value_parser,
        default_value_t = 65
    )]
    min_base_quality: u8,

    /// Maximum mismatches we allow in the read
    #[clap(
        long = "maximum_edit_distance_per_read",
        value_parser,
        default_value_t = 15
    )]
    max_edit_distance_per_read: i8,

    /// Mimimum mismatches we require in the read
    #[clap(
        long = "minimum_edit_distance_per_read",
        value_parser,
        default_value_t = 1
    )]
    min_edit_distance_per_read: i8,

    /// Mimimum average base quality of the read
    #[clap(
        long = "minimum_average_base_quality",
        value_parser,
        default_value_t = 25.
    )]
    min_average_base_quality: f32,

    /// Length of fragments to be considered in the analysis
    #[clap(long = "fragment_length_intervals", value_parser, action = clap::ArgAction::Append, multiple_values=true, default_values(&["100-150", "250-325"]))]
    fragment_length_intervals: Vec<String>,

    /// only use the overlap of the two reads
    #[clap(long = "only_overlaps", value_parser, default_value_t = false , action = clap::ArgAction::SetTrue)]
    only_overlaps: bool,

    /// only analyse mismatches if the read pair agree (does not restrict to only overlap)
    #[clap(long = "strict_overlaps", value_parser, default_value_t = false , action = clap::ArgAction::SetTrue)]
    strict_overlap: bool,

    /// overwrite previous results
    #[clap(long, value_parser, default_value_t = false , action = clap::ArgAction::SetTrue)]
    overwrite: bool,

    /// Bams to analyse
    #[clap(value_hint = ValueHint::FilePath, action = clap::ArgAction::Append, multiple_values=true)]
    bams: Vec<PathBuf>,
}

fn main() {
    let cli = Options::parse();

    SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .env()
        .init()
        .unwrap();

    // ensure output folder exists
    fs::create_dir_all(&cli.output_folder).expect("Couldnt open output folder");

    // make this an option with lapper being default
    let lapper = true;
    //create white list object
    let bed;
    if lapper {
        bed = match cli.whitelist_file {
            Some(file) => {
                debug!(
                    "Parsing bed file {} as white list to lapper object",
                    file.display()
                );
                Some(bamreader::filter::region::BedObject::lapper_from_bed(file))
            }
            None => {
                debug!("No white list bed file detected, including all reads in analysis");
                None
            }
        };
    } else {
        bed = match cli.whitelist_file {
            Some(file) => {
                debug!(
                    "Parsing bed file {} as white list to lapper object",
                    file.display()
                );
                Some(bamreader::filter::region::BedObject::interval_tree_from_bed(file))
            }
            None => {
                debug!("No white list bed file detected, including all reads in analysis");
                None
            }
        };
    }

    let mut gnomad = match cli.germline_file {
        Some(file) => Some(
            bamreader::filter::germline::GermlineResource::load_echtvars_file(
                file.to_str().unwrap(),
            ),
        ),
        None => None,
    };

    let mut fragment_length_intervals: Vec<RangeInclusive<i64>> = Vec::new();

    for ivl in cli.fragment_length_intervals {
        let ivl: Vec<&str> = ivl.split('-').collect();
        if ivl.len() != 2 {
            panic!("We require the fragment length intervals to be in the shape of <x>-<y>")
        } else {
            let start = ivl[0].parse::<i64>().expect("Could not convert start value to numeric");
            let end = ivl[1].parse::<i64>().expect("Could not convert end value to numeric");
            fragment_length_intervals.push(start..=end);
        }
    }

    for bam_path in cli.bams {
        //prepare file names
        let base = bam_path.file_stem().unwrap();
        let mut bam = bam::Reader::from_path(&bam_path).unwrap();
        let vcf_file = cli
            .output_folder
            .join(format!("{}_bamsites.vcf.gz", base.to_str().unwrap()));

        if !vcf_file.exists() || cli.overwrite {
            // we create the file here, so that we can easier run the analysis in parallel
            File::create(&vcf_file).expect("Could not create vcf file");

            // Cast i8 CLI args to u8 for find_mismatches
            let min_ed: u8 = if cli.min_edit_distance_per_read < 0 { 0 } else { cli.min_edit_distance_per_read as u8 };
            let max_ed: u8 = if cli.max_edit_distance_per_read < 0 { 0 } else { cli.max_edit_distance_per_read as u8 };

            let mut mismatches = mismatchfinder::bamreader::find_mismatches(
                &mut bam,
                &bed,
                &fragment_length_intervals,
                min_ed,
                max_ed,
                cli.min_mapping_quality,
                cli.min_average_base_quality,
                cli.min_base_quality,
                cli.only_overlaps,
                cli.strict_overlap,
            );
            info!("Found {} mismatches ", mismatches.len());

            // If we have a germline resource, filter germline co-occurrence
            if let Some(ref mut g) = gnomad {
                // Build a counts view that the germline API expects (choose a metric; here: read-pair count)
                let mut co_counts: BTreeMap<Mismatch, usize> = mismatches
                    .iter()
                    .map(|(k, v)| (k.clone(), v.rp_names.len()))
                    .collect();

                // Filter in-place
                g.filter_germline_cooccurance(&mut co_counts);

                // Keep only variants that survived in co_counts
                mismatches.retain(|k, _| co_counts.contains_key(k));

                info!("Found {} somatic mismatches", mismatches.len());
            }

            // If your write_vcf takes (mismatches, file, include_small_vars, something_else) adjust call accordingly.
            match output::write_vcf(&mismatches, vcf_file, true, true) {
                Err(_) => error!("Could not write mismatch file"),
                Ok(_) => {}
            }
        }
    }
}

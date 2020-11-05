extern crate bincode;
extern crate bio;
extern crate bitvector;
extern crate clap;
extern crate indicatif;
extern crate itertools;
extern crate libradicl;
extern crate num_format;
extern crate pretty_env_logger;
extern crate quickersort;
extern crate rust_htslib;
extern crate serde;
extern crate sprs;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

pub mod mapping;
pub mod fragmentation;

pub const MIL: usize = 1_000_000;
pub const TMIL: usize = 10_000_000;
pub const HMIL: usize = 100_000_000;

pub const FRAG_DIST: i64 = 10;
pub const MATE_MIN_DISTANCE: i64 = 20;
pub const MATE_MAX_DISTANCE: i64 = 662;
pub const MIN_MAPQ: u8 = 30;
pub const CB_LENGTH: usize = 16;
pub const TN5_LEFT_OFFSET: i64 = 4;
pub const TN5_RIGHT_OFFSET: i64 = 5;
pub const IS_WTL_FWD: bool = true;
pub const NUM_SUPPORT_CB: u32 = 5;
pub const MIN_FEAT_COUNT: u64 = 5;
pub const WINDOW_SIZE: i64 = 500;
pub const PILEUP_THRESHOLD: u16 = 15;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("flash")
        .version("0.1.0")
        .author("Avi Srivastava, Tim Stuart, Bingjie Zhang, Rahul Satija")
        .about("A set of fast helper functions for Cut&Tag/ATAC data analysis.")
        .subcommand(
            SubCommand::with_name("filter")
                .about("A subcommand to filter BAM and generate BED.")
                .arg(
                    Arg::with_name("ibam")
                        .long("ibam")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BAM file"),
                )
                .arg(
                    Arg::with_name("obed")
                        .long("obed")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output bed file"),
                )
                .arg(
                    Arg::with_name("tenx")
                        .long("tenx")
                        .help("use tag CB from 10x generated BAM."),
                )
                .arg(
                    Arg::with_name("stats")
                        .long("stats")
                        .help("Don't write the output BED, just produce stats."),
                )
                .arg(
                    Arg::with_name("mitostr")
                        .long("mitostr")
                        .short("m")
                        .takes_value(true)
                        .default_value("chrM")
                        .help("String to identify the mitochondrial chromosome"),
                ),
        )
        .subcommand(
            SubCommand::with_name("correct")
                .about("A subcommand to sequence correct the cb sequences.")
                .arg(
                    Arg::with_name("whitelist")
                        .long("whitelist")
                        .short("w")
                        .takes_value(true)
                        .required(true)
                        .help("path to the list of known whitelist CB."),
                )
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with CB sequences."),
                ),
        )
        .subcommand(
            SubCommand::with_name("sort")
                .about("A subcommand to sort the file by a chromosome names.")
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with CB sequences."),
                ),
        )
        .subcommand(
            SubCommand::with_name("group")
                .about("A subcommand to group the file by (chr, start, end, CB)")
                .arg(
                    Arg::with_name("allcb")
                        .long("allcb")
                        .help("report all cb instead of count."),
                )
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with CB sequences."),
                ),
        )
        .subcommand(
            SubCommand::with_name("callpeak")
                .about("A subcommand to call peaks from a grouped bed file")
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the grouped BED file."),
                ),
        )
        .subcommand(
            SubCommand::with_name("count")
                .about("A subcommand to generate peak v cell count matrix")
                .arg(
                    Arg::with_name("bam")
                        .long("bam")
                        .short("b")
                        .takes_value(true)
                        .help("path to the bam file with the chromosome names."),
                )
                .arg(
                    Arg::with_name("pbed")
                        .long("pbed")
                        .short("p")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with the annotated peaks."),
                )
                .arg(
                    Arg::with_name("cbed")
                        .long("cbed")
                        .short("c")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file fragments and CB."),
                ),
        )
        .subcommand(
            SubCommand::with_name("text")
                .about("A subcommand to convert binary bed to text.")
                .arg(
                    Arg::with_name("cbtext")
                        .long("cbtext")
                        .help("writes the last column as CB sequence."),
                )
                .arg(
                    Arg::with_name("bam")
                        .long("bam")
                        .short("b")
                        .takes_value(true)
                        .help("path to the bam file with the chromosome names."),
                )
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with CB sequences."),
                ),
        )
        .subcommand(
            SubCommand::with_name("stats")
                .about("A subcommand to summary stats of the binary bed.")
                .arg(
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file."),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    if let Some(sub_m) = matches.subcommand_matches("filter") {
        filter::filter(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("correct") {
        barcode::correct(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("sort") {
        sort::sort(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("group") {
        group::dedup(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("callpeak") {
        peak::callpeak(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("count") {
        count::count(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("text") {
        text::convert(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("stats") {
        stats::stats(&sub_m)?
    }

    Ok(())
}
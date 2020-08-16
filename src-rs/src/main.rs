extern crate bincode;
extern crate bio;
extern crate clap;
extern crate indicatif;
extern crate itertools;
extern crate libradicl;
extern crate num_format;
extern crate pretty_env_logger;
extern crate rust_htslib;
extern crate serde;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod barcode;
mod count;
mod filter;
mod fragments;
mod fstats;
mod group;
mod peak;
mod sort;
mod text;

pub const MIL: usize = 1_000_000;
pub const TMIL: usize = 10_000_000;
pub const HMIL: usize = 100_000_000;

pub const FRAG_DIST: i64 = 10;
pub const MATE_MIN_DISTANCE: i64 = 20;
pub const MATE_MAX_DISTANCE: i64 = 5_000;
pub const MIN_MAPQ: u8 = 30;
pub const CB_LENGTH: usize = 16;
pub const TN5_LEFT_OFFSET: i64 = 4;
pub const TN5_RIGHT_OFFSET: i64 = 5;
pub const IS_TENX: bool = true;

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
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with peaks and frequency."),
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
                    Arg::with_name("ibed")
                        .long("ibed")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BED file with CB sequences."),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    match matches.subcommand_matches("filter") {
        Some(sub_m) => filter::filter(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("correct") {
        Some(sub_m) => barcode::correct(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("sort") {
        Some(sub_m) => sort::sort(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("group") {
        Some(sub_m) => group::dedup(&sub_m)?,
        None => (),
    }

    match matches.subcommand_matches("callpeak") {
        Some(sub_m) => peak::callpeak(&sub_m)?,
        None => (),
    }

    match matches.subcommand_matches("count") {
        Some(sub_m) => count::count(&sub_m)?,
        None => (),
    }

    match matches.subcommand_matches("text") {
        Some(sub_m) => text::convert(&sub_m)?,
        None => (),
    };

    Ok(())
}

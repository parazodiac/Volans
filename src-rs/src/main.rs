extern crate bincode;
extern crate bio;
extern crate clap;
extern crate flate2;
extern crate itertools;
extern crate libradicl;
extern crate pretty_env_logger;
extern crate rust_htslib;
extern crate serde;
extern crate indicatif;
extern crate num_format;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod barcode;
mod filter;
mod fragments;
mod fstats;
mod sort;
mod text;

pub const MIL: usize = 1_000_000;
pub const TMIL: usize = 10_000_000;
pub const HMIL: usize = 100_000_000;

pub const MATE_MIN_DISTANCE: i64 = 20;
pub const MATE_MAX_DISTANCE: i64 = 5_000;
pub const MIN_MAPQ: u8 = 30;
pub const CB_LENGTH: usize = 16;
pub const TN5_LEFT_OFFSET: i64 = 4;
pub const TN5_RIGHT_OFFSET: i64 = 5;
pub const CB_ORIENT_FW: bool = false;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("flash")
        .version("0.1.0")
        .author("Avi Srivastava")
        .about("A set of fast helper functions for Cut&Tag/ATAC data.")
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
            SubCommand::with_name("correct")
                .about("A subcommand to sequence correct the cb sequences.")
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
            SubCommand::with_name("text")
                .about("A subcommand to convert binary bed to text.")
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
        .get_matches();
    pretty_env_logger::init_timed();

    match matches.subcommand_matches("sort") {
        Some(sub_m) => sort::sort(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("text") {
        Some(sub_m) => text::convert(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("correct") {
        Some(sub_m) => barcode::correct(&sub_m)?,
        None => (),
    };

    match matches.subcommand_matches("filter") {
        Some(sub_m) => filter::filter(&sub_m)?,
        None => (),
    };

    Ok(())
}

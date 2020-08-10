extern crate bio;
extern crate clap;
extern crate flate2;
extern crate itertools;
extern crate libradicl;
extern crate pretty_env_logger;
extern crate rust_htslib;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod bam;
mod barcode;
mod fragments;

pub const MIL: usize = 1_000_000;
pub const MATE_DISTANCE: i64 = 5_000;
pub const MIN_MAPQ: u8 = 30;
pub const CB_LENGTH: usize = 16;
pub const TN5_LEFT_OFFSET: i64 = 4;
pub const TN5_RIGHT_OFFSET: i64 = 5;
pub const CB_ORIENT_FW: bool = false;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("flash")
        .version("0.1.0")
        .author("Avi Srivastava")
        .about("A set of fast helper functions for Cut&Tag data.")
        .subcommand(
            SubCommand::with_name("cbcorrect")
                .about("A barcode sequence correction method.")
                .arg(
                    Arg::with_name("barcode")
                        .long("barcode")
                        .short("b")
                        .takes_value(true)
                        .required(true)
                        .multiple(true)
                        .help("path to the cellular barcode file"),
                ),
        )
        .subcommand(
            SubCommand::with_name("fragments")
                .about("A barcode sequence correction method.")
                .arg(
                    Arg::with_name("bam")
                        .long("bam")
                        .short("b")
                        .takes_value(true)
                        .required(true)
                        .help("path to the BAM file"),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    match matches.subcommand_matches("cbcorrect") {
        Some(sub_m) => {
            let ret = barcode::correct(&sub_m);
            return ret;
        }
        None => (),
    };

    let tfile = "/home/srivastavaa/avi/bingzie/flash/workflow/output/test.bed";
    match matches.subcommand_matches("fragments") {
        Some(sub_m) => bam::fragments(&sub_m, &tfile)?,
        None => (),
    };

    Ok(())
}

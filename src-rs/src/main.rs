extern crate bio;
extern crate clap;
extern crate flate2;
extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod barcode;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("signac-kit")
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
                        .help("path to the cellular barcode file"),
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

    Ok(())
}

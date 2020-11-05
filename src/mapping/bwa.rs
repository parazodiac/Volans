use clap::ArgMatches;
use std::error::Error;

use crate::carina::*;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let ibam_path = carina::file::file_path_from_clap(sub_m, "ibam");

    Ok(())
}
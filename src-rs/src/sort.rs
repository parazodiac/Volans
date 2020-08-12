use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use clap::ArgMatches;

pub fn sort(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED files: {:?}", bed_file_path);

    let input_bed = BufReader::new(File::open(bed_file_path).expect("Can't open BED file"));
    for record in input_bed.lines() {
        let line = record.expect("Error reading record.");
        let toks: Vec<&str> = line.split("\t").collect();
        if toks[1] == "chr1" {
            let b = 2;
        } else {
            continue;
        }
    }

    Ok(())
}

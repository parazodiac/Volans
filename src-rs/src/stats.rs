use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

use crate::fragments::FragmentFile;
use clap::ArgMatches;
use num_format::{Locale, ToFormattedString};

pub fn stats(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);
    let input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let mut num_lines = 0;
    for _ in FragmentFile::new(input_bed).into_iter() {
        if num_lines % crate::TMIL == 0 {
            print!("\rDone processing {}0M reads", num_lines / crate::TMIL);
            std::io::stdout().flush().expect("Can't flush output");
        }

        num_lines += 1;
    }

    println!();
    info!(
        "Saw total {} lines in the BED file.",
        (num_lines).to_formatted_string(&Locale::en)
    );
    Ok(())
}

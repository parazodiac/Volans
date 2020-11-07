use clap::ArgMatches;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use carina::barcode::*;

use crate::quantify::fragments::Fragment;
use num_format::{Locale, ToFormattedString};
use std::collections::HashSet;

pub fn correct(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let wtl_file_path = Path::new(
        sub_m
            .value_of("whitelist")
            .expect("can't find whitelist flag"),
    )
    .canonicalize()
    .expect("can't find absolute path of input whitelist file");
    info!("Found whitelist CB file: {:?}", wtl_file_path);

    info!("Importing Whitelist barcode");
    let mut wtl_file = BufReader::new(File::open(wtl_file_path).expect("Unable to open file"));

    let mut wtl_strings = String::new();
    wtl_file
        .read_to_string(&mut wtl_strings)
        .expect("Unable to read the file");

    let mut wtl_barcodes: HashSet<u64> = HashSet::new();
    for cb_str in wtl_strings.split_terminator('\n') {
        let cb_bytes = match crate::configs::IS_WTL_FWD {
            true => cb_str.as_bytes().to_owned(),
            false => bio::alphabets::dna::revcomp(cb_str.as_bytes()).to_owned(),
        };

        let cb_id = cb_string_to_u64(&cb_bytes)?;
        wtl_barcodes.insert(cb_id);
    }

    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);
    let mut input_bed =
        BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let correct_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".corrected.bed";
    info!("Creating CB corrected BED file: {:?}", correct_file_path);
    let mut output_bed =
        BufWriter::new(File::create(correct_file_path).expect("Can't create BED file"));

    let mut mem_block = [0; 28];
    let mut num_lines = 0;
    let mut num_corrected = 0;
    while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
        num_lines += 1;
        if num_lines % crate::configs::TMIL == 0 {
            print!("\rDone processing {}0M reads", num_lines / crate::configs::TMIL);
            std::io::stdout().flush().expect("Can't flush output");
        }

        if wtl_barcodes.contains(&frag.cb) {
            frag.write(&mut output_bed, "binary")?;
            num_corrected += 1;
        }
    }

    println!();
    info!(
        "Total Fragments Passed {} out of {} ({:.2}%)",
        num_corrected.to_formatted_string(&Locale::en),
        num_lines.to_formatted_string(&Locale::en),
        num_corrected as f32 * 100.0 / num_lines as f32,
    );

    Ok(())
}

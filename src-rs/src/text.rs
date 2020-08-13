use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fragments::Fragment;
use clap::ArgMatches;

pub fn convert(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);

    let text_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".text.bed";
    info!("Creating text BED file: {:?}", text_file_path);

    let mut input_bed = BufReader::new(File::open(bed_file_path).expect("Can't open BED file"));
    let mut output_bed =
        BufWriter::new(File::create(text_file_path).expect("Can't open output BED file"));

    let mut num_lines = 0;
    let mut mem_block = [0; 28];
    while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
        num_lines += 1;
        if num_lines % crate::TMIL == 0 {
            print!("\rDone processing {}0M reads", num_lines / crate::TMIL);
            std::io::stdout().flush().expect("Can't flush output");
        }

        frag.write(&mut output_bed, "text")?;
    }

    Ok(())
}
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::collections::HashSet;

use crate::fragments::Fragment;
use clap::ArgMatches;

use indicatif::{ProgressBar, ProgressStyle};

pub fn sort(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);

    let sorted_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".sort.bed";
    info!("Creating sorted BED file: {:?}", sorted_file_path);
    let mut output_bed = BufWriter::new(File::create(sorted_file_path).expect("Can't create BED file"));

    let mut mem_block = [0; 28];
    let mut chr_names: HashSet<u32> = HashSet::new();
    {
        info!("Finding unique chromosome names in top 100M lines");
        let mut num_lines = 0;
        let mut input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));
        while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
            num_lines += 1;
            chr_names.insert(frag.chr);

            if num_lines > crate::HMIL { break; }
        }
        info!("Found Total {:?} unique chromosome names", chr_names.len());
    }
    
    let mut chr_names: Vec<u32> = chr_names.into_iter().collect();
    chr_names.sort();

    let pbar = ProgressBar::new(chr_names.len() as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    info!("Sorting files");
    for chr_name in chr_names {
        let mut input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));
        while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
            if frag.chr == chr_name {
                frag.write(&mut output_bed, "binary")?;
            }
        }
        pbar.inc(1);
    }

    Ok(())
}

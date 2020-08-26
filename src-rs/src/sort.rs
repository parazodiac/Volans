use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter};
use std::path::Path;

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
        + ".sorted.bed";

    let mut num_lines = 0;
    let mut mem_block = [0; 28];
    let mut chr_names: HashSet<u32> = HashSet::new();
    {
        info!("Finding unique chromosome names in the full file");
        let mut input_bed =
            BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));
        while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
            num_lines += 1;
            chr_names.insert(frag.chr);
        }
        info!("Found Total {:?} unique chromosome names", chr_names.len());
    }

    let mut chr_names: Vec<u32> = chr_names.into_iter().collect();
    chr_names.sort();
    {
        let mut file_handles: Vec<BufWriter<_>> = chr_names
            .iter()
            .map(|x| {
                BufWriter::new(
                    File::create(sorted_file_path.clone() + &x.to_string())
                        .expect("Can't create BED file"),
                )
            })
            .collect();

        let mut chr_names_to_idx = HashMap::new();
        for (idx, chr_name) in chr_names.iter().enumerate() {
            chr_names_to_idx.insert(chr_name, idx);
        }

        let pbar = ProgressBar::new(num_lines);
        pbar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}% ({eta})",
                )
                .progress_chars("╢▌▌░╟"),
        );
        pbar.set_draw_delta(num_lines / 100);

        info!("Sorting files");
        let mut input_bed = BufReader::new(File::open(bed_file_path).expect("Can't open BED file"));
        while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
            pbar.inc(1);
            match chr_names_to_idx.get(&frag.chr) {
                Some(fh_idx) => frag.write(&mut file_handles[*fh_idx], "binary")?,
                None => unreachable!(),
            };
        }

        pbar.finish_with_message("Segregation Complete");
    } // closing all files.

    info!("Merging files into {}", sorted_file_path);
    let file_names: Vec<_> = chr_names
        .iter()
        .map(|x| sorted_file_path.clone() + &x.to_string())
        .collect();

    // renaming first file
    std::fs::rename(file_names.first().unwrap(), sorted_file_path.clone())?;
    let file = OpenOptions::new()
        .append(true)
        .open(sorted_file_path)
        .expect("Can't create BED file");
    let mut master_fh = BufWriter::new(file);

    for file_name in file_names.iter().skip(1) {
        let mut file_name = BufReader::new(File::open(file_name)?);

        let mut mem_block = [0; 28];
        let mut frags: Vec<Fragment> = Fragment::read(&mut file_name, &mut mem_block)
            .into_iter()
            .collect();
        quickersort::sort_by(&mut frags, &|a, b| a.start().cmp(&b.start()));

        for frag in frags {
            frag.write(&mut master_fh, "binary")?;
        }
    }

    info!("Deleting temporary files.");
    for file_name in file_names.into_iter().skip(1) {
        std::fs::remove_file(file_name)?;
    }

    Ok(())
}

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fragments::{Fragment, FragmentFile};
use clap::ArgMatches;
use itertools::Itertools;

pub fn dedup(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);
    let input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let grouped_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".grouped.bed";
    info!("Creating grouped BED file: {:?}", grouped_file_path);
    let mut output_bed =
        BufWriter::new(File::create(grouped_file_path).expect("Can't create BED file"));

    let mut joint_class: HashMap<Fragment, u16> = HashMap::with_capacity(1_000);
    for (chr, chr_group) in FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr)
        .into_iter()
    {
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        for frag in chr_group {
            let stat = joint_class.entry(frag).or_insert(0);
            *stat += 1;
        }

        for (mut frag, freq) in joint_class.drain() {
            frag.cb = freq as u64;
            frag.write(&mut output_bed, "binary")?;
        }
    }

    println!();
    Ok(())
}

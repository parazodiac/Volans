use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fragments::{Fragment, FragmentFile};
use bio::data_structures::interval_tree::IntervalTree;
use clap::ArgMatches;
use itertools::Itertools;

pub fn count(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input BED file");
    info!("Found BED file: {:?}", bed_file_path);
    let input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let cbed_file_path = Path::new(sub_m.value_of("cbed").expect("can't find CB BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input cb BED file");
    info!("Found BED file: {:?}", cbed_file_path);
    let cb_input_bed =
        BufReader::new(File::open(cbed_file_path.clone()).expect("Can't open CB BED file"));

    let mtx_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".counts.mtx";
    info!("Creating output MTX file: {:?}", mtx_file_path);
    let mut output_mtx =
        BufWriter::new(File::create(mtx_file_path).expect("Can't create MTX file"));

    let frag_group = FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr);
    let mut frag_iter = frag_group.into_iter();

    let frag_cb_group = FragmentFile::new(cb_input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr);
    let mut frag_cb_iter = frag_cb_group.into_iter();

    while let (Some(chr_frag), Some(chr_cb_frag)) = (frag_iter.next(), frag_cb_iter.next()) {
        let (chr, chr_group) = chr_frag;
        let (chr_cb, chr_cb_group) = chr_cb_frag;

        assert_eq!(chr, chr_cb);
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        let mut itree = IntervalTree::new();
        for (frag_idx, frag) in chr_group.into_iter().enumerate() {
            assert!(
                frag.start < frag.end,
                format!("{}, {}", frag.start, frag.end)
            );
            itree.insert(frag.start..frag.end, (frag_idx, frag.cb))
        }

        let frags: Vec<Fragment> = chr_cb_group.collect();
    }

    println!();
    Ok(())
}

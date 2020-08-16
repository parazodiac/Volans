use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fragments::{Fragment, FragmentFile};
use bio::data_structures::interval_tree::IntervalTree;
use clap::ArgMatches;
use itertools::Itertools;

pub fn callpeak(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);
    let input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let peak_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".peaks.bed";
    info!("Creating peak BED file: {:?}", peak_file_path);
    let mut output_bed =
        BufWriter::new(File::create(peak_file_path).expect("Can't create BED file"));

    let mut total_peaks = 0;
    for (chr, chr_group) in FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr)
        .into_iter()
    {
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        let mut frags: Vec<Fragment> = chr_group.collect();
        frags.sort_unstable_by(|a, b| a.cb.cmp(&b.cb).reverse());

        //let mut itree = IntervalTree::new();
        //for (fidx, frag) in frags.iter().enumerate() {
        //    assert!(
        //        frag.start < frag.end,
        //        format!("{}, {}", frag.start, frag.end)
        //    );
        //    itree.insert(frag.start..frag.end, fidx)
        //}

        //let mut peak_bed: Vec<Fragment> = Vec::with_capacity(frags.len() / 10);
        //let mut tasks: HashSet<usize> = (0..frags.len()).into_iter().collect();
        //for (fidx, frag) in frags.into_iter().enumerate() {
        //    if tasks.len() == 0 {
        //        break;
        //    }
        //    if !tasks.contains(&fidx) {
        //        continue;
        //    }

        //    let start = frag.start;
        //    let end = frag.end;
        //    let mut count = frag.cb;

        //    tasks.remove(&fidx);
        //    let start_search = std::cmp::max(0, start as i64 - crate::FRAG_DIST) as u64;
        //    for intv in itree.find(start_search..end + crate::FRAG_DIST as u64) {
        //        let intv_frag_idx = intv.data().0;
        //        if !tasks.contains(&intv_frag_idx) {
        //            continue;
        //        }

        //        let intv_frag_start = intv.interval().start;
        //        let intv_frag_end = intv.interval().end;

        //        if ((intv_frag_end as i64 - end as i64).abs() <= crate::FRAG_DIST
        //            && intv_frag_start >= start)
        //            || ((intv_frag_start as i64 - start as i64).abs() <= crate::FRAG_DIST
        //                && intv_frag_end <= end)
        //        {
        //            count += intv.data().1;
        //            tasks.remove(&intv_frag_idx);
        //        }
        //    } // end for itree.find()

        //    peak_bed.push(Fragment {
        //        chr: frag.chr,
        //        start: start,
        //        end: end,
        //        cb: count,
        //    });
        //} // end for frags.into_iter()

        //for frag in peak_bed {
        //    frag.write(&mut output_bed, "binary")?;
        //    total_peaks += 1;
        //}

        break;
    }

    println!();
    info!("Found total {} peaks", total_peaks);
    Ok(())
}

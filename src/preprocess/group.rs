use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::ops::Range;
use std::path::Path;

use crate::quantify::fragments::{Fragment, FragmentFile};
use clap::ArgMatches;
use itertools::Itertools;
use num_format::{Locale, ToFormattedString};

pub fn dedup(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let report_all_cb = match sub_m.occurrences_of("allcb") {
        0 => false,
        _ => true,
    };

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

    let mut total_frag = 0;
    let mut total_group = 0;
    let mut total_classes = 0;

    let mut joint_class = HashMap::with_capacity(500);
    for (chr, chr_group) in FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr)
        .into_iter()
    {
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        for frag in chr_group {
            let val = joint_class
                .entry(Range {
                    start: frag.start,
                    end: frag.end,
                })
                .or_insert_with(Vec::new);
            val.push(frag.cb);
        }

        for (range, mut cbs) in joint_class.drain() {
            total_frag += cbs.len();
            cbs.sort_unstable();
            cbs.dedup();
            total_classes += 1;
            total_group += cbs.len();

            if report_all_cb {
                for cb in cbs {
                    let frag = Fragment {
                        start: range.start,
                        end: range.end,
                        cb,
                        chr,
                    };
                    frag.write(&mut output_bed, "binary")?;
                }
            } else {
                let frag = Fragment {
                    start: range.start,
                    end: range.end,
                    cb: cbs.len() as u64,
                    chr,
                };
                frag.write(&mut output_bed, "binary")?;
            }
        }
    }

    println!();
    info!(
        "Saw total {} fragment and grouped into {} classes w/ {} deduplicated fragments.",
        (total_frag).to_formatted_string(&Locale::en),
        (total_classes).to_formatted_string(&Locale::en),
        (total_group).to_formatted_string(&Locale::en)
    );
    Ok(())
}

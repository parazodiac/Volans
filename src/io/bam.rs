use std::error::Error;
use std::io::Write;

use clap::ArgMatches;

use rust_htslib::bam;
use rust_htslib::bam::{Read, Record};

use itertools::Itertools;

use carina::barcode::cb_string_to_u64;

use crate::fragments::counting_stats;
use crate::fragments::filter;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bam_file_path = carina::file::file_path_from_clap(sub_m, "ibam")?;
    let mut obed_file = carina::file::bufwriter_from_clap(sub_m, "obed")?;

    let mut input_bam = bam::Reader::from_path(bam_file_path).expect("Can't open BAM file");
    let bam_header = input_bam.header().clone();
    input_bam.set_threads(4).unwrap();

    let (is_tenx, cb_extractor): (bool, fn(&Record) -> u64) = match sub_m.occurrences_of("tenx") {
        0 => (false, |aln: &Record| -> u64 {
            let qname = aln.qname();
            cb_string_to_u64(&qname[(qname.len() - crate::configs::CB_LENGTH)..])
                .expect("can't convert cb string to u64")
        }),
        _ => (true, |aln: &Record| -> u64 {
            cb_string_to_u64(&aln.aux(b"CB").unwrap().string()[..crate::configs::CB_LENGTH])
                .expect("can't convert cb string to u64")
        }),
    };

    let just_stats = match sub_m.occurrences_of("stats") {
        0 => false,
        _ => true,
    };

    let mito_string = sub_m
        .value_of("mitostr")
        .expect("can't find string to identify mitochondrial chromosome");
    let mito_tid = bam_header
        .tid(mito_string.as_bytes())
        .expect("Can't find mito string in the BAM file");
    info!(
        "Using {} as Mitochondrial Chromosome with {} as id.",
        mito_string, mito_tid
    );

    let mut counter = counting_stats::FragStats {
        ..Default::default()
    };
    for (_, read_group) in input_bam
        .records()
        .map(|res| res.unwrap())
        .group_by(|rec| rec.qname().to_owned())
        .into_iter()
    {
        counter.total_reads += 1;
        if counter.total_reads % crate::configs::MIL == 0 {
            print!(
                "\rDone processing {}M reads",
                counter.total_reads / crate::configs::MIL
            );
            std::io::stdout().flush().expect("Can't flush output");
        }

        let alignments: Vec<Record> = read_group.collect();
        match filter::callback(
            alignments,
            &mut counter,
            just_stats,
            is_tenx,
            mito_tid,
            cb_extractor,
        ) {
            Some(frag) => frag.write(&mut obed_file, "binary")?,
            None => continue,
        };
    }

    println!("{}", counter);
    Ok(())
}

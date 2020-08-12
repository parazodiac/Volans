use std::convert::TryInto;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use clap::ArgMatches;

use rust_htslib::bam;
use rust_htslib::bam::{Read, Record};

use crate::fragments::Fragment;
use itertools::Itertools;

use crate::fragments;
use crate::fstats::FragStats;
use crate::{MATE_MAX_DISTANCE, MATE_MIN_DISTANCE, MIN_MAPQ};

fn is_unmapped(alignments: &Vec<Record>, counter: &mut FragStats) -> bool {
    let mut no_map_count = 0;
    for alignment in alignments {
        if alignment.is_unmapped() {
            no_map_count += 1;
        }
    }

    if no_map_count > 0 {
        counter.unmap_skip += 1;
        if no_map_count != alignments.len() {
            counter.unmap_orphan += 1;
        }
        return true;
    }

    false
}

fn is_multi_mapping(alignments: &Vec<Record>, counter: &mut FragStats) -> bool {
    if alignments.len() > 2
        || alignments.first().unwrap().aux("XA".as_bytes()).is_some()
        || alignments.last().unwrap().aux("XA".as_bytes()).is_some()
    {
        counter.mm_reads += 1;
        return true;
    }

    false
}

fn is_mitochondrial(alignments: &Vec<Record>, counter: &mut FragStats, mito_tid: u32) -> bool {
    if alignments.first().expect("no alignments found").tid()
        == mito_tid
            .try_into()
            .expect("bam header inconsist for mito chromosome")
    {
        counter.mito_skip += 1;
        return true;
    }
    false
}

fn is_high_quality(alignments: &Vec<Record>, counter: &mut FragStats) -> bool {
    let mut min_quality = u8::MAX;
    for alignment in alignments {
        min_quality = std::cmp::min(min_quality, alignment.mapq());
    }

    if min_quality < MIN_MAPQ {
        counter.mapq_skip += 1;
        return false;
    }

    true
}

fn is_chimeric(alignments: &Vec<Record>, counter: &mut FragStats) -> bool {
    assert_eq!(alignments.len(), 2);
    let mut aln = alignments.first().unwrap();
    let mut maln = alignments.last().unwrap();
    if aln.is_reverse() {
        std::mem::swap(&mut aln, &mut maln);
    }

    // both reads mapped to different chromosome
    if aln.tid() != maln.tid() {
        counter.chimeric_tids += 1;
        return true;
    }
    // if first read is reverse or mate is not reverse
    // Note: we already swapped the first and second
    if aln.is_reverse() || !maln.is_reverse() {
        counter.chimeric_strand += 1;
        return true;
    }
    // if mate-pairs mapped too far
    if (aln.pos() - maln.pos()).abs() > MATE_MAX_DISTANCE {
        counter.chimeric_max_distance += 1;
        return true;
    }
    // if first is after second
    if fragments::soft_clip_pos(aln) + MATE_MIN_DISTANCE > fragments::soft_clip_pos(maln) {
        counter.chimeric_min_distance += 1;
        return true;
    }

    false
}

pub fn filter(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bam_file_path = Path::new(sub_m.value_of("ibam").expect("can't find BAM flag"))
        .canonicalize()
        .expect("can't find absolute path of input file");

    info!("Found BAM files: {:?}", bam_file_path);
    let mut input_bam = bam::Reader::from_path(bam_file_path).expect("Can't open BAM file");
    input_bam.set_threads(6).unwrap();
    let bam_header = input_bam.header().clone();

    let bed_file_path = sub_m.value_of("obed").expect("can't find BED flag");
    let mut obed_file = BufWriter::new(File::create(bed_file_path)?);
    let bed_file_path = Path::new(bed_file_path)
        .canonicalize()
        .expect("cann't resolved absolute bed file path");
    info!("Created BED file: {:?}", bed_file_path);

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

    let mut counter = FragStats {
        ..Default::default()
    };
    for (_, read_group) in input_bam
        .records()
        .map(|res| res.unwrap())
        .group_by(|rec| rec.qname().to_owned())
        .into_iter()
    {
        counter.total_reads += 1;
        if counter.total_reads % crate::MIL == 0 {
            print!(
                "\rDone processing {}M reads",
                counter.total_reads / crate::MIL
            );
            std::io::stdout().flush().expect("Can't flush output");
        }

        let alignments: Vec<Record> = read_group.collect();
        if is_unmapped(&alignments, &mut counter) {
            continue;
        }
        if is_multi_mapping(&alignments, &mut counter) {
            continue;
        }
        if is_mitochondrial(&alignments, &mut counter, mito_tid) {
            continue;
        }
        if is_high_quality(&alignments, &mut counter) {
            continue;
        }
        if is_chimeric(&alignments, &mut counter) {
            continue;
        }

        let aln = alignments.first().unwrap();
        let maln = alignments.last().unwrap();
        let frag = match alignments.first().unwrap().is_reverse() {
            true => Fragment::new(maln, aln),
            false => Fragment::new(aln, maln),
        };
        frag.write(&mut obed_file)?;
    }

    println!("{}", counter);
    Ok(())
}

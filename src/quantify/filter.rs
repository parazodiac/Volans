use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use clap::ArgMatches;

use rust_htslib::bam;
use rust_htslib::bam::{Read, Record};

use crate::quantify::fragments::Fragment;
use itertools::Itertools;

use crate::quantify::fragments;
use crate::quantify::fstats::FragStats;
use crate::configs::{MATE_MAX_DISTANCE, MATE_MIN_DISTANCE, MIN_MAPQ};

fn is_unmapped(alignments: &[Record], counter: &mut FragStats) -> bool {
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

fn is_multi_mapping(alignments: &[Record], counter: &mut FragStats) -> bool {
    if alignments.len() > 2
        || alignments.first().unwrap().aux(b"XA").is_some()
        || alignments.last().unwrap().aux(b"XA").is_some()
    {
        counter.mm_reads += 1;
        return true;
    }

    false
}

fn is_mitochondrial(alignments: &[Record], counter: &mut FragStats, mito_tid: u32) -> bool {
    if alignments.first().expect("no alignments found").tid() == mito_tid as i32 {
        counter.mito_skip += 1;
        return true;
    }
    false
}

fn is_high_quality(alignments: &[Record], counter: &mut FragStats) -> bool {
    let mut min_quality = u8::MAX;
    for alignment in alignments {
        min_quality = std::cmp::min(min_quality, alignment.mapq());
    }

    if min_quality < MIN_MAPQ {
        counter.mapq_skip += 1;
        return true;
    }

    false
}

fn is_chimeric(alignments: &[Record], counter: &mut FragStats) -> bool {
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

fn has_barcode_tag(alignments: &[Record], counter: &mut FragStats) -> bool {
    for aln in alignments {
        if aln.aux(b"CB").is_none() {
            counter.cb_skip += 1;
            return false;
        }
    }

    true
}

pub fn filter(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bam_file_path = Path::new(sub_m.value_of("ibam").expect("can't find BAM flag"))
        .canonicalize()
        .expect("can't find absolute path of input file");

    info!("Found BAM files: {:?}", bam_file_path);
    let mut input_bam = bam::Reader::from_path(bam_file_path).expect("Can't open BAM file");
    input_bam.set_threads(4).unwrap();
    let bam_header = input_bam.header().clone();

    let bed_file_path = sub_m.value_of("obed").expect("can't find BED flag");
    let mut obed_file = BufWriter::new(File::create(bed_file_path)?);
    let bed_file_path = Path::new(bed_file_path)
        .canonicalize()
        .expect("can't resolved absolute bed file path");
    info!("Created BED file: {:?}", bed_file_path);

    let (is_tenx, cb_extractor): (bool, fn(&Record) -> u64) = match sub_m.occurrences_of("tenx") {
        0 => (false, |aln: &Record| -> u64 {
            let qname = aln.qname();
            fragments::cb_string_to_u64(&qname[(qname.len() - crate::configs::CB_LENGTH)..])
                .expect("can't convert cb string to u64")
        }),
        _ => (true, |aln: &Record| -> u64 {
            fragments::cb_string_to_u64(&aln.aux(b"CB").unwrap().string()[..crate::configs::CB_LENGTH])
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
        if counter.total_reads % crate::configs::MIL == 0 {
            print!(
                "\rDone processing {}M reads",
                counter.total_reads / crate::configs::MIL
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

        if is_tenx && !has_barcode_tag(&alignments, &mut counter) {
            continue;
        }

        if just_stats {
            continue;
        }

        let aln = alignments.first().unwrap();
        let maln = alignments.last().unwrap();
        let frag = match alignments.first().unwrap().is_reverse() {
            true => Fragment::new(maln, aln, cb_extractor),
            false => Fragment::new(aln, maln, cb_extractor),
        };
        frag.write(&mut obed_file, "binary")?;
    }

    println!("{}", counter);
    Ok(())
}

use std::convert::TryInto;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use clap::ArgMatches;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{Read, Record};

use crate::fragments::Fragment;
use itertools::Itertools;

use crate::{CB_LENGTH, MATE_DISTANCE, MIN_MAPQ, TN5_LEFT_OFFSET, TN5_RIGHT_OFFSET};

fn soft_clip_pos(aln: &Record) -> i64 {
    let mut softclip_offset = 0;
    for cigar in aln.cigar().iter() {
        match cigar {
            Cigar::SoftClip(val) => softclip_offset += val,
            _ => (),
        }
    }

    match aln.is_reverse() {
        true => aln.pos() - softclip_offset as i64,
        false => aln.cigar().end_pos() + softclip_offset as i64,
    }
}

fn is_chimeric(aln: &Record, maln: &Record) -> bool {
    // both reads mapped to different chromosome
    if aln.tid() != maln.tid() {
        return true;
    }
    // if first read is reverse or mate is not reverse
    // Note: we already swapped the first and second
    if aln.is_reverse() || !maln.is_reverse() {
        return true;
    }
    // if mate-pairs mapped too far
    if (aln.pos() - maln.pos()).abs() > MATE_DISTANCE {
        return true;
    }
    // if first is after second
    if soft_clip_pos(aln) > soft_clip_pos(maln) {
        return true;
    }

    false
}

pub fn fragments(sub_m: &ArgMatches, tfile: &str) -> Result<(), Box<dyn Error>> {
    let bam_file_path = Path::new(sub_m.value_of("bam").expect("can't find BAM flag"))
        .canonicalize()
        .expect("can't find absolute path of input file");

    info!("Found BAM files: {:?}", bam_file_path);

    let mut input_bam = bam::Reader::from_path(bam_file_path).expect("Can't open BAM file");
    input_bam.set_threads(6).unwrap();

    let mito_string = "chrM";
    let bam_header = input_bam.header().clone();
    let mito_tid = bam_header.tid(mito_string.as_bytes()).unwrap();
    info!(
        "Using {} as Mitochondrial Chromosome with {} as id.",
        mito_string, mito_tid
    );

    let mut mm_reads = 0;
    let mut mapq_skip = 0;
    let mut read_count = 0;
    let mut mapq_orphan = 0;
    let mut unmap_skip = 0;
    let mut unmap_orphan = 0;
    let mut chimeric_skip = 0;
    let mut mito_skip = 0;

    let mut file = BufWriter::new(File::create(tfile)?);

    for (qname, read_group) in input_bam
        .records()
        .map(|res| res.unwrap())
        .group_by(|rec| rec.qname().to_owned())
        .into_iter()
    {
        read_count += 1;
        if read_count % crate::MIL == 0 {
            print!("\rDone processing {}M reads", read_count / crate::MIL);
            std::io::stdout().flush().expect("Can't flush output");
        }

        let alignments: Vec<Record> = read_group.collect();
        if alignments.len() > 2 {
            mm_reads += 1;
            continue;
        }

        let mut first_mate_pair = alignments.first().unwrap();
        let mut second_mate_pair = alignments.last().unwrap();

        let is_one_unmapped = first_mate_pair.is_unmapped() || second_mate_pair.is_unmapped();
        if is_one_unmapped {
            unmap_skip += 1;
            let is_one_mapped = !first_mate_pair.is_unmapped() || !second_mate_pair.is_unmapped();
            if is_one_mapped {
                unmap_orphan += 1;
            }
            continue;
        }

        if first_mate_pair.tid()
            == mito_tid
                .try_into()
                .expect("bam header inconsistency for mito chromosome")
        {
            mito_skip += 1;
            continue;
        }

        // assume the first one is always in forward strand if not then swap
        if first_mate_pair.is_reverse() {
            std::mem::swap(&mut first_mate_pair, &mut second_mate_pair);
        }

        if is_chimeric(first_mate_pair, second_mate_pair) {
            chimeric_skip += 1;
            continue;
        }

        let mapq_first = first_mate_pair.mapq();
        let mapq_second = first_mate_pair.mapq();
        if std::cmp::min(mapq_first, mapq_second) < MIN_MAPQ {
            mapq_skip += 1;
            if mapq_first >= MIN_MAPQ || mapq_second >= MIN_MAPQ {
                mapq_orphan += 1;
            }
            continue;
        }

        let frag = Fragment::new(
            first_mate_pair.tid() as u32,
            (soft_clip_pos(first_mate_pair) + TN5_LEFT_OFFSET) as u64,
            (soft_clip_pos(second_mate_pair) - TN5_RIGHT_OFFSET) as u64,
            &qname[(qname.len() - CB_LENGTH)..],
        );
        frag.write(&mut file)?;
    }

    println!();
    let total_skipped = mm_reads + mapq_skip + unmap_skip + chimeric_skip + mito_skip;

    info!("STATS: Total Reads: {}", read_count);
    info!(
        "STATS: Total MultiMapping Reads: {}({:.02}%)",
        mm_reads,
        mm_reads as f32 * 100.0 / read_count as f32
    );
    info!(
        "STATS: Total Chimeric Reads: {}({:.02}%)",
        chimeric_skip,
        chimeric_skip as f32 * 100.0 / read_count as f32
    );
    info!(
        "STATS: Total Mitochondrial Reads: {}({:.02}%)",
        mito_skip,
        mito_skip as f32 * 100.0 / read_count as f32
    );
    info!(
        "STATS: Total MAPQ skip: {}({:.02}%) with {}({:.02}%) HighQ orphan",
        mapq_skip,
        mapq_skip as f32 * 100.0 / read_count as f32,
        mapq_orphan,
        mapq_orphan as f32 * 100.0 / read_count as f32,
    );
    info!(
        "STATS: Total Unmapped skip: {}({:.02}%) with {}({:.02}%) HighQ orphan",
        unmap_skip,
        unmap_skip as f32 * 100.0 / read_count as f32,
        unmap_orphan,
        unmap_orphan as f32 * 100.0 / read_count as f32,
    );
    info!(
        "STATS: Total Reads skipped: {}({:.02}%)",
        total_skipped,
        total_skipped as f32 * 100.0 / read_count as f32,
    );

    Ok(())
}

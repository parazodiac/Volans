use rust_htslib::bam::Record;

use crate::fragments::counting_stats::FragStats;
use crate::fragments::schema::{soft_clip_pos, Fragment};

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
    if soft_clip_pos(aln) + MATE_MIN_DISTANCE > soft_clip_pos(maln) {
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

pub fn callback(
    alignments: Vec<Record>,
    mut counter: &mut FragStats,
    stats_only: bool,
    is_tenx: bool,
    mito_tid: u32,
    cb_extractor: fn(&Record) -> u64,
) -> Option<Fragment> {
    if is_unmapped(&alignments, &mut counter)
        || is_multi_mapping(&alignments, &mut counter)
        || is_mitochondrial(&alignments, &mut counter, mito_tid)
        || is_high_quality(&alignments, &mut counter)
        || is_chimeric(&alignments, &mut counter)
        || (is_tenx && !has_barcode_tag(&alignments, &mut counter))
        || stats_only
    {
        return None;
    }

    let aln = alignments.first().unwrap();
    let maln = alignments.last().unwrap();
    let frag = match alignments.first().unwrap().is_reverse() {
        true => Fragment::new(maln, aln, cb_extractor),
        false => Fragment::new(aln, maln, cb_extractor),
    };

    Some(frag)
}

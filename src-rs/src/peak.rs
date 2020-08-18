use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fragments::{Feature, Fragment, FragmentFile, Interval};
use clap::ArgMatches;
use itertools::Itertools;

use bitvector::*;
use num_format::{Locale, ToFormattedString};

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

// fn process_feature_group(chr: u32, features: &Vec<&Feature>) -> Vec<Feature> {
//     let mut peaks: Vec<Feature> = Vec::new();

//     if chr == 3 {
//         let start = features.first().unwrap().start;
//         let end = features.last().unwrap().end;
//     //         // let start = feat.start;
//     //         // let end = feat.end;

//         for feat in features {
//             if feat.start <= 6_780_455 && feat.end >= 6_772_512 {
//                 peaks.push(*feat.clone());
//             }
//         }
//     }

//     // for feat in features {
//     //    peaks.push(Feature {
//     //        start: feat.start, end: feat.end, count: feat.count
//     //    });
//     // }

//     peaks
// }

fn process_pileup(feat_pile_up: &[u16]) -> Option<f32> {
    let feature_len = feat_pile_up.len();
    let one_deviation_distance = feature_len / 8;

    let total_mass: u32 = feat_pile_up.iter().map(|x| *x as u32).sum();
    let mut extrema_mass: u32 = feat_pile_up
        .iter()
        .take(2 * one_deviation_distance)
        .map(|x| *x as u32)
        .sum();
    extrema_mass += feat_pile_up
        .iter()
        .skip(6 * one_deviation_distance)
        .map(|x| *x as u32)
        .sum::<u32>();

    let prob: f32 = extrema_mass as f32 / total_mass as f32;
    if prob >= 0.75 {
        Some(prob)
    } else {
        None
    }
}

fn process_feature_group(features: &Vec<&Feature>) -> Option<Vec<Feature>> {
    let region_start = features.first().unwrap().start as i64;
    let region_end = features.iter().rev().map(|x| x.end).max().unwrap() as i64;
    let region_size: usize = std::cmp::max(0, region_end - region_start + 1) as usize;

    assert_ne!(region_size, 0);
    let mut pile_up: Vec<u16> = vec![0; region_size];
    for feat in features {
        (feat.start..feat.end).into_iter().for_each(|x| {
            assert!(x >= region_start as u32);

            pile_up[x as usize - region_start as usize] += feat.count as u16;
        });
    }

    let mut bitvec = BitVector::new(region_size);
    let mut sorted_feature_indices: Vec<usize> = (0..region_size).collect();
    sorted_feature_indices.sort_unstable_by(|a, b| pile_up[*a].cmp(&pile_up[*b]).reverse());

    let mut peaks: Vec<Feature> = Vec::new();
    for feat_idx in sorted_feature_indices {
        if bitvec.contains(feat_idx) {
            continue;
        }
        bitvec.insert(feat_idx);

        let start = std::cmp::max(0, feat_idx as i64 - (crate::WINDOW_SIZE / 2)) as usize;
        let end = std::cmp::min(
            pile_up.len(),
            (feat_idx as i64 + (crate::WINDOW_SIZE / 2)) as usize,
        );

        let peak_prob: Option<f32> = process_pileup(&pile_up[start..end]);
        match peak_prob {
            None => (),
            Some(prob) => {
                for i in start..end + 1 {
                    bitvec.insert(i);
                }
                peaks.push(Feature {
                    start: region_start as u32 + start as u32,
                    end: region_start as u32 + end as u32,
                    count: (prob * 100.0) as u32,
                });
            }
        };
    } // end for

    if peaks.len() == 0 {
        return None;
    }
    Some(peaks)
}

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

    let mut total_classes = 0;
    let mut total_groups = 0;
    let mut total_peaks = 0;
    let mut noise = 0;
    let mut features = Vec::with_capacity(500);
    for (chr, chr_group) in FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr)
        .into_iter()
    {
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        chr_group.for_each(|frag| {
            let feature = Feature {
                start: frag.start as u32,
                end: frag.end as u32,
                count: frag.cb as u32,
            };

            features.push(feature);
        });

        features.sort_unstable_by(|a, b| a.start.cmp(&b.start));

        let num_regions: usize;
        let mut features_group: Vec<Vec<&Feature>> = Vec::new();
        {
            // grouping the features into region
            let mut regions: Vec<Interval> = Vec::new();
            let first_feat = features.first().unwrap();
            let mut cur_region = Interval {
                start: first_feat.start,
                end: first_feat.end,
            };
            let mut cur_fgroup = vec![first_feat];

            for feature in features.iter().skip(1) {
                if feature.start >= cur_region.end {
                    regions.push(cur_region);
                    cur_region = Interval {
                        start: feature.start,
                        end: feature.end,
                    };

                    features_group.push(cur_fgroup);
                    cur_fgroup = vec![feature];

                    continue;
                }

                cur_fgroup.push(feature);
                let position_diff: i64 = feature.end as i64 - cur_region.end as i64;
                if position_diff > 0 {
                    cur_region.end = feature.end;
                }
            }

            if regions.len() > 0 && (*regions.last().unwrap() != cur_region) {
                regions.push(cur_region);
                features_group.push(cur_fgroup);
            }

            num_regions = regions.len();
        }

        for i in 0..num_regions {
            // if i > 1 { break; }
            let features = &features_group[i];
            let num_supporting_barcodes = features.iter().map(|x| x.count).sum::<u32>();

            total_groups += num_supporting_barcodes;
            if num_supporting_barcodes < crate::NUM_SUPPORT_CB {
                noise += num_supporting_barcodes;
                continue;
            }

            total_classes += 1;
            let peaks = match process_feature_group(features) {
                Some(peaks) => peaks,
                None => continue,
            };

            total_peaks += peaks.len();
            for peak in peaks {
                let frag = Fragment {
                    chr: chr,
                    start: peak.start as u64,
                    end: peak.end as u64,
                    cb: peak.count as u64,
                };

                frag.write(&mut output_bed, "text")?;
            }
        }

        features.clear();
        // break;
    }

    println!();
    info!(
        "Removed {} noisy groups (<5 CB evidence) out of {}({:.2}%)",
        (noise).to_formatted_string(&Locale::en),
        (total_groups).to_formatted_string(&Locale::en),
        (noise as f32 * 100.0 / total_groups as f32)
    );

    assert!(noise <= total_groups);
    total_groups -= noise;
    info!(
        "Found total {} classes out of {} HighQ groups. ({:.02}% reduction).",
        (total_classes).to_formatted_string(&Locale::en),
        (total_groups).to_formatted_string(&Locale::en),
        100.0 - (total_classes as f32 * 100.0 / total_groups as f32)
    );

    info!(
        "Found total {} peaks out of {} HighQ classes. ({:.02}% reduction).",
        (total_peaks).to_formatted_string(&Locale::en),
        (total_classes).to_formatted_string(&Locale::en),
        100.0 - (total_peaks as f32 * 100.0 / total_classes as f32)
    );

    Ok(())
}

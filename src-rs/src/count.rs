use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::ops::Range;
use std::path::Path;

use crate::fragments::FragmentFile;
use bio::data_structures::interval_tree::{Entry, IntervalTree};
use clap::ArgMatches;
use itertools::Itertools;
use sprs::TriMat;

use rust_htslib::bam;
use crate::rust_htslib::bam::Read;
use num_format::{Locale, ToFormattedString};

pub fn count(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed_file_path = Path::new(sub_m.value_of("pbed").expect("can't find peak BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input peak BED file");
    info!("Found BED file: {:?}", bed_file_path);
    let input_bed =
        BufReader::new(File::open(bed_file_path.clone()).expect("Can't open peak BED file"));

    let cbed_file_path = Path::new(sub_m.value_of("cbed").expect("can't find CB BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input cb BED file");
    info!("Found BED file: {:?}", cbed_file_path);
    let cb_input_bed =
        BufReader::new(File::open(cbed_file_path.clone()).expect("Can't open CB BED file"));

    let frag_group = FragmentFile::new(input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr);
    let mut frag_iter = frag_group.into_iter();

    let frag_cb_group = FragmentFile::new(cb_input_bed)
        .map(|maybe_frag| maybe_frag)
        .group_by(|frag| frag.chr);
    let mut frag_cb_iter = frag_cb_group.into_iter();

    let mut total_reads = 0;
    let mut skipped_chr_indices = Vec::new();
    let mut counts: Vec<HashMap<Range<u32>, HashMap<u64, u16>>> = Vec::new();
    while let (Some(chr_frag), Some(chr_cb_frag)) = (frag_iter.next(), frag_cb_iter.next()) {
        let (chr, chr_group) = chr_frag;
        let (chr_cb, chr_cb_group) = chr_cb_frag;
        if chr != chr_cb {
            warn!(
                "Early breaking after {} as some chromosomes does not appear in both files",
                chr
            );
            break;
        }

        skipped_chr_indices.push(chr);
        print!("\rWorking on Chromosome: {}", chr);
        std::io::stdout().flush().expect("Can't flush output");

        let peaks: Vec<Range<u32>> = chr_group
            .map(|x| Range {
                start: x.start as u32,
                end: x.end as u32,
            })
            .collect();

        // making interval tree for the peaks to find the overlap
        let mut peaks_tree = IntervalTree::new();
        for (peak_idx, peak) in peaks.iter().enumerate() {
            assert!(
                peak.start < peak.end,
                format!("{}, {}", peak.start, peak.end)
            );
            peaks_tree.insert(peak, peak_idx);
        }

        counts.push(HashMap::new());
        let count_matrix = counts.last_mut().unwrap();

        for (range, class) in &chr_cb_group.group_by(|frag| Range {
            start: frag.start as u32,
            end: frag.end as u32,
        }) {
            let overlapping_peaks: Vec<Entry<u32, usize>> = peaks_tree.find(range).collect();
            if overlapping_peaks.len() != 1 {
                continue;
            }

            let intv = overlapping_peaks.first().unwrap().interval();
            let peak = Range {
                start: intv.start,
                end: intv.end,
            };

            let cb_group: Vec<u64> = class.map(|x| x.cb).collect();
            for cb in cb_group {
                let val = count_matrix
                    .entry(peak.clone())
                    .or_insert(HashMap::new())
                    .entry(cb)
                    .or_insert(0);
                *val += 1;
                total_reads += 1;
            }
        } // end-for
    } // end-while

    println!();
    info!("Found total {} reads in the matrix", 
        total_reads.to_formatted_string(&Locale::en)
    );

    let mut col_names: HashMap<u64, usize> = HashMap::new();
    counts.iter().for_each(|matrix| {
        matrix.iter().for_each(|(_, row_data)| {
            row_data.iter().for_each(|(col_id, _)| {
                if !col_names.contains_key(col_id) {
                    let name_len = col_names.len();
                    col_names.insert(*col_id, name_len);
                }
            });
        });
    });

    let mut row_ids = Vec::new();
    let mut col_ids = Vec::new();
    let mut data = Vec::new();

    let bam_header = match sub_m.value_of("bam") {
        Some(path) => {
            let bam_file_path = Path::new(path).canonicalize()
                .expect("can't find absolute path of input BAM file");
            info!("Found BAM file: {:?}", bam_file_path);

            let input_bam = bam::Reader::from_path(bam_file_path)
                .expect("Can't open BAM file");
            Some(input_bam.header().clone())
        },
        None => None
    };

    let mut row_names: HashMap<String, usize> = HashMap::new();
    for (chr_idx, count) in counts.into_iter().enumerate() {
        let skip_chr_id = skipped_chr_indices[chr_idx];
        let chr_name = skip_chr_id.to_string();
        for (row_range, row_data) in count {
            if row_data.len() < 10 {
                continue;
            }

            let mut row_name = match &bam_header {
                Some(header) => std::str::from_utf8(header.tid2name(skip_chr_id)).unwrap().to_string(),
                None => chr_name.clone(),
            };
            row_name.push_str("_");
            row_name.push_str(&row_range.start.to_string());
            row_name.push_str(":");
            row_name.push_str(&row_range.end.to_string());

            let row_id = match row_names.contains_key(&row_name) {
                true => *row_names.get(&row_name).unwrap(),
                false => {
                    let name_len = row_names.len();
                    row_names.insert(row_name, name_len);
                    name_len
                }
            };

            for (col_name, col_data) in row_data {
                let col_id = col_names.get(&col_name).unwrap();
                row_ids.push(row_id);
                col_ids.push(*col_id);
                data.push(col_data);
            }
        }
    }

    let matrix =
        TriMat::from_triplets((row_names.len(), col_names.len()), row_ids, col_ids, data).to_csr();

    let mtx_file_path = bed_file_path
        .parent()
        .unwrap()
        .join("counts.mtx");

    info!("Creating output MTX file: {:?}", mtx_file_path);
    sprs::io::write_matrix_market(mtx_file_path, &matrix)?;

    {
        let rows_file_path = bed_file_path
            .parent()
            .unwrap()
            .join("counts_rows.txt");
        let mut file = BufWriter::new(File::create(rows_file_path)?);
        let mut sorted_row_names = vec![String::new(); row_names.len()];
        row_names.into_iter().for_each(|(k,v)| {
            sorted_row_names[v] = k;
        });

        for row_name in sorted_row_names {
            file.write((row_name + "\n").as_bytes())?;
        }
    }

    {
        let cols_file_path = bed_file_path
            .parent()
            .unwrap()
            .join("counts_cols.txt");
        let mut file = BufWriter::new(File::create(cols_file_path)?);
        let mut sorted_col_names = vec![String::new(); col_names.len()];
        col_names.into_iter().for_each(|(k,v)| {
            sorted_col_names[v] = crate::fragments::u64_to_cb_string(k).unwrap();
        });

        for col_name in sorted_col_names {
            file.write((col_name + "\n").as_bytes())?;
        }
    }

    Ok(())
}

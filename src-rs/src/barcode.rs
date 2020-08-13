use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::io::Write;
use std::path::Path;
use std::io::Read;
use clap::ArgMatches;

use crate::fragments::{Fragment, cb_string_to_u64};
use std::collections::HashSet;
use flate2::read::MultiGzDecoder;
use num_format::{Locale, ToFormattedString};

pub fn correct(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {    
    info!("Importing Whitelist barcode");
    let mut wtl_file = MultiGzDecoder::new(
        File::open("./ext/737K-cratac-v1.txt.gz").expect("Unable to open file"),
    );

    let mut wtl_strings = String::new();
    wtl_file
        .read_to_string(&mut wtl_strings)
        .expect("Unable to read the file");

    let mut wtl_barcodes: HashSet<u64> = HashSet::new();
    for cb_str in wtl_strings.split_terminator("\n") {
        let cb_bytes = match crate::CB_ORIENT_FW {
            true => cb_str.as_bytes().to_owned(),
            false => bio::alphabets::dna::revcomp(cb_str.as_bytes()).to_owned(),
        };

        let cb_id = cb_string_to_u64(&cb_bytes)?;
        wtl_barcodes.insert(cb_id);
    }
   
    let bed_file_path = Path::new(sub_m.value_of("ibed").expect("can't find BED flag"))
        .canonicalize()
        .expect("can't find absolute path of input bed file");
    info!("Found BED file: {:?}", bed_file_path);
    let mut input_bed = BufReader::new(File::open(bed_file_path.clone()).expect("Can't open BED file"));

    let correct_file_path = bed_file_path
        .parent()
        .unwrap()
        .join(bed_file_path.file_stem().unwrap())
        .to_str()
        .unwrap()
        .to_owned()
        + ".correct.bed";
    info!("Creating CB corrected BED file: {:?}", correct_file_path);
    let mut output_bed = BufWriter::new(File::create(correct_file_path).expect("Can't create BED file"));

    let mut mem_block = [0; 28];
    let mut num_lines = 0;
    let mut num_filtered = 0;
    while let Ok(frag) = Fragment::read(&mut input_bed, &mut mem_block) {
        num_lines += 1;
        if num_lines % crate::TMIL == 0 {
            print!("\rDone processing {}0M reads", num_lines / crate::TMIL);
            std::io::stdout().flush().expect("Can't flush output");
        }

        
        if wtl_barcodes.contains(&frag.cb) | wtl_barcodes.contains(&!frag.cb) {
            frag.write(&mut output_bed, "binary")?;
            num_filtered += 1;
        }
    }

    println!();
    info!("Filtered {} out of {} ({:.2}%)", (num_lines-num_filtered).to_formatted_string(&Locale::en), 
        (num_lines).to_formatted_string(&Locale::en),
        (num_lines-num_filtered) as f32 * 100.0 / num_lines as f32,
    );

    Ok(())
}

// pub fn correct(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
//     let mut fastq_file_paths = Vec::new();
//     if let Some(in_v) = sub_m.values_of("barcode") {
//         for in_file in in_v {
//             let in_path = Path::new(in_file)
//                 .canonicalize()
//                 .expect("can't find absolute path of input file");

//             fastq_file_paths.push(in_path);
//         }
//     }
//     info!("Found FASTQ files: {:?}", fastq_file_paths);

//     let wtl_barcodes: HashMap<u64, u64>;
//     {
//         info!("Importing Whitelist barcode");
//         let mut wtl_file = MultiGzDecoder::new(
//             File::open("./ext/737K-cratac-v1.txt.gz").expect("Unable to open file"),
//         );

//         let mut wtl_strings = String::new();
//         wtl_file
//             .read_to_string(&mut wtl_strings)
//             .expect("Unable to read the file");

//         let wtl_cb_ids: Vec<u64> = wtl_strings.split_terminator("\n")
//             .map(|cb_str| cb_string_to_u64(cb_str.as_bytes()).expect("can't convert string to barcode"))
//             .collect();
    
//         wtl_barcodes = libradicl::utils::generate_permitlist_map(&wtl_cb_ids, crate::CB_LENGTH)?;
//         for cb in wtl_cb_ids {
//             wtl_barcodes.insert(cb, cb);
//         }
//     }

//     //info!("{:?}", wtl_barcodes.len());

//     let mut reads_with_wtl_cb = 0;
//     let mut reads_with_correction = 0;
//     let mut total_reads = 0;
//     //let mut cb_freq_counter = HashMap::new();
//     for fastq_file_path in fastq_file_paths {
//         info!("Working on: {:?}", fastq_file_path);

//         let fq_file = fastq::Reader::new(MultiGzDecoder::new(
//             File::open(fastq_file_path).expect("can't open file"),
//         ));

//         for result in fq_file.records() {
//             total_reads += 1;
//             if total_reads % crate::MIL == 0 {
//                 print!(
//                     "\r Done processing {:?} Million reads",
//                     total_reads / crate::MIL
//                 );
//                 std::io::stdout().flush().expect("Can't flush output");
//             }

//             let record = result.expect("Error during fastq record parsing");
//             let cb: u64 = match crate::CB_ORIENT_FW {
//                 true => cb_string_to_u64(record.seq())?,
//                 false => cb_string_to_u64(&bio::alphabets::dna::revcomp(record.seq()))?,
//             };

//             if wtl_barcodes.contains(&cb) { 
//                 reads_with_wtl_cb += 1; 
//             } else {
//                 let mut neighbors = HashSet::new();
//                 crate::libradicl::utils::get_all_one_edit_neighbors(cb, crate::CB_LENGTH, &mut neighbors)?;
//                 for neighbor in neighbors {
//                     if wtl_barcodes.contains(&neighbor) {
//                         reads_with_correction += 1;
//                         break;
//                     }
//                 }
//             }
//         }
//         println!();
//     }

//     info!(
//         "Found barcodes in fwd = {} strand with Exact {}({:.02}%) and corrected {}({:.02}%) sequences.",
//         crate::CB_ORIENT_FW,
//         reads_with_wtl_cb,
//         (reads_with_wtl_cb as f32 * 100.0) / total_reads as f32,
//         reads_with_correction,
//         (reads_with_correction as f32 * 100.0) / total_reads as f32,
//     );
//     Ok(())
// }

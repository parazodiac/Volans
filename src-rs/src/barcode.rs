use flate2::read::MultiGzDecoder;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::path::Path;

use clap::ArgMatches;
use std::collections::HashSet;

use bio::io::fastq;

pub fn correct(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let mut fastq_file_paths = Vec::new();
    if let Some(in_v) = sub_m.values_of("barcode") {
        for in_file in in_v {
            let in_path = Path::new(in_file)
                .canonicalize()
                .expect("can't find absolute path of input file");

            fastq_file_paths.push(in_path);
        }
    }
    info!("Found FASTQ files: {:?}", fastq_file_paths);

    info!("Importing Whitelist barcode");
    let mut wtl_file = File::open("./ext/737K-cratac-v1.txt").expect("Unable to open file");

    let mut wtl_strings = String::new();
    wtl_file
        .read_to_string(&mut wtl_strings)
        .expect("Unable to read the file");

    let wtl_barcodes = wtl_strings
        .split("\n")
        .map(|x| x.as_bytes())
        .collect::<HashSet<&[u8]>>();

    let mut total_reads = 0;
    let mut found = 0;
    for fastq_file_path in fastq_file_paths {
        info!("Working on: {:?}", fastq_file_path);

        let fq_file = MultiGzDecoder::new(File::open(fastq_file_path).expect("can't open file"));

        let fq_file = fastq::Reader::new(fq_file);
        for result in fq_file.records() {
            total_reads += 1;
            let record = result.expect("Error during fastq record parsing");
            let seq = record.seq();
            if wtl_barcodes.contains(seq) {
                found += 1;
            }

            if total_reads % 1_000_000 == 0 {
                print!(
                    "\r Done processing {:?} Million reads",
                    total_reads / 1_000_000
                );
                std::io::stdout().flush().unwrap();
            }
        }
        println!();
    }

    info!(
        "Found: {}/{}=({:.02}%)",
        found,
        total_reads,
        (found as f32 * 100.0) / total_reads as f32
    );
    Ok(())
}

use std::error::Error;
use std::fs::File;
use std::path::PathBuf;

use crate::carina::fastq::FastqFeeder3;
use crate::fragments::schema::Fragment;
use bwa::BwaAligner;
use carina::barcode::cb_string_to_u64;
use clap::ArgMatches;
use rust_htslib::bam::Record;

pub fn process_reads(
    fq_feeder: FastqFeeder3<File>,
    bwa: BwaAligner,
    opath: PathBuf,
) -> Result<(), Box<dyn Error>> {
    // open the file to write the fragments
    //let file = crate::carina::file::bufreader_from_filepath(opath)?;

    let mut read_count = 0;
    let mut mm_read_count = 0;
    let mut orphan_read_count = 0;
    let mut chimeric_read_count = 0;

    let cb_extractor = |aln: &Record| -> u64 {
        let qname = aln.seq().as_bytes();
        cb_string_to_u64(&qname[(qname.len() - crate::configs::CB_LENGTH)..])
            .expect("can't convert cb string to u64")
    };

    for record in fq_feeder {
        read_count += 1;

        let (r1_alns, r2_alns) = bwa.align_read_pair(
            b"jdoe",
            record.0.seq(),
            record.0.qual(),
            record.2.seq(),
            record.2.qual(),
        );

        if r1_alns.len() > 1 || r2_alns.len() > 1 {
            mm_read_count += 1;
            continue;
        }

        if r1_alns.len() == 0 || r2_alns.len() == 0 {
            orphan_read_count += 1;
            continue;
        }

        if r1_alns[0].tid() != r2_alns[0].tid() {
            chimeric_read_count += 1;
            continue;
        }

        let frag = Fragment::new(&r1_alns[0], &r2_alns[0], cb_extractor);
        println!("{:?}", frag);

        if read_count % 10_000 == 0 {
            break;
        }
    }

    info!("Total reads: {}", read_count);
    info!("Total orphan reads: {}", orphan_read_count);
    info!("Total multimapping reads: {}", mm_read_count);
    info!("Total chimeric reads: {}", chimeric_read_count);

    Ok(())
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let fastq_one_path = carina::file::file_path_from_clap(sub_m, "one")?;
    let fastq_second_path = carina::file::file_path_from_clap(sub_m, "two")?;
    let fastq_third_path = carina::file::file_path_from_clap(sub_m, "three")?;

    let fastq_paths = vec![fastq_one_path, fastq_second_path, fastq_third_path.clone()];
    let fq_feeder: FastqFeeder3<File> = FastqFeeder3::<File>::new(fastq_paths);

    let index_path = carina::file::file_path_from_clap(sub_m, "index")?;
    let bwa = BwaAligner::from_path(index_path)?;

    //let output_path = carina::file::file_path_from_clap(sub_m, "obed")?;

    process_reads(fq_feeder, bwa, fastq_third_path)?;
    Ok(())
}

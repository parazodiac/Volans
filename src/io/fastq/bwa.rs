use std::error::Error;
use std::fs::File;

use crate::carina::fastq::FastqFeeder3;
use bwa::BwaAligner;
use clap::ArgMatches;

pub fn process_reads(fq_feeder: FastqFeeder3<File>, bwa: BwaAligner) -> Result<(), Box<dyn Error>> {
    let mut count = 0;
    for record in fq_feeder {
        count += 1;

        let (r1_alns, _r2_alns) = bwa.align_read_pair(
            b"dummy",
            record.0.seq(),
            record.0.qual(),
            record.2.seq(),
            record.2.qual(),
        );

        println!(
            "r1 mapping -- tid: {}, pos: {}",
            r1_alns[0].tid(),
            r1_alns[0].pos()
        );

        if count % 10_000 == 0 {
            break;
        }
    }

    println!("{}", count);
    Ok(())
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let fastq_one_path = carina::file::file_path_from_clap(sub_m, "one")?;
    let fastq_second_path = carina::file::file_path_from_clap(sub_m, "two")?;
    let fastq_third_path = carina::file::file_path_from_clap(sub_m, "three")?;

    let fastq_paths = vec![fastq_one_path, fastq_second_path, fastq_third_path];
    let fq_feeder: FastqFeeder3<File> = FastqFeeder3::<File>::new(fastq_paths);

    let index_path = carina::file::file_path_from_clap(sub_m, "index")?;
    let bwa = BwaAligner::from_path(index_path)?;

    process_reads(fq_feeder, bwa)?;
    Ok(())
}

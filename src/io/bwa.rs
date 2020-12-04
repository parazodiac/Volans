use std::error::Error;
use std::fs::File;
use std::io::Write;

use crate::carina::fastq::FastqFeeder3;
use crate::fragments::schema::Fragment;
use bwa::BwaAligner;
use carina::barcode::cb_string_to_u64;
use clap::ArgMatches;

use crate::fragments::count_stats;
use crate::fragments::filter;

pub fn process_reads(
    fq_feeder: FastqFeeder3<File>,
    bwa: BwaAligner,
    mut obed_file: std::io::BufWriter<std::fs::File>,
) -> Result<(), Box<dyn Error>> {
    let mut counter = count_stats::FragStats {
        ..Default::default()
    };

    for record in fq_feeder {
        counter.total_reads += 1;
        if counter.total_reads % crate::configs::TKILO == 0 {
            print!(
                "\rDone processing {}0 K reads",
                counter.total_reads / crate::configs::TKILO
            );
            std::io::stdout().flush().expect("Can't flush output");
        }

        let (mut r1_alns, mut r2_alns) = bwa.align_read_pair(
            b"jdoe",
            record.0.seq(),
            record.0.qual(),
            record.2.seq(),
            record.2.qual(),
        );
        assert!(
            r1_alns.len() == r2_alns.len(),
            "uneven number of alignments"
        );

        r1_alns.append(&mut r2_alns);
        match filter::callback(&r1_alns, &mut counter, false, false, 10_000) {
            Some((aln, maln)) => {
                let cb_name = record.1.seq();
                let cb = cb_string_to_u64(&cb_name[(cb_name.len() - crate::configs::CB_LENGTH)..])
                    .expect("can't convert cb string to u64");

                let frag = Fragment::new_with_cb(aln, maln, cb);
                frag.write(&mut obed_file, "text")?;
            }
            None => continue,
        };
    }

    println!("{}", counter);
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

    let obed_file = carina::file::bufwriter_from_clap(sub_m, "obed")?;

    process_reads(fq_feeder, bwa, obed_file)?;
    Ok(())
}

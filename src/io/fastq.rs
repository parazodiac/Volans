use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use crate::carina::fastq::FastqFeeder3;
use bio::io::fastq;
use clap::ArgMatches;

pub fn process_reads(
    fq_feeder: FastqFeeder3<File>,
    out_first_file: std::io::BufWriter<std::fs::File>,
    out_second_file: std::io::BufWriter<std::fs::File>,
) -> Result<(), Box<dyn Error>> {
    let mut total_reads = 0;
    let mut file_1 = fastq::Writer::new(out_first_file);
    let mut file_2 = fastq::Writer::new(out_second_file);

    for record in fq_feeder {
        total_reads += 1;
        if total_reads % crate::configs::TKILO == 0 {
            print!(
                "\rDone processing {}0 K reads",
                total_reads / crate::configs::TKILO
            );
            std::io::stdout().flush().expect("Can't flush output");
        }

        let seq_str = std::str::from_utf8(&record.2.seq()).unwrap();
        let new_read_one_name = format!("{}:{}", seq_str, record.0.id());
        let new_read_two_name = format!("{}:{}", seq_str, record.0.id());

        let new_r1 =
            fastq::Record::with_attrs(&new_read_one_name, None, record.0.seq(), record.0.qual());
        let new_r2 =
            fastq::Record::with_attrs(&new_read_two_name, None, record.1.seq(), record.1.qual());

        file_1.write_record(&new_r1)?;
        file_2.write_record(&new_r2)?;
    }

    println!("{}", total_reads);
    Ok(())
}

pub fn bufwriter_from_clap_with_name(
    sub_m: &ArgMatches,
    clap_id: &str,
    name: &str,
) -> Result<BufWriter<File>, Box<dyn Error>> {
    let file_str = sub_m
        .value_of(clap_id)
        .expect(&format!("can't find the flag: {}", clap_id));

    let file_path = Path::new(file_str)
        .canonicalize()
        .expect(&format!("Can't find absolute path of file tag {}", clap_id));

    let parent = file_path.parent().unwrap();
    let new_file = parent.join(name);

    let file = BufWriter::new(File::create(new_file)?);
    let file_path = Path::new(file_str)
        .canonicalize()
        .expect(&format!("Can't find absolute path of the created file"));

    info!("Created File at {:?}", file_path);

    Ok(file)
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let fastq_one_path = carina::file::file_path_from_clap(sub_m, "one")?;
    let fastq_second_path = carina::file::file_path_from_clap(sub_m, "two")?;
    let fastq_barcode_path = carina::file::file_path_from_clap(sub_m, "barcode")?;

    let fastq_paths = vec![fastq_one_path, fastq_second_path, fastq_barcode_path];
    let fq_feeder: FastqFeeder3<File> = FastqFeeder3::<File>::new(fastq_paths);

    let out_first_file = bufwriter_from_clap_with_name(sub_m, "barcode", "barcoded.1.fastq")?;
    let out_second_file = bufwriter_from_clap_with_name(sub_m, "barcode", "barcoded.2.fastq")?;

    process_reads(fq_feeder, out_first_file, out_second_file)?;
    Ok(())
}

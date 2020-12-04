extern crate bincode;
extern crate bio;
extern crate bitvector;
extern crate clap;
extern crate indicatif;
extern crate itertools;
extern crate num_format;
extern crate pretty_env_logger;
extern crate quickersort;
extern crate rust_htslib;
extern crate serde;
extern crate sprs;

extern crate bwa;
extern crate carina;

#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

pub mod configs;
pub mod fragments;
pub mod io;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("volans")
        .version("0.1.0")
        .author("Avi Srivastava, Tim Stuart, Bingjie Zhang, Rahul Satija")
        .about("A set of fast helper functions for Cut&Tag/ATAC data analysis.")
        .subcommand(
            SubCommand::with_name("bwa")
                .about("A subcommand to map the fastq reads using bwa")
                .arg(
                    Arg::with_name("one")
                        .short("1")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input 1 fastq file"),
                )
                .arg(
                    Arg::with_name("two")
                        .short("2")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input 2 fastq file"),
                )
                .arg(
                    Arg::with_name("three")
                        .short("3")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input 3 fastq file"),
                )
                .arg(
                    Arg::with_name("index")
                        .long("index")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the fasta file with the bwa index"),
                )
                .arg(
                    Arg::with_name("obed")
                        .long("obed")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output bed file"),
                ),
        )
        .subcommand(
            SubCommand::with_name("extract")
                .about("A subcommand to covert 3 fastq files to 2 fastq files ")
                .arg(
                    Arg::with_name("one")
                        .short("1")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input 1 fastq file"),
                )
                .arg(
                    Arg::with_name("two")
                        .short("2")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input 2 fastq file"),
                )
                .arg(
                    Arg::with_name("barcode")
                        .short("b")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input fastq file w/ barcodes"),
                ),
        )
        //.subcommand(
        //    SubCommand::with_name("filter")
        //        .about("A subcommand to filter BAM and generate BED.")
        //        .arg(
        //            Arg::with_name("ibam")
        //                .long("ibam")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BAM file"),
        //        )
        //        .arg(
        //            Arg::with_name("obed")
        //                .long("obed")
        //                .short("o")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the output bed file"),
        //        )
        //        .arg(
        //            Arg::with_name("tenx")
        //                .long("tenx")
        //                .help("use tag CB from 10x generated BAM."),
        //        )
        //        .arg(
        //            Arg::with_name("stats")
        //                .long("stats")
        //                .help("Don't write the output BED, just produce stats."),
        //        )
        //        .arg(
        //            Arg::with_name("mitostr")
        //                .long("mitostr")
        //                .short("m")
        //                .takes_value(true)
        //                .default_value("chrM")
        //                .help("String to identify the mitochondrial chromosome"),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("correct")
        //        .about("A subcommand to sequence correct the cb sequences.")
        //        .arg(
        //            Arg::with_name("whitelist")
        //                .long("whitelist")
        //                .short("w")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the list of known whitelist CB."),
        //        )
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file with CB sequences."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("sort")
        //        .about("A subcommand to sort the file by a chromosome names.")
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file with CB sequences."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("group")
        //        .about("A subcommand to group the file by (chr, start, end, CB)")
        //        .arg(
        //            Arg::with_name("allcb")
        //                .long("allcb")
        //                .help("report all cb instead of count."),
        //        )
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file with CB sequences."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("callpeak")
        //        .about("A subcommand to call peaks from a grouped bed file")
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the grouped BED file."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("count")
        //        .about("A subcommand to generate peak v cell count matrix")
        //        .arg(
        //            Arg::with_name("bam")
        //                .long("bam")
        //                .short("b")
        //                .takes_value(true)
        //                .help("path to the bam file with the chromosome names."),
        //        )
        //        .arg(
        //            Arg::with_name("pbed")
        //                .long("pbed")
        //                .short("p")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file with the annotated peaks."),
        //        )
        //        .arg(
        //            Arg::with_name("cbed")
        //                .long("cbed")
        //                .short("c")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file fragments and CB."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("text")
        //        .about("A subcommand to convert binary bed to text.")
        //        .arg(
        //            Arg::with_name("cbtext")
        //                .long("cbtext")
        //                .help("writes the last column as CB sequence."),
        //        )
        //        .arg(
        //            Arg::with_name("bam")
        //                .long("bam")
        //                .short("b")
        //                .takes_value(true)
        //                .help("path to the bam file with the chromosome names."),
        //        )
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file with CB sequences."),
        //        ),
        //)
        //.subcommand(
        //    SubCommand::with_name("stats")
        //        .about("A subcommand to summary stats of the binary bed.")
        //        .arg(
        //            Arg::with_name("ibed")
        //                .long("ibed")
        //                .short("i")
        //                .takes_value(true)
        //                .required(true)
        //                .help("path to the BED file."),
        //        ),
        //)
        .get_matches();
    pretty_env_logger::init_timed();

    if let Some(sub_m) = matches.subcommand_matches("bwa") {
        io::bwa::callback(&sub_m)?
    }
    if let Some(sub_m) = matches.subcommand_matches("extract") {
        io::fastq::callback(&sub_m)?
    }
    //if let Some(sub_m) = matches.subcommand_matches("filter") {
    //    quantify::filter::filter(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("correct") {
    //    quantify::barcode::correct(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("sort") {
    //    quantify::sort::sort(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("group") {
    //    quantify::group::dedup(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("callpeak") {
    //    quantify::peak::callpeak(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("count") {
    //    quantify::count::count(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("text") {
    //    quantify::text::convert(&sub_m)?
    //}
    //if let Some(sub_m) = matches.subcommand_matches("stats") {
    //    quantify::stats::stats(&sub_m)?
    //}

    Ok(())
}

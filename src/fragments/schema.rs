use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, BufWriter};

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

use crate::configs::{CB_LENGTH, TN5_LEFT_OFFSET, TN5_RIGHT_OFFSET};
use serde::{Deserialize, Serialize};

use carina::barcode::u64_to_cb_string;

#[derive(Debug, Clone, Copy)]
pub struct Feature {
    pub start: u32,
    pub end: u32,
    pub count: u32,
}

#[derive(Serialize, Deserialize, Hash, PartialEq, Eq, Debug)]
pub struct Fragment {
    pub chr: u32,
    pub start: u64,
    pub end: u64,
    pub cb: u64,
}

impl Fragment {
    pub fn new(aln: &Record, maln: &Record, extract_cb: fn(&Record) -> u64) -> Fragment {
        assert_eq!(aln.is_reverse(), false);

        let chr = aln.tid() as u32;
        let start = std::cmp::max(0, soft_clip_pos(aln) + TN5_LEFT_OFFSET);
        let end = std::cmp::max(0, soft_clip_pos(maln) - TN5_RIGHT_OFFSET);
        let cb_id = extract_cb(aln);

        Fragment {
            chr: chr as u32,
            start: start as u64,
            end: end as u64,
            cb: cb_id,
        }
    }

    pub fn new_with_cb(aln: &Record, maln: &Record, cb: u64) -> Fragment {
        assert_eq!(aln.is_reverse(), false);

        let chr = aln.tid() as u32;
        let start = std::cmp::max(0, soft_clip_pos(aln) + TN5_LEFT_OFFSET);
        let end = std::cmp::max(0, soft_clip_pos(maln) - TN5_RIGHT_OFFSET);
        let cb_id = cb;

        Fragment {
            chr: chr as u32,
            start: start as u64,
            end: end as u64,
            cb: cb_id,
        }
    }

    pub fn write(
        &self,
        mut file: &mut BufWriter<File>,
        write_mode: &str,
    ) -> Result<(), Box<dyn Error>> {
        match write_mode {
            "text" => writeln!(
                &mut file,
                "{}\t{}\t{}\t{}",
                self.chr, self.start, self.end, self.cb
            )?,
            "cb_text" => writeln!(
                &mut file,
                "{}\t{}\t{}\t{}",
                self.chr,
                self.start,
                self.end,
                u64_to_cb_string(self.cb, CB_LENGTH)?
            )?,
            "binary" => {
                let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
                file.write_all(&encoded)?
            }
            _ => unreachable!(),
        }

        Ok(())
    }

    pub fn write_with_name(
        &self,
        mut file: &mut BufWriter<File>,
        write_mode: &str,
        name: &str,
    ) -> Result<(), Box<dyn Error>> {
        match write_mode {
            "text" => writeln!(
                &mut file,
                "{}\t{}\t{}\t{}",
                name, self.start, self.end, self.cb
            )?,
            "cb_text" => writeln!(
                &mut file,
                "{}\t{}\t{}\t{}",
                name,
                self.start,
                self.end,
                u64_to_cb_string(self.cb, CB_LENGTH)?
            )?,
            "binary" => {
                let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
                file.write_all(&encoded)?
            }
            _ => unreachable!(),
        }

        Ok(())
    }

    pub fn read(
        file: &mut BufReader<File>,
        mem_block: &mut [u8; 28],
    ) -> Result<Fragment, Box<bincode::ErrorKind>> {
        file.read_exact(mem_block)?;
        bincode::deserialize(&mem_block[..])
    }

    pub fn start(&self) -> u64 {
        self.start
    }
}

pub struct FragmentFile {
    buf: BufReader<File>,
    mem_block: [u8; 28],
}

impl FragmentFile {
    pub fn new(buf: BufReader<File>) -> FragmentFile {
        FragmentFile {
            buf,
            mem_block: [0; 28],
        }
    }
}

impl Iterator for FragmentFile {
    type Item = Fragment;

    fn next(&mut self) -> Option<Self::Item> {
        match Fragment::read(&mut self.buf, &mut self.mem_block) {
            Ok(frag) => Some(frag),
            Err(_) => None,
        }
    }
}

pub fn soft_clip_pos(aln: &Record) -> i64 {
    let mut softclip_offset = 0;
    for cigar in aln.cigar().iter() {
        if let Cigar::SoftClip(val) = cigar {
            softclip_offset += val
        }
    }

    match aln.is_reverse() {
        true => aln.cigar().end_pos() + softclip_offset as i64, // reverse strand
        false => aln.pos() - softclip_offset as i64,            // forward strand
    }
}

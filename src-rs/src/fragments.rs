use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io::{BufReader, BufWriter};

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

use crate::{CB_LENGTH, TN5_LEFT_OFFSET, TN5_RIGHT_OFFSET};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Fragment {
    pub chr: u32,
    pub start: u64,
    pub end: u64,
    pub cb: u64,
}

impl std::fmt::Display for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}\t{}", self.chr, self.start, self.end, self.cb)
    }
}

impl Fragment {
    pub fn new(aln: &Record, maln: &Record) -> Fragment {
        assert_eq!(aln.is_reverse(), false);

        let chr = aln.tid() as u32;
        let start = soft_clip_pos(aln) + TN5_LEFT_OFFSET;
        let end = soft_clip_pos(maln) - TN5_RIGHT_OFFSET;

        let qname = aln.qname();
        let cb_id = cb_string_to_u64(&qname[(qname.len() - CB_LENGTH)..])
            .expect("can't convert cb string to u64");

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
        if write_mode == "text" {
            write!(
                &mut file,
                "{}\t{}\t{}\t{}\n",
                self.chr, self.start, self.end, self.cb
            )?;
        } else if write_mode == "binary" {
            let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
            file.write_all(&encoded)?;
        } else {
            unreachable!();
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
}

pub fn soft_clip_pos(aln: &Record) -> i64 {
    let mut softclip_offset = 0;
    for cigar in aln.cigar().iter() {
        match cigar {
            Cigar::SoftClip(val) => softclip_offset += val,
            _ => (),
        }
    }

    match aln.is_reverse() {
        true => aln.cigar().end_pos() + softclip_offset as i64, // reverse strand
        false => aln.pos() - softclip_offset as i64,            // forward strand
    }
}

pub fn cb_string_to_u64(cb_str: &[u8]) -> Result<u64, Box<dyn Error>> {
    let mut cb_id: u64 = 0;
    for (idx, nt) in cb_str.iter().rev().enumerate() {
        let offset = idx * 2;
        match nt {
            65 | 78 => (),              // A | N 00
            67 => cb_id |= 1 << offset, // C 01
            71 => cb_id |= 2 << offset, // G 10
            84 => cb_id |= 3 << offset, // T 11
            _ => panic!("unknown nucleotide"),
        };
    }

    Ok(cb_id)
}

#[cfg(test)]
mod tests {
    use crate::fragments::*;

    #[test]
    fn test_cb_string_to_u64() {
        let cb_id = cb_string_to_u64("A".repeat(16).as_bytes()).unwrap();
        assert_eq!(cb_id, 0);

        let cb_id = cb_string_to_u64("T".repeat(16).as_bytes()).unwrap();
        assert_eq!(cb_id, u32::MAX as u64);
    }
}

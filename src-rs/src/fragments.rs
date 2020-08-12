use std::error::Error;
use std::io::BufWriter;
use std::io::Write;

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

use crate::{CB_LENGTH, TN5_LEFT_OFFSET, TN5_RIGHT_OFFSET};

pub struct Fragment {
    chr: u32,
    start: u64,
    end: u64,
    cb: u64,
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

    fn _as_bytes(&self) -> Option<Vec<u8>> {
        let mut frag_bytes = Vec::new();

        frag_bytes.append(&mut self.chr.to_le_bytes().to_vec());
        frag_bytes.append(&mut self.start.to_le_bytes().to_vec());
        frag_bytes.append(&mut self.end.to_le_bytes().to_vec());
        frag_bytes.append(&mut self.cb.to_le_bytes().to_vec());

        Some(frag_bytes)
    }

    pub fn write<T: std::io::Write>(
        &self,
        mut file: &mut BufWriter<T>,
    ) -> Result<(), Box<dyn Error>> {
        write!(
            &mut file,
            "{}\t{}\t{}\t{}\n",
            self.chr, self.start, self.end, self.cb
        )?;

        Ok(())
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

use std::error::Error;
use std::io::BufWriter;
use std::io::Write;

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
    pub fn new(chr: u32, start: u64, end: u64, cb_str: &[u8]) -> Fragment {
        let cb = cb_string_to_u64(cb_str).unwrap();
        Fragment {
            chr,
            start,
            end,
            cb,
        }
    }

    fn as_bytes(&self) -> Option<Vec<u8>> {
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

pub fn cb_string_to_u64(cb_str: &[u8]) -> Result<u64, Box<dyn Error>> {
    let mut cb_id = 0;
    for (idx, nt) in cb_str.iter().rev().enumerate() {
        match nt {
            65 | 78 => (),                                                              // A | N 00
            67 => cb_id += 2_u64.pow(idx as u32 * 2),                                   // C 01
            71 => cb_id += 2_u64.pow((idx as u32 * 2) + 1),                             // G 10
            84 => cb_id += 2_u64.pow(idx as u32 * 2) + 2_u64.pow((idx as u32 * 2) + 1), // T 11
            _ => panic!("unknown nucleotide"),
        };
    }

    Ok(cb_id)
}

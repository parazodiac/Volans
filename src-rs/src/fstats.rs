use num_format::{Locale, ToFormattedString};

#[derive(Default)]
pub struct FragStats {
    pub total_reads: usize,
    pub mm_reads: usize,
    pub mapq_skip: usize,
    pub chimeric_tids: usize,
    pub chimeric_strand: usize,
    pub chimeric_max_distance: usize,
    pub chimeric_min_distance: usize,
    pub mito_skip: usize,
    pub unmap_skip: usize,
    pub unmap_orphan: usize,
}

impl FragStats {
    fn num_skipped(&self) -> usize {
        return self.mm_reads
            + self.mapq_skip
            + self.chimeric_tids
            + self.chimeric_strand
            + self.chimeric_max_distance
            + self.chimeric_min_distance
            + self.mito_skip
            + self.unmap_skip;
    }

    fn num_chimeric(&self) -> usize {
        return self.chimeric_tids
            + self.chimeric_strand
            + self.chimeric_max_distance
            + self.chimeric_min_distance;
    }

    fn percent_total(&self, num: usize) -> f32 {
        assert!(self.total_reads != 0);
        num as f32 * 100.0 / self.total_reads as f32
    }
}

impl std::fmt::Display for FragStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let total_skipped = self.num_skipped();
        let total_chimeric = self.num_chimeric();

        let mut stats: String = String::new();
        stats += &format!(
            "\n\n\nSTATS: Total Reads: {}\n",
            (self.total_reads).to_formatted_string(&Locale::en)
        );
        stats += &format!(
            "STATS: Total Unmapped skip: {}({:.02}%) with {}({:.02}%) HighQ orphan.\n",
            (self.unmap_skip).to_formatted_string(&Locale::en),
            self.percent_total(self.unmap_skip),
            (self.unmap_orphan).to_formatted_string(&Locale::en),
            self.percent_total(self.unmap_orphan)
        );
        stats += &format!(
            "STATS: Total MultiMapping Reads: {}({:.02}%).\n",
            (self.mm_reads).to_formatted_string(&Locale::en),
            self.percent_total(self.mm_reads)
        );
        stats += &format!(
            "STATS: Total Mitochondrial Reads: {}({:.02}%).\n",
            (self.mito_skip).to_formatted_string(&Locale::en),
            self.percent_total(self.mito_skip)
        );
        stats += &format!(
            "STATS: Total MAPQ skip: {}({:.02}%).\n",
            (self.mapq_skip).to_formatted_string(&Locale::en),
            self.percent_total(self.mapq_skip)
        );
        stats += &format!(
            "STATS: Total Chimeric Reads: {}({:.02}%)\n",
            (total_chimeric).to_formatted_string(&Locale::en),
            self.percent_total(total_chimeric)
        );
        stats += &format!(
            "STATS: Total Reads skipped: {}({:.02}%)\n",
            (total_skipped).to_formatted_string(&Locale::en),
            self.percent_total(total_skipped),
        );

        write!(f, "{}", stats)
    }
}

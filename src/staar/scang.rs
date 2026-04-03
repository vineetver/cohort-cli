use super::masks::{AnnotatedVariant, MaskGroup};

/// SCANG parameters: data-adaptive variable-size scanning windows.
///
/// Instead of a single fixed window size, SCANG tests a range of sizes
/// (measured in variant count) and reports all windows that pass a
/// significance threshold adjusted for the total number of windows tested.
pub struct ScangParams {
    pub lmin: usize,
    pub lmax: usize,
    pub step: usize,
}

impl Default for ScangParams {
    fn default() -> Self {
        Self { lmin: 40, lmax: 300, step: 10 }
    }
}

/// Build SCANG window groups across multiple variant-count sizes.
///
/// Returns `(window_length_variants, groups)` for each tested size.
/// Each group is scored independently via the standard STAAR pipeline.
///
/// Windows slide by 1 variant position for maximum resolution.
/// The caller is responsible for multiple testing correction across
/// the total number of windows tested.
pub fn build_scang_windows(
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    chromosome: &str,
    params: &ScangParams,
) -> Vec<(u32, Vec<MaskGroup>)> {
    let n = chrom_indices.len();
    if n < params.lmin {
        return Vec::new();
    }

    let lmax = params.lmax.min(n);
    let mut result = Vec::new();

    let mut wsize = params.lmin;
    while wsize <= lmax {
        let mut groups = Vec::new();

        for start in 0..=(n - wsize) {
            let window_indices: Vec<usize> =
                chrom_indices[start..start + wsize].to_vec();
            let first_pos = variants[window_indices[0]].position;
            let last_pos = variants[*window_indices.last().unwrap()].position;

            groups.push(MaskGroup {
                name: format!("scang_L{}_{chromosome}:{first_pos}-{last_pos}", wsize),
                chromosome: chromosome.to_string(),
                start: first_pos,
                end: last_pos,
                variant_indices: window_indices,
            });
        }

        result.push((wsize as u32, groups));
        wsize += params.step;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn var_at(pos: u32) -> AnnotatedVariant {
        AnnotatedVariant {
            chromosome: "1".into(),
            position: pos,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
            maf: 0.005,
            gene_name: String::new(),
            region_type: String::new(),
            consequence: String::new(),
            cadd_phred: 10.0,
            revel: 0.0,
            annotation_weights: vec![0.5; 11],
            is_cage_promoter: false,
            is_cage_enhancer: false,
            is_ccre_promoter: false,
            is_ccre_enhancer: false,
        }
    }

    #[test]
    fn scang_windows_basic() {
        let variants: Vec<AnnotatedVariant> = (0..100).map(|i| var_at(i * 100)).collect();
        let indices: Vec<usize> = (0..100).collect();
        let params = ScangParams { lmin: 10, lmax: 30, step: 10 };
        let result = build_scang_windows(&variants, &indices, "1", &params);

        // 3 window sizes: 10, 20, 30
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, 10);
        assert_eq!(result[1].0, 20);
        assert_eq!(result[2].0, 30);

        // Window counts: 100-10+1=91, 100-20+1=81, 100-30+1=71
        assert_eq!(result[0].1.len(), 91);
        assert_eq!(result[1].1.len(), 81);
        assert_eq!(result[2].1.len(), 71);
    }

    #[test]
    fn scang_too_few_variants() {
        let variants: Vec<AnnotatedVariant> = (0..5).map(|i| var_at(i * 100)).collect();
        let indices: Vec<usize> = (0..5).collect();
        let params = ScangParams::default(); // lmin=40
        let result = build_scang_windows(&variants, &indices, "1", &params);
        assert!(result.is_empty());
    }

}

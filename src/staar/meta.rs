use std::collections::HashMap;
use std::fs::File;

use arrow::array::{Array, AsArray, Float64Array, Int32Array, ListArray};
use faer::Mat;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use super::masks::MaskGroup;
use super::score;
use super::sumstats::StudyMeta;
use super::GeneResult;
use crate::db::DuckEngine;
use crate::error::FavorError;
use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, FunctionalAnnotation, MetaVariant,
    RegulatoryFlags,
};

/// Load ALL segments for a (study, chromosome) by reading the parquet file directly.
pub fn load_all_segments(
    study: &StudyHandle,
    chrom: &str,
) -> Result<HashMap<i32, SegmentCov>, FavorError> {
    let seg_path = study.path.join(format!("chromosome={chrom}/segments.parquet"));
    if !seg_path.exists() { return Ok(HashMap::new()); }

    let file = File::open(&seg_path)
        .map_err(|e| FavorError::Analysis(format!("Cannot open {}: {e}", seg_path.display())))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| FavorError::Analysis(format!("Bad parquet {}: {e}", seg_path.display())))?
        .build()
        .map_err(|e| FavorError::Analysis(format!("Parquet reader error {}: {e}", seg_path.display())))?;

    let mut result = HashMap::new();
    for batch in reader {
        let batch = batch.map_err(|e| FavorError::Analysis(format!("Parquet read error: {e}")))?;
        let n = batch.num_rows();

        let seg_ids = batch.column_by_name("segment_id")
            .ok_or_else(|| FavorError::Analysis("Missing segment_id column".into()))?
            .as_any().downcast_ref::<Int32Array>()
            .ok_or_else(|| FavorError::Analysis("segment_id is not Int32".into()))?;
        let positions_col = batch.column_by_name("positions")
            .ok_or_else(|| FavorError::Analysis("Missing positions column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("positions is not List".into()))?;
        let refs_col = batch.column_by_name("refs")
            .ok_or_else(|| FavorError::Analysis("Missing refs column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("refs is not List".into()))?;
        let alts_col = batch.column_by_name("alts")
            .ok_or_else(|| FavorError::Analysis("Missing alts column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("alts is not List".into()))?;
        let cov_col = batch.column_by_name("cov_lower")
            .ok_or_else(|| FavorError::Analysis("Missing cov_lower column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("cov_lower is not List".into()))?;

        for row in 0..n {
            let seg_id = seg_ids.value(row);

            let pos_arr = positions_col.value(row);
            let pos_vals = pos_arr.as_any().downcast_ref::<Int32Array>()
                .ok_or_else(|| FavorError::Analysis("positions inner is not Int32".into()))?;
            let positions: Vec<u32> = (0..pos_vals.len()).map(|i| pos_vals.value(i) as u32).collect();

            let refs_arr = refs_col.value(row);
            let refs_str = refs_arr.as_string::<i32>();
            let refs: Vec<String> = (0..refs_str.len()).map(|i| refs_str.value(i).to_string()).collect();

            let alts_arr = alts_col.value(row);
            let alts_str = alts_arr.as_string::<i32>();
            let alts: Vec<String> = (0..alts_str.len()).map(|i| alts_str.value(i).to_string()).collect();

            let cov_arr = cov_col.value(row);
            let cov_vals = cov_arr.as_any().downcast_ref::<Float64Array>()
                .ok_or_else(|| FavorError::Analysis("cov_lower inner is not Float64".into()))?;
            let cov_lower: Vec<f64> = (0..cov_vals.len()).map(|i| cov_vals.value(i)).collect();

            let mut pos_index: HashMap<u32, Vec<usize>> = HashMap::new();
            for (i, &p) in positions.iter().enumerate() {
                pos_index.entry(p).or_default().push(i);
            }

            result.insert(seg_id, SegmentCov { refs, alts, cov_lower, pos_index });
        }
    }
    Ok(result)
}

pub struct StudyHandle {
    pub path: std::path::PathBuf,
    pub meta: StudyMeta,
}

/// Segment covariance read from one study.
pub struct SegmentCov {
    refs: Vec<String>,
    alts: Vec<String>,
    cov_lower: Vec<f64>,
    pos_index: HashMap<u32, Vec<usize>>,
}

impl SegmentCov {
    /// Extract sub-matrix for a subset of variants identified by (position, ref, alt).
    fn extract_submatrix(&self, keys: &[(u32, &str, &str)]) -> Mat<f64> {
        let m = keys.len();
        let mut mat = Mat::zeros(m, m);

        let mut key_to_local: Vec<Option<usize>> = Vec::with_capacity(m);
        for &(pos, ref_a, alt_a) in keys {
            let local = self.pos_index.get(&pos).and_then(|candidates| {
                candidates.iter().find(|&&i| self.refs[i] == ref_a && self.alts[i] == alt_a).copied()
            });
            key_to_local.push(local);
        }

        for i in 0..m {
            let Some(li) = key_to_local[i] else { continue };
            for j in 0..=i {
                let Some(lj) = key_to_local[j] else { continue };
                let (row, col) = if li >= lj { (li, lj) } else { (lj, li) };
                let idx = row * (row + 1) / 2 + col;
                if idx < self.cov_lower.len() {
                    mat[(i, j)] = self.cov_lower[idx];
                    mat[(j, i)] = self.cov_lower[idx];
                }
            }
        }
        mat
    }
}

/// Load and validate study directories.
pub fn load_studies(paths: &[std::path::PathBuf]) -> Result<Vec<StudyHandle>, FavorError> {
    let mut studies = Vec::with_capacity(paths.len());
    for path in paths {
        let meta_path = path.join("meta_staar.json");
        if !meta_path.exists() {
            return Err(FavorError::Input(format!(
                "Not a MetaSTAAR study directory: {}. Missing meta_staar.json. \
                 Run `favor staar --emit-sumstats` first.", path.display()
            )));
        }
        let content = std::fs::read_to_string(&meta_path)?;
        let meta: StudyMeta = serde_json::from_str(&content)
            .map_err(|e| FavorError::Input(format!("Invalid meta_staar.json in {}: {e}", path.display())))?;
        if meta.favor_meta_version != 1 {
            return Err(FavorError::Input(format!(
                "Unsupported meta version {} in {}. Expected 1.", meta.favor_meta_version, path.display()
            )));
        }
        studies.push(StudyHandle { path: path.clone(), meta });
    }

    if studies.len() < 2 {
        return Err(FavorError::Input("MetaSTAAR requires at least 2 studies.".into()));
    }

    let first_type = &studies[0].meta.trait_type;
    let first_seg = studies[0].meta.segment_size;
    for s in &studies[1..] {
        if s.meta.trait_type != *first_type {
            return Err(FavorError::Input(format!(
                "Trait type mismatch: {} vs {}. All studies must have the same trait type.",
                first_type, s.meta.trait_type
            )));
        }
        if s.meta.segment_size != first_seg {
            return Err(FavorError::Input(format!(
                "Segment size mismatch: {} vs {}. All studies must use the same segment size.",
                first_seg, s.meta.segment_size
            )));
        }
    }

    Ok(studies)
}

/// Build union variant catalog for one chromosome across all studies.
pub fn merge_chromosome(
    engine: &DuckEngine,
    studies: &[StudyHandle],
    chrom: &str,
    maf_cutoff: f64,
) -> Result<Vec<MetaVariant>, FavorError> {
    let mut union_parts = Vec::new();
    for (idx, study) in studies.iter().enumerate() {
        let var_path = study.path.join(format!("chromosome={chrom}/variants.parquet"));
        if !var_path.exists() { continue; }
        union_parts.push(format!(
            "SELECT {idx} AS study_idx, * FROM read_parquet('{}')",
            var_path.display()
        ));
    }
    if union_parts.is_empty() { return Ok(Vec::new()); }

    engine.execute(&format!(
        "CREATE OR REPLACE TEMP TABLE _study_variants AS {}",
        union_parts.join(" UNION ALL ")
    ))?;

    engine.execute(&format!(
        "CREATE OR REPLACE TEMP TABLE _meta_variants AS \
         SELECT \
             position, ref_allele, alt_allele, \
             SUM(u_stat) AS u_meta, \
             SUM(mac) AS mac_total, \
             SUM(n_obs) AS n_total, \
             FIRST(gene_name) FILTER (WHERE gene_name != '') AS gene_name, \
             FIRST(region_type) FILTER (WHERE region_type != '') AS region_type, \
             FIRST(consequence) FILTER (WHERE consequence != '') AS consequence, \
             FIRST(cadd_phred) AS cadd_phred, \
             FIRST(revel) AS revel, \
             BOOL_OR(is_cage_promoter) AS is_cage_promoter, \
             BOOL_OR(is_cage_enhancer) AS is_cage_enhancer, \
             BOOL_OR(is_ccre_promoter) AS is_ccre_promoter, \
             BOOL_OR(is_ccre_enhancer) AS is_ccre_enhancer, \
             FIRST(w_cadd) AS w_cadd, FIRST(w_linsight) AS w_linsight, \
             FIRST(w_fathmm_xf) AS w_fathmm_xf, \
             FIRST(w_apc_epi_active) AS w_apc_epi_active, \
             FIRST(w_apc_epi_repressed) AS w_apc_epi_repressed, \
             FIRST(w_apc_epi_transcription) AS w_apc_epi_transcription, \
             FIRST(w_apc_conservation) AS w_apc_conservation, \
             FIRST(w_apc_protein_function) AS w_apc_protein_function, \
             FIRST(w_apc_local_nd) AS w_apc_local_nd, \
             FIRST(w_apc_mutation_density) AS w_apc_mutation_density, \
             FIRST(w_apc_tf) AS w_apc_tf, \
             list(struct_pack(s := study_idx, seg := segment_id)) AS study_segs \
         FROM _study_variants \
         WHERE maf < {maf_cutoff} \
         GROUP BY position, ref_allele, alt_allele \
         ORDER BY position"
    ))?;

    let conn = engine.connection();
    let mut stmt = conn.prepare(
        "SELECT position, ref_allele, alt_allele, u_meta, \
         mac_total, n_total, gene_name, region_type, consequence, \
         cadd_phred, revel, \
         is_cage_promoter, is_cage_enhancer, is_ccre_promoter, is_ccre_enhancer, \
         w_cadd, w_linsight, w_fathmm_xf, \
         w_apc_epi_active, w_apc_epi_repressed, w_apc_epi_transcription, \
         w_apc_conservation, w_apc_protein_function, w_apc_local_nd, \
         w_apc_mutation_density, w_apc_tf, \
         study_segs \
         FROM _meta_variants ORDER BY position"
    ).map_err(|e| FavorError::Analysis(format!("{e}")))?;

    let mut rows = stmt.query([]).map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut result = Vec::new();

    let chrom_parsed: Chromosome = chrom.parse().unwrap_or(Chromosome::Autosome(1));

    while let Ok(Some(row)) = rows.next() {
        let mut weights = [0.0f64; 11];
        for (i, w) in weights.iter_mut().enumerate() {
            *w = row.get::<_, f64>(15 + i).unwrap_or(0.0);
        }

        let segs_str: String = row.get(26).unwrap_or_default();
        let study_segments = parse_study_segments(&segs_str);

        let mac_total: i64 = row.get(4).unwrap_or(0);
        let n_total: i64 = row.get(5).unwrap_or(0);
        let maf = if n_total > 0 { mac_total as f64 / n_total as f64 } else { 0.0 };

        result.push(MetaVariant {
            variant: AnnotatedVariant {
                chromosome: chrom_parsed,
                position: row.get::<_, i32>(0).unwrap_or(0) as u32,
                ref_allele: row.get(1).unwrap_or_default(),
                alt_allele: row.get(2).unwrap_or_default(),
                maf,
                gene_name: row.get(6).unwrap_or_default(),
                annotation: FunctionalAnnotation {
                    region_type: row.get(7).unwrap_or_default(),
                    consequence: row.get(8).unwrap_or_default(),
                    cadd_phred: row.get(9).unwrap_or(0.0),
                    revel: row.get(10).unwrap_or(0.0),
                    regulatory: RegulatoryFlags {
                        cage_promoter: row.get::<_, bool>(11).unwrap_or(false),
                        cage_enhancer: row.get::<_, bool>(12).unwrap_or(false),
                        ccre_promoter: row.get::<_, bool>(13).unwrap_or(false),
                        ccre_enhancer: row.get::<_, bool>(14).unwrap_or(false),
                    },
                    weights: AnnotationWeights(weights),
                },
            },
            u_meta: row.get(3).unwrap_or(0.0),
            mac_total,
            n_total,
            study_segments,
        });
    }

    engine.execute("DROP TABLE IF EXISTS _study_variants")?;
    engine.execute("DROP TABLE IF EXISTS _meta_variants")?;

    Ok(result)
}

/// Run meta-analysis score tests for one gene/mask group.
pub fn meta_score_gene(
    group: &MaskGroup,
    meta_variants: &[MetaVariant],
    studies: &[StudyHandle],
    _chrom: &str,
    segment_cache: &HashMap<(usize, i32), SegmentCov>,
) -> Option<GeneResult> {
    let indices: Vec<usize> = group.variant_indices.iter()
        .filter(|&&i| i < meta_variants.len())
        .copied()
        .collect();
    if indices.len() < 2 { return None; }

    let m = indices.len();

    let mut u = Mat::zeros(m, 1);
    for (local, &gi) in indices.iter().enumerate() {
        u[(local, 0)] = meta_variants[gi].u_meta;
    }

    let keys: Vec<(u32, &str, &str)> = indices.iter()
        .map(|&gi| (meta_variants[gi].variant.position, meta_variants[gi].variant.ref_allele.as_str(), meta_variants[gi].variant.alt_allele.as_str()))
        .collect();

    let mut cov = Mat::zeros(m, m);
    for study_idx in 0..studies.len() {
        let mut needed_segments: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for &gi in &indices {
            for &(sidx, seg_id) in &meta_variants[gi].study_segments {
                if sidx == study_idx {
                    needed_segments.insert(seg_id);
                }
            }
        }

        for seg_id in needed_segments {
            let cache_key = (study_idx, seg_id);
            if let Some(seg) = segment_cache.get(&cache_key) {
                let sub = seg.extract_submatrix(&keys);
                for i in 0..m {
                    for j in 0..m {
                        cov[(i, j)] += sub[(i, j)];
                    }
                }
            }
        }
    }

    let mafs: Vec<f64> = indices.iter().map(|&gi| {
        let mv = &meta_variants[gi];
        if mv.n_total > 0 { mv.mac_total as f64 / mv.n_total as f64 } else { 0.0 }
    }).collect();

    let ann_matrix: Vec<Vec<f64>> = (0..11).map(|ch| {
        indices.iter().map(|&gi| {
            meta_variants[gi].variant.annotation.weights.0[ch]
        }).collect()
    }).collect();

    let sr = score::run_staar_from_sumstats(&u, &cov, &ann_matrix, &mafs);

    let cmac: i64 = indices.iter().map(|&gi| meta_variants[gi].mac_total).sum();

    Some(GeneResult {
        ensembl_id: group.name.clone(),
        gene_symbol: group.name.clone(),
        chromosome: group.chromosome.clone(),
        start: group.start,
        end: group.end,
        n_variants: m as u32,
        cumulative_mac: cmac as u32,
        staar: sr,
    })
}

fn parse_study_segments(s: &str) -> Vec<(usize, i32)> {
    let mut result = Vec::new();
    for part in s.split('{') {
        let part = part.trim_matches(|c: char| c == '[' || c == ']' || c == ',' || c == ' ' || c == '}');
        if part.is_empty() { continue; }
        let mut study = None;
        let mut seg = None;
        for kv in part.split(',') {
            let kv = kv.trim();
            if let Some(val) = kv.strip_prefix("s:").or_else(|| kv.strip_prefix("s :")) {
                study = val.trim().parse().ok();
            } else if let Some(val) = kv.strip_prefix("seg:").or_else(|| kv.strip_prefix("seg :")) {
                seg = val.trim().parse().ok();
            }
        }
        if let (Some(s), Some(g)) = (study, seg) {
            result.push((s, g));
        }
    }
    result
}

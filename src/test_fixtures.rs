use std::path::Path;

use crate::staar::masks::AnnotatedVariant;

// Ground truth annotation weights from real OR11H1 stopgain on chr22.
// Order matches staar::weights::ANNOTATION_CHANNELS exactly.
const GROUND_TRUTH_WEIGHTS: [f64; 11] = [
    0.9957,  // cadd: 1 - 10^(-23.7/10)
    0.2149,  // linsight
    1.7630,  // fathmm_xf
    0.0126,  // apc_epigenetics_active
    0.3103,  // apc_epigenetics_repressed
    0.3234,  // apc_epigenetics_transcription
    6.3541,  // apc_conservation_v2
    2.9695,  // apc_protein_function_v3
    13.3931, // apc_local_nucleotide_diversity_v3
    13.3340, // apc_mutation_density
    3.1428,  // apc_transcription_factor
];

pub fn base_variant() -> AnnotatedVariant {
    AnnotatedVariant {
        chromosome: "22".into(),
        position: 15528164,
        ref_allele: "C".into(),
        alt_allele: "T".into(),
        maf: 0.0007,
        gene_name: "OR11H1".into(),
        region_type: "exonic".into(),
        consequence: "stopgain".into(),
        cadd_phred: 23.7,
        revel: 0.0,
        annotation_weights: GROUND_TRUTH_WEIGHTS.to_vec(),
        is_cage_promoter: false,
        is_cage_enhancer: false,
        is_ccre_promoter: false,
        is_ccre_enhancer: false,
    }
}

pub fn stopgain() -> AnnotatedVariant { base_variant() }

pub fn splice() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528200,
        region_type: "splicing".into(),
        consequence: String::new(),
        cadd_phred: 25.1,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn missense_high() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528300,
        consequence: "nonsynonymous SNV".into(),
        cadd_phred: 28.0,
        revel: 0.8,
        ..base_variant()
    }
}

pub fn missense_low() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528400,
        consequence: "nonsynonymous SNV".into(),
        cadd_phred: 8.0,
        revel: 0.2,
        ..base_variant()
    }
}

pub fn synonymous() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528500,
        consequence: "synonymous SNV".into(),
        cadd_phred: 6.4,
        ..base_variant()
    }
}

pub fn frameshift() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528600,
        consequence: "frameshift deletion".into(),
        cadd_phred: 20.9,
        ..base_variant()
    }
}

pub fn stoploss() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528700,
        consequence: "stoploss".into(),
        cadd_phred: 15.0,
        ..base_variant()
    }
}

pub fn upstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15527000,
        region_type: "upstream".into(),
        consequence: String::new(),
        cadd_phred: 5.6,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn downstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15530000,
        region_type: "downstream".into(),
        consequence: String::new(),
        cadd_phred: 4.2,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn upstream_downstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15531000,
        region_type: "upstream;downstream".into(),
        consequence: String::new(),
        cadd_phred: 3.0,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn utr3() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15529000,
        region_type: "UTR3".into(),
        consequence: String::new(),
        cadd_phred: 7.4,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn cage_promoter() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15526000,
        region_type: "intergenic".into(),
        consequence: String::new(),
        cadd_phred: 12.3,
        revel: 0.0,
        is_cage_promoter: true,
        ..base_variant()
    }
}

pub fn cage_enhancer() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15525000,
        region_type: "intergenic".into(),
        consequence: String::new(),
        cadd_phred: 8.1,
        revel: 0.0,
        is_cage_enhancer: true,
        ..base_variant()
    }
}

pub fn ccre_pls() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15524000,
        region_type: "intergenic".into(),
        consequence: String::new(),
        cadd_phred: 7.0,
        revel: 0.0,
        is_ccre_promoter: true,
        ..base_variant()
    }
}

pub fn ccre_els() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15523000,
        region_type: "intergenic".into(),
        consequence: String::new(),
        cadd_phred: 8.7,
        revel: 0.0,
        is_ccre_enhancer: true,
        ..base_variant()
    }
}

pub fn ncrna() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15532000,
        gene_name: "RF00004".into(),
        region_type: "ncRNA_exonic".into(),
        consequence: String::new(),
        cadd_phred: 9.9,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn intronic() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528800,
        gene_name: "POTEH".into(),
        region_type: "intronic".into(),
        consequence: String::new(),
        cadd_phred: 13.9,
        revel: 0.0,
        ..base_variant()
    }
}

pub fn all_variants() -> Vec<AnnotatedVariant> {
    vec![
        stopgain(), splice(), missense_high(), missense_low(),
        synonymous(), frameshift(), stoploss(), upstream(),
        downstream(), upstream_downstream(), utr3(),
        cage_promoter(), cage_enhancer(), ccre_pls(), ccre_els(),
        ncrna(), intronic(),
    ]
}

/// SQL that creates a FAVOR-schema-correct annotation parquet.
/// Struct nesting matches the real full-tier parquet exactly for every column STAAR touches.
fn annotation_rows_sql() -> String {
    let rows: Vec<String> = all_variants().iter().map(|v| {
        let cage: &str = if v.is_cage_promoter && v.is_cage_enhancer {
            "{'cage_enhancer': 'e1@chr22', 'cage_promoter': 'p1@chr22', 'cage_tc': NULL}"
        } else if v.is_cage_promoter {
            "{'cage_enhancer': NULL, 'cage_promoter': 'p1@chr22', 'cage_tc': NULL}"
        } else if v.is_cage_enhancer {
            "{'cage_enhancer': 'e1@chr22', 'cage_promoter': NULL, 'cage_tc': NULL}"
        } else {
            "NULL::STRUCT(cage_enhancer VARCHAR, cage_promoter VARCHAR, cage_tc VARCHAR)"
        };

        let ccre_ann = if v.is_ccre_promoter { "'PLS'" }
            else if v.is_ccre_enhancer { "'dELS'" }
            else { "NULL" };

        let consequence_sql = if v.consequence.is_empty() { "NULL".into() }
            else { format!("'{}'", v.consequence) };

        let revel_sql = if v.revel == 0.0 { "NULL".into() }
            else { format!("{}", v.revel) };

        format!(
            "SELECT \
                22::BIGINT AS chromosome, \
                {pos} AS position, \
                'T' AS ref_vcf, \
                'A' AS alt_vcf, \
                {{'region_type': '{rt}', 'genes': ['{gene}']::VARCHAR[], \
                  'consequence': {csq}, \
                  'transcripts': NULL::STRUCT(gene VARCHAR, transcript_id VARCHAR, \
                      \"location\" VARCHAR, hgvsc VARCHAR, hgvsp VARCHAR)[]}} AS gencode, \
                {{'cadd': {{'raw': 0.0, 'phred': {cadd}::FLOAT}}}} AS main, \
                {{'revel': {revel}::FLOAT}} AS dbnsfp, \
                {linsight}::FLOAT AS linsight, \
                {fathmm}::FLOAT AS fathmm_xf, \
                {cage} AS cage, \
                {{'ids': NULL::VARCHAR, 'accessions': NULL::VARCHAR, \
                  'annotations': {ccre_ann}::VARCHAR, 'count': 0::UTINYINT}} AS ccre, \
                {{'conservation_v2': {w6}::FLOAT, 'epigenetics': 0.0::FLOAT, \
                  'epigenetics_active': {w3}::FLOAT, 'epigenetics_repressed': {w4}::FLOAT, \
                  'epigenetics_transcription': {w5}::FLOAT, \
                  'local_nucleotide_diversity_v3': {w8}::FLOAT, \
                  'mappability': 0.0::FLOAT, 'micro_rna': 0.0::FLOAT, \
                  'mutation_density': {w9}::FLOAT, \
                  'protein_function_v3': {w7}::FLOAT, \
                  'proximity_to_coding_v2': 0.0::FLOAT, \
                  'proximity_to_tsstes': 0.0::FLOAT, \
                  'transcription_factor': {w10}::FLOAT}} AS apc",
            pos = v.position,
            rt = v.region_type,
            gene = v.gene_name,
            csq = consequence_sql,
            cadd = v.cadd_phred,
            revel = revel_sql,
            linsight = v.annotation_weights[1],
            fathmm = v.annotation_weights[2],
            cage = cage,
            ccre_ann = ccre_ann,
            w3 = v.annotation_weights[3],
            w4 = v.annotation_weights[4],
            w5 = v.annotation_weights[5],
            w6 = v.annotation_weights[6],
            w7 = v.annotation_weights[7],
            w8 = v.annotation_weights[8],
            w9 = v.annotation_weights[9],
            w10 = v.annotation_weights[10],
        )
    }).collect();

    rows.join(" UNION ALL ")
}

pub fn write_test_annotation_parquet(dir: &Path) -> std::path::PathBuf {
    let out = dir.join("test_annotations.parquet");
    let db = duckdb::Connection::open_in_memory().unwrap();
    let sql = format!(
        "COPY ({}) TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)",
        annotation_rows_sql(),
        out.display(),
    );
    db.execute_batch(&sql).unwrap();
    out
}

pub fn write_test_genotype_parquet(dir: &Path) -> std::path::PathBuf {
    let out = dir.join("test_genotypes.parquet");
    let db = duckdb::Connection::open_in_memory().unwrap();

    let rows: Vec<String> = all_variants().iter().map(|v| {
        format!(
            "SELECT '22' AS chromosome, {} AS position, 'T' AS ref, 'A' AS alt, \
             {}::FLOAT AS maf, [0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0]::FLOAT[] AS dosages",
            v.position, v.maf,
        )
    }).collect();

    let sql = format!(
        "COPY ({}) TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)",
        rows.join(" UNION ALL "),
        out.display(),
    );
    db.execute_batch(&sql).unwrap();
    out
}

/// The exact SQL STAAR uses to build the _rare table, parameterized for test paths.
pub fn staar_rare_sql(geno_path: &str, ann_path: &str) -> String {
    format!(
        "CREATE TEMP TABLE _rare AS
         SELECT
             g.chromosome::VARCHAR AS chrom,
             g.position AS pos,
             g.ref AS ref_allele,
             g.alt AS alt_allele,
             g.maf,
             COALESCE(a.gencode.genes[1], '') AS gene_name,
             COALESCE(a.gencode.region_type, '') AS region_type,
             COALESCE(a.gencode.consequence, '') AS consequence,
             COALESCE(a.main.cadd.phred, 0) AS cadd_phred,
             COALESCE(a.dbnsfp.revel, 0) AS revel,
             a.cage.cage_promoter IS NOT NULL AS is_cage_promoter,
             a.cage.cage_enhancer IS NOT NULL AS is_cage_enhancer,
             COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%PLS%', false) AS is_ccre_promoter,
             COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%ELS%', false) AS is_ccre_enhancer,
             CASE WHEN a.main.cadd.phred > 0 THEN 1.0 - POW(10.0, -a.main.cadd.phred / 10.0) ELSE 0.0 END AS w_cadd,
             COALESCE(a.linsight, 0)::DOUBLE AS w_linsight,
             COALESCE(a.fathmm_xf, 0)::DOUBLE AS w_fathmm_xf,
             COALESCE(a.apc.epigenetics_active, 0)::DOUBLE AS w_apc_epi_active,
             COALESCE(a.apc.epigenetics_repressed, 0)::DOUBLE AS w_apc_epi_repressed,
             COALESCE(a.apc.epigenetics_transcription, 0)::DOUBLE AS w_apc_epi_transcription,
             COALESCE(a.apc.conservation_v2, 0)::DOUBLE AS w_apc_conservation,
             COALESCE(a.apc.protein_function_v3, 0)::DOUBLE AS w_apc_protein_function,
             COALESCE(a.apc.local_nucleotide_diversity_v3, 0)::DOUBLE AS w_apc_local_nd,
             COALESCE(a.apc.mutation_density, 0)::DOUBLE AS w_apc_mutation_density,
             COALESCE(a.apc.transcription_factor, 0)::DOUBLE AS w_apc_tf
         FROM read_parquet('{geno}') g
         INNER JOIN read_parquet('{ann}') a
             ON g.chromosome::VARCHAR = a.chromosome::VARCHAR
             AND g.position = a.position AND g.ref = a.ref_vcf AND g.alt = a.alt_vcf
         WHERE g.maf > 0 AND g.maf < 0.01",
        geno = geno_path, ann = ann_path,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::masks;
    use crate::staar::MaskType;

    fn read_extracted_variants(db: &duckdb::Connection) -> Vec<AnnotatedVariant> {
        let mut stmt = db.prepare(
            "SELECT chrom, pos, ref_allele, alt_allele, maf, gene_name, region_type, consequence, \
             cadd_phred, revel, \
             is_cage_promoter, is_cage_enhancer, is_ccre_promoter, is_ccre_enhancer, \
             w_cadd, w_linsight, w_fathmm_xf, \
             w_apc_epi_active, w_apc_epi_repressed, w_apc_epi_transcription, \
             w_apc_conservation, w_apc_protein_function, w_apc_local_nd, \
             w_apc_mutation_density, w_apc_tf \
             FROM _rare ORDER BY pos"
        ).unwrap();
        let mut rows = stmt.query([]).unwrap();
        let mut variants = Vec::new();
        while let Ok(Some(row)) = rows.next() {
            let weights: Vec<f64> = (14..25).map(|c| row.get::<_, f64>(c).unwrap_or(0.0)).collect();
            variants.push(AnnotatedVariant {
                chromosome: row.get(0).unwrap_or_default(),
                position: row.get::<_, i32>(1).unwrap_or(0) as u32,
                ref_allele: row.get(2).unwrap_or_default(),
                alt_allele: row.get(3).unwrap_or_default(),
                maf: row.get(4).unwrap_or(0.0),
                gene_name: row.get(5).unwrap_or_default(),
                region_type: row.get(6).unwrap_or_default(),
                consequence: row.get(7).unwrap_or_default(),
                cadd_phred: row.get(8).unwrap_or(0.0),
                revel: row.get(9).unwrap_or(0.0),
                annotation_weights: weights,
                is_cage_promoter: row.get::<_, bool>(10).unwrap_or(false),
                is_cage_enhancer: row.get::<_, bool>(11).unwrap_or(false),
                is_ccre_promoter: row.get::<_, bool>(12).unwrap_or(false),
                is_ccre_enhancer: row.get::<_, bool>(13).unwrap_or(false),
            });
        }
        variants
    }

    fn setup_rare_table() -> (tempfile::TempDir, duckdb::Connection, Vec<AnnotatedVariant>) {
        let dir = tempfile::tempdir().unwrap();
        let ann = write_test_annotation_parquet(dir.path());
        let geno = write_test_genotype_parquet(dir.path());
        let db = duckdb::Connection::open_in_memory().unwrap();
        db.execute_batch(&staar_rare_sql(
            &geno.to_string_lossy(), &ann.to_string_lossy(),
        )).unwrap();
        let variants = read_extracted_variants(&db);
        (dir, db, variants)
    }

    #[test]
    fn extraction_preserves_all_variants() {
        let (_dir, _db, variants) = setup_rare_table();
        assert_eq!(variants.len(), all_variants().len());
    }

    #[test]
    fn extraction_stopgain_fields() {
        let (_dir, _db, variants) = setup_rare_table();
        let v = variants.iter().find(|v| v.position == 15528164).unwrap();
        assert_eq!(v.gene_name, "OR11H1");
        assert_eq!(v.region_type, "exonic");
        assert_eq!(v.consequence, "stopgain");
        assert!((v.cadd_phred - 23.7).abs() < 0.1);
        assert!(!v.is_cage_promoter);
        assert!(!v.is_cage_enhancer);
        assert!(!v.is_ccre_promoter);
        assert!(!v.is_ccre_enhancer);
    }

    #[test]
    fn extraction_splice_has_empty_consequence() {
        let (_dir, _db, variants) = setup_rare_table();
        let v = variants.iter().find(|v| v.region_type == "splicing").unwrap();
        assert_eq!(v.consequence, "");
        assert!((v.cadd_phred - 25.1).abs() < 0.1);
    }

    #[test]
    fn extraction_cage_flags() {
        let (_dir, _db, variants) = setup_rare_table();
        let prom = variants.iter().find(|v| v.position == 15526000).unwrap();
        assert!(prom.is_cage_promoter);
        assert!(!prom.is_cage_enhancer);

        let enh = variants.iter().find(|v| v.position == 15525000).unwrap();
        assert!(!enh.is_cage_promoter);
        assert!(enh.is_cage_enhancer);
    }

    #[test]
    fn extraction_ccre_flags() {
        let (_dir, _db, variants) = setup_rare_table();
        let pls = variants.iter().find(|v| v.position == 15524000).unwrap();
        assert!(pls.is_ccre_promoter);
        assert!(!pls.is_ccre_enhancer);

        let els = variants.iter().find(|v| v.position == 15523000).unwrap();
        assert!(!els.is_ccre_promoter);
        assert!(els.is_ccre_enhancer);
    }

    #[test]
    fn extraction_null_cage_is_false() {
        let (_dir, _db, variants) = setup_rare_table();
        let intronic = variants.iter().find(|v| v.region_type == "intronic").unwrap();
        assert!(!intronic.is_cage_promoter);
        assert!(!intronic.is_cage_enhancer);
    }

    #[test]
    fn extraction_annotation_weights_nonzero() {
        let (_dir, _db, variants) = setup_rare_table();
        let v = variants.iter().find(|v| v.consequence == "stopgain").unwrap();
        assert!(v.annotation_weights.len() == 11);
        assert!(v.annotation_weights[0] > 0.99); // cadd weight for phred=23.7
        assert!(v.annotation_weights[1] > 0.0);  // linsight
    }

    // -- mask predicates on extracted variants --

    #[test]
    fn coding_masks_on_extracted() {
        let (_dir, _db, variants) = setup_rare_table();
        let coding = masks::build_coding_masks(&variants, 1);

        let plof_genes: Vec<&str> = coding.iter()
            .filter(|(mt, _)| *mt == MaskType::PLof)
            .flat_map(|(_, groups)| groups.iter().map(|g| g.name.as_str()))
            .collect();
        assert!(plof_genes.contains(&"OR11H1"));

        let ptv = coding.iter().find(|(mt, _)| *mt == MaskType::Ptv).unwrap();
        for group in &ptv.1 {
            for &idx in &group.variant_indices {
                let csq = &variants[idx].consequence;
                assert!(csq == "stopgain" || csq.starts_with("frameshift"),
                    "ptv should not include: {csq}");
            }
        }
    }

    #[test]
    fn noncoding_masks_on_extracted() {
        let (_dir, _db, variants) = setup_rare_table();
        let noncoding = masks::build_noncoding_masks(&variants, 1);

        let has_mask = |mt: MaskType| -> bool {
            noncoding.iter().any(|(m, groups)| *m == mt && !groups.is_empty())
        };

        assert!(has_mask(MaskType::Upstream));
        assert!(has_mask(MaskType::Downstream));
        assert!(has_mask(MaskType::Utr));
        assert!(has_mask(MaskType::PromoterCage));
        assert!(has_mask(MaskType::PromoterDhs));
        assert!(has_mask(MaskType::EnhancerCage));
        assert!(has_mask(MaskType::EnhancerDhs));
        assert!(has_mask(MaskType::Ncrna));
    }

    #[test]
    fn intronic_in_no_mask() {
        let (_dir, _db, variants) = setup_rare_table();
        let intronic_idx = variants.iter().position(|v| v.region_type == "intronic").unwrap();

        let coding = masks::build_coding_masks(&variants, 1);
        let noncoding = masks::build_noncoding_masks(&variants, 1);

        for (_, groups) in coding.iter().chain(noncoding.iter()) {
            for group in groups {
                assert!(!group.variant_indices.contains(&intronic_idx),
                    "intronic variant should not appear in any mask");
            }
        }
    }

    #[test]
    fn sliding_windows_on_extracted() {
        let (_dir, _db, variants) = setup_rare_table();
        let indices: Vec<usize> = (0..variants.len()).collect();
        let windows = masks::build_sliding_windows(&variants, &indices, "22", 2000, 2000);
        assert!(!windows.is_empty());
        for w in &windows {
            assert!(w.variant_indices.len() >= 2);
        }
    }
}

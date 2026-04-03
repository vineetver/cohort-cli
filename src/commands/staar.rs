use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Builder, Int32Builder, StringBuilder, UInt32Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use faer::Mat;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use serde_json::json;

use crate::commands::dry_run;
use crate::config::Config;
use crate::db::engine::DuckEngine;
use crate::db::query::{query_scalar, query_strings};
use crate::variant_set::VariantSet;
use crate::error::FavorError;
use crate::output::Output;
use crate::resource::Resources;
use crate::staar::masks::{AnnotatedVariant, MaskGroup};
use crate::staar::{self, GeneResult, MaskCategory, MaskType};

const GB: u64 = 1024 * 1024 * 1024;
const MB: u64 = 1024 * 1024;

const WEIGHT_COL_START: usize = 14;
const WEIGHT_COL_COUNT: usize = 11;

struct StaarConfig {
    genotypes: PathBuf,
    phenotype: PathBuf,
    annotations: PathBuf,
    trait_names: Vec<String>,
    covariates: Vec<String>,
    mask_categories: Vec<MaskCategory>,
    maf_cutoff: f64,
    window_size: u32,
    individual: bool,
    spa: bool,
    ancestry_col: Option<String>,
    scang_params: staar::scang::ScangParams,
    known_loci: Option<PathBuf>,
    emit_sumstats: bool,
    output_dir: PathBuf,
}

/// Resource budgets computed from system detection.
struct StaarResources {
    resources: Resources,
    geno_budget: u64,
}

struct ScoringContext<'a> {
    null_model: &'a staar::null_model::NullModel,
    extra_nulls: &'a [staar::null_model::NullModel],
    use_spa: bool,
    ancestry: Option<&'a staar::ancestry::AncestryInfo>,
    n_samples: usize,
}

pub fn run(
    genotypes: PathBuf,
    phenotype: PathBuf,
    trait_names: Vec<String>,
    covariates: Vec<String>,
    annotations: Option<PathBuf>,
    masks: Vec<String>,
    maf_cutoff: f64,
    window_size: u32,
    individual: bool,
    spa: bool,
    ancestry_col: Option<String>,
    scang_lmin: usize,
    scang_lmax: usize,
    scang_step: usize,
    known_loci: Option<PathBuf>,
    emit_sumstats: bool,
    output_path: Option<PathBuf>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), FavorError> {
    let config = validate_and_parse(genotypes, phenotype, trait_names, covariates, annotations, masks, maf_cutoff, window_size, individual, spa, ancestry_col, scang_lmin, scang_lmax, scang_step, known_loci, emit_sumstats, output_path)?;

    if dry_run {
        return emit_dry_run(&config, out);
    }

    std::fs::create_dir_all(&config.output_dir)?;

    let res = setup_resources(out)?;
    let geno = extract_geno(&config, &res, out)?;
    let engine = join_with_annotations(&config, &geno, &res, out)?;

    let is_multi = config.trait_names.len() > 1;
    let primary_trait = &config.trait_names[0];

    let (y, mut x, trait_type, n) = load_phenotype(&engine, &config, &geno, primary_trait, out)?;

    if let Some(ref loci_path) = config.known_loci {
        let x_cond = load_known_loci(&engine, &geno, loci_path, n, out)?;
        x = augment_covariates(&x, &x_cond);
        out.status(&format!("  Conditional: {} known loci added as covariates", x_cond.ncols()));
    }

    let use_spa = config.spa && trait_type == staar::TraitType::Binary;
    if config.spa && trait_type == staar::TraitType::Continuous {
        out.warn("--spa ignored: saddlepoint approximation only applies to binary traits");
    }
    if use_spa {
        out.status("  SPA enabled: saddlepoint approximation for Burden and ACAT-V");
    }

    // AI-STAAR: load ancestry groups if specified
    let ancestry_info: Option<staar::ancestry::AncestryInfo> = if let Some(ref col) = config.ancestry_col {
        let groups = load_ancestry_groups(&engine, col, &geno, out)?;
        let n_pops = *groups.iter().max().unwrap_or(&0) + 1;
        out.status(&format!("  AI-STAAR: {n_pops} populations from column '{col}'"));
        Some(staar::ancestry::AncestryInfo { group: groups, n_pops })
    } else {
        None
    };

    let null_model = fit_null_model(&y, &x, trait_type, out)?;

    if config.emit_sumstats {
        let variants = read_variant_metadata(&engine)?;
        let meta = staar::sumstats::StudyMeta {
            favor_meta_version: 1,
            trait_type: format!("{trait_type:?}"),
            trait_name: config.trait_names[0].clone(),
            n_samples: n,
            sigma2: null_model.sigma2,
            maf_cutoff: config.maf_cutoff,
            covariates: config.covariates.clone(),
            segment_size: 500_000,
        };
        return staar::sumstats::emit(
            &engine, &geno, &variants, &null_model, n, &config.output_dir, &meta, out,
        );
    }

    // MultiSTAAR: fit additional null models for secondary traits
    let mut extra_nulls: Vec<staar::null_model::NullModel> = Vec::new();
    if is_multi {
        out.status(&format!("  MultiSTAAR: {} traits", config.trait_names.len()));
        for trait_name in &config.trait_names[1..] {
            let (y_k, _, tt_k, _) = load_phenotype(&engine, &config, &geno, trait_name, out)?;
            extra_nulls.push(fit_null_model(&y_k, &x, tt_k, out)?);
        }
    }

    let ctx = ScoringContext {
        null_model: &null_model,
        extra_nulls: &extra_nulls,
        use_spa,
        ancestry: ancestry_info.as_ref(),
        n_samples: n,
    };
    let (results, individual_pvals, variants) = run_chromosome_tests(
        &engine, &geno, &config, res.geno_budget, &ctx, out,
    )?;

    if config.individual && !individual_pvals.is_empty() {
        write_individual_results(&engine, &individual_pvals, &variants, &config, out)?;
    }

    write_results(&engine, &results, &config, &null_model, trait_type, n, out)?;
    Ok(())
}

fn validate_and_parse(
    genotypes: PathBuf,
    phenotype: PathBuf,
    trait_names: Vec<String>,
    covariates: Vec<String>,
    annotations: Option<PathBuf>,
    masks: Vec<String>,
    maf_cutoff: f64,
    window_size: u32,
    individual: bool,
    spa: bool,
    ancestry_col: Option<String>,
    scang_lmin: usize,
    scang_lmax: usize,
    scang_step: usize,
    known_loci: Option<PathBuf>,
    emit_sumstats: bool,
    output_path: Option<PathBuf>,
) -> Result<StaarConfig, FavorError> {
    if !genotypes.exists() {
        return Err(FavorError::Input(format!("File not found: {}", genotypes.display())));
    }
    if !phenotype.exists() {
        return Err(FavorError::Input(format!("File not found: {}", phenotype.display())));
    }
    if maf_cutoff <= 0.0 || maf_cutoff >= 0.5 {
        return Err(FavorError::Input(format!("MAF cutoff must be in (0, 0.5), got {maf_cutoff}")));
    }
    let annotations_path = annotations.ok_or_else(|| FavorError::Input(
        "STAAR requires --annotations <path> from `favor annotate`.".into(),
    ))?;
    if !annotations_path.exists() {
        return Err(FavorError::Input(format!(
            "Annotations not found: {}. Run `favor annotate` first.", annotations_path.display()
        )));
    }
    let ann_vs = VariantSet::open(&annotations_path).map_err(|_| FavorError::Input(format!(
        "{} is not a valid variant set. Run `favor annotate` to produce one.",
        annotations_path.display()
    )))?;
    ann_vs.require_annotated()?;
    let mask_categories: Vec<MaskCategory> = masks.iter().map(|s| {
        s.parse::<MaskCategory>().map_err(|_| FavorError::Input(format!(
            "Unknown mask '{s}'. Available: coding, noncoding, sliding-window, scang, custom"
        )))
    }).collect::<Result<_, _>>()?;

    let output_dir = output_path.unwrap_or_else(|| {
        let stem = genotypes.file_stem().unwrap_or_default().to_string_lossy();
        genotypes.with_file_name(format!("{stem}.staar"))
    });

    if let Some(ref loci) = known_loci {
        if !loci.exists() {
            return Err(FavorError::Input(format!("Known loci file not found: {}", loci.display())));
        }
    }

    if trait_names.is_empty() {
        return Err(FavorError::Input("At least one --trait-name is required.".into()));
    }

    Ok(StaarConfig {
        genotypes,
        phenotype,
        annotations: annotations_path,
        trait_names,
        covariates,
        mask_categories,
        maf_cutoff,
        window_size,
        individual,
        spa,
        ancestry_col,
        scang_params: staar::scang::ScangParams {
            lmin: scang_lmin,
            lmax: scang_lmax,
            step: scang_step,
        },
        known_loci,
        emit_sumstats,
        output_dir,
    })
}

fn emit_dry_run(config: &StaarConfig, out: &dyn Output) -> Result<(), FavorError> {
    let is_bgzf = config.genotypes.extension().map(|e| e == "gz" || e == "bgz").unwrap_or(false);
    let sample_names = staar::genotype::read_sample_names(&config.genotypes, is_bgzf)?;
    let n_samples = sample_names.len();

    let ann_rows = parquet_row_count(&config.annotations);
    let est_rare = (ann_rows as f64 * 0.02) as u64;

    let bytes_per_variant = (n_samples as u64) * 8;
    let chr1_geno = (est_rare as f64 * 0.10) as u64 * bytes_per_variant;
    let overhead = 4 * GB;
    let recommended = chr1_geno + overhead;

    let plan = dry_run::DryRunPlan {
        command: "staar".into(),
        inputs: json!({
            "genotypes": config.genotypes.to_string_lossy(),
            "genotype_size": dry_run::file_size(&config.genotypes),
            "annotations": config.annotations.to_string_lossy(),
            "annotation_rows": ann_rows,
            "n_samples": n_samples,
            "estimated_rare_variants": est_rare,
            "trait": config.trait_names[0],
            "maf_cutoff": config.maf_cutoff,
        }),
        memory: dry_run::MemoryEstimate {
            minimum: "4G".into(),
            recommended: dry_run::human_bytes(recommended.max(8 * GB)),
            minimum_bytes: 4 * GB,
            recommended_bytes: recommended.max(8 * GB),
        },
        output_path: config.output_dir.to_string_lossy().into(),
    };
    dry_run::emit(&plan, out);
    Ok(())
}

fn setup_resources(out: &dyn Output) -> Result<StaarResources, FavorError> {
    let config_resources = Config::load()
        .map(|c| c.resources)
        .unwrap_or_default();
    let resources = Resources::detect_with_config(&config_resources);

    rayon::ThreadPoolBuilder::new()
        .num_threads(resources.threads)
        .build_global()
        .ok();

    // Memory budget partitioning:
    //   - DuckDB 20% (capped 8G): sufficient for joining genotypes x annotations
    //     and temp tables. 20% is generous; DuckDB streams most joins. Cap at 8G
    //     because DuckDB spills to disk beyond its buffer anyway.
    //   - Rayon 128M/thread: covers per-thread stack + intermediate faer matrices
    //     during score tests. Measured empirically on genes with <=5K variants.
    //   - Remainder -> genotype buffer: the dominant consumer, holding f64 dosage
    //     matrices (n_samples x n_variants per chromosome/batch).
    // Adjust if: genotype matrices are sparse (could compress), or DuckDB joins
    // grow (e.g., noncoding masks with whole-genome scans).
    let total = resources.memory_bytes;
    let duckdb_budget = (total / 5).min(8 * GB);
    let rayon_overhead = 128 * MB * resources.threads as u64;
    let geno_budget = total.saturating_sub(duckdb_budget + rayon_overhead);

    out.status(&format!("STAAR: {} memory, {} threads ({})",
        resources.memory_human(), resources.threads, resources.environment()));
    out.status(&format!("  Budget: {:.1}G genotypes, {:.1}G DuckDB, {:.1}G rayon ({} threads)",
        geno_budget as f64 / GB as f64, duckdb_budget as f64 / GB as f64,
        rayon_overhead as f64 / GB as f64, resources.threads));

    Ok(StaarResources { resources, geno_budget })
}

fn extract_geno(
    config: &StaarConfig,
    res: &StaarResources,
    out: &dyn Output,
) -> Result<staar::genotype::GenotypeResult, FavorError> {
    out.status("Step 1/6: Extracting genotypes from VCF...");
    staar::genotype::extract_genotypes(
        &config.genotypes, &config.output_dir, res.geno_budget, out,
    )
}

fn join_with_annotations(
    config: &StaarConfig,
    geno: &staar::genotype::GenotypeResult,
    res: &StaarResources,
    out: &dyn Output,
) -> Result<DuckEngine, FavorError> {
    out.status("Step 2/6: Joining genotypes with annotations...");
    let engine = DuckEngine::new(&res.resources)?;
    engine.execute("SET preserve_insertion_order = false;")?;

    let ann_vs = VariantSet::open(&config.annotations)?;

    // Validate annotation columns exist before running a multi-hour job
    let ann_cols = query_strings(&engine, &format!(
        "SELECT column_name FROM (DESCRIBE SELECT * FROM {} LIMIT 0)",
        ann_vs.read_all()
    ))?;
    let required = ["gencode", "main", "cage", "apc"];
    for r in required {
        if !ann_cols.iter().any(|c| c == r) {
            let tier_hint = match ann_vs.kind() {
                Some(crate::variant_set::VariantSetKind::Annotated { tier }) => format!(" (tier: {tier})"),
                _ => String::new(),
            };
            return Err(FavorError::DataMissing(format!(
                "Annotation column '{r}' missing in {}{tier_hint}. \
                 STAAR requires favor-full annotations. Re-run: `favor annotate --full`.",
                config.annotations.display()
            )));
        }
    }

    let ann_read = ann_vs.read_all();
    let geno_glob = format!("{}/chromosome=*/data.parquet", geno.output_dir.display());
    let maf_cutoff = config.maf_cutoff;

    engine.execute(&format!(
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
         FROM read_parquet('{geno_glob}', hive_partitioning=true) g
         INNER JOIN ({ann_read}) a
             ON g.chromosome::VARCHAR = a.chromosome::VARCHAR
             AND g.position = a.position AND g.ref = a.ref_vcf AND g.alt = a.alt_vcf
         WHERE g.maf > 0 AND g.maf < {maf_cutoff}",
        geno_glob = geno_glob, ann_read = ann_read, maf_cutoff = maf_cutoff,
    ))?;

    let n_rare = query_scalar(&engine, "SELECT COUNT(*) FROM _rare")?;
    let n_genes = query_scalar(&engine, "SELECT COUNT(DISTINCT gene_name) FROM _rare WHERE gene_name != ''")?;
    out.status(&format!("  {} rare variants, {} genes", n_rare, n_genes));

    if n_rare == 0 {
        return Err(FavorError::Analysis(
            "No rare variants found after joining genotypes with annotations.".into(),
        ));
    }

    Ok(engine)
}

fn load_phenotype(
    engine: &DuckEngine,
    config: &StaarConfig,
    geno: &staar::genotype::GenotypeResult,
    trait_name: &str,
    out: &dyn Output,
) -> Result<(Mat<f64>, Mat<f64>, staar::TraitType, usize), FavorError> {
    out.status("Step 3/6: Loading phenotype...");

    // Create _pheno table if it doesn't exist yet (first trait)
    engine.execute(&format!(
        "CREATE TEMP TABLE IF NOT EXISTS _pheno AS SELECT * FROM read_csv_auto('{}')",
        config.phenotype.display(),
    ))?;

    let pheno_cols = query_strings(engine, "SELECT column_name FROM (DESCRIBE _pheno)")?;
    if !pheno_cols.contains(&trait_name.to_string()) {
        return Err(FavorError::Input(format!(
            "Trait '{}' not in phenotype. Available: {}", trait_name, pheno_cols.join(", ")
        )));
    }
    for cov in &config.covariates {
        if !pheno_cols.contains(cov) {
            return Err(FavorError::Input(format!("Covariate '{cov}' not in phenotype")));
        }
    }

    let trait_type = if query_scalar(engine, &format!(
        "SELECT COUNT(DISTINCT \"{trait_name}\") FROM _pheno"))? <= 2 {
        staar::TraitType::Binary
    } else {
        staar::TraitType::Continuous
    };
    out.status(&format!("  Trait '{trait_name}' -> {:?}", trait_type));

    let id_col = &pheno_cols[0];
    let cov_select = if config.covariates.is_empty() { String::new() }
        else { format!(", {}", config.covariates.iter().map(|c| format!("p.\"{c}\"")).collect::<Vec<_>>().join(", ")) };

    let sample_list = geno.sample_names.iter().map(|s| format!("'{s}'")).collect::<Vec<_>>().join(",");
    let pheno_sql = format!(
        "SELECT p.\"{trait_name}\" {cov_select} FROM _pheno p \
         WHERE p.\"{id_col}\" IN ({sample_list}) AND p.\"{trait_name}\" IS NOT NULL \
         ORDER BY p.\"{id_col}\""
    );

    let (y, x) = read_phenotype_matrix(engine, &pheno_sql, config.covariates.len())?;
    let n = y.nrows();
    out.status(&format!("  {} samples with phenotype + genotype", n));
    if n < 10 {
        return Err(FavorError::Analysis(format!("Only {n} samples. Need >= 10.")));
    }

    Ok((y, x, trait_type, n))
}

fn load_ancestry_groups(
    engine: &DuckEngine,
    col: &str,
    geno: &staar::genotype::GenotypeResult,
    out: &dyn Output,
) -> Result<Vec<usize>, FavorError> {
    let pheno_cols = query_strings(engine, "SELECT column_name FROM (DESCRIBE _pheno)")?;
    if !pheno_cols.contains(&col.to_string()) {
        return Err(FavorError::Input(format!(
            "Ancestry column '{col}' not in phenotype. Available: {}", pheno_cols.join(", ")
        )));
    }

    let id_col = &pheno_cols[0];
    let sample_list = geno.sample_names.iter().map(|s| format!("'{s}'")).collect::<Vec<_>>().join(",");

    // Get distinct population labels and map to 0-indexed integers
    let labels = query_strings(engine, &format!(
        "SELECT DISTINCT \"{col}\" FROM _pheno WHERE \"{id_col}\" IN ({sample_list}) ORDER BY \"{col}\""
    ))?;
    out.status(&format!("  Populations: {}", labels.join(", ")));

    let sql = format!(
        "SELECT \"{col}\" FROM _pheno WHERE \"{id_col}\" IN ({sample_list}) ORDER BY \"{id_col}\""
    );
    let conn = engine.connection();
    let mut stmt = conn.prepare(&sql).map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut rows = stmt.query([]).map_err(|e| FavorError::Analysis(format!("{e}")))?;

    let mut groups = Vec::new();
    while let Ok(Some(row)) = rows.next() {
        let label: String = row.get(0).unwrap_or_default();
        let idx = labels.iter().position(|l| l == &label).unwrap_or(0);
        groups.push(idx);
    }
    Ok(groups)
}

fn load_known_loci(
    engine: &DuckEngine,
    geno: &staar::genotype::GenotypeResult,
    loci_path: &Path,
    n_samples: usize,
    out: &dyn Output,
) -> Result<Mat<f64>, FavorError> {
    let content = std::fs::read_to_string(loci_path)?;
    let loci: Vec<(&str, i32)> = content.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let parts: Vec<&str> = l.split(':').collect();
            if parts.len() >= 2 {
                let chrom = parts[0].strip_prefix("chr").unwrap_or(parts[0]);
                parts[1].parse::<i32>().ok().map(|pos| (chrom, pos))
            } else { None }
        })
        .collect();

    if loci.is_empty() {
        return Err(FavorError::Input("Known loci file is empty or unparseable.".into()));
    }

    out.status(&format!("  Loading {} known loci for conditional analysis...", loci.len()));

    let extract_cols = staar::geno_load::dosage_columns(n_samples);

    let mut by_chrom: std::collections::HashMap<&str, Vec<(usize, i32)>> = std::collections::HashMap::new();
    for (col, &(chrom, pos)) in loci.iter().enumerate() {
        by_chrom.entry(chrom).or_default().push((col, pos));
    }

    let mut x_cond = Mat::zeros(n_samples, loci.len());

    for (chrom, chrom_loci) in &by_chrom {
        let geno_path = format!("{}/chromosome={chrom}/data.parquet", geno.output_dir.display());
        let positions = chrom_loci.iter().map(|(_, p)| p.to_string()).collect::<Vec<_>>().join(",");
        let sql = format!(
            "SELECT position, {extract_cols} FROM read_parquet('{geno_path}') WHERE position IN ({positions}) ORDER BY position"
        );
        let conn = engine.connection();
        let mut stmt = conn.prepare(&sql)
            .map_err(|e| FavorError::Analysis(format!("Known loci {chrom}: {e}")))?;
        let mut rows = stmt.query([])
            .map_err(|e| FavorError::Analysis(format!("Known loci {chrom}: {e}")))?;

        let pos_to_cols: std::collections::HashMap<i32, Vec<usize>> = {
            let mut m: std::collections::HashMap<i32, Vec<usize>> = std::collections::HashMap::new();
            for &(col, pos) in chrom_loci {
                m.entry(pos).or_default().push(col);
            }
            m
        };

        while let Ok(Some(row)) = rows.next() {
            let pos: i32 = row.get(0).unwrap_or(0);
            if let Some(cols) = pos_to_cols.get(&pos) {
                for &col in cols {
                    for i in 0..n_samples {
                        x_cond[(i, col)] = row.get::<_, f64>(i + 1).unwrap_or(0.0);
                    }
                }
            }
        }
    }
    Ok(x_cond)
}

fn augment_covariates(x: &Mat<f64>, x_cond: &Mat<f64>) -> Mat<f64> {
    let n = x.nrows();
    let k_orig = x.ncols();
    let k_cond = x_cond.ncols();
    let mut x_new = Mat::zeros(n, k_orig + k_cond);
    for i in 0..n {
        for j in 0..k_orig { x_new[(i, j)] = x[(i, j)]; }
        for j in 0..k_cond { x_new[(i, k_orig + j)] = x_cond[(i, j)]; }
    }
    x_new
}

fn fit_null_model(
    y: &Mat<f64>,
    x: &Mat<f64>,
    trait_type: staar::TraitType,
    out: &dyn Output,
) -> Result<staar::null_model::NullModel, FavorError> {
    out.status("Step 4/6: Fitting null model...");
    let null_model = match trait_type {
        staar::TraitType::Continuous => staar::null_model::fit_glm(y, x),
        staar::TraitType::Binary => staar::null_model::fit_logistic(y, x, 25),
    };
    out.status(&format!("  sigma2 = {:.4}", null_model.sigma2));
    Ok(null_model)
}

fn run_chromosome_tests(
    engine: &DuckEngine,
    geno: &staar::genotype::GenotypeResult,
    config: &StaarConfig,
    geno_budget: u64,
    ctx: &ScoringContext,
    out: &dyn Output,
) -> Result<(Vec<(MaskType, Vec<GeneResult>)>, Vec<(usize, f64)>, Vec<AnnotatedVariant>), FavorError> {
    out.status("Step 5/6: Running score tests (adaptive memory)...");

    let n = ctx.n_samples;
    let variants = read_variant_metadata(engine)?;
    out.status(&format!("  {} rare variants across chromosomes", variants.len()));

    let bytes_per_variant = (n as u64) * 8;
    let max_variants_in_ram = (geno_budget / bytes_per_variant).max(1000) as usize;
    out.status(&format!("  Genotype budget: {} variants in RAM ({} per variant x {} samples)",
        max_variants_in_ram,
        crate::commands::data::human_size(bytes_per_variant),
        n,
    ));

    let extract_cols = staar::geno_load::dosage_columns(n);

    let chromosomes = query_strings(engine,
        "SELECT DISTINCT chrom FROM _rare ORDER BY chrom")?;

    let mut gene_masks: Vec<(MaskType, Vec<MaskGroup>)> = Vec::new();
    let mut has_windows = false;
    let mut has_scang = false;
    for cat in &config.mask_categories {
        match cat {
            MaskCategory::Coding => gene_masks.extend(staar::masks::build_coding_masks(&variants, 2)),
            MaskCategory::Noncoding => gene_masks.extend(staar::masks::build_noncoding_masks(&variants, 2)),
            MaskCategory::SlidingWindow => has_windows = true,
            MaskCategory::Scang => has_scang = true,
            MaskCategory::Custom => out.warn("Custom BED: not yet implemented"),
        }
    }

    let mut all_results: Vec<(MaskType, Vec<GeneResult>)> = gene_masks.iter()
        .map(|(mt, _)| (mt.clone(), Vec::new()))
        .collect();

    let window_idx = if has_windows {
        all_results.push((MaskType::SlidingWindow, Vec::new()));
        Some(all_results.len() - 1)
    } else {
        None
    };

    let scang_idx = if has_scang {
        all_results.push((MaskType::Scang, Vec::new()));
        Some(all_results.len() - 1)
    } else {
        None
    };

    let mut individual_pvals: Vec<(usize, f64)> = Vec::new();

    for chrom in &chromosomes {
        let chrom_variant_indices: Vec<usize> = variants.iter().enumerate()
            .filter(|(_, v)| v.chromosome == *chrom)
            .map(|(i, _)| i)
            .collect();

        let n_chrom = chrom_variant_indices.len();
        if n_chrom == 0 { continue; }

        let geno_path = format!(
            "{}/chromosome={chrom}/data.parquet",
            geno.output_dir.display()
        );

        let window_groups: Vec<MaskGroup> = if has_windows {
            staar::masks::build_sliding_windows(
                &variants, &chrom_variant_indices, chrom,
                config.window_size, config.window_size / 2,
            )
        } else {
            Vec::new()
        };

        let scang_all: Vec<(u32, Vec<MaskGroup>)> = if has_scang {
            staar::scang::build_scang_windows(
                &variants, &chrom_variant_indices, chrom, &config.scang_params,
            )
        } else {
            Vec::new()
        };

        if n_chrom <= max_variants_in_ram {
            out.status(&format!("  chr{chrom}: {} variants (fits in RAM), loading...", n_chrom));

            let mut global_to_local: std::collections::HashMap<usize, usize> =
                std::collections::HashMap::with_capacity(n_chrom);
            for (local, &global) in chrom_variant_indices.iter().enumerate() {
                global_to_local.insert(global, local);
            }

            // Load entire chromosome without per-position IN-list — avoids
            // serializing thousands of positions into the SQL string.
            let (geno_flat, loaded_positions) = staar::geno_load::load_all(
                engine, &geno_path, n, &extract_cols,
            )?;

            // Re-map global indices: loaded parquet positions may be a superset
            // of chrom_variant_indices (genotype parquet has all extracted variants,
            // _rare has only those passing MAF + annotation join). Build position
            // lookup from loaded data.
            let pos_to_loaded: std::collections::HashMap<u32, usize> = loaded_positions
                .iter().enumerate().map(|(i, &p)| (p, i)).collect();
            // Remap global_to_local to point into loaded_positions order
            let mut global_to_local_loaded: std::collections::HashMap<usize, usize> =
                std::collections::HashMap::with_capacity(n_chrom);
            for &gi in &chrom_variant_indices {
                if let Some(&loaded_idx) = pos_to_loaded.get(&variants[gi].position) {
                    global_to_local_loaded.insert(gi, loaded_idx);
                }
            }
            let global_to_local = global_to_local_loaded;

            for (idx, (mask_type, groups)) in gene_masks.iter().enumerate() {
                let r = score_chrom_genes(
                    groups, chrom, &variants, &geno_flat, &global_to_local, ctx,
                );
                if !r.is_empty() {
                    out.status(&format!("    {}: {} groups on chr{chrom}",
                        mask_type.file_stem(), r.len()));
                    all_results[idx].1.extend(r);
                }
            }

            if let Some(wi) = window_idx {
                if !window_groups.is_empty() {
                    let r: Vec<GeneResult> = window_groups.par_iter()
                        .filter_map(|g| score_gene(g, &variants, &geno_flat, &global_to_local, ctx))
                        .collect();
                    if !r.is_empty() {
                        out.status(&format!("    sliding_window: {} windows on chr{chrom}", r.len()));
                        all_results[wi].1.extend(r);
                    }
                }
            }

            if let Some(si) = scang_idx {
                for (wsize, groups) in &scang_all {
                    let r: Vec<GeneResult> = groups.par_iter()
                        .filter_map(|g| score_gene(g, &variants, &geno_flat, &global_to_local, ctx))
                        .collect();
                    if !r.is_empty() {
                        out.status(&format!("    scang L={wsize}: {} windows on chr{chrom}", r.len()));
                        all_results[si].1.extend(r);
                    }
                }
            }

            if config.individual {
                let g = staar::geno_load::to_mat(&geno_flat, n, n_chrom);
                let pvals = staar::score_test::individual_tests(&g, ctx.null_model, ctx.use_spa);
                for (local, pval) in pvals.into_iter().enumerate() {
                    individual_pvals.push((chrom_variant_indices[local], pval));
                }
            }
        } else {
            out.status(&format!("  chr{chrom}: {} variants (exceeds budget of {}), batching...",
                n_chrom, max_variants_in_ram));

            for (idx, (mask_type, groups)) in gene_masks.iter().enumerate() {
                score_batched(
                    groups, chrom, &variants, engine, &geno_path, &extract_cols,
                    max_variants_in_ram, mask_type, ctx, out,
                    &mut all_results[idx].1,
                )?;
            }

            if let Some(wi) = window_idx {
                if !window_groups.is_empty() {
                    score_batched(
                        &window_groups, chrom, &variants, engine, &geno_path, &extract_cols,
                        max_variants_in_ram, &MaskType::SlidingWindow, ctx, out,
                        &mut all_results[wi].1,
                    )?;
                }
            }

            if let Some(si) = scang_idx {
                for (wsize, groups) in &scang_all {
                    score_batched(
                        groups, chrom, &variants, engine, &geno_path, &extract_cols,
                        max_variants_in_ram, &MaskType::Scang, ctx, out,
                        &mut all_results[si].1,
                    )?;
                    out.status(&format!("    scang L={wsize}: {} windows batched on chr{chrom}", groups.len()));
                }
            }

            if config.individual {
                for chunk in chrom_variant_indices.chunks(max_variants_in_ram) {
                    let positions: Vec<u32> = chunk.iter().map(|&i| variants[i].position).collect();
                    let geno_flat = staar::geno_load::load(
                        engine, &geno_path, &positions, n, &extract_cols,
                    )?;
                    let g = staar::geno_load::to_mat(&geno_flat, n, positions.len());
                    let pvals = staar::score_test::individual_tests(&g, ctx.null_model, ctx.use_spa);
                    for (local, pval) in pvals.into_iter().enumerate() {
                        individual_pvals.push((chunk[local], pval));
                    }
                }
            }
        }
    }

    Ok((all_results, individual_pvals, variants))
}

fn score_batched(
    groups: &[MaskGroup], chrom: &str, variants: &[AnnotatedVariant],
    engine: &DuckEngine, geno_path: &str, extract_cols: &str,
    max_variants: usize, mask_type: &MaskType,
    ctx: &ScoringContext, out: &dyn Output,
    results: &mut Vec<GeneResult>,
) -> Result<(), FavorError> {
    let mut chrom_groups: Vec<&MaskGroup> = groups.iter()
        .filter(|g| g.chromosome == chrom).collect();
    if chrom_groups.is_empty() { return Ok(()); }

    chrom_groups.sort_by_key(|g|
        g.variant_indices.first().map(|&i| variants[i].position).unwrap_or(0)
    );

    let batches = pack_gene_batches(&chrom_groups, variants, max_variants);
    out.status(&format!("    {}: {} groups in {} batches on chr{chrom}",
        mask_type.file_stem(), chrom_groups.len(), batches.len()));

    run_batches(&batches, variants, engine, geno_path, extract_cols, chrom, ctx, out, results)
}

fn run_batches(
    batches: &[Vec<&MaskGroup>],
    variants: &[AnnotatedVariant],
    engine: &DuckEngine,
    geno_path: &str,
    extract_cols: &str,
    _chrom: &str,
    ctx: &ScoringContext,
    out: &dyn Output,
    results: &mut Vec<GeneResult>,
) -> Result<(), FavorError> {
    for (batch_idx, batch) in batches.iter().enumerate() {
        let mut positions: Vec<u32> = batch.iter()
            .flat_map(|g| g.variant_indices.iter().map(|&i| variants[i].position))
            .collect();
        positions.sort_unstable();
        positions.dedup();

        let geno_flat = staar::geno_load::load(
            engine, geno_path, &positions, ctx.n_samples, extract_cols,
        )?;

        let pos_to_local: std::collections::HashMap<u32, usize> = positions.iter()
            .enumerate().map(|(l, &p)| (p, l)).collect();

        let mut g2l: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for group in batch.iter() {
            for &gi in &group.variant_indices {
                if let std::collections::hash_map::Entry::Vacant(e) = g2l.entry(gi) {
                    if let Some(&local) = pos_to_local.get(&variants[gi].position) {
                        e.insert(local);
                    }
                }
            }
        }

        let batch_results: Vec<GeneResult> = batch.par_iter()
            .filter_map(|g| score_gene(g, variants, &geno_flat, &g2l, ctx))
            .collect();

        if !batch_results.is_empty() {
            out.status(&format!("      batch {}/{}: {} groups, {} variants",
                batch_idx + 1, batches.len(), batch_results.len(), positions.len()));
            results.extend(batch_results);
        }
    }
    Ok(())
}

fn write_individual_results(
    _engine: &DuckEngine,
    pvals: &[(usize, f64)],
    variants: &[AnnotatedVariant],
    config: &StaarConfig,
    out: &dyn Output,
) -> Result<(), FavorError> {
    let out_path = config.output_dir.join("individual.parquet");
    let n = pvals.len();

    // Sort by p-value for output ordering
    let mut sorted: Vec<(usize, f64)> = pvals.to_vec();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut b_chrom = StringBuilder::with_capacity(n, n * 2);
    let mut b_pos = Int32Builder::with_capacity(n);
    let mut b_ref = StringBuilder::with_capacity(n, n * 4);
    let mut b_alt = StringBuilder::with_capacity(n, n * 4);
    let mut b_maf = Float64Builder::with_capacity(n);
    let mut b_gene = StringBuilder::with_capacity(n, n * 8);
    let mut b_region = StringBuilder::with_capacity(n, n * 8);
    let mut b_consequence = StringBuilder::with_capacity(n, n * 8);
    let mut b_cadd = Float64Builder::with_capacity(n);
    let mut b_pvalue = Float64Builder::with_capacity(n);

    for &(idx, pval) in &sorted {
        let v = &variants[idx];
        b_chrom.append_value(&v.chromosome);
        b_pos.append_value(v.position as i32);
        b_ref.append_value(&v.ref_allele);
        b_alt.append_value(&v.alt_allele);
        b_maf.append_value(v.maf);
        b_gene.append_value(&v.gene_name);
        b_region.append_value(&v.region_type);
        b_consequence.append_value(&v.consequence);
        b_cadd.append_value(v.cadd_phred);
        b_pvalue.append_value(pval);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("maf", DataType::Float64, false),
        Field::new("gene_name", DataType::Utf8, false),
        Field::new("region_type", DataType::Utf8, false),
        Field::new("consequence", DataType::Utf8, false),
        Field::new("cadd_phred", DataType::Float64, false),
        Field::new("pvalue", DataType::Float64, false),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(b_chrom.finish()), Arc::new(b_pos.finish()),
        Arc::new(b_ref.finish()), Arc::new(b_alt.finish()),
        Arc::new(b_maf.finish()), Arc::new(b_gene.finish()),
        Arc::new(b_region.finish()), Arc::new(b_consequence.finish()),
        Arc::new(b_cadd.finish()), Arc::new(b_pvalue.finish()),
    ];
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(&out_path)
        .map_err(|e| FavorError::Resource(format!("Create individual.parquet: {e}")))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| FavorError::Resource(format!("Parquet writer: {e}")))?;
    writer.write(&batch).map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    let n_sig = pvals.iter().filter(|(_, p)| *p < 5e-8).count();
    out.success(&format!("  individual -> {} variants, {} genome-wide significant",
        pvals.len(), n_sig));
    Ok(())
}

fn write_results(
    engine: &DuckEngine,
    all_mask_results: &[(MaskType, Vec<GeneResult>)],
    config: &StaarConfig,
    null_model: &staar::null_model::NullModel,
    trait_type: staar::TraitType,
    n: usize,
    out: &dyn Output,
) -> Result<(), FavorError> {
    out.status("Step 6/6: Writing results...");
    let mut significant_genes: Vec<serde_json::Value> = Vec::new();

    let n_rare: i64 = query_scalar(engine, "SELECT COUNT(*) FROM _rare")?;

    let channels = staar::weights::ANNOTATION_CHANNELS;
    let n_channels = channels.len();

    for (mask_type, results) in all_mask_results {
        if results.is_empty() { continue; }
        let out_path = config.output_dir.join(format!("{}.parquet", mask_type.file_stem()));

        // Sort by STAAR-O p-value
        let mut sorted_results: Vec<&GeneResult> = results.iter().collect();
        sorted_results.sort_by(|a, b| a.staar.staar_o.partial_cmp(&b.staar.staar_o).unwrap_or(std::cmp::Ordering::Equal));

        let nr = sorted_results.len();
        let mut b_ensembl = StringBuilder::with_capacity(nr, nr * 16);
        let mut b_symbol = StringBuilder::with_capacity(nr, nr * 12);
        let mut b_chrom = StringBuilder::with_capacity(nr, nr * 2);
        let mut b_start = UInt32Builder::with_capacity(nr);
        let mut b_end = UInt32Builder::with_capacity(nr);
        let mut b_nvariants = UInt32Builder::with_capacity(nr);
        let mut b_cmac = UInt32Builder::with_capacity(nr);

        // 6 base + 6*n_channels per-annotation + 6 per-test omnibus + acat_o + staar_o = 6 + 6*11 + 6 + 2 = 80
        let n_pval_cols = 6 + 6 * n_channels + 6 + 2;
        let mut pval_builders: Vec<Float64Builder> = (0..n_pval_cols)
            .map(|_| Float64Builder::with_capacity(nr))
            .collect();

        for r in &sorted_results {
            let s = &r.staar;
            b_ensembl.append_value(&r.ensembl_id);
            b_symbol.append_value(&r.gene_symbol);
            b_chrom.append_value(&r.chromosome);
            b_start.append_value(r.start);
            b_end.append_value(r.end);
            b_nvariants.append_value(r.n_variants);
            b_cmac.append_value(r.cumulative_mac);

            let mut pi = 0;
            for p in [s.burden_1_25, s.burden_1_1, s.skat_1_25, s.skat_1_1, s.acat_v_1_25, s.acat_v_1_1] {
                pval_builders[pi].append_value(p);
                pi += 1;
            }
            for ann_p in &s.per_annotation {
                for &v in ann_p {
                    pval_builders[pi].append_value(v);
                    pi += 1;
                }
            }
            // Pad if fewer annotation channels than expected
            while pi < 6 + 6 * n_channels {
                pval_builders[pi].append_value(f64::NAN);
                pi += 1;
            }
            for p in [s.staar_b_1_25, s.staar_b_1_1, s.staar_s_1_25, s.staar_s_1_1, s.staar_a_1_25, s.staar_a_1_1, s.acat_o, s.staar_o] {
                pval_builders[pi].append_value(p);
                pi += 1;
            }
        }

        // Build schema
        let test_names = ["Burden(1,25)", "Burden(1,1)", "SKAT(1,25)", "SKAT(1,1)", "ACAT-V(1,25)", "ACAT-V(1,1)"];
        let mut fields = vec![
            Field::new("ensembl_id", DataType::Utf8, false),
            Field::new("gene_symbol", DataType::Utf8, false),
            Field::new("chromosome", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("n_variants", DataType::UInt32, false),
            Field::new("cMAC", DataType::UInt32, false),
        ];
        for test in &test_names {
            fields.push(Field::new(*test, DataType::Float64, true));
        }
        for ch in channels {
            for test in &test_names {
                fields.push(Field::new(format!("{test}-{ch}"), DataType::Float64, true));
            }
        }
        for name in ["STAAR-B(1,25)", "STAAR-B(1,1)", "STAAR-S(1,25)", "STAAR-S(1,1)", "STAAR-A(1,25)", "STAAR-A(1,1)", "ACAT-O", "STAAR-O"] {
            fields.push(Field::new(name, DataType::Float64, true));
        }
        let schema = Arc::new(Schema::new(fields));

        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(b_ensembl.finish()), Arc::new(b_symbol.finish()),
            Arc::new(b_chrom.finish()), Arc::new(b_start.finish()),
            Arc::new(b_end.finish()), Arc::new(b_nvariants.finish()),
            Arc::new(b_cmac.finish()),
        ];
        for b in &mut pval_builders {
            columns.push(Arc::new(b.finish()));
        }

        let batch = RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;
        let file = File::create(&out_path)
            .map_err(|e| FavorError::Resource(format!("Create {}: {e}", out_path.display())))?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| FavorError::Resource(format!("Parquet writer: {e}")))?;
        writer.write(&batch).map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
        writer.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

        let n_sig = results.iter().filter(|r| r.staar.staar_o < 2.5e-6).count();
        for r in results.iter().filter(|r| r.staar.staar_o < 2.5e-6) {
            significant_genes.push(json!({
                "gene": r.ensembl_id, "mask": mask_type.file_stem(),
                "STAAR-O": r.staar.staar_o, "n_variants": r.n_variants,
            }));
        }
        out.success(&format!("  {} -> {} genes, {} significant",
            mask_type.file_stem(), results.len(), n_sig));
    }

    let meta = json!({
        "favor_staar_version": 1,
        "traits": config.trait_names, "trait_type": format!("{:?}", trait_type),
        "n_samples": n, "n_rare_variants": n_rare, "maf_cutoff": config.maf_cutoff,
        "sigma2": null_model.sigma2, "significant_genes": significant_genes,
    });
    let _ = std::fs::write(config.output_dir.join("staar.meta.json"),
        serde_json::to_string_pretty(&meta).unwrap_or_default());

    match staar::summary::generate_report(all_mask_results, &config.trait_names, n, n_rare, &config.output_dir, "STAAR Rare Variant Association") {
        Ok(()) => out.success(&format!("  summary.html -> {}", config.output_dir.join("summary.html").display())),
        Err(e) => out.warn(&format!("  Summary report failed: {e}")),
    }

    out.success(&format!("STAAR complete -> {}", config.output_dir.display()));
    out.result_json(&meta);
    Ok(())
}


fn score_gene(
    group: &MaskGroup, variants: &[AnnotatedVariant],
    geno_flat: &[f64], global_to_local: &std::collections::HashMap<usize, usize>,
    ctx: &ScoringContext,
) -> Option<GeneResult> {
    let n_samples = ctx.n_samples;
    let local_indices: Vec<usize> = group.variant_indices.iter()
        .filter_map(|&gi| global_to_local.get(&gi).copied())
        .collect();

    if local_indices.len() < 2 { return None; }
    let m = local_indices.len();

    let g = Mat::from_fn(n_samples, m, |row, col| {
        geno_flat[local_indices[col] * n_samples + row]
    });

    let mafs: Vec<f64> = group.variant_indices.iter()
        .filter(|gi| global_to_local.contains_key(gi))
        .map(|&i| variants[i].maf)
        .collect();

    let n_channels = staar::weights::ANNOTATION_CHANNELS.len();
    let ann_matrix: Vec<Vec<f64>> = (0..n_channels)
        .map(|ch| group.variant_indices.iter()
            .filter(|gi| global_to_local.contains_key(gi))
            .map(|&i| variants[i].annotation_weights.get(ch).copied().unwrap_or(0.0))
            .collect())
        .collect();

    let sr = if let Some(anc) = ctx.ancestry {
        let pop_mafs: Vec<Vec<f64>> = (0..m).map(|col| {
            let li = local_indices[col];
            let base = li * n_samples;
            let mut pop_sums = vec![0.0; anc.n_pops];
            let mut pop_counts = vec![0usize; anc.n_pops];
            for i in 0..n_samples {
                let p = anc.group[i];
                pop_sums[p] += geno_flat[base + i];
                pop_counts[p] += 1;
            }
            (0..anc.n_pops).map(|p| {
                if pop_counts[p] > 0 { (pop_sums[p] / (2.0 * pop_counts[p] as f64)).clamp(1e-10, 0.499) }
                else { 0.0 }
            }).collect()
        }).collect();
        staar::ancestry::run_ai_staar(&g, &ann_matrix, &pop_mafs, anc, ctx.null_model, ctx.use_spa)
    } else if !ctx.extra_nulls.is_empty() {
        let mut all_nulls: Vec<&staar::null_model::NullModel> = vec![ctx.null_model];
        all_nulls.extend(ctx.extra_nulls.iter());
        let multi = staar::multi::run_multi_staar(&g, &ann_matrix, &mafs, &all_nulls, ctx.use_spa);
        let mut primary = multi.per_trait.into_iter().next()?;
        primary.staar_o = multi.multi_staar_o;
        primary.acat_o = multi.multi_burden;
        primary
    } else {
        staar::score_test::run_staar(&g, &ann_matrix, &mafs, ctx.null_model, ctx.use_spa)
    };

    let cmac: u32 = group.variant_indices.iter()
        .filter(|gi| global_to_local.contains_key(gi))
        .map(|&i| (2.0 * variants[i].maf * n_samples as f64).round() as u32)
        .sum();

    Some(GeneResult {
        ensembl_id: group.name.clone(),
        gene_symbol: group.name.clone(),
        chromosome: group.chromosome.clone(),
        start: group.start, end: group.end,
        n_variants: m as u32, cumulative_mac: cmac,
        staar: sr,
    })
}

fn score_chrom_genes(
    groups: &[MaskGroup], chrom: &str,
    variants: &[AnnotatedVariant],
    geno_flat: &[f64], global_to_local: &std::collections::HashMap<usize, usize>,
    ctx: &ScoringContext,
) -> Vec<GeneResult> {
    let chrom_groups: Vec<&MaskGroup> = groups.iter()
        .filter(|g| g.chromosome == chrom)
        .collect();
    if chrom_groups.is_empty() { return Vec::new(); }

    chrom_groups.par_iter().filter_map(|group| {
        score_gene(group, variants, geno_flat, global_to_local, ctx)
    }).collect()
}

fn pack_gene_batches<'a>(
    genes: &[&'a MaskGroup], variants: &[AnnotatedVariant], max_variants: usize,
) -> Vec<Vec<&'a MaskGroup>> {
    let mut batches: Vec<Vec<&'a MaskGroup>> = Vec::new();
    let mut current_batch: Vec<&'a MaskGroup> = Vec::new();
    let mut current_positions: std::collections::HashSet<u32> = std::collections::HashSet::new();

    for &gene in genes {
        // Compute what the union would be if we added this gene
        let new_positions: Vec<u32> = gene.variant_indices.iter()
            .map(|&i| variants[i].position)
            .filter(|p| !current_positions.contains(p))
            .collect();

        let union_size = current_positions.len() + new_positions.len();

        if !current_batch.is_empty() && union_size > max_variants {
            // This gene would overflow — flush current batch
            batches.push(std::mem::take(&mut current_batch));
            current_positions.clear();
            // Re-add this gene's positions to the fresh batch
            for &i in &gene.variant_indices {
                current_positions.insert(variants[i].position);
            }
            current_batch.push(gene);
        } else {
            // Gene fits — add it
            for p in new_positions {
                current_positions.insert(p);
            }
            current_batch.push(gene);
        }
    }

    if !current_batch.is_empty() {
        batches.push(current_batch);
    }
    batches
}


fn read_variant_metadata(engine: &DuckEngine) -> Result<Vec<AnnotatedVariant>, FavorError> {
    let conn = engine.connection();
    let mut stmt = conn.prepare(
        "SELECT chrom, pos, ref_allele, alt_allele, maf, gene_name, region_type, consequence, \
         cadd_phred, revel, \
         is_cage_promoter, is_cage_enhancer, is_ccre_promoter, is_ccre_enhancer, \
         w_cadd, w_linsight, w_fathmm_xf, \
         w_apc_epi_active, w_apc_epi_repressed, w_apc_epi_transcription, \
         w_apc_conservation, w_apc_protein_function, w_apc_local_nd, \
         w_apc_mutation_density, w_apc_tf \
         FROM _rare ORDER BY chrom, pos"
    ).map_err(|e| FavorError::Analysis(format!("{e}")))?;

    let mut rows = stmt.query([]).map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut variants = Vec::new();

    while let Ok(Some(row)) = rows.next() {
        let weights: Vec<f64> = (WEIGHT_COL_START..WEIGHT_COL_START + WEIGHT_COL_COUNT)
            .map(|c| row.get::<_, f64>(c).unwrap_or(0.0))
            .collect();
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
    Ok(variants)
}

fn read_phenotype_matrix(
    engine: &DuckEngine, sql: &str, n_cov: usize,
) -> Result<(Mat<f64>, Mat<f64>), FavorError> {
    let conn = engine.connection();
    let mut stmt = conn.prepare(sql).map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut rows = stmt.query([]).map_err(|e| FavorError::Analysis(format!("{e}")))?;

    let mut y_vec = Vec::new();
    let mut x_vecs: Vec<Vec<f64>> = vec![Vec::new(); n_cov];

    while let Ok(Some(row)) = rows.next() {
        y_vec.push(row.get::<_, f64>(0).unwrap_or(0.0));
        for j in 0..n_cov {
            x_vecs[j].push(row.get::<_, f64>(j + 1).unwrap_or(0.0));
        }
    }

    let n = y_vec.len();
    let mut y = Mat::zeros(n, 1);
    let mut x = Mat::zeros(n, 1 + n_cov);
    for i in 0..n {
        y[(i, 0)] = y_vec[i];
        x[(i, 0)] = 1.0; // intercept
        for j in 0..n_cov {
            x[(i, j + 1)] = x_vecs[j][i];
        }
    }
    Ok((y, x))
}

fn parquet_row_count(path: &Path) -> i64 {
    let file = match std::fs::File::open(path) { Ok(f) => f, Err(_) => return 0 };
    let reader = match parquet::file::reader::SerializedFileReader::new(file) { Ok(r) => r, Err(_) => return 0 };
    use parquet::file::reader::FileReader;
    reader.metadata().file_metadata().num_rows()
}

use std::path::Path;

use faer::Mat;

use crate::db::DuckEngine;
use crate::db::query_strings;
use crate::error::FavorError;
use crate::output::Output;
use crate::staar;
use crate::staar::genotype::GenotypeResult;

pub struct PhenotypeData {
    pub y: Mat<f64>,
    pub x: Mat<f64>,
    pub trait_type: staar::TraitType,
    pub n: usize,
}

pub fn load_phenotype(
    engine: &DuckEngine,
    phenotype: &Path,
    covariates: &[String],
    geno: &GenotypeResult,
    trait_name: &str,
    out: &dyn Output,
) -> Result<PhenotypeData, FavorError> {
    out.status("Step 3/6: Loading phenotype...");

    // Create _pheno table if it doesn't exist yet (first trait)
    engine.execute(&format!(
        "CREATE TEMP TABLE IF NOT EXISTS _pheno AS SELECT * FROM read_csv_auto('{}')",
        phenotype.display(),
    ))?;

    let pheno_cols = query_strings(engine, "SELECT column_name FROM (DESCRIBE _pheno)")?;
    if !pheno_cols.contains(&trait_name.to_string()) {
        return Err(FavorError::Input(format!(
            "Trait '{}' not in phenotype. Available: {}", trait_name, pheno_cols.join(", ")
        )));
    }
    for cov in covariates {
        if !pheno_cols.contains(cov) {
            return Err(FavorError::Input(format!("Covariate '{cov}' not in phenotype")));
        }
    }

    let trait_type = if crate::db::query_scalar(engine, &format!(
        "SELECT COUNT(DISTINCT \"{trait_name}\") FROM _pheno"))? <= 2 {
        staar::TraitType::Binary
    } else {
        staar::TraitType::Continuous
    };
    out.status(&format!("  Trait '{trait_name}' -> {:?}", trait_type));

    let id_col = &pheno_cols[0];
    let cov_select = if covariates.is_empty() { String::new() }
        else { format!(", {}", covariates.iter().map(|c| format!("p.\"{c}\"")).collect::<Vec<_>>().join(", ")) };

    let sample_list = geno.sample_names.iter().map(|s| format!("'{s}'")).collect::<Vec<_>>().join(",");
    let pheno_sql = format!(
        "SELECT p.\"{trait_name}\" {cov_select} FROM _pheno p \
         WHERE p.\"{id_col}\" IN ({sample_list}) AND p.\"{trait_name}\" IS NOT NULL \
         ORDER BY p.\"{id_col}\""
    );

    let (y, x) = read_phenotype_matrix(engine, &pheno_sql, covariates.len())?;
    let n = y.nrows();
    out.status(&format!("  {} samples with phenotype + genotype", n));
    if n < 10 {
        return Err(FavorError::Analysis(format!("Only {n} samples. Need >= 10.")));
    }

    Ok(PhenotypeData { y, x, trait_type, n })
}

pub fn load_ancestry_groups(
    engine: &DuckEngine,
    col: &str,
    geno: &GenotypeResult,
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

pub fn read_phenotype_matrix(
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

pub fn load_known_loci(
    engine: &DuckEngine,
    geno: &GenotypeResult,
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

    let extract_cols = staar::genotype::dosage_columns(n_samples);

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

pub fn augment_covariates(x: &Mat<f64>, x_cond: &Mat<f64>) -> Mat<f64> {
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

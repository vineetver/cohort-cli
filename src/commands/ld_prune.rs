//! `favor ld-prune` subcommand.
//!
//! Forward-selection LD pruning on conditional score-test p-values.
//! Mirrors STAARpipeline R/LD_pruning.R for the gaussian, unrelated,
//! single-trait path.

use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use serde_json::json;

use crate::commands;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::staar::ld_prune::{self, Candidate, KeptVariant, LdPruneParams};
use crate::staar::model::load_phenotype;
use crate::store::cohort::CohortId;
use crate::types::Chromosome;

const GB: u64 = 1024 * 1024 * 1024;

pub struct LdPruneArgs {
    pub cohort: String,
    pub phenotype: PathBuf,
    pub trait_name: String,
    pub covariates: Vec<String>,
    pub variants: PathBuf,
    pub maf_cutoff: f64,
    pub cond_p_thresh: f64,
    pub column_map: Vec<String>,
    pub output: Option<PathBuf>,
}

pub struct LdPruneConfig {
    pub cohort_id: CohortId,
    pub phenotype: PathBuf,
    pub trait_name: String,
    pub covariates: Vec<String>,
    pub variants: PathBuf,
    pub maf_cutoff: f64,
    pub cond_p_thresh: f64,
    pub column_map: HashMap<String, String>,
    pub output: PathBuf,
}

pub fn run(
    engine: &Engine,
    args: LdPruneArgs,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    let config = build_config(args)?;
    if dry_run {
        return emit_dry_run(&config, out);
    }
    run_ld_prune(engine, &config, out)
}

fn build_config(args: LdPruneArgs) -> Result<LdPruneConfig, CohortError> {
    if !args.phenotype.exists() {
        return Err(CohortError::Input(format!(
            "Phenotype file not found: '{}'",
            args.phenotype.display()
        )));
    }
    if !args.variants.exists() {
        return Err(CohortError::Input(format!(
            "Variants file not found: '{}'",
            args.variants.display()
        )));
    }
    if args.trait_name.trim().is_empty() {
        return Err(CohortError::Input("--trait-name is required".into()));
    }
    if args.maf_cutoff <= 0.0 || args.maf_cutoff >= 0.5 {
        return Err(CohortError::Input(format!(
            "MAF cutoff must be in (0, 0.5), got {}",
            args.maf_cutoff
        )));
    }
    if args.cond_p_thresh <= 0.0 || args.cond_p_thresh >= 1.0 {
        return Err(CohortError::Input(format!(
            "--cond-p-thresh must be in (0, 1), got {}",
            args.cond_p_thresh
        )));
    }

    let cohort_id = CohortId::new(args.cohort.trim().to_string());
    let output = args.output.unwrap_or_else(|| {
        PathBuf::from(format!("{}.ld_pruned.tsv", cohort_id.as_str()))
    });

    let mut column_map = HashMap::new();
    for entry in &args.column_map {
        let (k, v) = entry.split_once('=').ok_or_else(|| {
            CohortError::Input(format!(
                "Invalid --column-map entry '{entry}'. Expected key=value."
            ))
        })?;
        column_map.insert(k.trim().to_string(), v.trim().to_string());
    }

    Ok(LdPruneConfig {
        cohort_id,
        phenotype: args.phenotype,
        trait_name: args.trait_name,
        covariates: args.covariates,
        variants: args.variants,
        maf_cutoff: args.maf_cutoff,
        cond_p_thresh: args.cond_p_thresh,
        column_map,
        output,
    })
}

fn emit_dry_run(config: &LdPruneConfig, out: &dyn Output) -> Result<(), CohortError> {
    let n_candidates = std::fs::read_to_string(&config.variants)
        .map(|s| s.lines().filter(|l| !l.is_empty() && !l.starts_with('#')).count())
        .unwrap_or(0);
    let plan = commands::DryRunPlan {
        command: "ld-prune".into(),
        inputs: json!({
            "cohort_id": config.cohort_id.as_str(),
            "phenotype": config.phenotype.to_string_lossy(),
            "trait": config.trait_name,
            "covariates": config.covariates,
            "variants": config.variants.to_string_lossy(),
            "n_candidates": n_candidates,
            "maf_cutoff": config.maf_cutoff,
            "cond_p_thresh": config.cond_p_thresh,
        }),
        memory: commands::MemoryEstimate {
            minimum: "4G".into(),
            recommended: "8G".into(),
            minimum_bytes: 4 * GB,
            recommended_bytes: 8 * GB,
        },
        runtime: None,
        output_path: config.output.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

fn run_ld_prune(
    engine: &Engine,
    config: &LdPruneConfig,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let cohort = engine.cohort(&config.cohort_id);
    let store = cohort.load()?;

    let pheno = load_phenotype(
        engine.df(),
        &config.phenotype,
        &config.covariates,
        &store.geno,
        std::slice::from_ref(&config.trait_name),
        None,
        5,
        0,
        &config.column_map,
        out,
    )?;
    if !matches!(pheno.trait_type, crate::staar::TraitType::Continuous) {
        return Err(CohortError::Input(
            "ld-prune currently supports continuous traits only".into(),
        ));
    }

    let candidates_by_chrom = parse_candidates(&config.variants)?;
    if candidates_by_chrom.is_empty() {
        return Err(CohortError::Input(format!(
            "No variants parsed from '{}'",
            config.variants.display()
        )));
    }

    out.status(&format!(
        "ld-prune: {} candidate variants across {} chromosome(s)",
        candidates_by_chrom.values().map(|v| v.len()).sum::<usize>(),
        candidates_by_chrom.len()
    ));

    let params = LdPruneParams {
        maf_cutoff: config.maf_cutoff,
        cond_p_thresh: config.cond_p_thresh,
    };

    let mut kept_all: Vec<KeptVariant> = Vec::new();
    for (chrom, cands) in &candidates_by_chrom {
        let view = match cohort.chromosome(chrom) {
            Ok(v) => v,
            Err(e) => {
                out.warn(&format!("  chr{}: skipped ({e})", chrom.label()));
                continue;
            }
        };
        let kept = ld_prune::ld_prune_chromosome(
            &view,
            *chrom,
            &pheno.y,
            &pheno.x,
            &pheno.pheno_mask,
            cands,
            &params,
        )?;
        out.status(&format!(
            "  chr{}: {} / {} variants kept",
            chrom.label(),
            kept.len(),
            cands.len(),
        ));
        kept_all.extend(kept);
    }

    write_output(&config.output, &kept_all)?;
    out.success(&format!(
        "ld-prune: {} variants retained -> {}",
        kept_all.len(),
        config.output.display()
    ));
    out.result_json(&json!({
        "status": "ok",
        "output_path": config.output.to_string_lossy(),
        "n_kept": kept_all.len(),
    }));
    Ok(())
}

fn parse_candidates(
    path: &std::path::Path,
) -> Result<std::collections::BTreeMap<Chromosome, Vec<Candidate>>, CohortError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;

    let mut by_chrom: std::collections::BTreeMap<Chromosome, Vec<Candidate>> =
        std::collections::BTreeMap::new();

    for (lineno, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split(':').collect()
        };
        if parts.len() < 4 {
            return Err(CohortError::Input(format!(
                "variants file {}:{}: expected CHR{{:\\t}}POS{{:\\t}}REF{{:\\t}}ALT, got '{raw}'",
                path.display(),
                lineno + 1
            )));
        }
        let chrom: Chromosome = parts[0].parse().map_err(|e: String| {
            CohortError::Input(format!(
                "variants file {}:{}: {e}",
                path.display(),
                lineno + 1
            ))
        })?;
        let position: u32 = parts[1].parse().map_err(|e| {
            CohortError::Input(format!(
                "variants file {}:{}: bad position '{}': {e}",
                path.display(),
                lineno + 1,
                parts[1]
            ))
        })?;
        by_chrom.entry(chrom).or_default().push(Candidate {
            position,
            ref_allele: parts[2].to_string(),
            alt_allele: parts[3].to_string(),
        });
    }

    Ok(by_chrom)
}

fn write_output(path: &std::path::Path, kept: &[KeptVariant]) -> Result<(), CohortError> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", parent.display()))
            })?;
        }
    }
    let mut f = std::fs::File::create(path)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", path.display())))?;
    writeln!(f, "CHR\tPOS\tREF\tALT\tentry_log10p")
        .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    for v in kept {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            v.chromosome.label(),
            v.position,
            v.ref_allele,
            v.alt_allele,
            format_log10p(v.entry_log10p),
        )
        .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    }
    Ok(())
}

fn format_log10p(lp: f64) -> String {
    if lp.is_infinite() {
        "Inf".into()
    } else {
        format!("{:.6}", lp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_candidates_colon_and_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("candidates.txt");
        std::fs::write(
            &path,
            "# header\n1:100:A:T\n1\t200\tG\tC\n2:300:C:A\n\n",
        )
        .unwrap();
        let map = parse_candidates(&path).unwrap();
        assert_eq!(map.get(&Chromosome::Autosome(1)).unwrap().len(), 2);
        assert_eq!(map.get(&Chromosome::Autosome(2)).unwrap().len(), 1);
    }

    #[test]
    fn parse_candidates_rejects_short_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.txt");
        std::fs::write(&path, "1:100:A\n").unwrap();
        assert!(parse_candidates(&path).is_err());
    }

    #[test]
    fn build_config_rejects_bad_maf() {
        let dir = tempfile::tempdir().unwrap();
        let pheno = dir.path().join("pheno.tsv");
        let vars = dir.path().join("vars.tsv");
        std::fs::write(&pheno, "id\ttrait\n").unwrap();
        std::fs::write(&vars, "1:100:A:T\n").unwrap();
        let args = LdPruneArgs {
            cohort: "c1".into(),
            phenotype: pheno.clone(),
            trait_name: "trait".into(),
            covariates: vec![],
            variants: vars.clone(),
            maf_cutoff: 0.8,
            cond_p_thresh: 1e-4,
            column_map: vec![],
            output: None,
        };
        assert!(matches!(build_config(args), Err(CohortError::Input(_))));
    }
}

use std::path::PathBuf;

use serde_json::json;

use crate::config::{Config, DirProbe, Tier};
use crate::error::FavorError;
use crate::mode::OutputMode;
use crate::output::Output;
use crate::packs::Pack;

pub fn run(
    path: Option<PathBuf>,
    force: bool,
    out: &dyn Output,
    mode: &OutputMode,
) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let project_dir = match path {
        Some(p) => {
            std::fs::create_dir_all(&p)?;
            p.canonicalize().map_err(|e| {
                FavorError::Input(format!("Cannot resolve '{}': {e}", p.display()))
            })?
        }
        None => std::env::current_dir()?,
    };

    let claude_path = project_dir.join("CLAUDE.md");
    let codex_dir = project_dir.join(".codex");
    let codex_path = codex_dir.join("instructions.md");

    let is_refresh = claude_path.exists();

    if is_refresh && !force && !mode.is_machine() {
        out.warn(
            "CLAUDE.md already exists. Use --force to overwrite, or re-run to refresh after installing new packs.",
        );
        return Ok(());
    }

    let probe = DirProbe::scan(&config.root_dir());
    let content = render(&config, &probe);

    std::fs::write(&claude_path, &content)?;
    std::fs::create_dir_all(&codex_dir)?;
    std::fs::write(&codex_path, &content)?;

    let action = if is_refresh { "Refreshed" } else { "Initialized" };
    out.success(&format!("{action} agent context in {}", project_dir.display()));
    out.status("  CLAUDE.md (Claude Code)");
    out.status("  .codex/instructions.md (Codex)");
    out.status(&format!(
        "Open your coding agent in {} and ask a research question.",
        project_dir.display()
    ));

    out.result_json(&json!({
        "project_dir": project_dir.to_string_lossy(),
        "files": [claude_path.to_string_lossy(), codex_path.to_string_lossy()],
        "tier": config.data.tier.as_str(),
        "packs": config.data.packs,
    }));

    Ok(())
}

fn render(config: &Config, probe: &DirProbe) -> String {
    let pack_list = build_pack_list(config, probe);

    let environment = config
        .resources
        .environment
        .map(|e| e.as_str().to_string())
        .unwrap_or_else(|| "auto-detect".to_string());

    let memory_budget = config
        .resources
        .memory_budget
        .clone()
        .unwrap_or_else(|| "auto-detect".to_string());

    let chrom_count = match config.data.tier {
        Tier::Base => probe.base_chroms,
        Tier::Full => probe.full_chroms,
    };

    TEMPLATE
        .replace("{root_dir}", &config.data.root_dir)
        .replace("{tier}", config.data.tier.as_str())
        .replace("{tier_size}", config.data.tier.size_human())
        .replace("{chrom_count}", &chrom_count.to_string())
        .replace("{environment}", &environment)
        .replace("{memory_budget}", &memory_budget)
        .replace("{pack_list}", &pack_list)
}

fn build_pack_list(config: &Config, probe: &DirProbe) -> String {
    let mut lines = Vec::new();
    for pack in Pack::all() {
        if pack.always_installed {
            continue;
        }
        let installed = config.data.packs.contains(&pack.id.to_string())
            || pack
                .tables
                .iter()
                .any(|t| probe.tissue_tables.contains(&t.to_string()));
        if installed {
            lines.push(format!(
                "- **{}** ({}): {}",
                pack.id, pack.size_human, pack.description
            ));
        }
    }
    if lines.is_empty() {
        lines.push(
            "- No optional packs. Run `favor data pull --pack <id>` then `favor init --force`."
                .to_string(),
        );
    }
    lines.join("\n")
}

const TEMPLATE: &str = r#"# FAVOR Genomic Analysis

Always pass `--format json` to every favor command. Use `--dry-run` before heavy computation.

## System

- Data root: {root_dir}
- Tier: {tier} ({tier_size}, {chrom_count}/24 chromosomes)
- Environment: {environment}
- Memory budget: {memory_budget}

## Installed packs

{pack_list}

## Pipeline

```
variant file -> favor ingest -> favor annotate -> favor enrich --tissue X
                                                -> favor staar (rare-variant association)
                                                -> favor interpret (variant-to-gene)
```

## Commands

| Command | What it does |
|---------|-------------|
| `favor ingest <file>` | Normalize VCF/TSV/CSV to parquet with variant ID |
| `favor annotate <file>` | Add CADD, ClinVar, gnomAD, REVEL, aPC scores, regulatory marks |
| `favor enrich <file> --tissue <name>` | Add tissue eQTL, ChromBPNet, enhancer-gene links |
| `favor staar --genotypes <vcf> --phenotype <tsv> --trait-name <col> --annotations <parquet>` | Rare-variant burden test (STAAR-O) |
| `favor schema [table]` | Show column names and types |
| `favor manifest` | List installed data and available commands |

## STAAR usage

```bash
favor staar --dry-run --format json \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2 \
  --annotations annotated.parquet \
  --masks coding
```

Use `--dry-run` first to check memory. On HPC: `srun --mem=64G -c 8 favor staar ...`

## Interpreting results

- **CADD phred > 20**: top 1% most deleterious variants genome-wide
- **REVEL > 0.5**: likely pathogenic missense
- **STAAR-O < 2.5e-6**: genome-wide significant (Bonferroni for ~20K genes)
- **pLoF mask significant**: protein-destroying variants cause the trait (strongest evidence)
- **Synonymous mask significant**: suspicious — likely confounding, not biology

STAAR runs 6 tests (Burden, SKAT, ACAT-V x two beta weights) across 11 annotation channels:
- Burden: all variants push trait same direction
- SKAT: mixed-effect variants
- ACAT-V: one or two variants drive everything

## Querying parquet output

Use DuckDB:

```sql
-- Damaging variants
SELECT vid, gencode.genes[1] AS gene, main.cadd.phred
FROM 'annotated.parquet' WHERE main.cadd.phred > 20;

-- Significant genes
SELECT * FROM 'staar_results/coding_pLoF_missense.parquet'
WHERE "STAAR-O" < 2.5e-6 ORDER BY "STAAR-O";

-- Tissue enrichment
SELECT a.vid, e.* FROM 'enriched/annotated.parquet' a
JOIN 'enriched/eqtl.parquet' e ON a.vid = e.vid
WHERE e.tissue_name LIKE '%Liver%';
```

## Output conventions

- Parquet files: main results
- `.meta.json`: parameters and counts
- stdout (--format json): structured result summary
- stderr: progress/status messages
- Exit codes: 0=ok, 1=input, 2=data missing, 3=resource, 4=analysis
"#;

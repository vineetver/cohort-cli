# Setup Guide

> [Back to README](../README.md)

## Install

```bash
curl -fsSL https://raw.githubusercontent.com/vineetver/favor-cli/master/install.sh | sh
```

From source:

```bash
cargo install --git https://github.com/vineetver/favor-cli
```

The binary installs to `~/.local/bin/favor` by default. Override with `FAVOR_INSTALL_DIR`.

## Configuration

`favor setup` writes `~/.cohort/config.toml` and downloads annotation data.

```bash
# base tier: smaller, enough for annotation and STAAR
favor setup --root /data/favor --tier base

# full tier: more annotation columns (REVEL, CAGE, additional aPC weights)
favor setup --root /data/favor --tier full

# with optional enrichment packs
favor setup --root /data/favor --tier full --packs eqtl,regulatory

# HPC with memory budget
favor setup --root /data/favor --tier full --environment hpc --memory-budget 64GB
```

The config file:

```toml
[data]
tier = "full"
root_dir = "/data/favor"
packs = ["eqtl", "regulatory"]

[resources]
environment = "hpc"
memory_budget = "64GB"
```

## Annotation tiers

| Tier | Size | What you get | What you can run |
|------|------|-------------|-----------------|
| **base** | ~200 GB | CADD, cCRE flags, gene assignments, aPC weights | ingest, annotate, enrich, STAAR |
| **full** | ~508 GB | base + REVEL, CAGE, additional aPC channels | all commands, richer annotation |

Both tiers support STAAR. Full tier adds more annotation weight channels for a richer omnibus test, but base tier has everything needed to run the pipeline end-to-end.

## Enrichment packs

Packs add tissue-specific data for `favor enrich`. Download them with `favor data pull --pack <id>`.

| Pack | Size | Source | Contents |
|------|------|--------|----------|
| reference | 40 MB | FAVOR | Gene index, cCRE regions, tissue vocabulary |
| rollups | 49 MB | FAVOR | Gene-level and region summaries |
| variant-index | 155 GB | FAVOR | Variant-to-region junction table |
| eqtl | 3 GB | GTEx v10 | Bulk eQTL, sQTL, apaQTL across 50+ tissues |
| eqtl-catalogue | 2 GB | EBI | 127 cell types from eQTL Catalogue |
| sc-eqtl | 48 GB | OneK1K, DICE, PsychENCODE | Single-cell eQTL |
| regulatory | 18 GB | ENCODE SCREEN v4, Roadmap | Chromatin states, accessibility, validated enhancers |
| enhancer-gene | 12 GB | ABC, EPIraction, rE2G, EpiMap | Enhancer-gene link predictions |
| tissue-scores | 5 GB | ChromBPNet, TLand, ENTEx | Tissue-specific variant effect scores |
| pgs | 75 GB | PGS Catalog | Polygenic risk scores for 2000+ phenotypes |
| genotypes | 3 GB | TOPMed | Reference LD genotype matrix |

The first three (reference, rollups, variant-index) are always installed during setup. The rest are optional.

```bash
# pull everything configured during setup
favor data pull

# add a pack later
favor data pull --pack eqtl

# check what's installed
favor data status

# verify checksums
favor data verify --checksums
```

## Data root layout

After setup and data pull, your data root looks like this:

```
/data/favor/
  full/
    chromosome=1/sorted.parquet
    chromosome=2/sorted.parquet
    ...
    chromosome=Y/sorted.parquet
  tissue/
    reference/
    rollups/
    variant_in_region/
    variant_eqtl/          # if eqtl pack installed
    variant_sqtl/
    ...
```

This directory is shared across all projects. Point multiple studies at the same data root. It's read-only after setup.

## Project directory structure

Each study gets its own working directory. Keep input files, intermediate outputs, and results together:

```
my_gwas_study/
  raw/
    cohort_chr1.vcf.gz     # input VCFs
    cohort_chr2.vcf.gz
    ...
    phenotypes.tsv          # sample_id, trait, covariates
  
  variants.ingested/        # favor ingest output
  variants.annotated/       # favor annotate output
  variants.enriched/        # favor enrich output (optional)
  
  staar_results/            # favor staar output
    coding_pLoF_missense.parquet
    coding_synonymous.parquet
    noncoding_promoter.parquet
    staar_run.meta.json
  
  .cohort/                  # auto-created project store
    cohorts/my_cohort/
      manifest.json
      samples.txt
      chromosome=*/
        sparse_g.bin
        variants.parquet
        membership.parquet
    cache/score_cache/      # reused across reruns
```

## Store resolution

The `.cohort/` store holds genotype matrices and caches. FAVOR CLI finds it by checking (in order):

1. `--store-path /path/to/.cohort` (explicit CLI flag)
2. `FAVOR_STORE=/path/to/.cohort` (environment variable)
3. Walk up from cwd looking for an existing `.cohort/` directory
4. Fall back to `<cwd>/.cohort/` (created on first use)

For HPC jobs, set `FAVOR_STORE` in your job script so the store is found regardless of the working directory.

## HPC setup

### Never build on login nodes

```bash
# build from source on a compute node
srun -p <partition> -c 8 --mem=16G cargo build --release
```

### SLURM job script

```bash
#!/bin/bash
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --tmp=100G

export TMPDIR=$TMPDIR   # SLURM sets this to local scratch
export FAVOR_STORE=/path/to/project/.cohort

favor staar \
  --genotypes /data/vcfs/cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2,PC3,PC4,PC5 \
  --annotations /path/to/variants.annotated \
  --masks coding,noncoding \
  --output results/
```

FAVOR CLI auto-detects SLURM allocations (`SLURM_MEM_PER_NODE`, `SLURM_CPUS_PER_TASK`) and sizes its memory pools accordingly. No need to pass `--threads` inside a SLURM job.

### Temp directory

Large intermediate files go to the temp directory. FAVOR CLI checks (in order): config `temp_dir`, `TMPDIR`, `SCRATCH`, `LOCAL_SCRATCH`, `/scratch`, system default.

On HPC, make sure this points to local SSD or fast scratch, not a network filesystem:

```bash
export TMPDIR=/scratch/$USER/tmp
```

### Memory

FAVOR CLI adapts to available memory. Use `--dry-run` to see what a command needs before committing resources:

```bash
favor staar --dry-run --format json \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --annotations variants.annotated
# prints recommended memory, estimated runtime, cache status
```

For kinship-aware analysis, the dense REML working set scales as `5 * n_samples^2 * 8 bytes`. The budget defaults to 16 GB but can be raised:

```bash
export FAVOR_KINSHIP_MEM_GB=64
```

## Environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `FAVOR_STORE` | Override store root path | walk up for `.cohort/`, then `<cwd>/.cohort/` |
| `FAVOR_MACHINE` | Force JSON output mode | unset (auto-detect tty) |
| `FAVOR_THREADS` | Override thread count | SLURM/cgroup/nproc auto-detect |
| `FAVOR_KINSHIP_MEM_GB` | Dense REML memory cap | 16 |
| `FAVOR_INSTALL_DIR` | Install location for the binary | `~/.local/bin` |

## Keeping things clean

**Annotation data** (`/data/favor/`) is immutable after download. Share it across projects.

**Project stores** (`.cohort/`) accumulate caches. Clean up with:

```bash
# remove orphaned caches from deleted cohorts
favor store gc

# nuclear option: delete the whole store and rebuild
rm -rf .cohort/
```

**Intermediate outputs** (`.ingested/`, `.annotated/`, `.enriched/`) can be deleted after STAAR runs if you don't need them. The genotype store in `.cohort/cohorts/` is the durable artifact. Rebuilding from VCF is the expensive step; everything else is cheap.

**Score caches** under `.cohort/cache/score_cache/` are keyed by cohort + trait + covariates + kinship. Changing masks or MAF cutoffs reuses the cache. Changing the trait or covariates builds a new one. Old caches are safe to delete.

## Uninstall

```bash
favor uninstall
# or just: rm ~/.local/bin/favor && rm -rf ~/.cohort/
```

This removes the binary and config. Your data root and project stores are untouched.

## Pipeline walkthrough

### Sites-only annotation (no genotypes)

```bash
favor ingest variants.vcf.gz
favor annotate variants.ingested --full
# query results with any parquet tool:
# SELECT * FROM 'variants.annotated/chromosome=1/data.parquet' WHERE cadd_phred > 20
```

### Full STAAR analysis

```bash
# 1. ingest variant sites
favor ingest cohort.vcf.gz --output sites.ingested

# 2. annotate with full tier
favor annotate sites.ingested --full

# 3. run STAAR (builds genotype store on first run, reuses on reruns)
favor staar \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2 \
  --annotations sites.annotated \
  --masks coding,noncoding \
  --output ldl_results/

# 4. second trait reuses the genotype store, only rebuilds score cache
favor staar \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name HDL \
  --covariates age,sex,PC1,PC2 \
  --annotations sites.annotated \
  --masks coding \
  --output hdl_results/
```

### Using a pre-built cohort

If you already built the cohort store (via `favor ingest --cohort-id`), skip the VCF:

```bash
favor staar \
  --cohort my_cohort \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2 \
  --masks coding
```

### Tissue enrichment

```bash
favor enrich variants.annotated --tissue Liver --output liver_enriched/
# join tissue data by vid:
# SELECT a.*, e.* FROM 'liver_enriched/annotated.parquet' a
# JOIN 'liver_enriched/eqtl.parquet' e ON a.vid = e.vid
```

### Cross-study meta-analysis

```bash
# each study exports summary statistics
favor staar --emit-sumstats ... --output study1/
favor staar --emit-sumstats ... --output study2/

# combine
favor meta-staar --studies study1/,study2/ --masks coding --output meta_results/
```

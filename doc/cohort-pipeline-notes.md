# FAVOR Pipeline Notes

Internal-only notes. This file is intentionally ignored by git.

## Why this exists

This is a plain-English map of the current pipeline as implemented in the Rust code today.
The goal is not to explain the statistics in depth. The goal is to answer:

- What are the stages?
- What goes in and out of each stage?
- What gets cached?
- Where are the memory / I/O hotspots?
- Which paths are actually wired today vs only partially implemented?

## One-screen view

```text
raw variants / cohort VCF
    |
    | favor ingest
    v
ingested variant set
    |
    | favor annotate
    v
annotated variant set
    |
    | favor staar
    |   Step 1: build or reuse genotype store
    |   Step 2: fit null model
    |   Step 2b: build or reuse score cache
    |   Step 3: run mask tests
    |   Step 4: write results
    v
STAAR result parquet files
```

## Simple mental model

- `ingest` is format cleanup and normalization.
- `annotate` is a big join against FAVOR annotation data.
- `staar` is not one giant monolith. It is a 3-layer system:
  - Layer 1: build a reusable sparse genotype store.
  - Layer 2: build a reusable per-phenotype score cache.
  - Layer 3: run many mask tests cheaply by slicing cached arrays.

If you only remember one thing, remember this:

```text
expensive once per VCF          expensive once per phenotype          cheap many times
----------------------         -------------------------------       ----------------
build genotype store     ->    build score cache               ->    test masks / windows
```

That split is the main architectural defense against spaghetti and repeated work.

## Current pipeline, stage by stage

### 1. `favor ingest`

Purpose:
- Accept VCF / TSV / CSV / parquet inputs.
- Normalize them into a standard variant-set directory.

Input:
- Raw variant files or a directory of them.

Output:
- `<name>.ingested/`
- `meta.json`
- `chromosome=*/data.parquet`

What it does:
- Detects input format.
- Expands directories into matching files.
- For VCF:
  - normalizes contigs
  - splits multi-allelic sites into biallelic records
  - writes a variant-set directory

Performance meaning:
- This is mostly parse + write work.
- It is not the main rare-variant math bottleneck.
- Memory estimate is generic and modest compared with STAAR.

### 2. `favor annotate`

Purpose:
- Join ingested variants against the FAVOR annotation database.

Input:
- `<name>.ingested/`

Output:
- `<name>.annotated/`
- `meta.json`
- `chromosome=*/data.parquet`

What it does:
- Opens the ingested variant set.
- Validates required join columns exist: chromosome, position, ref, alt.
- Joins per chromosome against FAVOR parquet.
- Produces annotated variant rows.

Performance meaning:
- This is large parquet join work.
- Main cost is I/O + query engine join cost.
- It is much more about data movement than dense math.

Important STAAR detail:
- STAAR requires full-tier annotations, not base-tier only.
- The code checks that the 11 STAAR annotation channels are present.

### 3. `favor staar`

Purpose:
- Run rare-variant association tests using the annotated variants plus a multi-sample genotype VCF and phenotype table.

Inputs:
- `--genotypes <multi-sample VCF>`
- `--phenotype <TSV>`
- `--trait-name <name>`
- `--covariates ...`
- `--annotations <annotated variant set>`

Main outputs:
- `<vcf>.staar/store/` or chosen output dir store
- `<vcf>.staar/results/<trait>/`
- optional `<vcf>.staar/sumstats/<trait>/` when `--emit-sumstats`

### 3A. STAAR Step 1: build or reuse genotype store

Wire:

```text
multi-sample VCF + annotated variants
    |
    | extract_genotypes()
    v
packed genotype parquet by chromosome
    |
    | join with annotations
    v
annotated rare-variant table
    |
    | build sparse store
    v
store/
  manifest.json
  samples.txt
  chromosome=1/
    sparse_g.bin
    variants.parquet
    membership.parquet
  ...
```

What this means in simple words:

- The tool reads the cohort VCF once and converts each variant into dosages across samples.
- It joins those cohort variants with the annotated FAVOR rows.
- Then it rewrites the result into a custom sparse on-disk structure optimized for rare variants.

Why sparse matters:

- Rare variants have very few carriers.
- Instead of storing a huge sample x variant matrix full of zeros, `sparse_g.bin` stores only the carriers.
- This changes work from "touch every sample for every variant" to "touch only carrier entries".

Key files:

- `sparse_g.bin`
  - one carrier list per variant
  - memory-mapped at runtime
  - O(1) access by `variant_vcf`
- `variants.parquet`
  - one aligned metadata row per variant
- `membership.parquet`
  - variant-to-gene mapping

Cache key:
- Based on the genotype VCF fingerprint plus annotation fingerprint.
- If either changes, the store rebuilds.

Performance meaning:

- This is the big up-front cost.
- Heavy on:
  - VCF parsing
  - parquet writing
  - annotation join
  - full-store materialization
- It is worth caching because later runs can skip it entirely.

### 3B. STAAR Step 2: fit null model

Purpose:
- Fit the background model before testing genes/windows.

Input:
- phenotype column
- covariates
- optional kinship inputs

Output:
- in-memory model objects, not a user-facing dataset

Simple meaning:
- This estimates what the phenotype looks like before asking whether a gene/window adds signal.

Performance meaning:
- This can be memory-heavy when kinship is used.
- Without kinship, this is much simpler.
- With kinship, dense matrix burden can become important.

### 3C. STAAR Step 2b: build or reuse score cache

Wire:

```text
store/ + null model
    |
    | for each chromosome
    | for each gene
    v
score_cache/<cache_key>/chromosome=*/scores.bin
```

What gets cached:

- `U`: one score value per variant
- `K`: one covariance matrix per gene, when gene is not too large

Simple meaning:
- The expensive per-phenotype math is done once.
- Later mask runs mostly slice these cached results instead of rereading genotype carriers.

Why this is a big deal:

- Changing mask definition or MAF cutoff does not force a full recomputation.
- It usually only repeats the cheap testing layer.

Current fallback:
- Very large genes store `U` only.
- Their `K` matrix is recomputed on demand from carrier data.

Performance meaning:
- This is the bridge between "expensive phenotype-specific setup" and "cheap repeated testing".
- Good target for reuse across experiments.

### 3D. STAAR Step 3: run mask tests

Current supported mask categories in the code:

- coding
- noncoding
- sliding-window
- scang

Simple meaning:
- A "mask" is just a rule for selecting a subset of variants together.
- Example: all rare coding variants in one gene, or all rare variants in a sliding window.

Important implementation detail:

- Standard gene masks usually run from cached `U/K`.
- SPA and AI-STAAR reload carrier lists for the selected variants.
- Very large genes also reload carriers because `K` was not cached.

So the testing layer is not always purely cache-only. It is:

```text
fast path:
  slice cached U/K -> test

slower fallback path:
  load carrier lists from sparse_g.bin -> recompute subset stats -> test
```

### 3E. STAAR Step 4: write results

Outputs:
- one parquet file per mask type
- under `results/<trait>/`

Examples:
- `coding_*.parquet`
- `noncoding_*.parquet`
- `sliding_window.parquet`
- `scang.parquet`

Contents:
- gene/window identifier
- chromosome/start/end
- number of variants
- cumulative MAC
- many p-values

## How the main STAAR variants differ

### Plain STAAR

What it is:
- Standard rare-variant test battery on a selected group of variants.

Current execution:
- Usually uses cached `U/K` for gene masks.
- Window masks use chromosome-wide cached `U` plus assembled `K`.

### AI-STAAR

What it is in simple words:
- Same rare-variant idea, but genotype contributions are reweighted by ancestry/population group.

Extra input:
- `--ancestry-col` in phenotype file

Current execution:
- For a selected gene/window, carrier lists are loaded from `sparse_g.bin`.
- The code builds dense matrices for the selected variants and runs multiple weighted STAAR runs, then combines them.

Performance meaning:
- More expensive than plain cached STAAR.
- More repeated per-mask work.
- Less able to stay entirely in the cheap "slice cached U/K" path.

### MultiSTAAR

What it should mean:
- Joint rare-variant testing across multiple traits.

Current repo reality:
- There is a `src/staar/multi.rs` implementation of the combine logic.
- But the main `favor staar` pipeline currently only uses `trait_names[0]`.
- I do not see the main pipeline actually dispatching MultiSTAAR end-to-end today.

Practical takeaway:
- Treat MultiSTAAR as partially implemented / not wired into the main path yet.
- If you are cleaning architecture, this is a good place to clarify intended behavior before optimizing.

### SCANG

What it is in simple words:
- Instead of predefining one fixed mask per gene, scan many windows of different sizes and test them.

Current execution:
- Builds many window groups per chromosome.
- For non-AI path, uses cached chromosome `U` and assembled `K`.
- For AI-STAAR path, reloads carriers for each window.

Performance meaning:
- Can explode the number of groups tested.
- Even with cheap per-window math, total runtime can grow fast because there are many windows.

### MetaSTAAR

Purpose:
- Cross-study meta-analysis from summary statistics, not raw per-study VCFs.

Wire:

```text
study 1: cohort staar --emit-sumstats
study 2: cohort staar --emit-sumstats
study 3: cohort staar --emit-sumstats
    |
    v
meta_staar study directories
    |
    | favor meta-staar
    v
merged meta-analysis result parquet files
```

Per-study sumstats output:
- `meta_staar.json`
- chromosome-level segment parquet files

Simple meaning:
- Each study exports enough summary information that a later stage can merge studies without rerunning raw genotype analysis together.

Performance meaning:
- Good for scaling across cohorts.
- Shifts work from raw genotype access to summary-stat merging.

## Answering the practical questions

### Does STAAR take batches and do vectorized math, or does it run one-by-one?

Short answer:
- Both, depending on the layer.

More concrete answer:

- VCF extraction is batched during parquet writing.
- Sparse genotype loading supports batch access with internal sorting for more sequential mmap access.
- Score-cache building works gene-by-gene, but inside each gene the math is vectorized over the selected variants.
- Standard gene-mask testing is not "variant one by one"; it slices arrays for the whole selected group.
- SCANG tests many windows, each window as a group.
- AI-STAAR and SPA paths reload carriers for each selected group and do more per-group work.

So the real picture is:

```text
outer loop: per chromosome / per gene / per window
inner work: vector/matrix math over the selected variant group
```

It is not one-variant-at-a-time statistics, but it is also not one giant whole-genome dense matrix multiply.

### Where is memory pressure highest?

Roughly:

1. VCF extraction and genotype staging for large cohorts
2. Kinship-aware null model fitting
3. Building large `K` matrices for genes with many variants
4. AI-STAAR dense reconstruction for selected groups

The code explicitly avoids caching `K` for very large genes, which is a sign that this matrix is a real memory concern.

### Where is I/O pressure highest?

Roughly:

1. initial VCF parse
2. annotate join
3. store build writes
4. score cache write
5. repeated carrier reloads in AI-STAAR / SPA / large-gene fallback / some window paths

### What actually gets reused?

- genotype store is reused across runs with same VCF + annotations
- score cache is reused across runs with same store + trait + covariates + known loci + kinship setup
- changing mask definitions or MAF cutoff usually reuses existing store and cache

That is the main anti-spaghetti structure already present in the code.

## Current architecture strengths

- Clear 3-layer split in STAAR
- Sparse genotype representation is a strong fit for rare variants
- Stable `variant_vcf` indexing avoids repeated joins in hot paths
- Store cache and score cache are both content-keyed
- Batch carrier loading sorts requests for better sequential mmap behavior

## Current architectural rough edges

- `favor staar` CLI exposes multiple trait names, but the main pipeline only uses the first trait today
- MultiSTAAR logic exists but is not clearly wired into the main end-to-end path
- `--no-store` is parsed in CLI dispatch but ignored in `main.rs`
- some "fast cached" paths fall back to carrier reloads depending on options
- AI-STAAR path reconstructs dense matrices for selected groups, which may weaken the sparse-first architecture
- SCANG can create a very large number of windows, so control-flow complexity and runtime can grow quickly

## If the goal is "make it faster / cleaner", look here first

### 1. Make the fast path more explicit

Separate these cases clearly in docs and code:

- cache-only testing
- sparse carrier reload testing
- dense reconstruction testing

That alone makes optimization work much easier to reason about.

### 2. Decide what MultiSTAAR should be

Right now the API suggests support, but the main pipeline still behaves as single-trait-first.
Before optimizing, decide:

- should `favor staar` truly support multi-trait end-to-end?
- or should MultiSTAAR become a separate explicit command/path?

### 3. Reduce repeated carrier reloads

Likely wins:

- reuse loaded carrier subsets across multiple masks for the same gene
- reuse per-window carrier material when scanning overlapping windows
- avoid dense reconstruction unless a method truly requires it

### 4. Make performance mode visible

For each run, it would help to log counts like:

- genes tested from cached `U/K`
- genes requiring sparse reload
- genes requiring large-gene `K` recomputation
- windows tested
- windows tested under AI-STAAR path

That would make hotspots obvious.

### 5. Keep boundaries strict

The cleanest shape is:

- ingest = normalization only
- annotate = join only
- store build = genotype + annotation materialization only
- score cache = phenotype-specific preprocessing only
- tests = mask/window selection plus final statistics only

Any cross-layer leakage should be treated as suspicious.

## Suggested simple vocabulary

If statistical terms are getting in the way, use these translations:

- null model = baseline phenotype model
- mask = a group of variants tested together
- U = per-variant score summary
- K = per-group covariance summary
- MAF cutoff = only keep variants rarer than this threshold
- cMAC = approximate total minor-allele count across the tested group
- SCANG = scan many windows instead of only gene boundaries
- MetaSTAAR = combine study summaries instead of combining raw genotypes

## Bottom line

The current pipeline is not random spaghetti. It already has a strong reusable shape:

```text
normalize variants
  -> annotate variants
  -> build sparse genotype store once
  -> build phenotype-specific score cache once
  -> run many tests cheaply from cached summaries
```

The main cleanup opportunities are not in the broad architecture. They are in:

- making the real execution modes explicit
- clarifying partially wired paths like MultiSTAAR
- reducing repeated carrier reloads
- instrumenting which groups hit the cheap cached path vs expensive fallback paths

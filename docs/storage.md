# Genotype Store

COHORT stores rare-variant genotypes as a sparse matrix over `(sample_id, variant_vcf) -> dosage`.
The main idea is simple:

- every variant in a chromosome gets a dense internal index, `variant_vcf`
- metadata, masks, weights, and score-cache entries all line up to that same index
- carrier data is stored once per variant, not once per sample

`variant_vcf` is an internal array index. Genomic positions stay 1-based, VCF-style.

> [Back to README](../README.md) · [Performance](performance.md)

## Core Model

For one chromosome with `N` variants and `S` samples:

```text
variant_vcf: [0, 1, 2, ..., N-1]
```

Everything is aligned to it:

```text
position[v]
ref_allele[v]
alt_allele[v]
vid[v]                 "chr-pos-ref-alt"
maf[v]
region_type[v]
consequence[v]
weights[channel][v]
G[s, v]                dosage for carriers only
```

That gives the storage and query model:

```text
gene("BRCA2")      -> [variant_vcf]
maf < 0.01         -> mask over variant_vcf
is_plof()          -> mask over variant_vcf
carriers(v)        -> sparse carrier list for variant_vcf v
```

Most analysis reduces to:

```text
mask     = compile(gene, maf_cutoff, annotation predicate)
weights  = weights[mask]
carriers = G[:, mask]
result   = score(carriers, weights)
```

## Invariants

`variant_vcf` is dense, 0-based, and fixed for the lifetime of a store.

1. `variant_vcf` is assigned at build time, sorted by `(position, ref_allele, alt_allele)`.
2. If the variant universe changes, the store is rebuilt and dependent caches are invalidated.
3. `variant_vcf = 42` always means the same physical variant inside one store build.
4. Metadata, masks, and score-cache entries all share this coordinate space.

This is why the current store is fast for gene- and mask-driven work, and also why incremental ingest is still a gap.

## On-Disk Layout

```text
.genotype_store/
  manifest.json
  samples.txt
  chromosome={chr}/
    sparse_g.bin
    variants.parquet
    membership.parquet
```

### `manifest.json`

Store-level metadata:

- schema version
- content fingerprint
- sample count
- variant count
- chromosome list
- build timestamp

The content key is derived from the VCF fingerprint and the annotation fingerprint. If either input changes, the current store is treated as stale.

### `samples.txt`

One sample ID per line, in VCF column order. Line number equals the sample index used in `sparse_g.bin`.

### `sparse_g.bin`

This is the sparse genotype matrix.

- header: counts, flags, offsets location
- carrier blocks: one variable-length block per `variant_vcf`
- offset table: byte offset for every variant

Carrier entry format:

| Mode | When | Size | Layout |
|---|---|---|---|
| narrow | `S <= 65,535` | 3 bytes | `u16 sample_idx` + `u8 dosage` |
| wide | `S > 65,535` | 5 bytes | `u32 sample_idx` + `u8 dosage` |

Only non-reference carriers are stored.

To read one variant:

1. read `offsets[v]`
2. jump to that carrier block
3. read `n_carriers`
4. read `n_carriers` entries

Point lookup by `variant_vcf` is O(1).

### `variants.parquet`

One row per `variant_vcf`.

Current columns:

- identity: `variant_vcf`, `position`, `ref_allele`, `alt_allele`, `vid`
- annotations: `maf`, `region_type`, `consequence`, `cadd_phred`, `revel`
- regulatory flags: CAGE and cCRE booleans
- 11 STAAR weight channels

The weight columns are still `Float64` today. That is likely heavier than necessary.

### `membership.parquet`

Many-to-many gene membership:

| Column | Type |
|---|---|
| `variant_vcf` | UInt32 |
| `gene_name` | Utf8 |

One row per `(variant, gene)` pair, sorted by `(gene_name, variant_vcf)`.

This is why gene queries are cheap without storing `gene_name` redundantly in `variants.parquet`.

## Query Patterns Today

This is the current shape of the store, not a wishlist.

### Fast paths

| Query | Today | Why it is fast |
|---|---|---|
| variant by `variant_vcf` | O(1) | direct offset lookup in `sparse_g.bin` |
| all variants in one gene | cheap | `membership.parquet` already stores gene -> variant set |
| gene + MAF + consequence mask | cheap | all metadata is already aligned to `variant_vcf` |
| burden / SKAT inputs for one gene | cheap | sparse carriers plus aligned weights |
| rescore with a different mask | cheap after cache build | score cache slices by `variant_vcf` |

### Works, but not at the ideal shape

| Query | Today | Current limitation |
|---|---|---|
| variant by `vid` | works | builds an in-memory `vid -> variant_vcf` map at load |
| region query, SNVs | works | binary search on start positions, no interval index |
| multi-region query | works | repeated range resolution, not a batched region primitive |
| PRS-style scan over selected variants | works | scan + sum path, not a specialized store primitive |

### Real gaps

| Query | Today | Gap |
|---|---|---|
| sample -> carried variants | poor | no sample-side index |
| trio / family lookups | poor | same root problem as sample queries |
| long indel / SV overlap query | incorrect | no `end_position`, no interval sidecar |
| add samples or variants incrementally | rebuild | one immutable file per chromosome |
| VCF/BCF round-trip export | missing | headers and raw VCF structure are not retained |
| object-store-backed reads | missing | loaders assume local files and mmap |

The query behavior used to live in `queries.md`. It is here now because it only makes sense in the context of the actual store layout.

## Why Sparse Carriers Matter

For rare variants, storing carriers instead of a dense sample-by-variant matrix changes both storage and compute.

Example: 200K samples, 20 variants, MAC=5

| Representation | Stored entries |
|---|---|
| dense matrix | 4,000,000 |
| sparse carriers | about 100 |

That changes the main kernels:

| Operation | Dense | Sparse |
|---|---|---|
| `G'r` | `O(S * m)` | `O(total_MAC)` |
| `G'G` | `O(S * m^2)` | `O(total_MAC * compound_hets)` |
| `G'X` | `O(S * m * k)` | `O(total_MAC * k)` |

This is the main reason COHORT is built around sparse carrier lists in the first place.

## Score Cache

Per-phenotype score statistics are cached under:

```text
.genotype_store/
  score_cache/
    {cache_key}/
      chromosome={chr}/
        scores.bin
```

Three layers:

| Layer | What | Typical cost |
|---|---|---|
| 1 | build store | hours |
| 2 | fit null model, compute U and K | minutes |
| 3 | apply masks, run test battery | seconds |

Current `scores.bin` layout:

- header: magic, version, counts, `sigma2`
- section 1: vid-keyed U values
- section 2: per-gene K blocks

Today the cache works, but it is still a custom parsed binary format. That is why score-cache layout is now tracked as a separate storage issue.

## Build and Load

### Build

```text
VCF + FAVOR annotations
  -> extract genotypes
  -> join with annotations
  -> deduplicate by (position, ref, alt)
  -> write sparse_g.bin
  -> write membership.parquet
  -> write variants.parquet
```

### Load

`VariantIndex` load:

1. read `variants.parquet`
2. read `membership.parquet`
3. build `vid -> variant_vcf`
4. materialize 11 weight vectors

`SparseG` load:

1. mmap `sparse_g.bin`
2. validate header
3. load offset table
4. serve carrier lookups by direct offset jump

Batch carrier loads sort requested `variant_vcf` values first so reads stay mostly sequential.

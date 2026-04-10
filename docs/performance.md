# Performance

Where time and memory go in COHORT today.

> [Back to README](../README.md) · [Genotype Store](storage.md)

## Hot Paths

For a warm-cache rare-variant run, the expensive parts are usually:

| Stage | Why it costs time |
|---|---|
| score-cache miss on large genes | building `G'G` gets expensive as masks grow |
| SKAT eigendecomposition on the largest masks | cubic in mask size |
| first-touch carrier reads on cold data | page faults on initial access |

If a run is slow after the store is already built, it is usually one of the first two.

## Memory

Current memory picture:

| Source | Notes |
|---|---|
| DataFusion operators | already tracked by the DataFusion memory pool |
| stats-kernel scratch | still too allocation-heavy in hot paths |
| score cache | in-memory structures can grow beyond what the pool sees |
| `sparse_g.bin` mmap | intentional, kernel-managed, not pool-managed |
| `variants.parquet` metadata | small per chromosome |

`pool.used()` is not the same thing as RSS today. Mmap-backed storage and some internal allocations still sit outside that number.

## Threads

Current thread story:

| Pool | Today |
|---|---|
| tokio runtime | separate runtime budget |
| rayon | separate global pool |
| per-chromosome work | naturally parallel |

Two independent pools means peak runnable threads can exceed physical cores. That is why `v0.5.0` exists.

## Current Milestones

### `v0.5.0 - memory and thread pool overhaul`

This milestone is about making memory accounting and thread ownership less fuzzy.

Main themes:

- one compute handle instead of multiple competing pools
- less scratch allocation in hot kernels
- bounded score-cache memory
- better machine-visible pool stats

Milestone: [v0.5.0](https://github.com/vineetver/favor-cli/milestone/5)

### `v0.6.0 - storage and query engine`

This is the store-side follow-on:

- sample-side index
- fragment-backed ingest
- interval-aware region queries
- object-store reads
- cleaner store query commands

Milestone: [v0.6.0](https://github.com/vineetver/favor-cli/milestone/6)

## How To Read The Numbers

The docs should describe current behavior, not ideal behavior.

When we quote timings or memory in the future, they should come from the bench harness and be tied to:

- fixture
- sample count
- variant count
- memory budget
- node type

If those details are missing, treat the number as anecdotal.

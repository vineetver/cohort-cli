# aGDS + STAARpipeline vs Our Custom Genotype Store

Simple version:

- original STAARpipeline keeps going back to the source data
- our implementation builds once and reuses it

| Original aGDS + STAARpipeline | Our custom implementation | Why ours is better | Practical impact |
|---|---|---|---|
| For each gene or mask, it filters the source data again, reads the subset again, and rebuilds the test inputs again. | Builds the store once with genotype, metadata, and gene membership already lined up. | Stops treating every test like a fresh extraction job. | Less repeated work. |
| For each selected subset, it reads dosage into a full in-memory matrix even though most entries are empty for rare variants. | Stores only samples that actually have each rare variant. | Avoids storing and rereading mostly empty data. | `200K x 20` can shrink from `4,000,000` values to about `100` useful entries. |
| When related masks share many of the same variants, it still redoes much of the same preparation work for each mask instead of reusing it. | Reuses one aligned variant index plus cached per-gene score pieces. | Shared work is reused. | Faster repeated testing. |
| Genotype reads, annotation reads, and mask-building all happen together in the same repeated per-test path. | Keeps genotype, metadata, and gene membership separate but aligned. | Can reuse metadata and sometimes avoid rereading genotype. | Less decoding and less data movement. |
| As it moves from one test subset to the next, it repeatedly rebuilds large in-memory genotype tables. | Keeps only useful genotype entries. | Avoids filling memory with mostly empty data. | Lower memory use. |
| The workflow is correct, but it is organized around repeatedly pulling subsets out of the source data rather than around repeated rare-variant testing. | Built around repeated gene and mask testing from the start. | Storage matches the workload better. | Better scaling. |

Original flow:

- filter again
- read genotype again
- rebuild inputs again
- run the test
- repeat

Our flow:

- build the store once
- store only the samples that actually have each rare variant
- keep metadata and gene membership aligned
- reuse cached score pieces

Main message:

**original STAARpipeline = repeatedly rereading and reconstructing test inputs**

**our implementation = build once, reuse many times**

We are still improving the architecture. GitHub currently has no issues assigned to milestone `v0.6.0`, but the repo already lists these next steps:

- sample-side index
- fragment-backed ingest
- interval-aware region queries
- object-store reads
- cleaner store query commands

# FAVOR CLI

## Session Voice

Speak like caveman in session. Short words. Direct words. No polished helper voice. No corporate tone. No AI theater.

Code and docs stay normal English. Only session chatter is caveman.

## Coding Principles

Use these as mental models, not cosplay. "House" means the lens to apply when choosing a design. "Commandment" means the concrete rule to follow in code review and edits.

### Houses

- House of Structure: shape data so edge cases disappear. Good: one canonical variant/list/store type with typed fields. Bad: parallel structs and stringly-typed flags.
- House of Latency: respect hardware and IO cost. Favor streaming, mmap, sort-merge joins, bounded buffers, and reuse over extra allocation or pointer-chasing.
- House of Simplicity: keep pure logic separate from orchestration and file/system effects. Parsing, scoring math, and mask logic should be testable without plumbing.
- House of Resurrection: failures must be explicit, contained, and recoverable. Return typed errors, write atomically, probe before rebuild, never hide corruption.

### Commandments

- Reduce code first. Default move is delete, merge, or reuse before adding.
  Bad: new helper, new trait, new module, same behavior
  Good: remove duplication, extend existing path, ship fewer lines
- Every line must earn its keep.
- No AI slop comments, banner comments, TODO litter, or commented-out code.
- Comments only for invariants, statistical correctness, file format details, or non-obvious why.
- Parse, do not vaguely validate. Turn raw input into a type that carries the proof.
  Bad: `fn send(email: &str) -> bool`
  Good: `fn send(email: Email) -> Result<()>`
- Make invalid states unrepresentable.
  Bad: `{ is_loading: bool, is_success: bool, is_error: bool }`
  Good: `enum RequestState { Loading, Success(Data), Error(String) }`
- Keep a single source of truth.
  Bad: store both `items: Vec<Item>` and `item_count: usize`
  Good: store `items` and derive `items.len()`
- Flatten control flow.
  Bad: four nested `if` blocks before the real work
  Good: early returns for bad cases, then one straight happy path
- Separate query from command when practical.
  Bad: `get_user()` also updates `last_access`
  Good: `get_user()` reads, `log_access()` writes
- Prefer composition over inheritance-shaped lies.
  Bad: `trait Bird { fn fly(&self); }` and `Penguin` panics in `fly()`
  Good: separate `Flyable` and `Swimmable`
- Inject dependencies, do not reach for them. Pass `&Engine`, `&Store`, `&dyn Output`, `&Resources` into the function that needs them. No globals, no `lazy_static`, no `thread_local!` registries, no constructors that read env vars or open files at module load. The argument list is the dependency contract — if it's hard to read, the function is doing too much.
  Bad:
  ```rust
  // hidden dependency on a global, untestable without env mutation
  static STORE: OnceLock<Store> = OnceLock::new();
  pub fn load_cohort(id: &CohortId) -> Result<CohortHandle<'static>> {
      STORE.get_or_init(|| Store::open(StoreConfig::resolve(None).unwrap()).unwrap())
           .cohort(id)
  }
  ```
  Good:
  ```rust
  // dependencies are explicit; tests pass a tempdir Store, prod passes the real one
  pub fn load_cohort<'a>(store: &'a Store, id: &CohortId) -> CohortHandle<'a> {
      store.cohort(id)
  }
  ```
  Why this matters here: `Engine` is the one composition root (`src/runtime.rs`). Every command gets `engine: &Engine` and pulls what it needs (`engine.store()`, `engine.df()`, `engine.resources()`). Subsystems take borrowed handles, never construct them. That is what makes `stage_ensure_store` swappable between `Fresh` and `Existing`, what lets the parquet round-trip test run without a tempdir Store, and what keeps `score_cache` keyed by `(cohort_id, key)` instead of a process-global path.
- Do not allocate in the hot path.
  Bad: `let mut buf = Vec::new();` inside every loop iteration
  Good: allocate once, `clear()` and reuse
- Use existing hooks before adding new modules, traits, or layers.
  Bad: build a plugin system for one config file
  Good: write the five-line loader that solves today's job
- Output goes through the output layer. Do not bare `println!`.
- Errors must say what failed and what user should do next.
- Refactors do not get to change p-values.
- Three explicit lines beat one premature abstraction.

## Current CLI Stage

FAVOR CLI is in the middle of the storage-engine migration.

- `src/store/` is real and active.
- Store root resolution exists now: CLI arg, then `FAVOR_STORE`, then walk up for `.cohort`, then fallback to `<cwd>/.cohort`.
- Core store pieces landed: backend, layout, manifest, scratch pool, cohort storage, list storage, annotation module, cache module.
- STAAR cohort storage is the current canonical sparse layout:

```text
.cohort/
  cohorts/<id>/...     # target design

legacy staar store today still materializes as:
<output>/store/
  manifest.json
  samples.txt
  chromosome={chr}/
    sparse_g.bin
    variants.parquet
    membership.parquet
```

- `src/store/cohort/` holds the sparse genotype store code.
- `src/store/list/` now owns variant-set containers used by ingest and annotate paths.
- `src/store/annotation/` owns annotation and tissue readers.
- `src/store/cache/score_cache.rs` exists, but parts of STAAR still think in legacy `store_dir/score_cache/...` terms.
- The migration is not finished. Some commands use new store modules; some call sites still carry old path/layout assumptions.

## What The CLI Does Right Now

- `ingest` builds chromosome-partitioned variant sets.
- `annotate` joins ingested variants against FAVOR base or full annotation tiers.
- `enrich` overlays tissue/regulatory tables.
- `staar` runs rare-variant association on annotated sets and uses the sparse cohort store.
- `meta-staar` exists for cross-study meta analysis.
- `inspect` reads artifacts and schemas.
- `interpret` is still a stub.

## Working Rules

- Machine mode is first-class. Prefer `--format json`.
- Use `--dry-run` before expensive resource creation or rebuilds.
- Respect HPC reality: do not assume login-node builds or unlimited `/tmp`.
- When editing storage code, match the migration direction in `docs/storage-engine-design.md`: one store handle, one resolver, subsystem ownership, fewer ad hoc paths.

# Contributing

## Build

```bash
cargo check                                      # iterate fast (no DuckDB linking)
srun -p test -c 8 --mem=16G cargo build --release # full build on HPC
```

First build takes ~10 min (DuckDB). `cargo check` takes seconds.

## Rules

- Zero warnings. Always.
- Every error message tells the user what command fixes it.
- No `println!`. All output goes through the `Output` trait.
- No hardcoded sizes. Derive from `resources.memory_bytes`.
- If it doesn't fit in memory, batch it. Never crash.
- Types over strings. Enums over magic values.
- Commands are thin. Logic lives in modules.
- Early returns. Flat control flow. No deep nesting.

## Add a command

1. Create `src/commands/<name>.rs`
2. Add variant to `Command` in `src/cli.rs`
3. Route in `src/main.rs`
4. Support `--dry-run` and `--format json`
5. Use `Resources::detect_with_config()` for memory (except setup, which uses `detect()`)

## Add a data pack

One entry in `PACKS` in `src/packs.rs`. Everything else picks it up.

## Test

```bash
cargo test
cargo check  # must be zero warnings
```

# z-fastq

Fast Zig FASTQ I/O library and CLI (work in progress).

## Status

Planning phase. See `plan/v0.0.1-plan.md` for the active milestone (plain parser + count). Public **v0.1.0** is a later target.

## Build (once implemented)

```bash
# Requires zig-0.16.0 at ./zig-0.16.0/zig and ./zig wrapper
./zig build
./zig build test --summary all
./zig build -Doptimize=ReleaseFast
```

## Data corpus (local only)

`data/` is gitignored. Fetch fixtures on each machine:

```bash
bash scripts/download_corpus.sh --benchmark
bash scripts/build_demo_extracts.sh
```

See `data/README.md` (local file after first download).

## License

TBD

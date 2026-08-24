# Benchmark tools

Status: **Active** (last updated: 2026-08-24)

`tools/` prepares the external commands used to check and compare z-fastq. It never installs system packages or changes your shell. Versions live in [`versions.sh`](versions.sh); build and runtime checks live in [`install.sh`](install.sh).

## Use

```bash
tools/install.sh <name>          # install or reuse one tool
tools/install.sh peers           # prepare every enabled peer
tools/install.sh all             # prepare the Python venv, adapters, and peers
tools/install.sh --check all     # check without downloading or building
tools/install.sh --list          # show local state
```

Benchmarks call commands through `tools/bin/`. These are relative links to the checked current installs under `tools/.local/`. A new runtime is tested before the links switch, so a failed update does not replace a working command.

Downloads and builds use `/tmp` and are removed after each tool. Only the runtime and a small receipt remain. `tools/bin/`, `tools/.local/`, and `tools/venv/` are generated and ignored. The root `.venv` is not used. Build tools come from `PATH` only when a build is required; prepared commands do not require their old build toolchains.

Every retained native executable is stripped, and `--check` verifies this. BBTools contains Java class files and shell scripts, so native executable stripping does not apply to it. Fasten builds only its sampler from the published crate and included Cargo lock. SeqFu builds with a temporary pinned Nim compiler and locked packages. Build files are discarded afterward.

## Current tools

Sizes are the rounded retained sizes from the current Linux x86-64 installation. They exclude temporary source and build trees.

| Tool                             |              Version |    Size | Use here                                                      |
| -------------------------------- | -------------------: | ------: | ------------------------------------------------------------- |
| seqtk                            |             1.5-r133 |  84 KiB | Count reference, compatible sampling, and interleave delivery |
| FastQValidator                   |               0.1.1a | 176 KiB | Historical validation behavior only                           |
| Fasten                           |                0.9.0 | 468 KiB | Plain interleaved probability-sampling peer                   |
| Needletail and Helicase adapters |      0.7.3 and 0.2.0 | 840 KiB | One-worker parser engine checks for count and aggregate stats |
| fqtools                          | 2.3 with HTSlib 1.24 | 892 KiB | Independent count and qualified validation cases              |
| SeqFu                            |               1.27.1 | 1.7 MiB | Count, stats, check, interleave, and deinterleave             |
| IRMA Core                        |               0.10.1 | 2.4 MiB | Exact sampling, interleave, and deinterleave                  |
| fq                               |               0.12.0 | 2.6 MiB | Validation and paired probability-sampling reference          |
| fastp                            |                1.3.6 | 2.9 MiB | Broader short-read QC and preprocessing context               |
| Rasusa                           |                5.1.0 | 5.4 MiB | Paired probability and broader sampling peer                  |
| fqkit                            |               0.4.14 | 5.8 MiB | Interleave, deinterleave, stats, sampling, and sharding       |
| SeqKit                           |               2.13.0 |  19 MiB | Broad stats, sampling, pairing, and conversion peer           |
| BBTools                          |                40.02 |  22 MiB | Interleaved sampling reference and read-pair repair           |

The complete generated tool area is about 64 MiB, plus a 96 KiB tools-only Python environment.

Fqtools needs [`patches/fqtools-gzfile-casts.patch`](patches/fqtools-gzfile-casts.patch) because its current source passes zlib `gzFile` handles with the wrong pointer depth. The patch changes fqtools, not zlib. Keeping this small patch visible is safer than rewriting downloaded source silently.

Installation checks prove that the retained commands start and complete representative work. The benchmark suites must still prove semantic parity for each measured job before reporting performance.

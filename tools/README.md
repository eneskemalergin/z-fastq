# Development tools

This directory contains optional maintainer tooling for checking external FASTQ commands and running local comparisons. End users do not need it to build or run z-fastq.

## Use

```bash
tools/install.sh <name>
tools/install.sh peers
tools/install.sh all
tools/install.sh --check all
tools/install.sh --list
```

The scripts prepare the tools needed by local comparison workflows. They do not change the shell or install system packages.

## Current tools

The following sizes are rounded retained executable sizes from the current Linux x86-64 installation. They show the scale of the comparison set; they exclude temporary source and build trees.

| Tool                             |              Version |    Size | Use here                                                      |
| -------------------------------- | -------------------: | ------: | ------------------------------------------------------------- |
| `seqtk`                          |             1.5-r133 |  84 KiB | Count reference, compatible sampling, and interleave delivery |
| `FastQValidator`                 |               0.1.1a | 176 KiB | Historical validation behavior                                |
| `Fasten`                         |                0.9.0 | 468 KiB | Plain interleaved probability sampling                        |
| Needletail and Helicase adapters |      0.7.3 and 0.2.0 | 840 KiB | One-worker parser checks for count and aggregate stats        |
| `fqtools`                        | 2.3 with HTSlib 1.24 | 892 KiB | Independent count and validation cases                        |
| `SeqFu`                          |               1.27.1 | 1.7 MiB | Count, stats, check, interleave, and deinterleave             |
| `IRMA Core`                      |               0.10.1 | 2.4 MiB | Exact sampling, interleave, and deinterleave                  |
| `fq`                             |               0.12.0 | 2.6 MiB | Validation and paired sampling reference                      |
| `fastp`                          |                1.3.6 | 2.9 MiB | Broader short-read QC and preprocessing context               |
| `Rasusa`                         |                5.1.0 | 5.4 MiB | Paired sampling and broader sampling peer                     |
| `fqkit`                          |               0.4.14 | 5.8 MiB | Interleave, deinterleave, stats, sampling, and sharding       |
| `SeqKit`                         |               2.13.0 |  19 MiB | Broad stats, sampling, pairing, and conversion peer           |
| `BBTools`                        |                40.02 |  23 MiB | Sampling, pair validation, interleave, and deinterleave       |

The table is descriptive local comparison context, not a universal ranking.

## Tracked contents

- `install.sh` prepares and checks external command installations.
- `versions.sh` owns the selected tool versions.
- `patches/` contains small source fixes required by pinned external recipes.
- `wrappers/` contains adapters for tools that need a common comparison interface.

Generated downloads, builds, environments, and command links are local state. They are not source files, release artifacts, or part of the z-fastq CLI contract.

## Scope of the comparisons

The tools support correctness checks and local development comparisons for overlapping FASTQ operations. A comparison must establish compatible output and failure behavior before timing anything.

Local measurements are development evidence. They are not a permanent public ranking, and this document does not define z-fastq's user-facing performance claims.

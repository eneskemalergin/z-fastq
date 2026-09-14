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

The rows below list the pinned comparison tools, the retained Linux x86-64 size, and their roles. The sizes were measured from the current local installation on 2026-09-14. Compiled-tool values are exact byte lengths of the stripped executable linked from `tools/bin/`. The BBTools value is the sum of regular files in its retained runtime tree because BBTools is a Java class tree rather than one executable. Versions are defined by `tools/versions.sh`.

| Tool                 |              Version | Size (bytes) | Use here                                                      |
| -------------------- | -------------------: | -----------: | ------------------------------------------------------------- |
| `seqtk`              |             1.5-r133 |       77,600 | Count reference, compatible sampling, and interleave delivery |
| `FastQValidator`     |               0.1.1a |      169,120 | Descriptive validation peer                                   |
| `Fasten`             |                0.9.0 |      470,336 | Plain interleaved probability sampling                        |
| `Needletail adapter` |                0.7.3 |      425,896 | One-worker parser checks for count and aggregate stats        |
| `Helicase adapter`   |                0.2.0 |      424,080 | One-worker parser checks for count and aggregate stats        |
| `fqtools`            | 2.3 with HTSlib 1.24 |      905,088 | Independent count and validation cases                        |
| `SeqFu`              |               1.27.1 |    1,673,280 | Count, stats, check, interleave, and deinterleave             |
| `IRMA Core`          |               0.10.1 |    2,485,768 | Exact sampling, interleave, and deinterleave                  |
| `fq`                 |               0.12.0 |    2,740,832 | Validation and paired sampling reference                      |
| `fastp`              |                1.3.7 |    3,007,424 | Short-read QC and preprocessing peer                          |
| `Rasusa`             |                5.1.0 |    5,640,352 | Paired sampling and broader sampling peer                     |
| `fqkit`              |               0.4.14 |    6,007,128 | Interleave, deinterleave, stats, sampling, and sharding       |
| `SeqKit`             |               2.13.0 |   20,076,696 | Broad stats, sampling, pairing, and conversion peer           |
| `BBTools`            |                40.02 |   19,359,140 | Sampling, pair validation, interleave, and deinterleave       |

The table describes local comparison inputs; it does not rank the tools.

## Tracked contents

- `install.sh` prepares and checks external command installations.
- `versions.sh` owns the selected tool versions.
- `patches/` contains small source fixes required by pinned external recipes.
- `wrappers/` contains adapters for tools that need a common comparison interface.

## Scope of the comparisons

The tools support correctness checks and local development comparisons for overlapping FASTQ operations. A comparison must establish compatible output and failure behavior before timing anything.

Local measurements are development evidence. They are not a permanent public ranking, and this document does not define z-fastq's user-facing performance claims.

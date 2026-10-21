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

The table lists the pinned comparison tools and the retained Linux x86-64 size of each pinned version. Compiled-tool values are exact byte lengths of the stripped executable linked from `tools/bin/`. The BBTools value is the sum of regular files in its retained runtime tree because BBTools is a Java class tree rather than one executable. Versions are defined by `tools/versions.sh`.

| Tool                 |              Version | Size (bytes) |
| -------------------- | -------------------: | -----------: |
| `seqtk`              |             1.5-r133 |       77,552 |
| `FastQValidator`     |               0.1.1a |      169,072 |
| `Fasten`             |                0.9.0 |      470,288 |
| `Needletail adapter` |                0.7.3 |      425,912 |
| `Helicase adapter`   |                0.2.0 |      428,944 |
| `fqtools`            | 2.3 with HTSlib 1.24 |      905,048 |
| `SeqFu`              |               1.29.0 |    2,227,064 |
| `IRMA Core`          |               0.11.0 |    2,601,752 |
| `fq`                 |               0.12.0 |    2,740,832 |
| `fastp`              |                1.4.0 |    3,007,424 |
| `Rasusa`             |                5.1.0 |    5,640,352 |
| `fqkit`              |               0.4.14 |    6,007,912 |
| `SeqKit`             |               2.14.0 |   20,473,816 |
| `BBTools`            |                40.05 |   27,434,020 |

What each tool is used for here:

- `seqtk`: count reference, compatible sampling, and interleave delivery.
- `FastQValidator`: descriptive validation peer.
- `Fasten`: plain interleaved probability sampling.
- Needletail and Helicase adapters: one-worker parser checks for count and aggregate stats.
- `fqtools`: independent count and validation cases.
- `SeqFu`: count, stats, check, interleave, and deinterleave.
- `IRMA Core`: exact sampling, interleave, and deinterleave.
- `fq`: validation and paired sampling reference.
- `fastp`: short-read QC and preprocessing peer.
- `Rasusa`: paired sampling and broader sampling peer.
- `fqkit`: interleave, deinterleave, stats, sampling, and sharding.
- `SeqKit`: broad stats, sampling, pairing, and conversion peer.
- `BBTools`: sampling, pair validation, interleave, and deinterleave.

This list describes local comparison inputs; it does not rank the tools.

## Report Python and zebrac

`tools/venv` is also the report environment. `tools/install.sh` (or `tools/install.sh venv`) installs pinned matplotlib 3.11.2, pandas 3.0.6, and tabulate 0.10.0 into that venv. Bench scripts call `tools/venv/bin/python` and do not use a second virtualenv.

The benchmark runner is [zebrac](https://github.com/eneskemalergin/zebrac). Keep the Linux binary at `tools/zebrac`. `tools/install.sh --check zebrac` and `tools/install.sh --check all` verify that path; they do not download or build zebrac.

```bash
tools/install.sh venv
tools/install.sh --check venv
tools/install.sh --check zebrac
```

## Tracked contents

- `install.sh` prepares and checks external command installations, the report Python packages in `tools/venv`, and the presence of `tools/zebrac`.
- `versions.sh` pins the selected tool versions.
- `patches/` contains small source fixes required by pinned external recipes.
- `wrappers/` contains adapters for tools that need a common comparison interface.

## What the comparisons cover

The tools support correctness checks and local development comparisons for overlapping FASTQ operations. A comparison must establish compatible output and failure behavior before timing anything.

Local measurements are development evidence. They are not a permanent public ranking, and this document does not define z-fastq's user-facing performance claims.

<!-- markdownlint-disable MD013 MD033 MD041 -->

# Benchmarks

How z-fastq compares with the tools people use for the same FASTQ jobs, at a glance. Each feature report below holds the detail: every input, format, and size, the edge cases, scaling, and the commands.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="summary/01-jobs-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="summary/01-jobs-light.png">
    <img src="summary/01-jobs-light.png" alt="Each job: the fastest complete and partial tools as multiples of z-fastq's time and memory" width="100%">
  </picture>
</p>

| Job                      | z-fastq time | z-fastq memory | Complete tool | Its time | Its memory | Partial tool | Its time | Its memory |
| ------------------------ | -----------: | -------------: | ------------- | -------: | ---------: | ------------ | -------: | ---------: |
| Count records            |       147 ms |       0.81 MiB | Needletail    |    1.37x |      2.57x | Helicase     |    1.22x |      3.19x |
| Aggregate statistics     |       280 ms |       0.81 MiB | Needletail    |    3.54x |      2.58x | SeqFu        |    1.70x |      6.41x |
| Validate a file          |       216 ms |       0.79 MiB | fq            |    4.08x |      5.29x | fqtools      |   16.14x |      4.11x |
| Validate a pair          |       435 ms |       1.06 MiB | fq            |    6.74x |     25.45x | -            |        - |          - |
| Sample by fraction       |       237 ms |       0.98 MiB | IRMA Core     |    5.05x |      4.79x | Rasusa       |    0.99x |      4.61x |
| Sample by count          |       421 ms |       1.34 MiB | IRMA Core     |    2.78x |      3.50x | Rasusa       |    1.00x |      7.78x |
| Sample pairs by fraction |       438 ms |       1.29 MiB | IRMA Core     |    3.71x |      3.35x | Rasusa       |    0.77x |      3.57x |
| Interleave a pair        |       471 ms |       1.06 MiB | IRMA Core     |    3.06x |      3.81x | fqkit        |    1.12x |     13.10x |
| Deinterleave a pair      |       782 ms |       1.06 MiB | IRMA Core     |    4.80x |      3.84x | fqkit        |    0.94x |     12.31x |
| Trim a pair              |       3.41 s |      14.68 MiB | -             |        - |          - | fastp        |   10.97x |     75.88x |
| QC tables for a pair     |       1.91 s |      10.17 MiB | -             |        - |          - | FastQC-Rust  |    2.13x |      2.83x |

- A *complete* tool does the whole job z-fastq does, checks included; a *partial* one does part of it (fewer checks, or no reads written). Partial tools are shown, never ranked with complete ones.
- Times and memory: z-fastq's own, then each tool's as a multiple of z-fastq's; geometric means over the 1 GiB plain inputs, every tool on one CPU.

## How the numbers are taken

Every tool's output is checked against an independent truth before any timing; a tool is timed only on the inputs it got right. Each case runs every tool in one batch, on one pinned CPU each (threaded variants on as many physical cores as threads), with the input in the page cache. Wall time comes from the runner, peak memory and CPU time from the kernel for the tool alone. Peak RSS can lag the true peak by up to about 0.5 MiB, so smaller memory differences are inconclusive.

## Broken input

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="summary/02-broken-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="summary/02-broken-light.png">
    <img src="summary/02-broken-light.png" alt="Of the malformed inputs each tool was given, how many it rejected" width="100%">
  </picture>
</p>

| Tool           | Given | Rejected | Passed on, exit 0 | Crashed |
| -------------- | ----: | -------: | ----------------: | ------: |
| z-fastq        |    54 |       54 |                 0 |       0 |
| IRMA Core      |    21 |       21 |                 0 |       0 |
| Needletail     |    16 |       16 |                 0 |       0 |
| fqtools        |    28 |       22 |                 6 |       0 |
| fastQValidator |     9 |        7 |                 2 |       0 |
| Rasusa         |     9 |        7 |                 2 |       0 |
| fq             |    18 |       13 |                 5 |       0 |
| FastQC         |     2 |        1 |                 1 |       0 |
| FastQC-Rust    |     2 |        1 |                 1 |       0 |
| Helicase       |    16 |        8 |                 4 |       4 |
| SeqFu          |    28 |       11 |                17 |       0 |
| SeqKit         |    18 |        7 |                11 |       0 |
| Fasten         |     9 |        2 |                 7 |       0 |
| fqkit          |    21 |        2 |                19 |       0 |
| fastp          |     6 |        0 |                 6 |       0 |
| seqtk          |    31 |        0 |                31 |       0 |

A tool that passes a broken record on with exit 0 leaves a script no way to tell; z-fastq stops with the record, byte offset, and rule it broke.

## One core against many

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="summary/03-threads-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="summary/03-threads-light.png">
    <img src="summary/03-threads-light.png" alt="Wall time of the other tools by thread count against z-fastq on one core" width="100%">
  </picture>
</p>

| Job                  | Tool        | 1 thread |       2 |      4 |      8 |     16 |
| -------------------- | ----------- | -------: | ------: | -----: | -----: | -----: |
| Trim a pair          | z-fastq     |   3.41 s |       - |      - |      - |      - |
| Trim a pair          | fastp       |  37.40 s | 17.51 s | 9.14 s | 5.06 s | 3.55 s |
| QC tables for a pair | z-fastq     |   1.91 s |       - |      - |      - |      - |
| QC tables for a pair | FastQC-Rust |   4.07 s |  2.23 s | 1.86 s | 1.69 s | 2.32 s |
| QC tables for a pair | FastQC      |  22.10 s | 17.64 s | 8.82 s | 5.48 s | 5.46 s |

Wall time on the 1 GiB plain inputs; z-fastq runs on one core. The process report gives CPU time at each count.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="summary/04-threads-memory-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="summary/04-threads-memory-light.png">
    <img src="summary/04-threads-memory-light.png" alt="The threaded runs on time against memory, relative to z-fastq on one core" width="100%">
  </picture>
</p>

| Job                  | Tool        |   1 thread |          2 |          4 |          8 |         16 |
| -------------------- | ----------- | ---------: | ---------: | ---------: | ---------: | ---------: |
| Trim a pair          | z-fastq     |   14.7 MiB |          - |          - |          - |          - |
| Trim a pair          | fastp       | 1113.6 MiB | 1150.6 MiB | 1186.1 MiB | 1249.5 MiB | 1355.2 MiB |
| QC tables for a pair | z-fastq     |   10.2 MiB |          - |          - |          - |          - |
| QC tables for a pair | FastQC-Rust |   28.8 MiB |   44.3 MiB |   63.5 MiB |   63.2 MiB |   60.7 MiB |
| QC tables for a pair | FastQC      |  251.0 MiB |  347.0 MiB |  379.2 MiB |  702.3 MiB |  717.3 MiB |

Peak memory of the same runs. Threads buy the other tools wall time, never memory: at the thread count where a tool first matches z-fastq's time, it still holds many times z-fastq's memory.

## Where z-fastq is not fastest

- **Sample by fraction:** Rasusa (partial) takes 0.99x z-fastq's time, with 4.6x its memory.
- **Sample pairs by fraction:** Rasusa (partial) takes 0.77x z-fastq's time, with 3.6x its memory.
- **Deinterleave a pair:** fqkit (partial) takes 0.94x z-fastq's time, with 12.3x its memory.
- **QC tables for a pair:** FastQC-Rust at 4 threads takes 0.97x z-fastq's time on one core, for 2.4x its CPU time.

## Feature reports

- [count](count/README.md): Counting records; timed 2026-10-28, z-fastq 0.1.0, source at 246bbb1.
- [stats](stats/README.md): Aggregate statistics; timed 2026-10-28, z-fastq 0.1.0, source at 246bbb1.
- [check](check/README.md): Validating FASTQ files; timed 2026-10-28, z-fastq 0.1.0, source at 246bbb1.
- [sample](sample/README.md): Sampling records; timed 2026-10-28, z-fastq 0.1.0, source at 246bbb1.
- [layout](layout/README.md): Converting paired layouts; timed 2026-10-28, z-fastq 0.1.0, source at 246bbb1.
- [process](process/README.md): Processing paired reads; timed 2026-10-29, z-fastq 0.1.0, source at 246bbb1.

## Limits

- One CPU per tool in every comparison except the threads figure; a tool built for many cores may do better with them.
- The real-run groups are uneven: three short fixed-length runs, one variable-length run, one synthetic shape per job.
- Copying a pair (the process report's first part) is the cost floor of `process`, not a job of its own, so it is not in the job figure.

## Reproduce

```bash
bench/process/run.sh release    # one feature: checks, timing, its README (each feature has its own run.sh)
bench/process/run.sh summary    # this page, from every feature's report tables
```

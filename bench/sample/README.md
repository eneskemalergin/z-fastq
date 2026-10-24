<!-- markdownlint-disable MD013 MD033 MD041 -->

# Sampling records
> <sub>`z-fastq sample` against 7 other samplers, on 7 real sequencing runs and 5 synthetic shapes, plain and gzip. A tool is timed only where its sample was right; time and memory come from the same run, on one CPU.</sub>

- `z-fastq`: *0.1.0 at 90e46b4*
- `machine`: *Ryzen 9 3950X, Linux 7.2.11; one CPU per run*
- `runner`: *zebrac 0.6.2 with zexec*
- `data`:
  - **Real**: *7 public sequencing runs, 1 GiB each, pooled in the bars by the kind of reads they hold*
    - `Short, fixed length`: *NovaSeq, 150 bp; NovaSeq, Casava headers; BGISEQ, 100 bp; miRNA, 36 bp*
    - `Short, variable length`: *10x mates, 28 and 150 bp*
    - `Long reads`: *PromethION, 30 to 153 kb; Revio HiFi, 4 to 21 kb*
  - **Synthetic**: *5 shapes generated from real reads, each changing one thing; one group in the bars at 1 GiB, and the scaling lines from 1 MiB to 1 GiB*
    - `dense`: *uniform 150 bp, bare plus line*
    - `endings`: *100 bp, CRLF and LF lines*
    - `variable`: *50 to 150 bp, annotated plus line*
    - `longread`: *nanopore headers, rare 300 kb reads*
    - `long`: *300 kb sequence lines*
- `date`: *2026-10-24*

## Findings

- **Right answers:** z-fastq returned the right sample on all 84 checked inputs and rejected every malformed edge case. 6 of the 7 other tools print a wrong sample or crash on at least one edge case.
- **Fair samples:** every tool passes the fairness tests except IRMA Core (fraction).
- **Matches seqtk:** z-fastq kept the same records as seqtk on 45 of 45 screened inputs (fraction 36 of 36, exact count 9 of 9).
- **Fraction speed:** the fastest other tool takes 0.69x z-fastq's time on plain input (Rasusa) and 1.01x on gzip (Rasusa); IRMA Core, the fastest other full-check tool, takes 3.28x z-fastq's time on plain input and 2.94x on gzip.
- **Exact-count speed:** the fastest other tool takes 0.90x z-fastq's time on plain input (Rasusa); IRMA Core, the fastest other full-check tool, takes 2.36x z-fastq's time on plain input.
- **Memory against K:** at K = 100,000, z-fastq peaks at 1.27 MiB and the other tools at 3.4x to 15.3x that; at K = 1,000,000, z-fastq peaks at 3.77 MiB and the other tools at 1.2x to 9.5x that.
- **Where z-fastq is not fastest, fraction:** Rasusa on plain synthetic (0.64x); seqtk on plain synthetic (0.90x); Rasusa on plain short, variable length (0.68x); Rasusa on gzip short, variable length (0.89x); Rasusa on plain long reads (0.76x); seqtk on plain long reads (0.89x); Rasusa on plain short, fixed length (0.73x); Rasusa on gzip short, fixed length (0.98x).
- **Where z-fastq is not fastest, exact count:** Rasusa on plain synthetic (0.86x); Rasusa on plain long reads (0.76x); seqtk on plain long reads (0.95x).

## Compared tools
> <sub>*complete* tools check records as much as z-fastq; *partial* tools check less, so they do less work per record. Modes: `yes` when the tool was right on every input of the mode, `partly` on some, `-` when the tool has no such mode. Stdin and gzip: measured on a 1 MiB file with the timed command (a two-pass exact count cannot read stdin). Seed: the same output twice with one seed, and different output with another.</sub>

| Tool      | Class    | Fraction | Exact count | Pairs, two files | Pairs, interleaved | Stdin | Gzip | Seed |
| --------- | -------- | -------- | ----------- | ---------------- | ------------------ | ----- | ---- | ---- |
| z-fastq   | subject  | yes      | yes         | yes              | yes                | yes   | yes  | yes  |
| Rasusa    | partial  | yes      | yes         | yes              | -                  | no    | yes  | yes  |
| seqtk     | partial  | partly   | partly      | -                | -                  | yes   | yes  | yes  |
| fq        | partial  | yes      | yes         | yes              | -                  | no    | yes  | yes  |
| fqkit     | partial  | -        | yes         | -                | -                  | no    | yes  | yes  |
| SeqKit    | partial  | yes      | yes         | -                | -                  | yes   | yes  | yes  |
| IRMA Core | complete | yes      | yes         | partly           | -                  | no    | yes  | yes  |
| Fasten    | partial  | -        | -           | -                | yes                | yes   | no   | no   |

### What each tool checks and changes
> <sub>Each kind of malformed input from the edge cases: `rejects`, or **passes** when the tool writes a sample from it with exit 0, so it does less work than z-fastq. Plus line and line endings written: on the `variable` (annotated plus line) and `endings` (CRLF and LF) shapes; `as input` and `kept` copy the input. Unfair in: the modes where its sample fails the fairness tests below.</sub>

| Tool      | Length mismatch | Truncated record | Bad plus line | Bad header | Quality range | Blank line | Broken gzip | Plus line written | Line endings written | Unfair in |
| --------- | --------------- | ---------------- | ------------- | ---------- | ------------- | ---------- | ----------- | ----------------- | -------------------- | --------- |
| z-fastq   | rejects         | rejects          | rejects       | rejects    | rejects       | rejects    | rejects     | kept              | LF                   | -         |
| Rasusa    | rejects         | rejects          | rejects       | rejects    | **passes**    | rejects    | rejects     | bare +            | CRLF                 | -         |
| seqtk     | **passes**      | **passes**       | **passes**    | **passes** | **passes**    | **passes** | **passes**  | bare +            | LF                   | -         |
| fq        | **passes**      | rejects          | **passes**    | **passes** | **passes**    | rejects    | rejects     | kept              | LF                   | -         |
| fqkit     | **passes**      | **passes**       | **passes**    | **passes** | **passes**    | **passes** | rejects     | bare +            | mixed                | -         |
| SeqKit    | rejects         | rejects          | rejects       | rejects    | **passes**    | **passes** | rejects     | bare +            | LF                   | -         |
| IRMA Core | rejects         | rejects          | rejects       | rejects    | rejects       | rejects    | rejects     | bare +            | LF                   | Fraction  |
| Fasten    | **passes**      | **passes**       | **passes**    | **passes** | **passes**    | **passes** | rejects     | -                 | -                    | -         |

- `Rasusa 5.1.0`: *-1 is its one-pass fraction; without it, two passes and exactly round(pN)*
- `seqtk 1.5`: *exact counts need -2 for input order*
- `fq 0.12.0`: *writes only to named files*
- `fqkit 0.4.14`: *exact count only; 4 threads unless given -@ 1*
- `SeqKit 2.14.0`: *exact counts need sample2 -2 -w 0; 4 threads unless given -j 1*
- `IRMA Core 0.11.0`: *whole percents only; an exact count (-t) on gzip input comes out of input order, so exact count is timed on plain input*
- `Fasten 0.9.0`: *stdin only, plain only, no seed*

> [!IMPORTANT]
> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.

## Matches seqtk
> <sub>z-fastq promises seqtk's selection on inputs with LF endings, nonempty reads, bare plus lines, printable headers, a seed that fits signed 64-bit, and a value seqtk reads as a fraction. Equivalent commands: `z-fastq sample --fraction P --seed S IN` and `seqtk sample -s S IN P`; `z-fastq sample --count K --seed S IN` and `seqtk sample -2 -s S IN K`.</sub>

- `fraction`: *the same records on 36 of 36 inputs*
- `exact count`: *the same records on 9 of 9 inputs*

> [!NOTE]
> seqtk without `-2` keeps the same records but writes them in reservoir order, reads `1.0` as a count of 1, and reads `0` as a fraction.

## Right answers
> <sub>A sample is right when every kept record is an input record, unchanged in header, sequence, and quality, in input order, with mates kept together, and the count is the one asked for (exact count) or passes a two-sided binomial test at 1e-6 (fraction); the subset checker, which shares no code with z-fastq, decides. The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>

z-fastq returned the right sample on all 84 inputs it was given. Every other tool was right on every input it was given, except: seqtk gave a wrong sample on 3 of its 69 inputs (PromethION, 30 to 153 kb) and is not timed on them; IRMA Core gave a wrong sample on 2 of its 79 inputs (paired: 100 bp mate pairs) and is not timed on them.

> [!CAUTION]
> Rasusa, seqtk, fq, fqkit, SeqKit and Fasten print a wrong sample and exit 0 on some malformed inputs, so a script cannot tell the sample is wrong.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/edges-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/edges-light.svg">
  <img src="results/figures/edges-light.svg" alt="What each tool does on small edge-case inputs" width="100%">
</picture>

> [!NOTE]
> A cell is the tool's verdict when asked for every record (K = 1,000, above every edge file's count, or p = 1): `accept` (exit 0, every record back unchanged and in order), `reject S003` (with the rule it named), `reject`, `wrong: accept` (exit 0 on a malformed input), `wrong: N kept` or `wrong: not a subset` (exit 0 on a valid input with another output: records lost, changed, or reordered), `refused` (rejected a valid input), or `crash`. "tool policy" rows have no single right answer.

<details><summary>Inputs where a tool gave a wrong sample</summary>

| Tool      | Input                                           | Exit | Printed             |
| --------- | ----------------------------------------------- | ---: | ------------------- |
| IRMA Core | bench/data/synthetic/paired/medium_R1.fastq     |    1 | reject              |
| seqtk     | bench/data/real/ont-err14758576/medium.fastq    |    0 | wrong: not a subset |
| seqtk     | bench/data/real/ont-err14758576/medium.fastq.gz |    0 | wrong: not a subset |

</details>

## Is the sample fair?
> <sub>200 seeds per tool and mode on files whose headers are the record index: 79,100 short reads (single-end modes) and 14,413 pairs (pair modes), p = 0.1 and K = 1,000; first and last record on a 200-record file with K = 20 or p = 0.1. Count bias: the mean kept count's distance from N p in standard errors (exact: every run kept K). Spread: the standard deviation over the binomial's (`fixed`: the same count every run). Worst tenth: the largest distance of a tenth of the file's kept share from its expected share, in standard errors.</sub>

| Mode                 | Tool      | Count bias (z) | Spread | Worst tenth (z) | First kept | Last kept | Same twice | Seeds differ | Mates | Fair |
| -------------------- | --------- | -------------: | -----: | --------------: | ---------: | --------: | ---------- | ------------ | ----- | ---- |
| Fraction             | z-fastq   |           -2.1 |   0.97 |             1.8 |         15 |        26 | yes        | yes          | -     | yes  |
| Fraction             | Rasusa    |           -0.4 |   1.02 |             1.5 |         24 |        16 | yes        | yes          | -     | yes  |
| Fraction             | seqtk     |           -2.1 |   0.97 |             1.8 |         15 |        26 | yes        | yes          | -     | yes  |
| Fraction             | fq        |           +0.7 |   0.99 |             3.2 |         20 |        18 | yes        | yes          | -     | yes  |
| Fraction             | SeqKit    |           +1.0 |   1.04 |             1.6 |         24 |        19 | yes        | yes          | -     | yes  |
| Fraction             | IRMA Core |           +0.0 |  fixed |            15.9 |         18 |        20 | yes        | yes          | -     | no   |
| Exact count          | z-fastq   |          exact |      - |             1.4 |         18 |        26 | yes        | yes          | -     | yes  |
| Exact count          | Rasusa    |          exact |      - |             1.4 |         19 |        16 | yes        | yes          | -     | yes  |
| Exact count          | seqtk     |          exact |      - |             1.4 |         18 |        26 | yes        | yes          | -     | yes  |
| Exact count          | fq        |          exact |      - |             2.5 |         26 |        23 | yes        | yes          | -     | yes  |
| Exact count          | fqkit     |          exact |      - |             1.3 |         14 |        17 | yes        | yes          | -     | yes  |
| Exact count          | SeqKit    |          exact |      - |             1.1 |         20 |        15 | yes        | yes          | -     | yes  |
| Exact count          | IRMA Core |          exact |      - |             1.9 |         18 |        20 | yes        | yes          | -     | yes  |
| Paired fraction      | z-fastq   |           -2.1 |   0.96 |             3.1 |          - |         - | yes        | yes          | yes   | yes  |
| Paired fraction      | Rasusa    |           -1.0 |   1.03 |             2.2 |          - |         - | yes        | yes          | yes   | yes  |
| Paired fraction      | fq        |           -1.1 |   1.06 |             2.1 |          - |         - | yes        | yes          | yes   | yes  |
| Paired fraction      | IRMA Core |           -0.1 |  fixed |             2.7 |          - |         - | yes        | yes          | yes   | yes  |
| Paired exact count   | z-fastq   |          exact |      - |             2.0 |          - |         - | yes        | yes          | yes   | yes  |
| Paired exact count   | Rasusa    |          exact |      - |             1.5 |          - |         - | yes        | yes          | yes   | yes  |
| Paired exact count   | fq        |          exact |      - |             2.6 |          - |         - | yes        | yes          | yes   | yes  |
| Paired exact count   | IRMA Core |          exact |      - |             1.5 |          - |         - | yes        | yes          | yes   | yes  |
| Interleaved fraction | z-fastq   |           -2.1 |   0.96 |             3.1 |          - |         - | yes        | yes          | yes   | yes  |
| Interleaved fraction | Fasten    |           +0.7 |   0.98 |             2.5 |          - |         - | no         | yes          | yes   | yes  |

A tool fails at a count bias beyond 4.9, a spread outside 0.5 to 1.5, a tenth beyond 5.3 (each two-sided 1e-6), the first or last record kept outside 3 to 43 times of 200, or the same kept set from seeds 11 and 12. At these sizes a tenth's share is tested to about 0.7% (fraction) and 2% (exact) of its value. Same twice is reproducibility, not fairness.

## Fraction sampling, p = 0.1
> <sub>The real runs and synthetic shapes at the medium tier, about 1 GiB of FASTQ each, pooled into kinds of input. Every tool runs on one CPU with the file in the page cache.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-overview-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-overview-light.svg">
  <img src="results/figures/fraction-overview-light.svg" alt="Time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes 415.4 ms and 0.89 MiB per plain input, 1762.8 ms and 0.89 MiB per gzip input (geometric means over 12 inputs each); the table gives every tool as a multiple of that.

| Tool      | Plain time | Plain memory | Gzip time | Gzip memory | Passes malformed |
| --------- | ---------: | -----------: | --------: | ----------: | ---------------- |
| z-fastq   |      1.00x |        1.00x |     1.00x |       1.00x | none             |
| Rasusa    |      0.69x |        4.95x |     1.01x |       5.42x | 2 of 9           |
| seqtk     |      0.96x |        3.22x |     1.08x |       3.22x | 9 of 9           |
| fq        |      1.10x |        4.70x |     1.75x |       4.83x | 5 of 9           |
| SeqKit    |      1.53x |       27.73x |     2.66x |      25.79x | 3 of 9           |
| IRMA Core |      3.28x |        5.25x |     2.94x |       5.16x | none             |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-speed-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-speed-light.svg">
  <img src="results/figures/fraction-speed-light.svg" alt="Wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-memory-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-memory-light.svg">
  <img src="results/figures/fraction-memory-light.svg" alt="Peak memory by kind of input" width="100%">
</picture>

<details><summary>Minor page faults by kind of input</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-faults-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-faults-light.svg">
  <img src="results/figures/fraction-faults-light.svg" alt="Minor page faults by kind of input" width="100%">
</picture>

</details>

### Scaling
> <sub>Each synthetic shape from 1 MiB to 1 GiB of FASTQ: throughput and peak memory against input size. A shape changes one thing a sampler pays for, so its line shows what that one thing costs each tool.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-scaling-speed-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-scaling-speed-light.svg">
  <img src="results/figures/fraction-scaling-speed-light.svg" alt="Throughput by input size" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/fraction-scaling-memory-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/fraction-scaling-memory-light.svg">
  <img src="results/figures/fraction-scaling-memory-light.svg" alt="Peak memory by input size" width="100%">
</picture>

At 1 GiB, z-fastq samples 1.5 to 5.3 GiB/s on plain input and 0.4 to 0.7 GiB/s on gzip input. Faster there: Rasusa on plain long (1.74x); Rasusa on plain variable (1.52x); Rasusa on plain longread (1.45x); Rasusa on gzip dense (1.01x); Rasusa on plain endings (1.49x); Rasusa on plain dense (1.66x); seqtk on plain long (1.45x); seqtk on plain longread (1.25x); seqtk on plain endings (1.16x); seqtk on plain dense (1.08x), throughput as a multiple of z-fastq.
z-fastq peaks at 0.76 to 1.77 MiB across sizes and shapes. Memory that grows with the input: seqtk by 1.4x on plain longread; fq by 1.3x on gzip longread; SeqKit by 1.5x on plain dense; IRMA Core by 1.4x on plain longread.

## Exact count, K = 100,000 (1,000 for long reads)
> <sub>The real runs and synthetic shapes at the medium tier, about 1 GiB of FASTQ each, pooled into kinds of input. Every tool runs on one CPU with the file in the page cache.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/exact-overview-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/exact-overview-light.svg">
  <img src="results/figures/exact-overview-light.svg" alt="Time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes 518.0 ms and 1.28 MiB per plain input (geometric mean over 12 inputs); the table gives every tool as a multiple of that.

| Tool      | Plain time | Plain memory | Passes malformed |
| --------- | ---------: | -----------: | ---------------- |
| z-fastq   |      1.00x |        1.00x | none             |
| Rasusa    |      0.90x |        8.18x | 2 of 9           |
| seqtk     |      1.43x |        3.08x | 9 of 9           |
| fq        |      1.17x |        3.51x | 5 of 9           |
| fqkit     |      1.38x |       27.96x | 7 of 9           |
| SeqKit    |      2.24x |       18.78x | 3 of 9           |
| IRMA Core |      2.36x |        3.67x | none             |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/exact-speed-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/exact-speed-light.svg">
  <img src="results/figures/exact-speed-light.svg" alt="Wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/exact-memory-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/exact-memory-light.svg">
  <img src="results/figures/exact-memory-light.svg" alt="Peak memory by kind of input" width="100%">
</picture>

<details><summary>Minor page faults by kind of input</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/exact-faults-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/exact-faults-light.svg">
  <img src="results/figures/exact-faults-light.svg" alt="Minor page faults by kind of input" width="100%">
</picture>

</details>

Plain input only: with gzip, most exact-count tools decompress the file twice, so the time is mostly zlib; fraction sampling shows the gzip cost. K is 100,000 reads for short-read inputs and 1,000 for long-read inputs (the Long reads group, `longread`, and `long`); the `long` shape holds 1,789 records, so its 1,000 is close to a copy.

### Time and memory against K

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/exact-k-lines-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/exact-k-lines-light.svg">
  <img src="results/figures/exact-k-lines-light.svg" alt="Time and peak memory against K" width="100%">
</picture>

|         K | z-fastq peak | Lowest other     | Highest other  |
| --------: | -----------: | ---------------- | -------------- |
|     1,000 |     0.77 MiB | 3.5x (seqtk)     | 37.9x (SeqKit) |
|    10,000 |     1.02 MiB | 2.9x (seqtk)     | 18.7x (SeqKit) |
|   100,000 |     1.27 MiB | 3.4x (IRMA Core) | 15.3x (Rasusa) |
| 1,000,000 |     3.77 MiB | 1.2x (IRMA Core) | 9.5x (SeqKit)  |

## Pairs
> <sub>Each tool's time and peak memory as a multiple of z-fastq's, geometric mean over the 4 paired runs and the `paired` shape at 1 GiB. Plain input: paired fraction and paired exact count read two mate files, interleaved fraction one interleaved file on stdin.</sub>

| Mode                 | Tool      | Plain time | Plain memory | Passes malformed |
| -------------------- | --------- | ---------: | -----------: | ---------------- |
| Paired fraction      | z-fastq   |      1.00x |        1.00x | none             |
| Paired fraction      | Rasusa    |      0.57x |        4.50x | 2 of 9           |
| Paired fraction      | fq        |      0.90x |        3.97x | 5 of 9           |
| Paired fraction      | IRMA Core |      2.54x |        4.32x | none             |
| Paired exact count   | z-fastq   |      1.00x |        1.00x | none             |
| Paired exact count   | Rasusa    |      0.74x |        9.65x | 2 of 9           |
| Paired exact count   | fq        |      0.73x |        2.86x | 5 of 9           |
| Paired exact count   | IRMA Core |      1.78x |        2.93x | none             |
| Interleaved fraction | z-fastq   |      1.00x |        1.00x | none             |
| Interleaved fraction | Fasten    |      2.59x |        3.12x | 7 of 9           |

## Reproducibility

<details><summary>How it was measured</summary>

Every tool runs through `zexec -r` inside one zebrac batch per case, held to one CPU core by a systemd scope with a CPU quota and a 4G memory cap, after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. The kernel reports peak memory in 256 KiB steps, so values closer than that are equal. Runs per tier, each after one warmup: tiny 31, xsmall 21, small 11, medium 3; sizes 1 MiB to 1 GiB. Generated by `bench/sample/run.sh release` on 2026-10-24.

</details>

<details><summary>Identity</summary>

- `z-fastq`: *0.1.0 at 90e46b4*
- `zebrac`: *0.6.2*
- `zexec`: *12448 bytes, static, no libc*
- `kernel`: *7.2.11-200.fc44.x86_64*
- `cpu`: *AMD Ryzen 9 3950X 16-Core Processor*
- `load at batch start`: *2.00 to 13.41*

</details>

<details><summary>Commands</summary>

- `z-fastq`, fraction: *zig-out/bin/z-fastq sample --fraction 0.1 --seed 11 IN*
- `Rasusa`, fraction: *tools/bin/rasusa reads -f 0.1 -1 -s 11 IN*
- `seqtk`, fraction: *tools/bin/seqtk sample -s 11 IN 0.1*
- `fq`, fraction: *tools/bin/fq subsample -p 0.1 -s 11 --r1-dst /dev/null IN*
- `SeqKit`, fraction: *tools/bin/seqkit sample -p 0.1 -s 11 -j 1 IN*
- `IRMA Core`, fraction: *tools/bin/irma-core sampler -p 10 -s 11 IN*
- `z-fastq`, exact: *zig-out/bin/z-fastq sample --count 100000 --seed 11 IN*
- `Rasusa`, exact: *tools/bin/rasusa reads -n 100000 -s 11 IN*
- `seqtk`, exact: *tools/bin/seqtk sample -2 -s 11 IN 100000*
- `fq`, exact: *tools/bin/fq subsample -n 100000 -s 11 --r1-dst /dev/null IN*
- `fqkit`, exact: *tools/bin/fqkit sample -n 100000 -s 11 -2 -@ 1 IN*
- `SeqKit`, exact: *tools/bin/seqkit sample2 -n 100000 -2 -w 0 -s 11 -j 1 IN*
- `IRMA Core`, exact: *tools/bin/irma-core sampler -t 100000 -s 11 IN*
- `z-fastq`, exact-k: *zig-out/bin/z-fastq sample --count 1000 --seed 11 IN*
- `Rasusa`, exact-k: *tools/bin/rasusa reads -n 1000 -s 11 IN*
- `seqtk`, exact-k: *tools/bin/seqtk sample -2 -s 11 IN 1000*
- `fq`, exact-k: *tools/bin/fq subsample -n 1000 -s 11 --r1-dst /dev/null IN*
- `fqkit`, exact-k: *tools/bin/fqkit sample -n 1000 -s 11 -2 -@ 1 IN*
- `SeqKit`, exact-k: *tools/bin/seqkit sample2 -n 1000 -2 -w 0 -s 11 -j 1 IN*
- `IRMA Core`, exact-k: *tools/bin/irma-core sampler -t 1000 -s 11 IN*
- `z-fastq`, paired-fraction: *zig-out/bin/z-fastq sample --paired --fraction 0.1 --seed 11 IN*
- `Rasusa`, paired-fraction: *tools/bin/rasusa reads -f 0.1 -1 -s 11 -o /dev/null -o /dev/null IN*
- `fq`, paired-fraction: *tools/bin/fq subsample -p 0.1 -s 11 --r1-dst /dev/null --r2-dst /dev/null IN*
- `IRMA Core`, paired-fraction: *tools/bin/irma-core sampler -p 10 -s 11 IN*
- `z-fastq`, paired-exact: *zig-out/bin/z-fastq sample --paired --count 100000 --seed 11 IN*
- `Rasusa`, paired-exact: *tools/bin/rasusa reads -n 100000 -s 11 -o /dev/null -o /dev/null IN*
- `fq`, paired-exact: *tools/bin/fq subsample -n 100000 -s 11 --r1-dst /dev/null --r2-dst /dev/null IN*
- `IRMA Core`, paired-exact: *tools/bin/irma-core sampler -t 100000 -s 11 IN*
- `z-fastq`, interleaved-fraction: *zig-out/bin/z-fastq sample --interleaved --fraction 0.1 --seed 11 IN*
- `Fasten`, interleaved-fraction: *tools/bin/fasten_sample -p -f 0.1*

</details>

<details><summary>Quick-check cases for peer and ab</summary>

Picked from these results: per group and format, the case where z-fastq leads its fastest full-check rival by the least. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.

- `fraction: Long reads, plain`: *PromethION, 30 to 153 kb*; rival IRMA Core at 2.76x z-fastq's time
- `fraction: Short, fixed length, plain`: *BGISEQ, 100 bp*; rival IRMA Core at 3.27x z-fastq's time
- `fraction: Short, variable length, plain`: *10x mates, 28 and 150 bp*; rival IRMA Core at 3.52x z-fastq's time
- `fraction: Synthetic, plain`: *long: 300 kb sequence lines*; rival IRMA Core at 2.35x z-fastq's time
- `fraction: Long reads, gzip`: *PromethION, 30 to 153 kb*; rival IRMA Core at 2.83x z-fastq's time
- `fraction: Short, fixed length, gzip`: *miRNA, 36 bp*; rival IRMA Core at 2.27x z-fastq's time
- `fraction: Short, variable length, gzip`: *10x mates, 28 and 150 bp*; rival IRMA Core at 2.78x z-fastq's time
- `fraction: Synthetic, gzip`: *endings: 100 bp, CRLF and LF lines*; rival IRMA Core at 2.80x z-fastq's time
- `exact: Long reads, plain`: *PromethION, 30 to 153 kb*; rival IRMA Core at 1.72x z-fastq's time
- `exact: Short, fixed length, plain`: *NovaSeq, Casava headers*; rival IRMA Core at 2.40x z-fastq's time
- `exact: Short, variable length, plain`: *10x mates, 28 and 150 bp*; rival IRMA Core at 2.90x z-fastq's time
- `exact: Synthetic, plain`: *long: 300 kb sequence lines*; rival IRMA Core at 1.51x z-fastq's time

</details>

<details><summary>Hardware counters (user space, plain medium-tier inputs, geometric means)</summary>

| Tool      | Instructions per byte | Instructions per cycle | Cache misses per MiB | Branch misses per MiB |
| --------- | --------------------: | ---------------------: | -------------------: | --------------------: |
| z-fastq   |                  2.76 |                   3.01 |                 2285 |                   430 |
| Rasusa    |                  1.09 |                   2.51 |                  881 |                   779 |
| seqtk     |                  2.00 |                   2.49 |                 4109 |                  2066 |
| fq        |                  2.90 |                   3.41 |                  487 |                  1602 |
| SeqKit    |                  3.59 |                   2.36 |                 1962 |                  2960 |
| IRMA Core |                  9.31 |                   3.14 |                 2337 |                  5297 |

</details>

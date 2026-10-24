<!-- markdownlint-disable MD013 MD033 MD041 -->

# Counting records
> <sub>`z-fastq count` against 5 other counters, on 11 real sequencing runs and 5 synthetic shapes, plain and gzip. A tool is timed only where its count was right; time and memory come from the same run, on one CPU.</sub>

- `z-fastq`: *0.1.0 at 9590e98*
- `machine`: *Ryzen 9 3950X, Linux 7.2.11; one CPU per run*
- `runner`: *zebrac 0.6.2 with zexec*
- `data`:
  - **Real**: *11 public sequencing runs, 100 MiB each, pooled in the bars by the kind of reads they hold*
    - `Short, fixed length`: *MiSeq, 151 bp; NovaSeq, 150 bp; NovaSeq, Casava headers; BGISEQ, 100 bp; miRNA, 36 bp*
    - `Short, variable length`: *MiSeq, 35 to 151 bp; 10x mates, 28 and 150 bp*
    - `Long reads`: *GridION, 0.1 to 92 kb; PromethION, 30 to 153 kb; Sequel II, 0.4 to 30 kb; Revio HiFi, 4 to 21 kb*
  - **Synthetic**: *5 shapes generated from real reads, each changing one thing; one group in the bars at 100 MiB, and the scaling lines from 1 MiB to 1 GiB*
    - `dense`: *uniform 150 bp, bare plus line*
    - `endings`: *100 bp, CRLF and LF lines*
    - `variable`: *50 to 150 bp, annotated plus line*
    - `longread`: *nanopore headers, rare 300 kb reads*
    - `long`: *300 kb sequence lines*
- `date`: *2026-10-24*

## Findings

- **Right answers:** z-fastq returned the right count on all 218 checked inputs and rejected every malformed edge case. 4 of the 5 other tools print a wrong count or crash on at least one.
- **Speed:** Needletail, the fastest other full-check tool, takes 1.29x z-fastq's time on plain input and 1.15x on gzip.
- **Memory:** every other tool peaks at 3.7x to 29x z-fastq's 0.56 MiB.
- **Where z-fastq is not fastest:** Needletail on plain long reads (0.94x).

## Compared tools
> <sub>*complete* tools check as much as z-fastq; *partial* tools check less, so they do less work. The bold count is how many of three malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, and `-` instead of `+`.</sub>

- `z-fastq 0.1.0` **(subject)**; *checks every record's structure*; **3 of 3 rejected**
- `Needletail 0.7.3` **(complete)**; *checks every record's structure; built for this CPU*; **3 of 3 rejected**
- `Helicase 0.2.0` **(partial)**; *misses a truncated last record and a bad plus line*; **1 of 3 rejected**
- `seqtk 1.5` **(partial)**; *no structure check*; **0 of 3 rejected**
- `SeqFu 1.29.0` **(partial)**; *accepts a truncated last record; starts 8 threads by default, held to one CPU here like every tool*; **2 of 3 rejected**
- `fqtools 2.3` **(complete)**; *checks every record's structure*; **3 of 3 rejected**

> [!IMPORTANT]
> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.

## Right answers
> <sub>Each count is compared with the count recorded when the corpus was cut, without z-fastq. The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>

z-fastq returned the right count on all 218 inputs: every corpus file and synthetic shape on disk, 1 MiB to 1 GiB, plain and gzip. Every other tool was right on every input it was given, except: fqtools missed the count on 2 of its 32 inputs (endings: 100 bp, CRLF and LF lines) and is not timed on them.

> [!CAUTION]
> Helicase, seqtk, SeqFu and fqtools print a wrong count and exit 0 on some malformed inputs, so a script cannot tell the count is wrong. Helicase crashes on some malformed gzip input.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/edges-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/edges-light.svg">
  <img src="results/figures/edges-light.svg" alt="What each tool does on small edge-case inputs" width="100%">
</picture>

> [!NOTE]
> A cell is the count printed, `reject` (exited with an error), `wrong: N` (printed another count, exit 0), `refused` (rejected a valid input), or `crash`. "tool policy" rows have no single right answer.

<details><summary>Inputs where a tool missed the count</summary>

| Tool    | Input                                       | Exit | Printed |
| ------- | ------------------------------------------- | ---: | ------- |
| fqtools | bench/data/synthetic/endings/small.fastq    |    1 | nothing |
| fqtools | bench/data/synthetic/endings/small.fastq.gz |    1 | nothing |

</details>

## Speed and memory
> <sub>The small-tier shapes, about 100 MiB of FASTQ each, pooled into kinds of input. Every tool runs on one CPU with the file in the page cache.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-overview-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-overview-light.svg">
  <img src="results/figures/file-overview-light.svg" alt="Time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes 16.6 ms and 0.56 MiB per plain shape, 146.1 ms and 0.56 MiB per gzip shape (geometric means over 16 shapes each); the table gives every tool as a multiple of that.

| Tool       | Plain time | Plain memory | Gzip time | Gzip memory |
| ---------- | ---------: | -----------: | --------: | ----------: |
| z-fastq    |      1.00x |        1.00x |     1.00x |       1.00x |
| Needletail |      1.29x |        3.75x |     1.15x |       4.15x |
| Helicase   |      1.17x |        4.56x |     1.06x |       4.71x |
| seqtk      |      2.06x |        4.55x |     1.24x |       4.54x |
| SeqFu      |      2.34x |        8.65x |     2.22x |      28.55x |
| fqtools    |     19.66x |        5.77x |     3.12x |       6.64x |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-speed-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-speed-light.svg">
  <img src="results/figures/file-speed-light.svg" alt="Wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-memory-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-memory-light.svg">
  <img src="results/figures/file-memory-light.svg" alt="Peak memory by kind of input" width="100%">
</picture>

<details><summary>Minor page faults by kind of input</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-faults-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-faults-light.svg">
  <img src="results/figures/file-faults-light.svg" alt="Minor page faults by kind of input" width="100%">
</picture>

</details>

## Scaling
> <sub>Each synthetic shape from 1 MiB to 1 GiB of FASTQ: throughput and peak memory against input size. A shape changes one thing a counter pays for, so its line shows what that one thing costs each tool.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-scaling-speed-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-scaling-speed-light.svg">
  <img src="results/figures/file-scaling-speed-light.svg" alt="Throughput by input size" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-scaling-memory-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-scaling-memory-light.svg">
  <img src="results/figures/file-scaling-memory-light.svg" alt="Peak memory by input size" width="100%">
</picture>

At 1 GiB, z-fastq counts 4.4 to 8.3 GiB/s on plain input and 0.4 to 1.0 GiB/s on gzip input. Faster there: Needletail on plain longread (1.18x); Helicase on plain longread (1.03x), throughput as a multiple of z-fastq.
z-fastq peaks at 0.56 MiB at every size and shape. Memory that grows with the input: Needletail by 1.5x on plain longread; Helicase by 1.3x on plain longread; seqtk by 1.4x on plain longread.

> [!NOTE]
> fqtools is timed at the small tier only: it is slow enough that sweeping sizes with it would take most of the run. Every tool is held to one CPU, so a tool that starts threads (SeqFu starts 8 by default) gains nothing from them here; its CPU time equals its wall time. Lines rise with size because each tool's fixed start cost spreads over more data; SeqFu has the largest start cost, so its line rises the most.

## Other input forms
> <sub>The same counting, with the input on stdin or as many files in one call.</sub>

### Stdin

The same file on stdin.

Each tool's time and peak memory reading the file from a pipe, as a multiple of its own with the path (geometric mean over the shared shapes):

| Tool       | Plain time | Plain memory | Gzip time | Gzip memory |
| ---------- | ---------: | -----------: | --------: | ----------: |
| z-fastq    |      1.06x |        1.00x |     1.04x |       1.00x |
| Needletail |      1.09x |        0.93x |     1.04x |       1.00x |
| Helicase   |      1.05x |        0.92x |     1.03x |       0.88x |
| seqtk      |      1.07x |        1.04x |     1.04x |       1.06x |
| SeqFu      |      1.07x |        0.98x |     0.99x |       0.99x |

### Many files

Every xsmall corpus file in one call.

Counting all 11 xsmall corpus files (110 MiB) in one call, z-fastq takes 18.5 ms on plain input and 186.2 ms on gzip, at 0.56 MiB peak. SeqFu, the only other tool that counts many files in one call, takes 2.7x and 2.8x as long with 10x and 138x the memory.

## Reproducibility

<details><summary>How it was measured</summary>

Every tool runs through `zexec -r` inside one zebrac batch per case, pinned to CPU 3 under a systemd scope with a CPU quota and a 4G memory cap, after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. The kernel reports peak memory in 256 KiB steps, so values closer than that are equal. Runs per tier, each after one warmup: tiny 31, xsmall 21, small 11, medium 5; sizes 1 MiB to 1 GiB. Generated by `bench/count/run.sh release` on 2026-10-24.

</details>

<details><summary>Identity</summary>

- `z-fastq`: *0.1.0 at 9590e98*
- `zebrac`: *0.6.2*
- `zexec`: *12448 bytes, static, no libc*
- `kernel`: *7.2.11-200.fc44.x86_64*
- `cpu`: *AMD Ryzen 9 3950X 16-Core Processor*
- `load at batch start`: *0.53 to 7.25*

</details>

<details><summary>Commands</summary>

- `z-fastq`: *zig-out/bin/z-fastq count IN*
- `Needletail`: *bench/cache/adapters/native/needletail-adapter count IN*
- `Helicase`: *bench/cache/adapters/native/helicase-adapter count IN*
- `seqtk`: *tools/bin/seqtk size IN*
- `SeqFu`: *tools/bin/seqfu count IN*
- `fqtools`: *tools/bin/fqtools count IN*

</details>

<details><summary>Quick-check cases for peer and ab</summary>

Picked from these results: per group and format, the case where z-fastq leads its fastest full-check rival by the least. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.

- `Long reads, plain`: *PromethION, 30 to 153 kb*; rival Needletail at 0.86x z-fastq's time
- `Short, fixed length, plain`: *NovaSeq, Casava headers*; rival Needletail at 1.38x z-fastq's time
- `Short, variable length, plain`: *MiSeq, 35 to 151 bp*; rival Needletail at 1.23x z-fastq's time
- `Synthetic, plain`: *longread: nanopore headers, rare 300 kb reads*; rival Needletail at 0.97x z-fastq's time
- `Long reads, gzip`: *PromethION, 30 to 153 kb*; rival Needletail at 1.04x z-fastq's time
- `Short, fixed length, gzip`: *NovaSeq, Casava headers*; rival Needletail at 1.16x z-fastq's time
- `Short, variable length, gzip`: *MiSeq, 35 to 151 bp*; rival Needletail at 1.15x z-fastq's time
- `Synthetic, gzip`: *longread: nanopore headers, rare 300 kb reads*; rival Needletail at 1.06x z-fastq's time

</details>

<details><summary>Hardware counters (user space, plain small-tier shapes, geometric means)</summary>

| Tool       | Instructions per byte | Instructions per cycle | Cache misses per MiB | Branch misses per MiB |
| ---------- | --------------------: | ---------------------: | -------------------: | --------------------: |
| z-fastq    |                  0.74 |                   4.33 |                  817 |                   182 |
| Needletail |                  0.61 |                   2.48 |                  555 |                   579 |
| Helicase   |                  0.76 |                   3.18 |                 1243 |                   384 |
| seqtk      |                  1.40 |                   2.44 |                 3209 |                  1632 |
| SeqFu      |                  1.79 |                   2.32 |                  768 |                  1644 |
| fqtools    |                 36.69 |                   3.62 |                 1117 |                  2852 |

</details>

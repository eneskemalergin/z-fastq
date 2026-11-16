<!-- markdownlint-disable MD013 MD033 MD041 -->

# Counting records
> <sub>`z-fastq count` against 5 other counters, on 7 real sequencing runs and 5 synthetic shapes, plain and gzip. A tool is timed only where its count was right; time and memory come from the same run, on one CPU.</sub>

- `z-fastq`: *0.3.0*
- `machine`: *Ryzen 9 3950X, Linux; one CPU per run*
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

## Findings

- **Right answers:** z-fastq returned the right count on all 218 checked inputs and rejected every malformed edge case. 4 of the 5 other tools print a wrong count or crash on at least one.
- **Speed:** Needletail, the fastest other full-check tool, takes 1.40x z-fastq's time on plain input and 1.39x on gzip.
- **Memory:** every other tool peaks at 4.1x to 29x z-fastq's 0.52 MiB.
- **Where z-fastq is not fastest:** on no kind of input; it is the fastest tool on every kind of input (geometric mean over the group), plain and gzip.
- **Wrong on valid input, exit 0:** two gzip members: Helicase (printed 2).

## Compared tools
> <sub>*complete* tools check as much as z-fastq; *partial* tools check less, so they do less work. The bold count is how many of three malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, and `-` instead of `+`.</sub>

- `z-fastq 0.3.0` **(subject)**; *checks every record's structure*; **3 of 3 rejected**
- `Needletail 0.7.3` **(complete)**; *checks every record's structure; built for this CPU*; **3 of 3 rejected**
- `Helicase 0.2.0` **(partial)**; *misses a truncated last record and a bad plus line*; **1 of 3 rejected**
- `SeqFu 1.29.0` **(partial)**; *accepts a truncated last record; starts 8 threads by default, held to one CPU here like every tool*; **2 of 3 rejected**
- `seqtk 1.5` **(partial)**; *no structure check*; **0 of 3 rejected**
- `fqtools 2.3` **(complete)**; *checks every record's structure*; **3 of 3 rejected**

> [!IMPORTANT]
> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.

## Right answers
> <sub>Each count is compared with the count recorded when the corpus was cut, without z-fastq. The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>

z-fastq returned the right count on all 218 inputs: every corpus file and synthetic shape on disk, 1 MiB to 1 GiB, plain and gzip. Every other tool was right on every input it was timed on except these, which are not timed where they missed the count:

| Tool    | Wrong |  Of | Inputs  |
| ------- | ----: | --: | ------- |
| fqtools |     2 |  24 | endings |

> [!CAUTION]
> On the edge cases below, Helicase, SeqFu, seqtk and fqtools print a wrong count and exit 0 on some malformed inputs, so a script cannot tell the count is wrong. Helicase crashes on some malformed gzip input.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/edges-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/edges-light.png">
  <img src="results/figures/edges-light.png" alt="What each tool does on small edge-case inputs" width="100%">
</picture>

> [!NOTE]
> A cell is the count printed, `reject` (exited with an error), `wrong: N` (printed another count, exit 0), `refused` (rejected a valid input), or `crash`. "tool policy" rows have no single right answer.

<details><summary>Inputs where a tool missed the count</summary>

| Tool    | Input                                        | Exit | Printed |
| ------- | -------------------------------------------- | ---: | ------- |
| fqtools | bench/data/synthetic/endings/medium.fastq    |    1 | nothing |
| fqtools | bench/data/synthetic/endings/medium.fastq.gz |    1 | nothing |

</details>

## Speed and memory
> <sub>The real runs and synthetic shapes at the medium tier, about 1 GiB of FASTQ each, pooled into kinds of input. Every tool runs on one CPU, with the file in the page cache.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-overview-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-overview-light.png">
  <img src="results/figures/file-overview-light.png" alt="Time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes 144.3 ms and 0.52 MiB per plain input, 1176.6 ms and 0.52 MiB per gzip input (geometric means over 12 inputs each); the table gives every other tool as a multiple of that.
Geometric-mean process CPU time is 135.7 ms for plain input and 1090.7 ms for gzip. Elapsed time also includes launcher overhead and time waiting to run or for I/O.

| Tool       | Plain time | Plain memory | Gzip time | Gzip memory |
| ---------- | ---------: | -----------: | --------: | ----------: |
| z-fastq    |     144 ms |     0.52 MiB |    1.18 s |    0.52 MiB |
| Needletail |      1.40x |        4.09x |     1.39x |       4.59x |
| Helicase   |      1.27x |        5.09x |     1.27x |       5.33x |
| SeqFu      |      2.16x |        9.49x |     2.52x |      28.98x |
| seqtk      |      2.36x |        5.11x |     1.53x |       5.04x |
| fqtools    |     23.66x |        6.29x |     4.06x |       7.30x |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-speed-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-speed-light.png">
  <img src="results/figures/file-speed-light.png" alt="Wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-memory-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-memory-light.png">
  <img src="results/figures/file-memory-light.png" alt="Peak memory by kind of input" width="100%">
</picture>

<details><summary>Minor page faults by kind of input</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-faults-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-faults-light.png">
  <img src="results/figures/file-faults-light.png" alt="Minor page faults by kind of input" width="100%">
</picture>

</details>

## Scaling
> <sub>Each synthetic shape from 1 MiB to 1 GiB of FASTQ: throughput and peak memory against input size. A shape changes one thing a counter pays for, so its line shows what that one thing costs each tool.</sub>


<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-scaling-speed-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-scaling-speed-light.png">
  <img src="results/figures/file-scaling-speed-light.png" alt="Throughput by input size" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/file-scaling-memory-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/file-scaling-memory-light.png">
  <img src="results/figures/file-scaling-memory-light.png" alt="Peak memory by input size" width="100%">
</picture>

At 1 GiB, z-fastq counts 5.5 to 9.5 GiB/s on plain input and 0.5 to 1.0 GiB/s on gzip input; no other tool is faster there.
z-fastq peaks at 0.52 MiB at every size of the dense, endings, long, longread and variable shapes. Memory that grows with the input: Needletail by 1.4x on gzip longread; Helicase by 1.3x on plain longread; seqtk by 1.4x on plain longread.

> [!NOTE]
> fqtools is timed at the medium tier only: it is slow enough that sweeping sizes with it would take most of the run. Every tool is held to one CPU, so a tool that starts threads (SeqFu starts 8 by default) cannot use extra CPU cores here. Elapsed time can still exceed process CPU time because of scheduling and I/O waits. Lines rise with size because each tool's fixed start cost spreads over more data; SeqFu has the largest start cost, so its line rises the most.

## Other input forms
> <sub>The same counting, with the input on stdin or as many files in one call.</sub>

### Stdin

The same file on stdin.

Each tool's time and peak memory reading the file from a pipe, as a multiple of its own with the path (geometric mean over the shared shapes):

| Tool       | Plain time | Plain memory | Gzip time | Gzip memory |
| ---------- | ---------: | -----------: | --------: | ----------: |
| z-fastq    |      0.96x |        1.00x |     0.93x |       1.00x |
| Needletail |      0.99x |        1.00x |     0.94x |       0.98x |
| Helicase   |      0.97x |        0.89x |     0.94x |       0.85x |
| SeqFu      |      0.98x |        1.01x |     0.94x |       0.99x |
| seqtk      |      0.97x |        0.99x |     0.92x |       1.01x |

### Many files

Every xsmall corpus file in one call.

Counting all 11 xsmall corpus files (110 MiB) in one call, z-fastq takes 15.2 ms on plain input and 130.1 ms on gzip, at 0.52 MiB peak. SeqFu, the only other tool that counts many files in one call, takes 2.8x and 3.2x as long with 10x and 150x the memory.

## Reproducibility

<details><summary>How it was measured</summary>

Every tool runs through `zexec -r` inside one zebrac batch per case, held to one CPU core by a systemd scope with a CPU quota and a 4G memory cap, after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. Peak RSS is the kernel's peak counter that `wait4` reports, which can lag the true peak by up to about 0.5 MiB, so memory differences smaller than that are inconclusive. Runs per tier, each after one warmup: tiny 31, xsmall 21, small 11, medium 3; sizes 1 MiB to 1 GiB. Timed by `bench/count/run.sh release`.

</details>

<details><summary>Identity</summary>

- `z-fastq`: *0.3.0*
- `zebrac`: *0.6.2*
- `zexec`: *15000 bytes, static, no libc*
- `cpu`: *AMD Ryzen 9 3950X 16-Core Processor*
- `load at batch start`: *1.51 to 3.46*

</details>

<details><summary>Commands</summary>

- `z-fastq`: *zig-out/bin/z-fastq count IN*
- `Needletail`: *bench/cache/adapters/native/needletail-adapter count IN*
- `Helicase`: *bench/cache/adapters/native/helicase-adapter count IN*
- `SeqFu`: *tools/bin/seqfu count IN*
- `seqtk`: *tools/bin/seqtk size IN*
- `fqtools`: *tools/bin/fqtools count IN*

</details>

<details><summary>Quick-check cases for peer and ab</summary>

Picked from these results: per group and format, the case where z-fastq leads its fastest full-check rival by the least. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.

- `Long reads, plain`: *PromethION, 30 to 153 kb*; rival Needletail at 0.99x z-fastq's time
- `Short, fixed length, plain`: *NovaSeq, 150 bp*; rival Needletail at 1.41x z-fastq's time
- `Short, variable length, plain`: *10x mates, 28 and 150 bp*; rival Needletail at 1.52x z-fastq's time
- `Synthetic, plain`: *longread: nanopore headers, rare 300 kb reads*; rival Needletail at 1.06x z-fastq's time
- `Long reads, gzip`: *PromethION, 30 to 153 kb*; rival Needletail at 1.22x z-fastq's time
- `Short, fixed length, gzip`: *BGISEQ, 100 bp*; rival Needletail at 1.42x z-fastq's time
- `Short, variable length, gzip`: *10x mates, 28 and 150 bp*; rival Needletail at 1.45x z-fastq's time
- `Synthetic, gzip`: *longread: nanopore headers, rare 300 kb reads*; rival Needletail at 1.26x z-fastq's time

</details>

<details><summary>Hardware counters (user space, plain medium-tier inputs, geometric means, tools at one thread)</summary>

| Tool       | Instructions per byte | Instructions per cycle | Cache misses per MiB | Branch misses per MiB |
| ---------- | --------------------: | ---------------------: | -------------------: | --------------------: |
| z-fastq    |                  0.46 |                   3.49 |                 1120 |                   113 |
| Needletail |                  0.68 |                   2.57 |                  426 |                   390 |
| Helicase   |                  0.80 |                   3.28 |                 1118 |                   231 |
| SeqFu      |                  1.35 |                   2.42 |                  440 |                  1251 |
| seqtk      |                  1.55 |                   2.49 |                 3178 |                  1463 |
| fqtools    |                 36.57 |                   3.51 |                  919 |                  2731 |

</details>

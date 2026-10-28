<!-- markdownlint-disable MD013 MD033 MD041 -->

# Processing paired reads
> <sub>`z-fastq process` against 3 other processors, on 1 real sequencing run and 0 synthetic shapes, plain and gzip. A tool is timed only where its outputs were right; time and memory come from the same run, on one CPU, a threaded variant on as many cores as its threads.</sub>

- `z-fastq`: *0.1.0 at b199bc6*
- `machine`: *Ryzen 9 3950X, Linux 7.2.11; one CPU per run, threaded variants on as many cores as threads*
- `runner`: *zebrac 0.6.2 with zexec*
- `data`:
  - **Real**: *1 public sequencing runs, 1 GiB each, pooled in the bars by the kind of reads they hold*
    - `Short, fixed length`: *NovaSeq, 150 bp*
  - **Synthetic**: *0 shapes generated from real reads, each changing one thing; one group in the bars at 1 GiB, and the scaling lines from 1 MiB to 10 MiB*
- `date`: *2026-10-28*

## Findings

- **Right answers:** z-fastq wrote the right copies, fastp's records, and FastQC's tables on all 0 checked inputs and rejected every broken pair and record. 3 of the 3 other tools write output from a broken pair or record with exit 0.
- **Copying into plain files:** fastp takes 18.2x z-fastq's time on plain input and 7.9x on gzip, with 484x its memory; at its fastest, 8 threads, 4.77x z-fastq's time for 17.7x its CPU time.
- **Copying into gzip files:** fastp takes 7.0x z-fastq's time on plain input and 5.5x on gzip, with 401x its memory; at its fastest, 8 threads, 1.26x z-fastq's time for 6.7x its CPU time.
- **Trimming, TruSeq settings:** fastp takes 9.8x z-fastq's time on plain input and 8.1x on gzip, with 29x its memory; at its fastest, 8 threads, 1.32x z-fastq's time for 8.7x its CPU time.
- **Trimming, Nextera settings:** fastp takes 8.6x z-fastq's time on plain input and 7.0x on gzip, with 29x its memory; at its fastest, 8 threads, 1.17x z-fastq's time for 7.4x its CPU time.
- **QC tables:** FastQC-Rust takes 2.0x z-fastq's time on plain input and 1.6x on gzip, with 3x its memory; at its fastest, 8 threads, 0.73x z-fastq's time for 3.9x its CPU time (faster than z-fastq from 4 threads); FastQC takes 10.7x z-fastq's time on plain input and 9.2x on gzip, with 17x its memory; at its fastest, 8 threads, 2.65x z-fastq's time for 13.5x its CPU time.

## Compared tools
> <sub>*partial* tools do less than z-fastq: fastp checks records less strictly, and the FastQC tools write no reads and check no pairs. A job is `yes` when the tool was right on every input of it, `partly` on some, `-` when the tool has no such job. Threads: z-fastq runs on one; every other tool runs at 1 thread and again at 2, 4, 8, and 16 on as many physical cores, at the 1 GiB inputs, in the same batches as z-fastq. Malformed rejected: of the broken pairs and records in the edge cases, how many the tool rejects.</sub>

| Tool        | Class   | Copy | Gzip copy | TruSeq trim | Nextera trim | QC tables | Threads      | Malformed rejected |
| ----------- | ------- | ---- | --------- | ----------- | ------------ | --------- | ------------ | -----------------: |
| z-fastq     | subject | -    | -         | -           | -            | -         | 1            |             8 of 8 |
| fastp       | partial | yes  | yes       | yes         | yes          | -         | 1 to 16 (-w) |             0 of 8 |
| FastQC-Rust | partial | -    | -         | -           | -            | yes       | 1 to 16 (-t) |             1 of 8 |
| FastQC      | partial | -    | -         | -           | -            | yes       | 1 to 16 (-t) |             1 of 8 |

- `fastp 1.4.0`: *writes its JSON and HTML report and measures duplication on every run; compression level 4*
- `FastQC-Rust 1.1.0`: *QC tables only, no reads written and no pair checks; the QC tables' reference*
- `FastQC 0.13.0`: *QC tables only, no reads written and no pair checks; a Java program whose JVM starts its own threads*

> [!IMPORTANT]
> fastp computes its QC report and pair duplication on every run, so z-fastq's trimming lanes write the QC summary and measure pair duplication too. In the QC job z-fastq also writes both mates, which it cannot skip, while the FastQC tools write tables only.

## Right answers
> <sub>Copies must equal the corpus pairs byte for byte. Trimmed pairs must hold the records fastp 1.4.0 keeps with the same settings, and every field of z-fastq's QC summary must equal fastp's (its own schema and tool names aside); fastp itself is that reference. QC tables must equal FastQC-Rust's `fastqc_data.txt` without its version and file name lines, tied overrepresented rows in sorted order. The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>

z-fastq returned right outputs on all 0 inputs it was given. Every other tool was right on every input it was given.

> [!CAUTION]
> fastp, FastQC-Rust and FastQC print wrong output and exit 0 on some malformed inputs, so a script cannot tell the pairs are broken.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/edges-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/edges-light.png">
  <img src="results/figures/edges-light.png" alt="What each tool does on small edge-case inputs" width="100%">
</picture>

> [!NOTE]
> A cell is the tool's verdict: `accept` (exit 0, outputs right), `reject P001` (with the rule it named), `reject`, `wrong: accept` (exit 0 on a broken pair or record), `wrong: records differ` or `wrong: table R1` (exit 0, other output), `no table` (no table written for a mate), `refused` (rejected a valid input), `crash`, or `-` (the tool has no such job). "tool policy" rows have no single right answer.

## Copying
> <sub>Plain output: two mate files into two plain files. Gzip output: two mate files into two gzip files. The real runs and the synthetic pair at 1 GiB, plain and gzip input, in the page cache. z-fastq runs on one CPU; every other tool at 1 thread, and at 2, 4, 8, and 16 threads on as many physical cores, in the same batch.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/copying-overview-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/copying-overview-light.png">
  <img src="results/figures/copying-overview-light.png" alt="Copying: time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes, geometric means over the 1 GiB inputs: Plain output, 0.64 s and 2.3 MiB on plain input, 1.59 s and 2.3 MiB on gzip; Gzip output, 2.58 s and 2.8 MiB on plain input, 3.39 s and 2.8 MiB on gzip. The table gives every tool as a multiple of that.

| Job          | Tool       | Threads | Plain time | Plain CPU | Plain memory | Gzip time | Gzip CPU | Gzip memory |
| ------------ | ---------- | ------: | ---------: | --------: | -----------: | --------: | -------: | ----------: |
| Plain output | z-fastq    |       1 |      1.00x |     1.00x |        1.00x |     1.00x |    1.00x |       1.00x |
| Plain output | fastp      |       1 |     18.21x |    18.86x |      484.36x |     7.85x |    8.48x |     488.06x |
| Plain output | fastp -w 2 |       2 |     11.00x |    17.08x |      506.52x |     4.17x |    7.79x |     509.60x |
| Plain output | fastp -w 4 |       4 |      6.10x |    17.09x |      527.13x |     2.52x |    7.49x |     530.51x |
| Plain output | fastp -w 8 |       8 |      4.77x |    17.71x |      548.27x |     1.83x |    7.81x |     548.70x |
| Gzip output  | z-fastq    |       1 |      1.00x |     1.00x |        1.00x |     1.00x |    1.00x |       1.00x |
| Gzip output  | fastp      |       1 |      6.97x |     6.91x |      401.34x |     5.52x |    5.47x |     404.47x |
| Gzip output  | fastp -w 2 |       2 |      5.40x |     6.43x |      397.95x |     4.38x |    5.09x |     401.71x |
| Gzip output  | fastp -w 4 |       4 |      2.41x |     6.36x |      406.59x |     1.91x |    4.99x |     409.45x |
| Gzip output  | fastp -w 8 |       8 |      1.26x |     6.72x |      422.25x |     0.98x |    5.28x |     425.35x |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/copying-speed-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/copying-speed-light.png">
  <img src="results/figures/copying-speed-light.png" alt="Copying: wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/copying-memory-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/copying-memory-light.png">
  <img src="results/figures/copying-memory-light.png" alt="Copying: peak memory by kind of input" width="100%">
</picture>

<details><summary>Scaling</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/copying-scaling-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/copying-scaling-light.png">
  <img src="results/figures/copying-scaling-light.png" alt="Copying: throughput and peak memory by input size" width="100%">
</picture>

**Plain output.** At 10 MiB, z-fastq processes 0.9 to 0.9 GiB/s on plain input and 0.4 to 0.4 GiB/s on gzip input; no other tool is faster there.
z-fastq peaks at 2.31 MiB at every size and shape. No other tool grows by more than a quarter either.

**Gzip output.** At 10 MiB, z-fastq processes 0.2 to 0.2 GiB/s on plain input and 0.2 to 0.2 GiB/s on gzip input; no other tool is faster there.
z-fastq peaks at 2.81 MiB at every size and shape. No other tool grows by more than a quarter either.

</details>

## Trimming
> <sub>TruSeq: TruSeq adapters, quality and N filters, minimum length 15 (rnaFusion's settings). Nextera: Nextera adapters, poly-G and quality-tail trimming, filters, minimum length 30 (FRAMEFusion's settings). The real runs and the synthetic pair at 1 GiB, plain and gzip input, in the page cache. z-fastq runs on one CPU; every other tool at 1 thread, and at 2, 4, 8, and 16 threads on as many physical cores, in the same batch.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/trimming-overview-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/trimming-overview-light.png">
  <img src="results/figures/trimming-overview-light.png" alt="Trimming: time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes, geometric means over the 1 GiB inputs: TruSeq, 3.98 s and 38.5 MiB on plain input, 4.92 s and 38.4 MiB on gzip; Nextera, 4.04 s and 38.5 MiB on plain input, 4.76 s and 38.4 MiB on gzip. The table gives every tool as a multiple of that.

| Job     | Tool       | Threads | Plain time | Plain CPU | Plain memory | Gzip time | Gzip CPU | Gzip memory |
| ------- | ---------- | ------: | ---------: | --------: | -----------: | --------: | -------: | ----------: |
| TruSeq  | z-fastq    |       1 |      1.00x |     1.00x |        1.00x |     1.00x |    1.00x |       1.00x |
| TruSeq  | fastp      |       1 |      9.78x |     9.68x |       29.07x |     8.12x |    8.12x |      29.33x |
| TruSeq  | fastp -w 2 |       2 |      5.29x |     9.05x |       30.10x |     4.04x |    7.63x |      30.51x |
| TruSeq  | fastp -w 4 |       4 |      2.63x |     8.57x |       31.11x |     2.26x |    7.32x |      31.48x |
| TruSeq  | fastp -w 8 |       8 |      1.32x |     8.71x |       33.13x |     1.11x |    7.40x |      33.35x |
| Nextera | z-fastq    |       1 |      1.00x |     1.00x |        1.00x |     1.00x |    1.00x |       1.00x |
| Nextera | fastp      |       1 |      8.55x |     8.06x |       29.08x |     6.99x |    6.95x |      29.33x |
| Nextera | fastp -w 2 |       2 |      4.70x |     7.61x |       30.25x |     4.03x |    6.61x |      30.55x |
| Nextera | fastp -w 4 |       4 |      2.20x |     7.18x |       31.44x |     1.95x |    6.18x |      31.43x |
| Nextera | fastp -w 8 |       8 |      1.17x |     7.44x |       33.13x |     1.00x |    6.39x |      33.37x |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/trimming-speed-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/trimming-speed-light.png">
  <img src="results/figures/trimming-speed-light.png" alt="Trimming: wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/trimming-memory-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/trimming-memory-light.png">
  <img src="results/figures/trimming-memory-light.png" alt="Trimming: peak memory by kind of input" width="100%">
</picture>

**TruSeq.** Output per plain 1 GiB input, mean over 1 inputs: fastp 1067.5 MB.

**Nextera.** Output per plain 1 GiB input, mean over 1 inputs: fastp 1065.6 MB.

<details><summary>Scaling</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/trimming-scaling-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/trimming-scaling-light.png">
  <img src="results/figures/trimming-scaling-light.png" alt="Trimming: throughput and peak memory by input size" width="100%">
</picture>

**TruSeq.** At 10 MiB, z-fastq processes 0.1 to 0.1 GiB/s on plain input and 0.1 to 0.1 GiB/s on gzip input; no other tool is faster there.
z-fastq peaks at 5.06 to 8.31 MiB across sizes and shapes. No other tool grows by more than a quarter either.

**Nextera.** At 10 MiB, z-fastq processes 0.2 to 0.2 GiB/s on plain input and 0.1 to 0.1 GiB/s on gzip input; no other tool is faster there.
z-fastq peaks at 5.06 to 8.31 MiB across sizes and shapes. No other tool grows by more than a quarter either.

</details>

## QC tables
> <sub>FastQC tables: FastQC tables for both mates, no trimming. The real runs and the synthetic pair at 1 GiB, plain and gzip input, in the page cache. z-fastq runs on one CPU; every other tool at 1 thread, and at 2, 4, 8, and 16 threads on as many physical cores, in the same batch.</sub>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/qc-overview-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/qc-overview-light.png">
  <img src="results/figures/qc-overview-light.png" alt="QC tables: time and peak RSS relative to z-fastq" width="100%">
</picture>

z-fastq itself takes, geometric means over the 1 GiB inputs: FastQC tables, 2.00 s and 14.1 MiB on plain input, 2.97 s and 14.0 MiB on gzip. The table gives every tool as a multiple of that.

| Job           | Tool             | Threads | Plain time | Plain CPU | Plain memory | Gzip time | Gzip CPU | Gzip memory |
| ------------- | ---------------- | ------: | ---------: | --------: | -----------: | --------: | -------: | ----------: |
| FastQC tables | z-fastq          |       1 |      1.00x |     1.00x |        1.00x |     1.00x |    1.00x |       1.00x |
| FastQC tables | FastQC-Rust      |       1 |      2.01x |     2.58x |        3.06x |     1.65x |    1.95x |       3.17x |
| FastQC tables | FastQC-Rust -t 2 |       2 |      1.02x |     2.52x |        5.16x |     0.89x |    2.02x |       5.48x |
| FastQC tables | FastQC-Rust -t 4 |       4 |      0.84x |     2.77x |        6.94x |     0.70x |    2.04x |       7.33x |
| FastQC tables | FastQC-Rust -t 8 |       8 |      0.73x |     3.92x |        6.76x |     0.51x |    2.69x |       9.72x |
| FastQC tables | FastQC           |       1 |     10.70x |    13.17x |       17.45x |     9.18x |   10.09x |      17.32x |
| FastQC tables | FastQC -t 2      |       2 |      9.68x |    13.66x |       22.12x |     7.48x |   10.20x |      20.86x |
| FastQC tables | FastQC -t 4      |       4 |      4.62x |    13.05x |       28.91x |     3.43x |    9.76x |      28.58x |
| FastQC tables | FastQC -t 8      |       8 |      2.65x |    13.47x |       56.91x |     1.78x |    9.93x |      51.10x |

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/qc-speed-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/qc-speed-light.png">
  <img src="results/figures/qc-speed-light.png" alt="QC tables: wall time by kind of input" width="100%">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/qc-memory-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/qc-memory-light.png">
  <img src="results/figures/qc-memory-light.png" alt="QC tables: peak memory by kind of input" width="100%">
</picture>

**FastQC tables.** Output per plain 1 GiB input, mean over 1 inputs: FastQC-Rust 0.9 MB, FastQC 0.8 MB.

<details><summary>Scaling</summary>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/qc-scaling-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/qc-scaling-light.png">
  <img src="results/figures/qc-scaling-light.png" alt="QC tables: throughput and peak memory by input size" width="100%">
</picture>

**FastQC tables.** At 10 MiB, z-fastq processes 0.3 to 0.3 GiB/s on plain input and 0.0 to 0.0 GiB/s on gzip input. Faster on a single shape there, as throughput over z-fastq's: FastQC-Rust on gzip paired (1.09x).
z-fastq peaks at 2.81 to 9.33 MiB across sizes and shapes. No other tool grows by more than a quarter either.

</details>

## Reproducibility

<details><summary>How it was measured</summary>

Every tool runs through `zexec -r` inside one zebrac batch per case, held to one CPU core by a systemd scope with a CPU quota and a 4G memory cap; a batch with threaded variants holds as many physical cores as its largest lane, one logical CPU each, with each lane pinned to its own cores and z-fastq to one, after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. The kernel reports peak memory in 256 KiB steps, so values closer than that are equal. Runs per tier, each after one warmup: tiny 31, xsmall 21, small 11, medium 3; sizes 1 MiB to 1 GiB. Generated by `bench/process/run.sh release` on 2026-10-28.

</details>

<details><summary>Identity</summary>

- `z-fastq`: *0.1.0 at b199bc6*
- `zebrac`: *0.6.2*
- `zexec`: *15000 bytes, static, no libc*
- `kernel`: *7.2.11-200.fc44.x86_64*
- `cpu`: *AMD Ryzen 9 3950X 16-Core Processor*
- `load at batch start`: *1.66 to 12.68*

</details>

<details><summary>Commands</summary>

- `z-fastq`: *zig-out/bin/z-fastq process --out1 {stem}_1.fq --out2 {stem}_2.fq IN*
- `fastp`: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 1 -j {stem}.json -h {stem}.html*
- `fastp -w 2`: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 2 -j {stem}.json -h {stem}.html*
- `fastp -w 4`: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 4 -j {stem}.json -h {stem}.html*
- `fastp -w 8`: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 8 -j {stem}.json -h {stem}.html*
- `fastp -w 16`: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 16 -j {stem}.json -h {stem}.html*
- `FastQC-Rust`: *tools/bin/fastqc-rust --threads 1 -o {stem}.qc IN*
- `FastQC-Rust -t 2`: *tools/bin/fastqc-rust --threads 2 -o {stem}.qc IN*
- `FastQC-Rust -t 4`: *tools/bin/fastqc-rust --threads 4 -o {stem}.qc IN*
- `FastQC-Rust -t 8`: *tools/bin/fastqc-rust --threads 8 -o {stem}.qc IN*
- `FastQC-Rust -t 16`: *tools/bin/fastqc-rust --threads 16 -o {stem}.qc IN*
- `FastQC`: *tools/bin/fastqc --threads 1 --quiet --outdir {stem}.qc IN*
- `FastQC -t 2`: *tools/bin/fastqc --threads 2 --quiet --outdir {stem}.qc IN*
- `FastQC -t 4`: *tools/bin/fastqc --threads 4 --quiet --outdir {stem}.qc IN*
- `FastQC -t 8`: *tools/bin/fastqc --threads 8 --quiet --outdir {stem}.qc IN*
- `FastQC -t 16`: *tools/bin/fastqc --threads 16 --quiet --outdir {stem}.qc IN*
- `z-fastq`, copy-gzip: *zig-out/bin/z-fastq process --gzip --out1 {stem}_1.fq.gz --out2 {stem}_2.fq.gz IN*
- `fastp`, copy-gzip: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 1 -j {stem}.json -h {stem}.html*
- `fastp -w 2`, copy-gzip: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 2 -j {stem}.json -h {stem}.html*
- `fastp -w 4`, copy-gzip: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 4 -j {stem}.json -h {stem}.html*
- `fastp -w 8`, copy-gzip: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 8 -j {stem}.json -h {stem}.html*
- `fastp -w 16`, copy-gzip: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 16 -j {stem}.json -h {stem}.html*
- `z-fastq`, trim-truseq: *zig-out/bin/z-fastq process --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --low-quality 15 --max-low-quality-percent 40 --max-n 5 --min-length 15 --pair-duplication --qc-summary {stem}.json --out1 {stem}_1.fq --out2 {stem}_2.fq IN*
- `fastp`, trim-truseq: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_trim_poly_g -q 15 -u 40 -n 5 -l 15 -w 1 -j {stem}.json -h {stem}.html*
- `fastp -w 2`, trim-truseq: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_trim_poly_g -q 15 -u 40 -n 5 -l 15 -w 2 -j {stem}.json -h {stem}.html*
- `fastp -w 4`, trim-truseq: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_trim_poly_g -q 15 -u 40 -n 5 -l 15 -w 4 -j {stem}.json -h {stem}.html*
- `fastp -w 8`, trim-truseq: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_trim_poly_g -q 15 -u 40 -n 5 -l 15 -w 8 -j {stem}.json -h {stem}.html*
- `fastp -w 16`, trim-truseq: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_trim_poly_g -q 15 -u 40 -n 5 -l 15 -w 16 -j {stem}.json -h {stem}.html*
- `z-fastq`, trim-nextera: *zig-out/bin/z-fastq process --adapter1 CTGTCTCTTATACACATCT --adapter2 CTGTCTCTTATACACATCT --trim-poly-g 10 --trim-tail-window 4 --trim-tail-quality 20 --low-quality 15 --max-low-quality-percent 40 --max-n 5 --min-length 30 --pair-duplication --qc-summary {stem}.json --out1 {stem}_1.fq --out2 {stem}_2.fq IN*
- `fastp`, trim-nextera: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -q 15 -u 40 -n 5 -l 30 -w 1 -j {stem}.json -h {stem}.html*
- `fastp -w 2`, trim-nextera: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -q 15 -u 40 -n 5 -l 30 -w 2 -j {stem}.json -h {stem}.html*
- `fastp -w 4`, trim-nextera: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -q 15 -u 40 -n 5 -l 30 -w 4 -j {stem}.json -h {stem}.html*
- `fastp -w 8`, trim-nextera: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -q 15 -u 40 -n 5 -l 30 -w 8 -j {stem}.json -h {stem}.html*
- `fastp -w 16`, trim-nextera: *tools/bin/fastp -i IN -I IN -o {stem}_1.fq -O {stem}_2.fq --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -q 15 -u 40 -n 5 -l 30 -w 16 -j {stem}.json -h {stem}.html*
- `z-fastq`, qc-tables: *zig-out/bin/z-fastq process --qc-tables {stem}.qc --contaminants tools/.local/fastqc-current/runtime/Configuration/contaminant_list.txt --out1 {stem}_1.fq --out2 {stem}_2.fq IN*

</details>

<details><summary>Quick-check cases for peer and ab</summary>

Picked from these results: per family, group, and format, the case where z-fastq is closest to, or furthest behind, its fastest rival of any class. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.

- `copy-plain: Short, fixed length, plain`: *NovaSeq, 150 bp*; rival fastp -w 8 at 4.77x z-fastq's time
- `copy-plain: Short, fixed length, gzip`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.83x z-fastq's time
- `copy-gzip: Short, fixed length, plain`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.26x z-fastq's time
- `copy-gzip: Short, fixed length, gzip`: *NovaSeq, 150 bp*; rival fastp -w 8 at 0.98x z-fastq's time
- `trim-truseq: Short, fixed length, plain`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.31x z-fastq's time
- `trim-truseq: Short, fixed length, gzip`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.11x z-fastq's time
- `trim-nextera: Short, fixed length, plain`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.17x z-fastq's time
- `trim-nextera: Short, fixed length, gzip`: *NovaSeq, 150 bp*; rival fastp -w 8 at 1.00x z-fastq's time
- `qc-tables: Short, fixed length, plain`: *NovaSeq, 150 bp*; rival FastQC-Rust -t 8 at 0.72x z-fastq's time
- `qc-tables: Short, fixed length, gzip`: *NovaSeq, 150 bp*; rival FastQC-Rust -t 8 at 0.51x z-fastq's time

</details>

<details><summary>Hardware counters (user space, plain medium-tier inputs, geometric means)</summary>

| Tool       | Instructions per byte | Instructions per cycle | Cache misses per MiB | Branch misses per MiB |
| ---------- | --------------------: | ---------------------: | -------------------: | --------------------: |
| z-fastq    |                  2.97 |                   2.98 |                 4470 |                   693 |
| fastp      |                 78.71 |                   2.70 |                85014 |                 61274 |
| fastp -w 2 |                 68.19 |                   2.42 |               119004 |                 56853 |
| fastp -w 4 |                 62.97 |                   2.11 |               150697 |                 54668 |
| fastp -w 8 |                 60.40 |                   1.92 |               165647 |                 54194 |

</details>

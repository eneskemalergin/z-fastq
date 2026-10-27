# Processing paired reads

`process` reads two mate files in lockstep, validates every pair as `check --paired` does, and writes each accepted pair to two new files, plain or gzip. With trim options it removes adapters, mate overlaps, 3' poly-G runs, and low-quality 3' ends, and with filter options it drops pairs by quality, `N` count, and length, all by fastp 1.1.0's rules; without them it changes and drops no read. With `--report` it also writes a JSON report with per-mate statistics.

```bash
z-fastq process --gzip --out1 clean_R1.fastq.gz --out2 clean_R2.fastq.gz --report sample.json R1.fastq.gz R2.fastq.gz
```

## Outputs

- `--out1` and `--out2` are created as new files. An existing path, including one of the inputs, fails with `output path already exists` (status `3`) before anything is written, and outputs created by the same run are removed.
- Records are written with LF line endings, in input order, mate 1 to `--out1` and mate 2 to `--out2`.
- `--gzip` compresses both outputs. A name ending in `.gz` without `--gzip` is rejected (status `2`), so a plain file never carries a gzip name. The gzip header holds no time stamp, so the same inputs and options give the same bytes.
- Standard input may be one of the two inputs; outputs cannot be standard output.

## Stop after N pairs

```bash
z-fastq process --max-pairs 1000000 --out1 sub_R1.fastq --out2 sub_R2.fastq R1.fastq.gz R2.fastq.gz
```

`--max-pairs N` writes the first N complete pairs in input order and stops reading. Nothing after pair N is checked, so a missing mate or a damaged gzip trailer later in the files goes unseen; the report says `"input_complete": false`. Use `check --paired` to validate whole files.

## Trimming

Each trim is off unless its options are given.

```bash
z-fastq process --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --out1 trimmed_R1.fastq --out2 trimmed_R2.fastq R1.fastq.gz R2.fastq.gz
z-fastq process --adapter1 CTGTCTCTTATACACATCT --adapter2 CTGTCTCTTATACACATCT --trim-poly-g 10 \
  --trim-tail-window 4 --trim-tail-quality 20 --out1 trimmed_R1.fastq --out2 trimmed_R2.fastq R1.fastq.gz R2.fastq.gz
```

The first command writes the records fastp 1.1.0 writes with `--adapter_sequence A1 --adapter_sequence_r2 A2 --disable_trim_poly_g -Q -L`, and the second those it writes with `--adapter_sequence A --adapter_sequence_r2 A --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 -Q -L`. fastp 1.4.0 gives the same result except on mate 2 reads that hold the IUPAC letters `S`, `W`, or `D`, where its result depends on the CPU.

The steps run in this order on each pair:

1. Quality tail (`--trim-tail-window W --trim-tail-quality Q`): windows of `W` bases are checked from the 3' end toward the 5' end, up to the first window whose mean Phred quality is at least `Q`. If that is the window at the 3' end, nothing is cut; otherwise the read keeps everything up to and including that window's first base. The window that starts at the first base is never checked, so a read with no passing window keeps one base. Trailing `N` bases are then removed. A read of `W` bases or fewer is emptied. When either mate is empty after this step, the pair is not trimmed further.
2. Poly-G (`--trim-poly-g N`): bases are scanned from the 3' end. The scan stops when more than 5 bases other than `G` have been seen, or, from the `N`th base scanned on, when they exceed one per 8 bases scanned. If the scan covered at least `N` bases, the read is cut before the leftmost `G` it saw. A run shorter than `N` can be cut when the scan carries on past it: with `N` 10, a read ending in 9 `G` after other bases loses the 9 `G`.
3. Adapters (`--adapter1 A1 --adapter2 A2`): the reverse complement of mate 2 is compared with mate 1 at successive offsets: first with mate 2 starting at or after the start of mate 1, while more than 30 bases of mate 1 remain, then with mate 2 starting 0, 1, 2, and more bases before mate 1, while more than 30 bases of mate 2 remain. The first offset with at most 5 mismatches, and at most one per 5 bases of overlap, in its first 50 positions wins; `N` counts as a mismatch. When mate 2 starts before mate 1 there (the insert is shorter than the reads), both mates are cut to the overlap length. Otherwise each mate is searched from the left for its own adapter: a full copy, or a partial copy of at least 5 bases at the 3' end, with one mismatch allowed per 8 bases compared. The search starts up to 4 bases before the read, so a read that begins with the end of an adapter is emptied. If no copy is found, a copy with one inserted or deleted base is tried. As in fastp, that test compares the adapter with the start of the read, so it finds such a copy only at the read start; when it passes only on a shorter length near the 3' end, the read is cut where that length starts, whether or not an adapter is there.

- Every trim removes a 3' end, cutting sequence and quality together; the header and plus lines are written unchanged.
- Bases compare as written: a lowercase base never matches an adapter base or a poly-G `G`.
- Adapters are 4 to 128 bases of uppercase `A`, `C`, `G`, and `T`, and are given together. The window is 1 to 1000 bases and the quality 1 to 30, given together; the poly-G minimum is at least 1.
- A pair whose mate is trimmed to nothing is dropped, and so is an adapter dimer; see [Filtering](#filtering).

## Filtering

Each filter is off unless its option is given. A pair is written only when both mates pass.

```bash
z-fastq process --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --low-quality 15 --max-low-quality-percent 40 --max-n 5 --min-length 15 \
  --out1 clean_R1.fastq --out2 clean_R2.fastq R1.fastq.gz R2.fastq.gz
```

This writes the records fastp 1.1.0 writes with `--adapter_sequence A1 --adapter_sequence_r2 A2 --disable_trim_poly_g` and its default filters (`-q 15 -u 40 -n 5 -l 15`). fastp turns those filters on by default; z-fastq applies only the filters given.

Each read is checked after trimming, in this order; the first check it fails is its reason:

1. A read with no bases left: `too_short`.
2. `--low-quality Q --max-low-quality-percent P`: more than `P` percent of its bases have a Phred quality below `Q`: `low_quality`. With 151 bases and 40 percent, 60 low-quality bases pass and 61 fail.
3. `--max-n N`: more than `N` uppercase `N` bases: `too_many_n`. A lowercase `n` is not counted, as in fastp.
4. `--min-length N`: fewer than `N` bases: `too_short`.

- With adapters given, a pair whose overlap or adapter step cut either mate and left both mates at 2 bases or fewer is an adapter dimer: both mates fail as `adapter_dimer`.
- The empty-read and adapter-dimer drops apply whenever any trim or filter option is given, as in fastp, where `-Q -L` leaves them on.
- A dropped pair is dropped whole; neither mate is written anywhere.
- The pair's reason is the higher of its mates' reasons in this order: `adapter_dimer`, `low_quality`, `too_short`, `too_many_n`. fastp counts both reads of a dropped pair under that reason, so its `*_reads` counts are twice the pair counts in the report.
- Dropping pairs never changes the exit status. A malformed record still fails the run, also after dropped pairs, because every pair is validated before it is trimmed or filtered.
- Ranges, as in fastp: `--min-length` 0 to 4294967295, `--max-n` 0 to 50, `--low-quality` 0 to 93, `--max-low-quality-percent` 0 to 100; the last two are given together.

## Duplication

The report always carries FastQC's Sequence Duplication Levels and Overrepresented sequences for each mate, raw and retained. Two options add to them, and both require `--report` (status `2` otherwise):

- `--contaminants PATH` reads a contaminant list in FastQC's format (`name`, one or more tabs, `sequence`; `#` lines and blank lines are skipped) and names each overrepresented sequence's best match, by FastQC's rules. Pass the same `contaminant_list.txt` you give FastQC. A list over 1 MiB fails with status `4`, a sequence with a base other than A, C, G, or T fails with status `2`, and an unreadable file with status `3`, before any output is created. z-fastq does not ship a list.
- `--pair-duplication` adds fastp's duplication rate for read pairs, the `duplication.rate` of fastp's JSON. fastp keeps two bit arrays of 512 MiB each; z-fastq keeps only the bits the pairs set, about 10 bytes per pair, and never more than fastp's 1 GiB. It is off by default.

## The report

`--report PATH` writes a `z-fastq/process-v1` document after both outputs are finished; `--report -` writes it to standard output. The statistics and QC modules are computed only when a report is requested, so a run without `--report` only validates, trims, filters, and writes. The report carries:

- `options`: the effective options.
- `inputs`: each input's path, and its size and modification time in nanoseconds when it is a regular file (`null` for a pipe).
- `outputs`: the output paths.
- `status`: `ok` or `error`; `input_complete`: whether both inputs were read to the end.
- `pairs`: pairs `read`, `written`, and `dropped`. On a complete run, `read` is `written` plus `dropped`.
- `mates`: one statistics object per mate, with the same fields as `stats --json`, under `raw` for the reads as they were read and `retained` for the reads as written. Without a trim or filter option the two are the same.
- `trimming`: one object per mate with an entry per step, `quality_tail`, `poly_g`, `overlap`, and `adapter`: the reads that step shortened and the bases it removed, or `null` when the step is off. It counts every pair read, including pairs dropped later.
- `filtering`: `dropped_pairs`, the dropped pairs by the pair's reason, and `failed_reads`, one object per mate with the reads that failed by their own reason (a passing mate of a dropped pair is not counted there); each holds `too_many_n`, `too_short`, `low_quality`, and `adapter_dimer`. `null` when no trim or filter option is given.
- `inserts`: insert sizes from the mate overlap, as fastp counts them, for every pair whose mates both reach the overlap step (neither emptied by the quality tail cut), before any filter. `pairs` counts them; `overlapping` counts those with an overlap found, so `overlapping` over `pairs` is the share of overlapping pairs; `unknown` counts no overlap or an insert of 512 or more; `peak` is the most frequent size (the lowest on a tie, 0 when none); `median` is the smallest size reached by half of the known sizes (`null` when none); `histogram` has 512 counts for sizes 0 to 511. An overlap after mate 1's start gives `len1 + len2 - overlap`, otherwise the overlap length, so a mate shorter than the fragment gives fastp's shorter value. `null` without `--adapter1`/`--adapter2`, because only they run the overlap step. fastp counts only every W-th batch of pairs with W workers; these counts equal fastp with `-w 1`.
- `pair_duplication`: with `--pair-duplication`, the `pairs` read, the `duplicates` among them, and the `rate`, as fastp computes them: every pair as read, before trimming and filtering, hashed over all of mate 1 and then mate 2. A pair counts as a duplicate when every bit it maps to was already set, so a hash collision counts too, as in fastp. The values equal fastp with `-w 1`; fastp with more threads can differ by a few pairs. `null` without the option.
- `failed_input`, `failed_output`, and `error`: the failure, in the form `check --paired --json` uses, or `null`.
- `qc`: FastQC's core modules for each mate, under `raw` for the reads as they were read and `retained` for the reads as written: `per_base_quality`, `per_base_content`, `per_base_n`, `length_distribution`, `per_sequence_quality`, `per_sequence_gc` (101 values, percent 0 to 100), `terminal_runs` (reads ending in at least 10 `G` or 10 `A`), `adapter_content`, `duplication`, and `overrepresented`. The values follow FastQC 0.13.0's definitions, including its groups of cycles for long reads; a value FastQC prints as `NaN` is `null`. A mate with a read longer than 1024 bases reports its per-cycle modules, length distribution, and adapter content as `null`, with the reason in `unavailable`; duplication and overrepresented sequences stay available, since they look only at the first 50 bases.
  - `adapter_content`: FastQC's bundled adapter list (Illumina Universal, Illumina Small RNA 3' and 5', Nextera Transposase, PolyA, PolyG, each 12 bases) under `adapters`, and `rows`, one per group of positions, each with its `position` label and one percentage per adapter: the share of all reads whose first exact match of that adapter (case-insensitive) starts at or before the position, and for a group the mean over its positions. Positions run to the longest read minus 11; `rows` is empty when no read is longer than 12 bases. A 12-base read that matches counts even before any longer read, as in FastQC-Rust; Java FastQC 0.13.0 skips it.
  - `duplication`: FastQC's Sequence Duplication Levels. Each read's key is its first 50 bases in uppercase; the first 100,000 distinct keys are stored, and a later read with a new key counts only in `reads`. `stored_sequences` and `reads_at_limit` (the reads counted when the last key was stored) are FastQC's inputs to its correction for unstored reads. `total_deduplicated_percentage` and 16 `levels` (`1` to `9`, `>10`, `>50`, `>100`, `>500`, `>1k`, `>5k`, `>10k+`), each the percentage of all reads in keys at that level, equal FastQC's values; with no reads the levels are `null`.
  - `overrepresented`: the stored keys with more than 0.1% of all reads, most reads first, ties in the order the keys first appeared (FastQC orders ties by its hash table). Each row has the `sequence`, its `count`, its unrounded `percentage`, and `possible_source`: the best contaminant match as `name`, `identity` percent, and `length` in bases, or `null` without `--contaminants` or without a match.

The report holds no time of the run or host name, so a rerun gives the same bytes.

## Failure and completion

A failure stops the run. Pairs already written stay in the outputs; if mate 2 cannot be written after mate 1, `--out1` holds one record more than `--out2`. The report, when requested, is still written with `"status": "error"`.

A run is complete only when the status is `0`, and, when a report was requested, the report says `"status": "ok"`. Treat any other result as a failed sample. Exit statuses are the usual ones: `1` validation, `2` usage, `3` I/O or allocation, `4` limits.

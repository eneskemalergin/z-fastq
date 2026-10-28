# Report reference

`process` writes two JSON documents: the **report** (`--report`), which holds everything about a run, and the **QC summary** (`--qc-summary`), a smaller file shaped for MultiQC. This page describes every field of both. For what the measurements mean, read [Quality control](Quality-Control); this page is for looking a field up.

## Shared conventions

Both documents begin with the same three fields:

- `schema` names the document's shape: `z-fastq/process-v1` or `z-fastq/qc-summary-v1`.
- `tool` gives the `name` and `version` of the z-fastq build that wrote it.
- `byte_strings` names how paths and the command line are written, `escaped-bytes-v1`: a byte outside printable ASCII becomes `\xHH` and a backslash becomes `\\`, before the usual JSON escaping, so any file name survives the trip.

A value that cannot be computed, such as a mean over zero reads, is `null`. Version 0.1.0 froze the `stats` and `check` documents ([Automation](Automation#json)); these two are new in 0.2.0 and may still change before it is released.

## Report

The report records what was asked, what was read, what was done, and what was measured. It holds no time of day or host name, so the same run gives the same bytes.

<details>
<summary>An example report, shortened</summary>

```json
{
  "schema": "z-fastq/process-v1",
  "tool": { "name": "z-fastq", "version": "0.1.0" },
  "byte_strings": "escaped-bytes-v1",
  "options": { "gzip": false, "max_pairs": null, "adapter1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "min_length": 15, "pair_duplication": true },
  "inputs": [ { "path": "R1.fastq.gz", "size": 183592, "mtime_ns": 1792773758730546172 } ],
  "status": "ok",
  "input_complete": true,
  "pairs": { "read": 1662, "written": 1662, "dropped": 0 },
  "mates": { "raw": [ { "reads": 1662, "bases": 240771, "mean_length": 144.8682, "gc_fraction": 0.5051, "q30_fraction": 0.8906 } ] },
  "trimming": [ { "quality_tail": null, "poly_g": null, "overlap": { "reads": 0, "bases": 0 }, "adapter": { "reads": 4, "bases": 22 } } ],
  "filtering": { "dropped_pairs": { "too_many_n": 0, "too_short": 0, "low_quality": 0, "adapter_dimer": 0 } },
  "inserts": { "pairs": 1662, "overlapping": 616, "unknown": 1046, "peak": 200, "median": 195 },
  "pair_duplication": { "pairs": 1662, "duplicates": 0, "rate": 0 },
  "qc": { "raw": [ { "per_base_quality": [ { "base": "1", "mean": 29.5156, "median": 30 } ] } ] },
  "error": null
}
```

Lists and objects are cut to their first entries; a real report holds both mates, both populations, and every module.

</details>

### The run

- **`options`**: every option of the run, as given or defaulted, under its name without the leading dashes and with `_` for `-`, so `--min-length` becomes `min_length`. An option not given is `null` or `false`. The [cheat sheet](Command-Cheat-Sheet#processing) lists them.
- **`inputs`**: each input's `path`, with its `size` and modification time `mtime_ns` in nanoseconds, or `null` for a pipe.
- **`outputs`**: each output's `path`.
- **`status`**: `ok` or `error`, and **`input_complete`**: whether both inputs were read to the end.
- **`pairs`**: the pairs `read`, `written`, and `dropped`. On a complete run, `read` equals `written` plus `dropped`.

### The reads

- **`mates`**: per-mate statistics with the fields of `stats --json` ([Counting and statistics](Counting-and-Statistics#stats)), under `raw` for the reads as they were read and `retained` for the reads as written, each a list of mate 1 then mate 2. Without a trim or filter option the two are the same.

### What trimming and filtering did

- **`trimming`**: one object per mate, with one entry per step: `quality_tail`, `poly_g`, `overlap`, and `adapter`. Each gives the `reads` that step shortened and the `bases` it removed, or `null` when the step is off. Pairs dropped later are counted too.
- **`filtering`**: `dropped_pairs`, counted by the pair's reason, and `failed_reads`, one object per mate counted by each read's own reason; a passing mate of a dropped pair is not counted there. Both use `too_many_n`, `too_short`, `low_quality`, and `adapter_dimer`. `null` without a trim or filter option.
- **`inserts`**, with adapters only: `pairs` that reached the overlap step, `overlapping` pairs with an overlap found, `unknown` for no overlap or an insert of 512 or more, the most frequent size `peak` (the lowest on a tie, 0 when none), the `median` known size, and a `histogram` of 512 counts for sizes 0 to 511.
- **`pair_duplication`**, with `--pair-duplication` only: the `pairs` read, the `duplicates` among them, and their ratio, `rate`.

### QC modules

**`qc`** has `raw` and `retained`, each a list of mate 1 then mate 2. Each mate holds one entry per module; cycle positions and groups, such as `10-14`, count from 1.

| Module | Each entry holds |
| --- | --- |
| `per_base_quality` | `base`, `mean`, `median`, quartiles, 10th and 90th percentiles |
| `per_base_content` | `base`, and `g`, `a`, `t`, `c` in percent |
| `per_base_n` | `base`, and `n` in percent |
| `length_distribution` | `length`, `count` |
| `per_sequence_quality` | `quality`, `count` |
| `per_sequence_gc` | 101 counts, for GC 0 to 100% |
| `terminal_runs` | `min_length` (10), `g`, `a` |
| `adapter_content` | `adapters`, and `rows` of `position` and `values` |
| `duplication` | totals and 16 `levels` |
| `overrepresented` | `sequence`, `count`, `percentage`, `possible_source` |
| `per_tile_quality` | `positions`, and `tiles` of `tile` and `deviations` |

A few need more than a field name:

- **Percentiles** in `per_base_quality` (`lower_quartile`, `upper_quartile`, `percentile_10`, `percentile_90`, and `median`) need more than 100 reads at a cycle, and are `null` otherwise.
- **`duplication`** holds `reads`, `stored_sequences`, `reads_at_limit` (the reads counted when the last key was stored), `total_deduplicated_percentage`, and `levels` from `1` to `>10k+`, each with its `percentage`.
- **`possible_source`** is `null`, or the best contaminant match as `name`, `identity` in percent, and `length` in bases.
- **`per_tile_unavailable`** and **`unavailable`** give the reason when per-tile quality, or every per-cycle module because of a read longer than 1024 bases, is `null`.

### Files and failures

- **`qc_table_files`**, with `--qc-tables` only: one entry per QC table, raw mate 1 and 2, then retained mate 1 and 2 when reads change. Each has `population`, `mate`, `path`, `written`, and `unavailable`, the reason a table was not written.
- **`error`**: `null` on success, or the failure's `code`, `message`, `record_index`, `byte_offset`, and `line_in_record`, as in `check --paired --json`. **`failed_input`** and **`failed_output`** name the file involved.

## QC summary

The QC summary keeps the field names of fastp's JSON report, so tools that read that report, MultiQC among them, can read this file. Its `schema` is `z-fastq/qc-summary-v1`. It is indented with tabs, which MultiQC needs to recognize it.

<details>
<summary>An example QC summary, shortened</summary>

```json
{
  "schema": "z-fastq/qc-summary-v1",
  "summary": {
    "before_filtering": { "total_reads": 3324, "total_bases": 481385, "q20_rate": 0.9085, "q30_rate": 0.8297, "read1_mean_length": 144, "read2_mean_length": 144, "gc_content": 0.5063 }
  },
  "filtering_result": { "passed_filter_reads": 3324, "low_quality_reads": 0, "too_many_N_reads": 0, "adapter_dimer_reads": 0, "too_short_reads": 0, "too_long_reads": 0 },
  "duplication": { "rate": 0 }
}
```

</details>

| Section | Holds |
| --- | --- |
| `summary` | totals before and after filtering: reads, bases, Q20 and Q30, mean lengths, GC |
| `filtering_result` | reads that passed, and reads dropped by reason |
| `duplication` | the pair duplication rate |
| `insert_size` | the insert peak, unknown inserts, and histogram (adapters only) |
| `adapter_cutting` | adapter-trimmed reads and bases, and the adapters (adapters only) |
| `read1_before_filtering` and three more | per-mate totals, and quality and base content per cycle |
| `command` | the z-fastq command line |

Everything named `before_filtering` describes the raw reads, and everything named `after_filtering` the reads written. Totals cover both mates; dropped pairs count both their reads, so the filter counts are twice the report's pair counts. An overlap trim counts both mates as adapter-trimmed, and a match that starts before the read counts the adapter bases it matched, as fastp counts them.

<details>
<summary>Every field name</summary>

- `summary.before_filtering`, `summary.after_filtering`: `total_reads`, `total_bases`, `q20_bases`, `q30_bases`, `q20_rate`, `q30_rate`, `read1_mean_length`, `read2_mean_length` (whole bases), `gc_content` (G and C over all bases).
- `filtering_result`: `passed_filter_reads`, `low_quality_reads`, `too_many_N_reads`, `adapter_dimer_reads` (with adapters), `too_short_reads`, `too_long_reads` (always 0).
- `duplication`: `rate`.
- `insert_size`: `peak`, `unknown`, `histogram` (512 counts).
- `adapter_cutting`: `adapter_trimmed_reads`, `adapter_trimmed_bases`, `read1_adapter_sequence`, `read2_adapter_sequence`.
- `read1_before_filtering`, `read2_before_filtering`, `read1_after_filtering`, `read2_after_filtering`: `total_reads`, `total_bases`, `q20_bases`, `q30_bases`, `q40_bases`, `total_cycles`, `quality_curves.mean`, and `content_curves` with `A`, `T`, `C`, `G`, `N`, `GC`. A mate with a read longer than 1024 bases has only the four totals.
- `command`.

</details>

fastp's own report holds more: per-base `A`, `T`, `C`, and `G` quality curves, `kmer_count`, `overrepresented_sequences`, the adapter sequence counts, `polyx_trimming`, `fastp_version`, and `sequencing`. This file leaves those out, and writes `insert_size` only when adapters are given, where fastp writes it for every paired run.

## Related pages

- [Quality control](Quality-Control): what each measurement means.
- [Processing paired reads](Processing): the trims and filters these counts describe.
- [Automation](Automation#json): reading JSON safely in scripts.

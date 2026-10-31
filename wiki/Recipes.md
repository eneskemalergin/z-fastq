# Recipes

Short patterns for everyday jobs. Each keeps the exit status and output behavior described on the linked pages.

## Count several lanes

```bash
z-fastq count lane1.fastq.gz lane2.fastq.gz lane3.fastq.gz
```

One count per input, in order.

## Stop when the pairs are broken

```bash
z-fastq check --paired R1.fastq.gz R2.fastq.gz || { echo "pairs failed validation" >&2; exit 1; }
z-fastq stats --json R1.fastq.gz > R1.stats.json
```

## Make an interleaved gzip file

```bash
z-fastq interleave R1.fastq.gz R2.fastq.gz | gzip > interleaved.fastq.gz
```

## Take a paired subset as two files

```bash
z-fastq sample --paired --fraction 0.05 --seed 2025 R1.fastq.gz R2.fastq.gz \
  | z-fastq deinterleave --out1 sub_R1.fastq --out2 sub_R2.fastq -
```

`sample --paired` writes one interleaved stream; `deinterleave` splits it again.

## Trim, filter, and report for MultiQC

```bash
z-fastq process --gzip \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --low-quality 15 --max-low-quality-percent 40 --max-n 5 --min-length 15 \
  --pair-duplication --qc-tables qc --qc-summary sample.json \
  --out1 trimmed_R1.fastq.gz --out2 trimmed_R2.fastq.gz R1.fastq.gz R2.fastq.gz
multiqc qc sample.json
```

These are fastp's default filters; see [Processing](Processing) and [Quality control](Quality-Control).

## Replace a file only on success

Write to a temporary path and move it after a zero exit status; [Automation](Automation#partial-output) shows the pattern and why it matters.

## Related pages

- [Automation](Automation): exit statuses and partial output, which every recipe relies on.
- [Command cheat sheet](Command-Cheat-Sheet): every option the recipes use.

# Recipes

These are short workflow patterns. Each command still owns its documented exit status and output behavior.

## Count several lanes

```bash
z-fastq count lane1.fastq.gz lane2.fastq.gz lane3.fastq.gz
```

The output has one count per successful input, in argument order.

## Gate a workflow on validation

```bash
if z-fastq check --paired R1.fastq.gz R2.fastq.gz; then
    z-fastq stats --json R1.fastq.gz > R1.stats.json
else
    echo "paired input failed validation" >&2
    exit 1
fi
```

Use `--pair-names exact` when the default Illumina-style normalization is too permissive for your workflow.

## Produce an interleaved gzip stream

```bash
z-fastq interleave R1.fastq.gz R2.fastq.gz | gzip > reads.interleaved.fastq.gz
```

The external `gzip` process owns output compression. z-fastq validates each pair before writing it.

## Sample a pair without breaking synchronization

```bash
z-fastq sample --paired --fraction 0.05 --seed 2025 R1.fastq R2.fastq \
    > sample.interleaved.fastq
```

The result is interleaved R1-then-R2 output. It is not two separate output files.

## Use exact-count sampling safely

```bash
temporary_sample="sample.fastq.tmp"
if z-fastq sample --count 10000 --seed 2025 reads.fastq > "$temporary_sample"; then
    mv "$temporary_sample" sample.fastq
else
    rm -f "$temporary_sample"
    exit 1
fi
```

Exact-count mode needs a regular file and two passes. The temporary path prevents a late input failure from replacing the previous sample.

## Keep data and diagnostics separate

```bash
z-fastq stats --json reads.fastq > stats.json 2> stats.err
status=$?
if [ "$status" -ne 0 ]; then
    cat stats.err >&2
    exit "$status"
fi
```

This matters for JSON and for commands that write FASTQ to standard output.


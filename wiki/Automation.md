# Automation

This page is for scripts and workflow engines that need to distinguish a valid result from a usable-looking partial result.

## Exit statuses

- `0`: The command completed successfully.
- `1`: FASTQ or paired-read validation failed.
- `2`: The invocation is invalid.
- `3`: I/O, gzip, output, pipe, or unexpected allocation failure.
- `4`: A configured line, sample-selection, staging, or arithmetic limit failed.

For independent `count`, `stats`, and single-input `check` inputs, z-fastq can continue after one input fails. The final status is the highest nonzero class. Paired operations stop at the first failed pair operation.

## Separate data from diagnostics

Human summaries and FASTQ data use standard output. Errors use standard error. Capture them separately when a workflow needs both:

```bash
if z-fastq stats reads.fastq > reads.stats.txt 2> reads.stats.err; then
    echo "statistics succeeded"
else
    echo "statistics failed" >&2
    cat reads.stats.err >&2
    exit 1
fi
```

FASTQ-producing commands can write earlier records before a later failure. Check the exit status before treating redirected output as complete.

## Machine-readable output

The current pre-alpha build exposes `--json` for `stats` and `check`:

```bash
z-fastq stats --json reads.fastq.gz > stats.json
z-fastq check --json reads.fastq.gz > check.json
```

The document shape is still evolving and is not a public compatibility contract. If you automate against it during pre-alpha, pin the z-fastq version and validate upgrades deliberately. This page documents the stable workflow boundary around output streams and exit statuses, not an evolving machine-readable format.

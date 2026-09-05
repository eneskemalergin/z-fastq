# Paired reads

z-fastq treats R1 and R2 as one logical unit whenever it samples, checks, interleaves, or deinterleaves paired data. A pair is either one record from each input file or two consecutive records in an interleaved input.

## Validate two mate files

```bash
z-fastq check --paired R1.fastq.gz R2.fastq.gz
```

The command reads the two inputs in lockstep. It validates both records before comparing their names. If one input ends first, the command reports a pair-count mismatch.

## Validate an interleaved pair stream

```bash
z-fastq check --interleaved reads.interleaved.fastq
```

Records are grouped as `(record 0, record 1)`, `(record 2, record 3)`, and so on. The first record in each group is R1. An odd final record has no mate.

## Default name policy

The default `illumina` policy handles the common forms used by Illumina and SRA-style FASTQ files. It compares a normalized identifier while leaving descriptive fields alone.

The policy:

1. Takes the first whitespace-delimited token from each header.
2. Removes a terminal `/1` or `/2` marker from that token.
3. Checks a leading `1:` or `2:` marker in the next token when present.
4. Checks a terminal `/1` or `/2` marker in that next token when present.
5. Requires the normalized identifiers to match.

These headers match under the default policy:

```text
@read42/1
@read42/2
```

These headers also match:

```text
@read42 1:N:0:1
@read42 2:N:0:1
```

Multiple recognized markers in one header must agree. A marker present on only one side is rejected. If neither header has a recognized marker, equal normalized identifiers are accepted. Fields after the recognized mate token are opaque.

## Exact name policy

Use `exact` when the complete first token must match byte for byte:

```bash
z-fastq check --paired --pair-names exact R1.fastq R2.fastq
```

This policy does not remove `/1` or `/2` and does not interpret a second token.

## Pair errors

- `P001`: Names do not match under the selected policy.
- `P002`: The pair inputs contain different numbers of records.

Pair diagnostics identify the zero-based pair index. Long header values are bounded in diagnostic output so an unusual name cannot create an unbounded error message.

## Where the pair unit matters

- `check --paired` compares two files.
- `check --interleaved` compares adjacent records.
- `sample --paired` selects both mates or neither mate and writes an interleaved stream.
- `sample --interleaved` preserves the same pair layout.
- `interleave` writes each validated pair as R1 then R2.
- `deinterleave` sends alternating records to R1 and R2.

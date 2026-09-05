# Validation

`check` answers one question: does this input satisfy the FASTQ and semantic rules that z-fastq currently accepts?

```bash
z-fastq check reads.fastq.gz
```

A valid input produces no standard output and exits with status `0`. This is deliberate. It makes `check` useful as a gate in a shell workflow.

## Record rules

Each record has four logical lines:

1. A header beginning with `@`.
2. A sequence line.
3. A plus line beginning with `+`.
4. A quality line with the same byte length as the sequence.

Quality bytes use Phred+33 and must be in the inclusive ASCII range 33 to 126. Empty sequence and quality lines are valid when both are empty. z-fastq expects sequence and quality on one logical line; it does not join wrapped lines into one record.

## Choose a sequence alphabet

The default `iupac` policy accepts upper- and lower-case IUPAC nucleotide symbols:

```bash
z-fastq check --alphabet iupac reads.fastq
```

Use `acgtn` when the workflow accepts only A, C, G, T, and N:

```bash
z-fastq check --alphabet acgtn reads.fastq
```

The policy validates the input. It does not rewrite the sequence or change its case.

## First failure and error codes

z-fastq reports the first structural or semantic failure it encounters. The codes are stable labels for the current validation contract:

- `S001`: The plus line does not begin with `+`.
- `S002`: A sequence byte is outside the selected alphabet.
- `S003`: The header does not begin with `@`.
- `S004`: The record ends before all four lines are present.
- `S005`: Sequence and quality lengths differ.
- `S006`: A quality byte is outside the Phred+33 range.

Diagnostics identify a zero-based record index, a zero-based offset in the decompressed FASTQ stream, and a one-based line number when those locations exist.

## Validate pairs

For two mate files:

```bash
z-fastq check --paired R1.fastq.gz R2.fastq.gz
```

For an interleaved file:

```bash
z-fastq check --interleaved reads.interleaved.fastq.gz
```

Read [Paired reads](Paired-Reads) before choosing a pair-name policy. Structural and semantic record errors take precedence over pair-name errors at the same position.

## Validation is not repair

`check` does not skip bad records, rewrite invalid quality bytes, repair names, or create a cleaned output. If you need a transformed file, first decide what the correction means for your data, then use a tool that owns that transformation explicitly.

# Choosing a FASTQ tool

No single FASTQ program is the right fit for every job. I think about the choice in terms of the work a pipeline needs, not a general claim that one tool is always faster.

## Use z-fastq when

- You need counting, aggregate statistics, validation, sampling, or pair layout changes.
- Plain and gzip input should follow the same command path.
- You want explicit behavior for paired names, validation codes, exit classes, and output creation.
- You want a streaming command that keeps normal process state bounded.
- You want to embed the reader, writer, statistics, or I/O adapters in a Zig program.
- You are comfortable with the current Linux x86-64 boundary and pre-`0.1.0` interface.

z-fastq is intentionally narrower than a general sequence manipulation suite. Its value is the consistency of these core operations and the fact that their failure behavior is part of the documented interface.

## Use seqtk when

[seqtk](https://github.com/lh3/seqtk) is a good fit for a compact, classic FASTA/Q utility with common transformations such as conversion, subsampling, trimming, and subsequence extraction. Its README is command-first and concise. If your workflow already depends on its established flags and output, changing tools may not be worth it.

## Use SeqKit when

[SeqKit](https://github.com/shenwei356/seqkit) is a better fit for a broad FASTA/Q toolbox with many subcommands, more compression formats, configurable threading, sequence searches, conversions, and larger workflow coverage. Its documentation is correspondingly broader. z-fastq does not try to reproduce that surface.

## Use fastp when

[fastp](https://github.com/OpenGene/fastp) is the right category of tool when the job includes short-read QC and preprocessing such as adapter handling, trimming, filtering, deduplication, and reports. z-fastq does not replace those transformations.

## Why keep another focused tool?

The gap I care about is smaller than a full preprocessing suite and more explicit than a minimal record loop. A workflow often needs to answer simple questions before the expensive steps begin:

- Are these records structurally valid?
- Are the sequence and quality fields internally consistent?
- Do these two files still describe the same pairs?
- How many records and bases are present?
- Can I take a repeatable subset without loading the complete file?
- What exactly should a script do when a record fails halfway through?

Those questions are the reason for z-fastq. They are not a reason to claim that it replaces every tool above.


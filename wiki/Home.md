# z-fastq

z-fastq is a focused FASTQ toolkit for counting, summarizing, checking, sampling, and working with paired reads. It runs as one Linux x86-64 executable and also exposes a small Zig library.

- Current version: `0.0.15`
- Current status: polishing

I am building z-fastq for a practical reason: FASTQ work should remain usable on ordinary machines. A tool can finish a job quickly and still be a poor workflow component if every process consumes more memory than the job needs. I want the common paths to stream, keep their state bounded, and make their behavior easy to inspect when something goes wrong.

Use the sidebar to start with installation, move to a task page, check behavior boundaries, and finish with recipes or the Zig API.

## Start here

Follow [Getting started](Getting-Started) for a small end-to-end run. It creates a test file, checks it, counts it, reads it as gzip and stdin, and makes a reproducible sample.

## What are you here to do?

- [Install or build the binary](Installation).
- [Find the right command](Command-Cheat-Sheet).
- [Count reads or inspect aggregate quality](Counting-and-Statistics).
- [Reject malformed or semantically invalid records](Validation).
- [Understand R1 and R2 matching](Paired-Reads).
- [Sample records without losing pair synchronization](Sampling).
- [Convert between two-file and interleaved pairs](Interleaving).
- [Use gzip, pipes, and output streams correctly](Input-and-Output).
- [Integrate results into a script](Automation).
- [Compare z-fastq's scope with other FASTQ tools](Choosing-a-FASTQ-Tool).
- [Recover from a failure](Troubleshooting).
- [Use the Zig API](Zig-Library).

## What z-fastq does

The CLI currently has six commands:

- `count` counts records from one or more plain or gzip inputs.
- `stats` calculates read, base, GC, and Phred quality summaries.
- `check` validates records and, when requested, paired names.
- `sample` selects records or read pairs by fraction or exact count.
- `interleave` combines two validated mate streams into one stream.
- `deinterleave` splits one validated interleaved stream into two files.

The reader accepts plain FASTQ, gzip-compressed FASTQ, and explicit standard input. FASTQ-producing commands write plain four-line records with LF endings. The current CLI does not write gzip output.

## What it is not

z-fastq is not a trimmer, adapter remover, aligner, duplicate remover, or general FASTA/Q Swiss Army knife. It does not repair malformed records silently. It does not support every compression container. Those boundaries are intentional while the core tool is being polished.

Read [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool) if you are deciding whether z-fastq belongs in a larger preprocessing workflow.

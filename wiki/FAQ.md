# FAQ

## Does z-fastq repair malformed FASTQ?

No. It reports the first problem and stops; it never skips records or guesses a fix. See [Validation](Validation).

## Does it load whole files into memory?

No. Every command streams, and exact-count sampling keeps only a small entry per chosen record. See [Limits](Limits-and-Supported-Formats#memory).

## Does the `.gz` name decide gzip decoding?

No, the first two bytes do. See [Input and output](Input-and-Output).

## Why does `check` print nothing on success?

So it can be a plain pass or fail step in a script: the exit status is the answer.

## Which quality encoding does it read?

Phred+33 only, ASCII 33 to 126.

## Can I sample an exact count from a pipe?

No, exact counts read the file twice. Use `--fraction` on a pipe. See [Sampling](Sampling).

## Are paired mates ever sampled separately?

Never. Both mates are kept or neither.

## Does it write gzip?

Only `process --gzip`. Pipe the other commands' output to `gzip`.

## Can an output file be overwritten?

No. `deinterleave` and `process` create new files and fail on an existing path.

## Does it support wrapped FASTQ or FASTA?

Not today. See [Limits](Limits-and-Supported-Formats).

## Is z-fastq a replacement for fastp, FastQC, or SeqKit?

No. It follows fastp's trimming rules and FastQC's QC definitions for the parts it covers, and does much less than any of them. See [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool).

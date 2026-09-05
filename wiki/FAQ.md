# FAQ

## Does z-fastq repair malformed FASTQ?

No. `check` reports the first failure and stops that input. It does not skip records or write a guessed correction.

## Does it load a complete file into memory?

The normal count, stats, validation, and streaming pair paths do not. Exact-count sampling stores selection state and reads a regular file twice. It does not retain the complete FASTQ.

## Does a `.gz` suffix control gzip decoding?

No. z-fastq checks the first two bytes. A valid gzip stream is decoded even when its name does not end in `.gz`.

## Does z-fastq support wrapped FASTQ?

No. The current record model has one logical line for each header, sequence, plus, and quality field.

## Which quality encoding does it use?

Phred+33 only. Quality bytes must be ASCII 33 through 126. Phred+64 is not supported.

## Why does `check` print nothing on success?

It is designed to act as a validation gate. A zero exit status is the success result; diagnostics appear only when something fails.

## Can I use stdin for exact-count sampling?

No. Exact-count sampling requires two passes. Use a regular path or use fraction sampling for a one-pass pipe.

## Does paired sampling select R1 and R2 independently?

No. A pair is the selection unit. Both mates are written to the interleaved output or both are skipped.

## Does z-fastq write gzip output?

No. Pipe FASTQ output through `gzip` when you need compression.

## Can `deinterleave` overwrite an existing file?

No. It creates output paths exclusively to avoid silently truncating an existing result.

## Is z-fastq a replacement for fastp or SeqKit?

No. It covers a narrower set of streaming FASTQ tasks. See [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool).


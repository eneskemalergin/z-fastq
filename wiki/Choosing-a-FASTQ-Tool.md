# Choosing a FASTQ tool

No FASTQ program is right for every job, and z-fastq is not trying to be. I think about the choice in terms of what a workflow needs, not which tool wins a race.

## Why z-fastq exists

Before the expensive steps of an analysis, a workflow usually needs plain answers to small questions:

- Are these records valid, and do the two files still describe the same pairs?
- How many reads and bases are there, and how good are they?
- Can I take a repeatable subset without loading the file?
- What exactly happens when a record fails halfway through?

z-fastq is built for those questions: a small, fast tool whose failures are as well defined as its results. It is narrower than a general toolkit on purpose.

## Use z-fastq when

- You need counting, statistics, validation, sampling, or pair layout changes, on plain or gzip input alike.
- You want exact behavior for pair names, error codes, exit statuses, and output files.
- Memory should stay small and bounded on any input.
- You want adapter, poly-G, and quality trimming with fastp's rules, and QC that MultiQC can read, from one streaming pass.
- You want the reader, writer, or statistics inside a Zig program.

## Use another tool when

- **[fastp](https://github.com/OpenGene/fastp)** for the full short-read preprocessing set: deduplication, base correction, pair merging, UMI handling, automatic adapter detection, and its own HTML report. z-fastq follows fastp's rules for the trims and filters it has, but does none of these.
- **[FastQC](https://github.com/s-andrews/FastQC)** for its HTML report and interactive viewer, Kmer Content, and quality encodings other than Phred+33. z-fastq measures FastQC's modules and writes their data, not the report.
- **[SeqKit](https://github.com/shenwei356/seqkit)** for a broad FASTA and FASTQ toolbox: many subcommands, more compression formats, threads, searching, and conversion.
- **[seqtk](https://github.com/lh3/seqtk)** for a compact classic with conversion, subsampling, trimming, and subsequences. If a workflow already depends on its flags and output, switching may not be worth it.

[References](References) lists these and the other tools z-fastq is compared with, and how each comparison is made.

## Related pages

- [Performance](Performance): the speed and memory comparison with each of these tools.
- [Home](Home): what z-fastq does, command by command.
- [Limits and supported formats](Limits-and-Supported-Formats): what is not supported today.
- [References](References): the tools z-fastq follows and is compared with.

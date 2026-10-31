# References

z-fastq is written from scratch in Zig, but it does not invent its own rules where people already rely on established ones. This page lists the tools it follows, the tools it is compared with, and the public data it is checked on, with versions, sources, and citations. If you use z-fastq's trimming or QC in published work, please cite the tools whose rules it follows as well.

## Tools z-fastq follows

z-fastq follows these tools' rules and file formats for the parts listed. It does not run them, and it is not a replacement for any of them: each does much more than the part z-fastq follows.

- **[fastp](https://github.com/OpenGene/fastp) 1.1.0** [1, 2]: the trimming and filtering rules of `process`, insert sizes, the pair duplication rate, and the layout of the QC summary. Results are also checked against fastp 1.4.0.
- **[FastQC](https://github.com/s-andrews/FastQC) 0.13.0** [3]: the definitions of the QC modules, their statuses at the default limits, and the layout of the QC tables.
- **[FastQC-Rust](https://github.com/ewels/FastQC-Rust) 1.1.0**: the QC tables are checked byte for byte against its output.
- **[MultiQC](https://github.com/MultiQC/MultiQC) 1.33** [4]: reads the QC tables with its FastQC module and the QC summary with its fastp module.

z-fastq itself is built with [Zig](https://ziglang.org) 0.16.0 and [zipir](https://github.com/eneskemalergin/zipir) 0.2.1 for gzip.

## Tools z-fastq is compared with

These tools are run beside z-fastq to check that results agree and to measure speed and memory on the same inputs. The comparisons are development evidence, not a ranking: every tool here was built for its own set of jobs. The pinned versions live in `tools/versions.sh`.

| Tool                                                        | Version  | Compared on                         | Cite    |
| ----------------------------------------------------------- | -------- | ----------------------------------- | ------- |
| [BBTools](https://github.com/bbushnell/BBTools)             | 40.05    | sampling, pairs, layout             |         |
| [Fasten](https://github.com/lskatz/fasten)                  | 0.9.0    | interleaved sampling                | [5]     |
| [fastp](https://github.com/OpenGene/fastp)                  | 1.4.0    | processing                          | [1, 2]  |
| [FastQC](https://github.com/s-andrews/FastQC)               | 0.13.0   | QC modules                          | [3]     |
| [FastQC-Rust](https://github.com/ewels/FastQC-Rust)         | 1.1.0    | QC modules and tables               |         |
| [FastQValidator](https://github.com/statgen/fastQValidator) | 0.1.1a   | validation                          |         |
| [fq](https://github.com/stjude-rust-labs/fq)                | 0.12.0   | validation, paired sampling         |         |
| [fqkit](https://github.com/BioinfoToolbox/fqkit)            | 0.4.14   | stats, sampling, layout             |         |
| [fqtools](https://github.com/alastair-droop/fqtools)        | 2.3      | counting, validation                | [6]     |
| [Helicase](https://github.com/imartayan/helicase)           | 0.2.0    | counting, stats                     |         |
| [IRMA Core](https://github.com/CDCgov/irma-core)            | 0.11.0   | sampling, layout                    |         |
| [Needletail](https://github.com/onecodex/needletail)        | 0.7.3    | counting, stats                     |         |
| [Rasusa](https://github.com/mbhall88/rasusa)                | 5.1.0    | sampling                            | [7]     |
| [SeqFu](https://github.com/telatin/seqfu2)                  | 1.29.0   | counting, stats, validation, layout | [8]     |
| [SeqKit](https://github.com/shenwei356/seqkit)              | 2.14.0   | stats, sampling, pairs              | [9, 10] |
| [seqtk](https://github.com/lh3/seqtk)                       | 1.5-r133 | counting, sampling, layout          |         |

- "Layout" means interleaving and deinterleaving.
- Helicase and Needletail are Rust parsing libraries; a small wrapper around each counts records and bases so they can be compared as commands.
- A blank citation means the tool has no paper to cite; cite its repository and version.

## Data z-fastq is checked on

### Public sequencing runs

Eleven public runs from the European Nucleotide Archive cover short and long reads, single and paired layouts, fixed and variable lengths, and several instruments.

| Accession                                                         | Platform | Instrument    | Sample            | Layout | Read length       |
| ----------------------------------------------------------------- | -------- | ------------- | ----------------- | ------ | ----------------- |
| [DRR018439](https://www.ebi.ac.uk/ena/browser/view/DRR018439)     | Illumina | MiSeq         | _E. coli_         | paired | 35 to 151         |
| [DRR018447](https://www.ebi.ac.uk/ena/browser/view/DRR018447)     | Illumina | MiSeq         | _E. coli_         | paired | 151               |
| [DRR1001190](https://www.ebi.ac.uk/ena/browser/view/DRR1001190)   | Illumina | NovaSeq 6000  | feces metagenome  | paired | 150               |
| [DRR147610](https://www.ebi.ac.uk/ena/browser/view/DRR147610)     | Illumina | HiSeq 3000    | quail miRNA       | single | 36                |
| [ERR4686103](https://www.ebi.ac.uk/ena/browser/view/ERR4686103)   | Illumina | NovaSeq 6000  | human             | paired | variable          |
| [SRR19571648](https://www.ebi.ac.uk/ena/browser/view/SRR19571648) | Illumina | NovaSeq 6000  | mouse single-cell | paired | 28 and 150        |
| [ERR10448913](https://www.ebi.ac.uk/ena/browser/view/ERR10448913) | MGI      | BGISEQ-500    | mouse             | paired | 100               |
| [ERR10496769](https://www.ebi.ac.uk/ena/browser/view/ERR10496769) | ONT      | GridION       | _E. coli_         | single | 110 to 92,287     |
| [ERR14758576](https://www.ebi.ac.uk/ena/browser/view/ERR14758576) | ONT      | PromethION    | human             | single | 30,000 to 291,910 |
| [ERR10472470](https://www.ebi.ac.uk/ena/browser/view/ERR10472470) | PacBio   | Sequel II CCS | _K. oxytoca_      | single | 403 to 31,646     |
| [SRR29396856](https://www.ebi.ac.uk/ena/browser/view/SRR29396856) | PacBio   | Revio HiFi    | human             | single | about 16 kb       |

Each run is cut into nested prefixes of whole records, sized by plain bytes with both mates counted: tiny (1 MiB), xsmall (10 MiB), and small (100 MiB) for every run, medium (1 GiB) for seven of them, and large (10 GiB) for five. Each original download's size and MD5 checksum are pinned in `bench/lib/corpus/real-data.json`, and `bench/lib/corpus/make.sh` downloads and cuts them.

### Generated files

- `bench/lib/shapes/` generates synthetic reads in six shapes (dense, variable, endings, paired, long, and long-read) for timing across sizes.
- `tests/data/synthetic/` holds small hand-made files, valid and malformed, that the tests check record by record.

## Citations

1. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. _Bioinformatics_. 2018;34(17):i884–i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
2. Chen S. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. _iMeta_. 2023;2(2):e107. [doi:10.1002/imt2.107](https://doi.org/10.1002/imt2.107)
3. Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. [bioinformatics.babraham.ac.uk/projects/fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics_. 2016;32(19):3047–3048. [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)
5. Katz LS, Phan J, den Bakker HC. Fasten: a toolkit for streaming operations on fastq files. _Journal of Open Source Software_. 2024;9(94):6030. [doi:10.21105/joss.06030](https://doi.org/10.21105/joss.06030)
6. Droop AP. fqtools: an efficient software suite for modern FASTQ file manipulation. _Bioinformatics_. 2016;32(12):1883–1884. [doi:10.1093/bioinformatics/btw088](https://doi.org/10.1093/bioinformatics/btw088)
7. Hall M. Rasusa: Randomly subsample sequencing reads to a specified coverage. _Journal of Open Source Software_. 2022;7(69):3941. [doi:10.21105/joss.03941](https://doi.org/10.21105/joss.03941)
8. Telatin A, Fariselli P, Birolo G. SeqFu: a suite of utilities for the robust and reproducible manipulation of sequence files. _Bioengineering_. 2021;8(5):59. [doi:10.3390/bioengineering8050059](https://doi.org/10.3390/bioengineering8050059)
9. Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. _PLOS ONE_. 2016;11(10):e0163962. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)
10. Shen W, Sipos B, Zhao L. SeqKit2: a Swiss army knife for sequence and alignment processing. _iMeta_. 2024;3(3):e191. [doi:10.1002/imt2.191](https://doi.org/10.1002/imt2.191)

## Related pages

- [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool): when to use one of these tools instead.
- [Processing paired reads](Processing) and [Quality control](Quality-Control): the parts that follow fastp and FastQC.

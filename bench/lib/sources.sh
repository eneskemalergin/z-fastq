# Corpus sources every feature shares, sourced by bench.sh: one short
# description per source, the bar groups the figures pool sources into, and
# the source lists the feature run.sh files declare their inputs from.
# The real sources with a medium (1 GiB) tier; the others are too small to cut one.
MEDIUM_SOURCES=illumina-drr1001190,illumina-err4686103,mgi-err10448913,illumina-drr147610,illumina-srr19571648,ont-err14758576,pacbio-srr29396856
SHAPE_SOURCES=dense,endings,variable,longread,long
# The long-read sources and shapes; every other source holds short reads. A
# family parameter written short:A,long:B takes its value from this class.
LONG_SOURCES=ont-err10496769,ont-err14758576,pacbio-err10472470,pacbio-srr29396856,longread,long

describe illumina-drr018447 "MiSeq, 151 bp"
describe illumina-drr1001190 "NovaSeq, 150 bp"
describe illumina-err4686103 "NovaSeq, Casava headers"
describe mgi-err10448913 "BGISEQ, 100 bp"
describe illumina-drr018439 "MiSeq, 35 to 151 bp"
describe illumina-srr19571648 "10x mates, 28 and 150 bp"
describe illumina-drr147610 "miRNA, 36 bp"
describe ont-err10496769 "GridION, 0.1 to 92 kb"
describe ont-err14758576 "PromethION, 30 to 153 kb"
describe pacbio-err10472470 "Sequel II, 0.4 to 30 kb"
describe pacbio-srr29396856 "Revio HiFi, 4 to 21 kb"
describe dense "uniform 150 bp, bare plus line"
describe endings "100 bp, CRLF and LF lines"
describe variable "50 to 150 bp, annotated plus line"
describe long "300 kb sequence lines"
describe longread "nanopore headers, rare 300 kb reads"
describe paired "100 bp mate pairs"
# The real sources with R1 and R2 files at medium.
PAIRED_SOURCES=illumina-drr1001190,illumina-err4686103,mgi-err10448913,illumina-srr19571648
describe all "every xsmall corpus file"

# Bar groups: real runs pooled by the kind of reads they hold, then the
# synthetic shapes, each one a marked dot in the last group.
group "Short, fixed length" illumina-drr018447,illumina-drr1001190,illumina-err4686103,mgi-err10448913,illumina-drr147610
group "Short, variable length" illumina-drr018439,illumina-srr19571648
group "Long reads" ont-err10496769,ont-err14758576,pacbio-err10472470,pacbio-srr29396856
group Synthetic "$SHAPE_SOURCES,paired"

# The inputs the headline file family of every feature declares: the real
# sources and shapes at 1 GiB, and the shapes at the smaller scaling tiers.
corpus_inputs() {
    input real "$MEDIUM_SOURCES" medium interleaved
    input shape "$SHAPE_SOURCES" medium
    input shape "$SHAPE_SOURCES" tiny,xsmall,small
}
# The stdin family's inputs: short, very short, and long reads, and CRLF lines.
stdin_inputs() {
    input real illumina-drr1001190,illumina-drr147610,ont-err14758576 medium interleaved
    input shape endings medium
}

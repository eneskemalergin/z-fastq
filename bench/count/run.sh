#!/usr/bin/env bash
# Benchmark: records in a FASTQ file. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "records in a FASTQ file"
truth "the record count the corpus cutter wrote into each manifest (and each shape's state.json), derived without z-fastq"
headline small
scaling dense,endings,variable,longread,long
pick auto
remark scaling "Every tool is held to one CPU, so a tool that starts threads (SeqFu starts 8 by default) gains nothing from them here; its CPU time equals its wall time. Lines rise with size because each tool's fixed start cost spreads over more data; SeqFu has the largest start cost, so its line rises the most."

lane z-fastq subject "{zfastq} count {in}"
lane fqtools complete "{bin} count {in}" headline-only
lane needletail complete "{bin} count {in}" note="built for this CPU"
lane helicase-checked partial "{bin} count {in}" note="misses a truncated last record and a bad plus line"
lane seqfu partial "{bin} count {in}" note="accepts a truncated last record; starts 8 threads by default, held to one CPU here like every tool"
lane seqtk partial "{bin} size {in}" note="no structure check"

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
describe all "every xsmall corpus file"

# Bar groups: real runs pooled by the kind of reads they hold, then the
# synthetic shapes, each one a marked dot in the last group.
group "Short, fixed length" illumina-drr018447,illumina-drr1001190,illumina-err4686103,mgi-err10448913,illumina-drr147610
group "Short, variable length" illumina-drr018439,illumina-srr19571648
group "Long reads" ont-err10496769,ont-err14758576,pacbio-err10472470,pacbio-srr29396856
group Synthetic dense,endings,variable,longread,long

family file "one file by path"
lanes all
input real illumina-drr018447,illumina-drr1001190,illumina-err4686103,mgi-err10448913,illumina-drr018439,illumina-srr19571648,illumina-drr147610,ont-err10496769,ont-err14758576,pacbio-err10472470,pacbio-srr29396856 small interleaved
input shape dense,endings,variable,longread,long small
input shape dense,endings,variable,longread,long tiny,xsmall,medium

family stdin "the same file on stdin"
lanes z-fastq,needletail,helicase-checked,seqfu,seqtk stdin
input real illumina-drr1001190,illumina-drr147610,ont-err14758576 small interleaved
input shape endings small

family many-files "every xsmall corpus file in one call"
lanes z-fastq,seqfu
input many xsmall interleaved

edges count

# One count per line: SeqFu prints NAME COUNT SE, seqtk prints COUNT BASES.
fact() {
    case "$LANE" in
        seqfu) awk '{ print $2 }' ;;
        seqtk) awk '{ print $1 }' ;;
        *) cat ;;
    esac
}

bench_main "$@"

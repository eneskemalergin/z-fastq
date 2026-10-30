#!/usr/bin/env bash
# Benchmark: records in a FASTQ file. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "records in a FASTQ file"
truth "the record count the corpus cutter wrote into each manifest (and each shape's state.json), derived without z-fastq"
remark scaling "Every tool is held to one CPU, so a tool that starts threads (SeqFu starts 8 by default) cannot use extra CPU cores here. Elapsed time can still exceed process CPU time because of scheduling and I/O waits. Lines rise with size because each tool's fixed start cost spreads over more data; SeqFu has the largest start cost, so its line rises the most."

lane z-fastq subject "{zfastq} count {in}"
lane fqtools complete "{bin} count {in}" headline-only
lane needletail complete "{bin} count {in}" note="built for this CPU"
lane helicase-checked partial "{bin} count {in}" note="misses a truncated last record and a bad plus line"
lane seqfu partial "{bin} count {in}" note="accepts a truncated last record; starts 8 threads by default, held to one CPU here like every tool"
lane seqtk partial "{bin} size {in}" note="no structure check"

family file "one file by path"
lanes all
corpus_inputs

family stdin "the same file on stdin"
lanes z-fastq,needletail,helicase-checked,seqfu,seqtk stdin
stdin_inputs

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

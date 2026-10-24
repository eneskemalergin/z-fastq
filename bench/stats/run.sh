#!/usr/bin/env bash
# Benchmark: aggregate statistics of a FASTQ file. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "aggregate statistics of a FASTQ file"
truth "the read count the corpus cutter wrote into each manifest; every other field from Needletail, an independent parser that must agree with z-fastq exactly"
agree needletail reads

word title "Aggregate statistics"
word subject "z-fastq stats"
word others "tools"
word timed "its values were right"
word right "the right values"
word wrong "wrong values"
word missed "got a value wrong"
word unseen "the values are wrong"
word legend "wrong value, exit 0"
word verb "summarizes"
word past "summarized"
word payer "a tool"
word same "statistics"
word many "Summarizing"
word truth "Reads are compared with the count recorded when the corpus was cut, and every other field with Needletail on the same input, an independent parser that must agree with z-fastq exactly. A printed decimal may differ by half its last digit."
word cells "A cell is \`ok\` (every value printed is right), \`reject\` (exited with an error), \`wrong: FIELD\` (printed a wrong value, exit 0), \`refused\` (rejected a valid input), or \`crash\`."
word probes "The bold count is how many of four malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, \`-\` instead of \`+\`, and a quality byte below 33."

lane z-fastq subject "{zfastq} stats {in}"
lane needletail complete "{bin} stats {in}" note="built for this CPU"
lane helicase-checked partial "{bin} stats {in}" note="accepts a bad plus line"
lane seqtk partial "{bin} fqchk {in}" note="no structure check; per-position table"
lane seqkit partial "{bin} stats -a -T -j 1 {in}" note="no structure check; keeps every read length"
lane seqfu partial "{bin} stats --threads 1 {in}" note="lengths only: no base or quality counts; accepts a bad plus line and quality bytes"

family file "one file by path"
lanes all
corpus_inputs

family stdin "the same file on stdin"
lanes all stdin
stdin_inputs

family many-files "every xsmall corpus file in one call"
lanes z-fastq,seqkit,seqfu
input many xsmall interleaved

edges stats

# FIELD VALUE [TOLERANCE] lines, "--" before each file; a printed decimal may
# differ by half its last digit. Field names are z-fastq's; a percentage a
# peer prints is compared with the same percentage from z-fastq's counts.
fact() {
    case "$LANE" in
        seqtk)
            awk -F'\t' 'NR == 1 { print "--"; n = split($0, p, "; ")
                    for (i = 1; i <= n; i++) if (split(p[i], kv, ": ") == 2) f[kv[1]] = kv[2]
                    print "min_length", f["min_len"]; print "max_length", f["max_len"]; print "mean_length", f["avg_len"], 0.005 }
                $1 == "ALL" { print "bases", $2; print "a_pct", $3, 0.05; print "c_pct", $4, 0.05; print "g_pct", $5, 0.05
                    print "t_pct", $6, 0.05; print "non_acgt_pct", $7, 0.05; print "mean_quality", $8, 0.05; print "q20_pct", $11, 0.05 }'
            ;;
        seqkit)
            awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
                { print "--"; print "reads", $c["num_seqs"]; print "bases", $c["sum_len"]; print "min_length", $c["min_len"]
                  print "max_length", $c["max_len"]; print "mean_length", $c["avg_len"], 0.05; print "n", $c["sum_n"]
                  print "q20_pct", $c["Q20(%)"], 0.5; print "q30_pct", $c["Q30(%)"], 0.5; print "gc_all_pct", $c["GC(%)"], 0.005 }'
            ;;
        seqfu)
            awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
                { print "--"; print "reads", $c["#Seq"]; print "bases", $c["Total bp"]; print "min_length", $c["Min"]
                  print "max_length", $c["Max"]; print "mean_length", $c["Avg"], 0.005 }'
            ;;
        *)
            # z-fastq and the adapters print "name: value" blocks, "-" for no value.
            awk -F': ' 'function close_block(  b) {
                    b = v["bases"]; if (!open || b <= 0) return
                    printf "gc_all_pct %.6f\nq20_pct %.6f\nq30_pct %.6f\nnon_acgt_pct %.6f\n", 100 * (v["g"] + v["c"]) / b,
                        100 * v["q20_bases"] / b, 100 * v["q30_bases"] / b, 100 * (v["n"] + v["other_bases"]) / b
                    printf "a_pct %.6f\nc_pct %.6f\ng_pct %.6f\nt_pct %.6f\n", 100 * v["a"] / b, 100 * v["c"] / b, 100 * v["g"] / b, 100 * v["t"] / b }
                $1 == "input" { close_block(); print "--"; open = 1; delete v; next }
                $2 != "-" { print $1, $2; v[$1] = $2 }
                END { close_block() }'
            ;;
    esac
}

bench_main "$@"

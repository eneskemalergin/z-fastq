#!/usr/bin/env bash
# Benchmark: validating FASTQ files. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "validating FASTQ files"
truth "every corpus file is valid FASTQ, so every tool must accept it; the edge cases and the defect sweep hold the inputs a tool must reject"
expect accept

word title "Validating FASTQ files"
word subject "z-fastq check"
word others "validators"
word timed "it accepted the valid input"
word right "the right verdict"
word wrong "a wrong verdict"
word missed "gave a wrong verdict"
word unseen "the file is broken"
word legend "passed with exit 0"
word verb "validates"
word past "validated"
word payer "a validator"
word many "Validating"
word forms "Other layouts"
word forms_sub "Mates in two files, with names checked across them, and mates interleaved in one file."
word truth "Every corpus file is valid, so the right verdict on each is to accept it."
word cells "A cell is the tool's verdict, with the rule it named where it names one: \`accept\`, \`reject S003\`, \`reject\` (no rule named), \`wrong: accept\` (passed a malformed input with exit 0), \`refused\` (rejected a valid input), or \`crash\`."
word probes "The bold count is how many of four malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, \`-\` instead of \`+\`, and a quality byte below 33."

lane z-fastq subject "{zfastq} check {in}"
lane fq complete "{bin} lint {in}" note="IUPAC codes and empty reads are invalid under its rules"
lane fqtools partial "{bin} -dualr validate {in}" note="rejects CRLF; accepts a blank line"
lane fastqvalidator partial "{bin} --file {in} --minReadLen 1 --disableSeqIDCheck" note="minimum length and duplicate-name checks off; reads the whole file before exiting"

family single "one file by path"
lanes z-fastq,fq,fqtools,fastqvalidator
corpus_inputs

family paired "mates in two files, names checked across them"
lanes z-fastq,fq
command z-fastq "{zfastq} check --paired {in}"
input pair "$PAIRED_SOURCES,paired" medium

family interleaved "mates alternating in one file"
lanes z-fastq,fqtools
command z-fastq "{zfastq} check --interleaved {in}"
command fqtools "{bin} -dualr -i validate {in}"
input real "$PAIRED_SOURCES" medium interleaved
input shape paired medium interleaved

edges check

# The verdict: accept on exit 0, otherwise reject with the rule the tool named.
fact() {
    if [[ "$RC" == 0 ]]; then
        printf 'accept\n'
        return
    fi
    local code
    code="$(sed 's/\x1b\[[0-9;]*m//g' "$ERR" | grep -o -m1 -E ': \[?[SP]00[1-6]\]?:?' | grep -o -E '[SP]00[1-6]' | head -1)" || code=""
    printf 'reject%s\n' "${code:+ $code}"
}

# --- defect sweep ---
# Real reads with one injected defect each; per tool and kind of defect, how
# many it rejected and how many of those it placed in the broken record.
SWEEP_SOURCES="illumina-drr1001190/tiny_interleaved.fastq illumina-drr018439/tiny_interleaved.fastq ont-err14758576/tiny.fastq"
SWEEP_KINDS="S001 S002 S003 S004 S005 S006-low S006-high"
SWEEP_PER_KIND=40
SWEEP="$RESULTS/sweep.tsv"
declare -A KIND_LABEL=([S001]="plus line" [S002]="sequence byte" [S003]="header" [S004]="truncated" [S005]="lengths differ" [S006-low]="quality below 33" [S006-high]="quality above 126")

# Writes FILE with one KIND defect in a record SEED picks; prints the defect's
# first and last file line (the broken record).
inject() { # source kind seed file
    awk -v kind="$2" -v seed="$3" -v out="$4" '
        { line[NR] = $0 }
        END {
            srand(seed); n = int(NR / 4); r = 1 + int(rand() * n); h = 4 * (r - 1) + 1; s = h + 1; q = h + 3
            p = 1 + int(rand() * length(line[s]))
            if (kind == "S001") line[h + 2] = "-"
            else if (kind == "S002") line[s] = substr(line[s], 1, p - 1) "*" substr(line[s], p + 1)
            else if (kind == "S003") line[h] = "X" substr(line[h], 2)
            else if (kind == "S005") line[q] = substr(line[q], 2)
            else if (kind == "S006-low") line[q] = substr(line[q], 1, p - 1) "\037" substr(line[q], p + 1)
            else if (kind == "S006-high") line[q] = substr(line[q], 1, p - 1) "\177" substr(line[q], p + 1)
            last = (kind == "S004") ? s : NR
            for (i = 1; i <= last; i++) print line[i] > out
            print h, (kind == "S004") ? s + 1 : h + 3
        }' "$REAL/$1"
}

# The file line a lane's message (stdout and stderr) points at, or nothing.
reported_line() { # lane errfile
    local text
    text="$(sed 's/\x1b\[[0-9;]*m//g' "$2")"
    case "$1" in
        z-fastq) sed -n -E 's/.*\(record ([0-9]+), line ([0-9]+),.*/\1 \2/p' <<<"$text" | head -1 | awk '{ print 4 * $1 + $2 }' ;;
        fq) grep -o -m1 -E ':[0-9]+:[0-9]+: \[' <<<"$text" | cut -d: -f2 ;;
        fastqvalidator) grep -o -m1 -E 'Line [0-9]+' <<<"$text" | cut -d' ' -f2 ;;
        fqtools) grep -o -m1 -E '\[line [0-9]+\]' <<<"$text" | tr -dc '0-9' ;;
    esac
}

sweep_one() { # "kind index"
    local kind="${1% *}" i="${1#* }" sources src seed file first last name rc line
    read -r -a sources <<<"$SWEEP_SOURCES"
    src="${sources[$((i % ${#sources[@]}))]}" seed=$((i * 7919 + ${#kind}))
    file="$(mktemp "/var/tmp/bench-check-sweep.XXXXXX")"
    read -r first last < <(inject "$src" "$kind" "$seed" "$file")
    for name in "${LANES[@]}"; do
        rc=0
        $(lane_cmd "$name" "$file") >"$file.err" 2>&1 || rc=$?
        line="$(reported_line "$name" "$file.err")" || line=""
        printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$kind" "$i" "$rc" \
            "$([[ -z "$line" ]] && echo - || { ((line >= first && line <= last)) && echo 1 || echo 0; })"
    done
    rm -f -- "$file" "$file.err"
}

sweep_jobs() {
    local kind i
    for kind in $SWEEP_KINDS; do
        for ((i = 0; i < SWEEP_PER_KIND; i++)); do printf '%s %s\n' "$kind" "$i"; done
    done
}
feature_check() {
    local i
    cached_step sweep "$(lane_keys; printf '%s %s %s\n' "$SWEEP_KINDS" "$SWEEP_PER_KIND" "$SWEEP_SOURCES"
        for i in $SWEEP_SOURCES; do fileid "$REAL/$i"; done; declare -f inject reported_line sweep_one | cksum)" sweep_jobs sweep_one "$SWEEP"
}

# --- report sections ---
sweep_table() {
    local labels="" kind
    for kind in $SWEEP_KINDS; do labels+="${KIND_LABEL[$kind]// /_} "; done
    awk -F'\t' -v kinds="$SWEEP_KINDS" -v labels="$labels" '
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        { n[$1, $2]++; if ($4 > 0 && $4 < 128) caught[$1, $2]++; if ($5 != "-") told[$1, $2]++; if ($5 == 1) at[$1, $2]++; lanes[$1] = 1 }
        END {
            k = split(kinds, ks, " "); split(labels, ls, " ")
            line = "Tool"; align = "l"; for (i = 1; i <= k; i++) { t = ls[i]; gsub("_", " ", t); line = line "\t" t; align = align "\tr" }
            print line; print align
            m = 0; for (l in lanes) list[++m] = l
            for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= m; i++) {
                l = list[i]; line = label[l]
                for (j = 1; j <= k; j++) {
                    c = caught[l, ks[j]] + 0; a = at[l, ks[j]] + 0; t = n[l, ks[j]]
                    line = line "\t" sprintf("%d%%", 100 * c / t) (told[l, ks[j]] ? sprintf(" (%d%%)", 100 * a / t) : "")
                }
                print line
            }
        }' "$REPORT_DIR/style.tsv" "$SWEEP" | mdtable
}

report_samples() {
    local file name rc first last text
    file="$(mktemp "/var/tmp/bench-check-sample.XXXXXX")"
    read -r first last < <(inject illumina-drr1001190/tiny_interleaved.fastq S006-low 11 "$file")
    printf 'One file of real reads with a quality byte below 33 in the record on lines %d to %d, and what each tool printed.\n\n' "$first" "$last"
    for name in "${LANES[@]}"; do
        rc=0
        $(lane_cmd "$name" "$file") >"$file.err" 2>&1 || rc=$?
        text="$(sed 's/\x1b\[[0-9;]*m//g' "$file.err" | grep -v -E ' INFO |^WARNING: Using incubator|^java |^Executing|^Version|^Picked up' |
            grep -v -E '^\s*$' | head -3 | sed "s#$file#reads.fastq#g")"
        printf -- '- `%s` exited %d:\n\n```text\n%s\n```\n\n' "${REG_LABEL[$name]}" "$rc" "${text:-(nothing printed)}"
    done
    rm -f -- "$file" "$file.err"
}

feature_report() {
    [[ -f "$SWEEP" ]] || return 0
    printf '## Defect sweep\n> <sub>Real reads with one injected defect each, %d files per kind, at random records of three 1 MiB files (short fixed, short variable, and long reads). Each cell is the share the tool rejected; in brackets, the share where the line it reported is in the broken record, left out for a tool that reports no line.</sub>\n\n' "$SWEEP_PER_KIND"
    sweep_table
    printf '\n## How each tool reports a problem\n\n'
    details "Messages for one broken record"
    report_samples
    printf '</details>\n\n'
}

bench_main "$@"

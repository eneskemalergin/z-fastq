#!/usr/bin/env bash
# Benchmark: converting paired reads between two mate files and one interleaved file. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "converting paired reads between two mate files and one interleaved file"
truth "the corpus holds every pair three ways, cut from the same records: R1, R2, and the interleaved file; a right output equals them byte for byte"
expect accept
facts piped
# count and check already run z-fastq on every corpus file; every pair here is in a family.
verify off
SUBSET_CHECKER=1
# Speed tables count what each tool lets through, since the work is not equal.
THROUGH=1
rival any
order z-fastq,fqkit,seqtk,irma-core,seqfu,fqtools
scaling paired

word title "Converting paired layouts"
word subject "z-fastq interleave\` and \`deinterleave"
word others "converters"
word timed "its output was byte-exact"
word right "byte-exact output"
word wrong "wrong output"
word missed "gave wrong output"
word unseen "the pairs are broken"
word legend "wrong output, exit 0"
word verb "converts"
word past "converted"
word payer "a converter"
word same "conversion"
word many "Converting"
word truth "The corpus holds every pair as R1, R2, and the interleaved file, cut from the same records, so a right interleave output equals the interleaved file and right deinterleave outputs equal R1 and R2, byte for byte; when bytes differ, the subset checker names how."
word cells "A cell is the tool's verdict: \`accept\` (exit 0, output byte-exact), \`reject P001\` (with the rule it named), \`reject\`, \`wrong: accept\` (exit 0 on a broken pair or record), \`wrong: plus line changed\` and the like (exit 0 on a valid input, other bytes), \`refused\` (rejected a valid input), \`crash\`, or \`-\` (the tool has no such command)."

# Interleave: two mate files in, one interleaved stream out on stdout.
lane z-fastq subject "{zfastq} interleave {in}"
lane fqkit partial "{bin} merge -1 {in1} -2 {in2} -@ 1" outputs=_r1.fq,_r2.fq note="4 threads unless given -@ 1"
lane seqtk partial "{bin} mergepe {in}"
lane irma-core complete "{bin} xleave {in}" outputs=_1.fq,_2.fq note="its pair check rejects mates whose names carry both /1 and 1:N:0"
lane seqfu partial "{bin} ilv -1 {in1} -2 {in2}" outputs=_R1.fq,_R2.fq
lane fqtools partial "{bin} -I view {in}" outputs=1.fastq,2.fastq

family interleave "two mate files into one interleaved file"
lanes all
look full
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

# Deinterleave: one interleaved file in, two mate files out, real files for
# every tool (z-fastq writes only new paths), removed after each run.
family deinterleave "one interleaved file into two mate files"
lanes z-fastq,fqkit,irma-core,seqfu,fqtools
look full
memory 8G
command z-fastq "{zfastq} deinterleave --out1 {stem}_1.fq --out2 {stem}_2.fq {in}"
command fqkit "{bin} split -p {name} -o {dir} -@ 1 {in}"
command irma-core "{bin} xleave -1 {stem}_1.fq -2 {stem}_2.fq {in}"
command seqfu "{bin} dei -o {stem} {in}"
command fqtools "{bin} -i -F F view -o {stem}% {in}"
input real "$PAIRED_SOURCES" medium interleaved
input shape paired tiny,xsmall,small,medium interleaved

# z-fastq's deinterleave outputs, by its command above.
LANE_OUTS[z-fastq]=_1.fq,_2.fq
edges layout

# --- edge cases ---------------------------------------------------------------
# Small pairs, each as R1, R2, and the interleaved file. Rows run either the
# interleave command on the two mate files (i_) or the deinterleave command on
# the interleaved file (d_).
PAIR_EDGES="$EDGES/pairs"
feature_edge_list() {
    cat <<'EOF'
i_ok	accept	valid pair (interleave)
i_plus	accept	annotated plus line (interleave)
i_crlf	policy	CRLF lines (interleave)
i_names	reject P001	mate names differ (interleave)
i_short	reject P002	R2 one record short (interleave)
i_lenmis	reject S005	length mismatch in R1 (interleave)
i_trunc	reject S004	truncated R1 record (interleave)
i_dash	reject S001	`-` instead of `+` in R1 (interleave)
i_qlow	reject S006	quality byte below 33 in R1 (interleave)
d_ok	accept	valid pair (deinterleave)
d_plus	accept	annotated plus line (deinterleave)
d_crlf	policy	CRLF lines (deinterleave)
d_names	reject P001	mate names differ (deinterleave)
d_odd	reject P002	odd record count (deinterleave)
d_lenmis	reject S005	length mismatch (deinterleave)
d_trunc	reject S004	truncated last record (deinterleave)
d_dash	reject S001	`-` instead of `+` (deinterleave)
d_qlow	reject S006	quality byte below 33 (deinterleave)
EOF
}
edge_family() { case "$1" in i_*) echo interleave ;; *) echo deinterleave ;; esac; }
edge_input() {
    local base="$PAIR_EDGES/${1#?_}"
    case "$1" in i_*) printf '%s_R1.fastq %s_R2.fastq\n' "$base" "$base" ;; *) printf '%s_interleaved.fastq\n' "$base" ;; esac
}
# Writes NAME_R1.fastq, NAME_R2.fastq, and NAME_interleaved.fastq from two
# record lists (one record per argument, fields separated by |).
pair_files() { # name r1-records r2-records [interleaved records]
    local name="$PAIR_EDGES/$1" rec
    rec() { tr '|' '\n' <<<"$1"; }
    : >"${name}_R1.fastq"; for r in $2; do rec "$r" >>"${name}_R1.fastq"; done
    : >"${name}_R2.fastq"; for r in $3; do rec "$r" >>"${name}_R2.fastq"; done
    if [[ -n "${4:-}" ]]; then
        : >"${name}_interleaved.fastq"; for r in $4; do rec "$r" >>"${name}_interleaved.fastq"; done
    else
        paste -d '\n' <(paste - - - - <"${name}_R1.fastq") <(paste - - - - <"${name}_R2.fastq") | tr '\t' '\n' >"${name}_interleaved.fastq"
    fi
}
make_feature_edges() {
    [[ -f "$PAIR_EDGES/qlow_interleaved.fastq" ]] && return 0
    mkdir -p -- "$PAIR_EDGES"
    local a1='@p1/1|ACGT|+|IIII' a2='@p2/1|GGCA|+|IIII' b1='@p1/2|TTGA|+|IIII' b2='@p2/2|CCAT|+|IIII'
    pair_files ok "$a1 $a2" "$b1 $b2"
    pair_files plus '@p1/1|ACGT|+p1/1|IIII @p2/1|GGCA|+p2/1|IIII' '@p1/2|TTGA|+p1/2|IIII @p2/2|CCAT|+p2/2|IIII'
    pair_files names "$a1 $a2" "@q1/2|TTGA|+|IIII $b2"
    pair_files short "$a1 $a2" "$b1"
    pair_files odd "$a1 $a2" "$b1 $b2" "$a1 $b1 $a2"
    pair_files lenmis "@p1/1|ACGT|+|III $a2" "$b1 $b2"
    pair_files trunc "$a1 @p2/1|GGCA" "$b1 $b2" "$a1 $b1 $a2 $b2 @p3/1|ACGT"
    pair_files dash "@p1/1|ACGT|-|IIII $a2" "$b1 $b2"
    # A quality byte below 33 (0x1f) in R1's first record.
    pair_files qlow "@p1/1|ACGT|+|IIXI $a2" "$b1 $b2"
    local f
    for f in qlow_R1 qlow_interleaved; do sed -i '4s/X/\x1f/' "$PAIR_EDGES/$f.fastq"; done
    for f in R1 R2 interleaved; do sed 's/$/\r/' "$PAIR_EDGES/ok_$f.fastq" >"$PAIR_EDGES/crlf_$f.fastq"; done
}

# --- truth --------------------------------------------------------------------
# The plain corpus file a path stands for: the copy beside a .gz.
plain_path() { printf '%s\n' "${1%.gz}"; }
# How an output that is not byte-exact differs, from the subset checker's facts.
how_wrong() { # facts expected-count
    awk -v n="$2" '{ f[$1] = $2 } END {
        if (f["order"] != "yes") print "wrong: not a subset"
        else if (("pairs" in f) && f["pairs"] != "yes") print "wrong: mates out of step"
        else if (f["kept"] + 0 < n + 0) print "wrong: records dropped"
        else if (f["plus"] + 0 > 0) print "wrong: plus line changed"
        else if (f["cr"] + 0 > 0) print "wrong: line endings changed"
        else print "wrong: bytes differ" }' <<<"$1"
}
# The verdict, then the checker's facts. Interleave: the stream on stdin
# against the interleaved file; deinterleave: OUT1 and OUT2 against R1 and R2.
fact() {
    if [[ "${RC:-0}" != 0 ]]; then
        reject_code
        return
    fi
    local -a in
    read -r -a in <<<"$IN"
    local r1 r2 il facts same=1 n
    if [[ "${FAMILY:-}" == interleave ]]; then
        r1="$(plain_path "${in[0]}")" r2="$(plain_path "${in[1]}")" il="${r1%_R1.fastq}_interleaved.fastq"
        local fd pid status
        status="$(mktemp /var/tmp/bench-layout-cmp.XXXXXX)"
        exec {fd}> >(cmp -s - "$il"; printf '%s\n' "$?" >"$status")
        pid=$!
        # tee keeps feeding the checker after cmp stops at a first difference.
        facts="$(tee -p /dev/fd/$fd | "$SUBSET" "$r1" "$r2")"
        exec {fd}>&-
        wait "$pid" 2>/dev/null || true
        [[ "$(<"$status")" == 0 ]] || same=0
        rm -f -- "$status"
        # The pair count: the case's truth, or counted for a small edge file.
        n="${TRUTH:-}"
        [[ -n "$n" && "$n" != - && -z "${EDGE:-}" ]] || n="$(($(wc -l <"$r1") / 4))"
    else
        il="$(plain_path "${in[0]}")" r1="${il%_interleaved.fastq}_R1.fastq" r2="${il%_interleaved.fastq}_R2.fastq"
        cmp -s -- "$OUT1" "$r1" && cmp -s -- "$OUT2" "$r2" || same=0
        local a b
        a="$("$SUBSET" "$r1" <"$OUT1")" b="$("$SUBSET" "$r2" <"$OUT2")"
        facts="$(awk 'FNR == NR { a[$1] = $2; next } { b[$1] = $2 }
            END { print "kept", a["kept"]; print "order", (a["order"] == "yes" && b["order"] == "yes") ? "yes" : "no"
                  print "pairs", (a["kept"] == b["kept"] && a["fingerprint"] == b["fingerprint"]) ? "yes" : "no"
                  print "plus", a["plus"] + b["plus"]; print "cr", a["cr"] + b["cr"] }' <(printf '%s\n' "$a") <(printf '%s\n' "$b"))"
        # The case's truth counts interleaved records: two per pair.
        n="${TRUTH:-}"
        if [[ -n "$n" && "$n" != - && -z "${EDGE:-}" ]]; then n=$((n / 2)); else n="$(($(wc -l <"$r1") / 4))"; fi
    fi
    if ((same)); then printf 'accept\n'; else how_wrong "$facts" "$n"; fi
    printf '%s\nn %s\n' "$facts" "$n"
}

# --- report -------------------------------------------------------------------
mode_name() { case "$1" in interleave) echo Interleave ;; deinterleave) echo Deinterleave ;; *) echo "$1" ;; esac; }
tools_table() {
    local name fam line
    {
        printf 'Tool\tClass\tInterleave\tDeinterleave\tMalformed rejected\n'
        printf 'l\tl\tl\tl\tr\n'
        for name in "${LANES[@]}"; do
            line="${REG_LABEL[$name]}"$'\t'"${LANE_CLASS[$name]}"
            for fam in interleave deinterleave; do
                if [[ " ${FAM_LANES[$fam]} " == *" $name "* ]]; then line+=$'\t'"$(mode_cell "$name" "$fam")"; else line+=$'\t-'; fi
            done
            printf '%s\t%s\n' "$line" "$(rejects_all "$name")"
        done
    } | mdtable
}
layout_findings() {
    awk -F'\t' -v others="$((${#LANES[@]} - 1))" 'FILENAME ~ /checks.tsv$/ { if ($1 == "z-fastq") { if ($3 == "pass") zp++; else zf++ } next }
        FNR > 1 && $6 != "policy" { if ($4 == "z-fastq") { if ($6 != "right") zbad++ } else if ($3 ~ /^reject/ && ($6 == "wrong" || $6 == "crash")) bad[$4] = 1 }
        END { line = zf ? sprintf("- **Right answers:** z-fastq gave wrong output on %d of %d checked inputs", zf, zp + zf) : sprintf("- **Right answers:** z-fastq wrote byte-exact output on all %d checked inputs", zp)
              line = line (zbad ? "." : " and rejected every broken pair and record.")
              k = 0; for (b in bad) k++
              if (k) line = line sprintf(" %d of the %d other tools write output from a broken pair or record with exit 0.", k, others)
              print line }' "$REPORT_DIR/checks.tsv" "$REPORT_DIR/edges.tsv"
    # Who notices mates out of step.
    awk -F'\t' 'FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0 } next }
        ($2 == "i_names" || $2 == "i_short" || $2 == "d_names" || $2 == "d_odd") && $3 ~ /^wrong/ { pass[$1] = 1 }
        END { n = 0; for (l in pass) list[++n] = l
              for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
              for (i = 1; i <= n; i++) out = out (i == 1 ? "" : i == n ? " and " : ", ") label[list[i]]
              print "- **Broken pairs:** " (out == "" ? "every tool rejects mates whose names differ or whose counts do not match." : out " write output from mates whose names differ or whose counts do not match, with exit 0.") }' \
        "$REPORT_DIR/style.tsv" "$RESULTS/edges.tsv"
    valid_wrong_bullet
    ANY=1 findings interleave speed Interleave
    ANY=1 findings deinterleave speed Deinterleave
    findings interleave memory Interleave
    findings deinterleave memory Deinterleave
    findings interleave where Interleave
    findings deinterleave where Deinterleave
}
feature_readme() {
    family_checks
    readme_head
    printf '## Findings\n\n'
    layout_findings
    printf '\n## Compared tools\n> <sub>*complete* tools check pairs and records as much as z-fastq; *partial* tools check less, so they do less work per record. A command is `yes` when the tool was right on every input of it, `partly` on some, `-` when the tool has no such command. Malformed rejected: of the broken pairs and records in the edge cases, how many the tool rejects.</sub>\n\n'
    tools_table
    printf '\n### What each tool checks\n> <sub>Each kind of broken input from the edge cases, over both commands: `rejects`, or **passes** when the tool writes output from it with exit 0, so it does less work than z-fastq; `-` when the tool has no command for it.</sub>\n\n'
    kinds_table "Names differ=i_names,d_names|Mate missing=i_short,d_odd|Length mismatch=i_lenmis,d_lenmis|Truncated record=i_trunc,d_trunc|Bad plus line=i_dash,d_dash|Quality range=i_qlow,d_qlow"
    printf '\n'
    local name
    for name in "${LANES[@]}"; do [[ -z "${LANE_NOTE[$name]}" ]] || printf -- '- `%s %s`: *%s*\n' "${REG_LABEL[$name]}" "$(lane_version "$name")" "${LANE_NOTE[$name]}"; done
    printf '\n> [!IMPORTANT]\n> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.\n'
    readme_right
    readme_speed interleave "Interleave"
    readme_scaling interleave "###"
    readme_speed deinterleave "Deinterleave"
    printf 'Every tool writes its two mate files to disk under the same rules: new files per run, removed when the run ends, in an 8 GB memory cap so the written pages stay in the page cache.\n\n'
    readme_scaling deinterleave "###"
    readme_repro
}

bench_main "$@"

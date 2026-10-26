#!/usr/bin/env bash
# Benchmark: validating paired reads into two output files, plain or gzip. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "validating paired reads into two output files, plain or gzip"
truth "the corpus holds every pair as R1 and R2; right outputs equal them byte for byte, after gzip decoding for gzip outputs"
expect accept
facts piped
# count and check already run z-fastq on every corpus file; every pair here is in a family.
verify off
rival any
order z-fastq,fastp
scaling paired

word title "Processing paired reads"
word subject "z-fastq process"
word others "paired processors"
word timed "its outputs were byte-exact"
word right "byte-exact outputs"
word wrong "wrong output"
word missed "gave wrong output"
word unseen "the pairs are broken"
word legend "wrong output, exit 0"
word verb "processes"
word past "processed"
word payer "a processor"
word same "processing"
word many "Processing"
word truth "The corpus holds every pair as R1 and R2, cut from the same records, so right outputs equal R1 and R2 byte for byte, after gzip decoding for gzip outputs."
word cells "A cell is the tool's verdict: \`accept\` (exit 0, outputs byte-exact), \`reject P001\` (with the rule it named), \`reject\`, \`wrong: accept\` (exit 0 on a broken pair or record), \`wrong: bytes differ\` (exit 0 on a valid input, other bytes), \`refused\` (rejected a valid input), or \`crash\`."

# Both families write two real files per run, removed when the run ends; the
# gzip family names them .gz, which is how fastp chooses to compress.
OUTPUTS=outputs=_1.fq,_2.fq,_1.fq.gz,_2.fq.gz,.json,.html
lane z-fastq subject "{zfastq} process --out1 {stem}_1.fq --out2 {stem}_2.fq {in}" "$OUTPUTS"
lane fastp partial "{bin} -i {in1} -I {in2} -o {stem}_1.fq -O {stem}_2.fq -A -G -Q -L -w 1 -j {stem}.json -h {stem}.html" "$OUTPUTS" note="trimming and filters off, one worker (-w 1), compression level 4"

family plain "two mate files into two plain files"
lanes all
look full
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

family gzip "two mate files into two gzip files"
lanes all
look full
command z-fastq "{zfastq} process --gzip --out1 {stem}_1.fq.gz --out2 {stem}_2.fq.gz {in}"
command fastp "{bin} -i {in1} -I {in2} -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L -w 1 -j {stem}.json -h {stem}.html"
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

edges process plain

# --- edge cases ---------------------------------------------------------------
# Small pairs as R1 and R2, run through the plain family's command.
PROCESS_EDGES="$EDGES/process-pairs"
feature_edge_list() {
    cat <<'EOF'
p_ok	accept	valid pair
p_plus	accept	annotated plus line
p_crlf	policy	CRLF lines
p_names	reject P001	mate names differ
p_short	reject P002	R2 one record short
p_lenmis	reject S005	length mismatch in R1
p_trunc	reject S004	truncated R1 record
p_dash	reject S001	`-` instead of `+` in R1
p_qlow	reject S006	quality byte below 33 in R1
EOF
}
edge_family() { echo plain; }
edge_input() { printf '%s_R1.fastq %s_R2.fastq\n' "$PROCESS_EDGES/${1#p_}" "$PROCESS_EDGES/${1#p_}"; }
# Writes NAME_R1.fastq and NAME_R2.fastq from two record lists (one record per
# argument, fields separated by |).
pair_files() { # name r1-records r2-records
    local name="$PROCESS_EDGES/$1" r
    : >"${name}_R1.fastq"; for r in $2; do tr '|' '\n' <<<"$r" >>"${name}_R1.fastq"; done
    : >"${name}_R2.fastq"; for r in $3; do tr '|' '\n' <<<"$r" >>"${name}_R2.fastq"; done
}
make_feature_edges() {
    [[ -f "$PROCESS_EDGES/qlow_R2.fastq" ]] && return 0
    mkdir -p -- "$PROCESS_EDGES"
    local a1='@p1/1|ACGT|+|IIII' a2='@p2/1|GGCA|+|IIII' b1='@p1/2|TTGA|+|IIII' b2='@p2/2|CCAT|+|IIII'
    pair_files ok "$a1 $a2" "$b1 $b2"
    pair_files plus '@p1/1|ACGT|+p1/1|IIII @p2/1|GGCA|+p2/1|IIII' '@p1/2|TTGA|+p1/2|IIII @p2/2|CCAT|+p2/2|IIII'
    pair_files names "$a1 $a2" "@q1/2|TTGA|+|IIII $b2"
    pair_files short "$a1 $a2" "$b1"
    pair_files lenmis "@p1/1|ACGT|+|III $a2" "$b1 $b2"
    pair_files trunc "$a1 @p2/1|GGCA" "$b1 $b2"
    pair_files dash "@p1/1|ACGT|-|IIII $a2" "$b1 $b2"
    # A quality byte below 33 (0x1f) in R1's first record.
    pair_files qlow "@p1/1|ACGT|+|IIXI $a2" "$b1 $b2"
    sed -i '4s/X/\x1f/' "$PROCESS_EDGES/qlow_R1.fastq"
    local f
    for f in R1 R2; do sed 's/$/\r/' "$PROCESS_EDGES/ok_$f.fastq" >"$PROCESS_EDGES/crlf_$f.fastq"; done
}

# --- truth --------------------------------------------------------------------
# The plain corpus file a path stands for: the copy beside a .gz.
plain_path() { printf '%s\n' "${1%.gz}"; }
# The verdict: both outputs, decoded when they are gzip, against R1 and R2.
fact() {
    if [[ "${RC:-0}" != 0 ]]; then
        reject_code
        return
    fi
    local -a in
    read -r -a in <<<"$IN"
    local r1 r2 o1="$OUT1" o2="$OUT2"
    r1="$(plain_path "${in[0]}")" r2="$(plain_path "${in[1]}")"
    [[ -f "$o1" ]] || o1="$OUT1.gz" o2="$OUT2.gz"
    if [[ "$o1" == *.gz ]]; then
        if cmp -s <(gzip -dc -- "$o1") "$r1" && cmp -s <(gzip -dc -- "$o2") "$r2"; then printf 'accept\n'; else printf 'wrong: bytes differ\n'; fi
    elif cmp -s -- "$o1" "$r1" && cmp -s -- "$o2" "$r2"; then
        printf 'accept\n'
    else
        printf 'wrong: bytes differ\n'
    fi
}

bench_main "$@"

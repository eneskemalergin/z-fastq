#!/usr/bin/env bash
# Benchmark: processing paired reads: copying into plain or gzip files, adapter and quality trimming, and QC tables, on one CPU, then against the peers on 2 to 16 threads. Modes: init, check, release, report, peer, ab REV.
BENCH_JOBS="${BENCH_JOBS:-8}"
BENCH_MEM="${BENCH_MEM:-16G}"
source "$(dirname -- "$0")/../lib/bench.sh"

job "processing paired reads: copying, trimming, and QC tables"
execution bounded
truth "copies equal the corpus pairs byte for byte; trimmed pairs and the QC summary equal fastp 1.4.0's with the same settings; QC tables equal FastQC-Rust's"
expect accept
facts piped
# count and check already run z-fastq on every corpus file; every pair here is in a family.
verify off
rival any
order "z-fastq,fastp$(printf ',fastp-w%s' 2 4 8 16),fastqc-rust$(printf ',fastqc-rust-t%s' 2 4 8 16),fastqc$(printf ',fastqc-t%s' 2 4 8 16)"
scaling paired

word title "Processing paired reads"
word subject "z-fastq process"
word others "tools"
word timed "its outputs were right"
word right "right outputs"
word wrong "wrong output"
word missed "gave wrong output"
word unseen "the pairs are broken"
word legend "wrong output, exit 0"
word verb "processes"
word past "processed"
word payer "a processor"
word same "processing"
word many "Processing"
word truth "Copies must equal the corpus pairs byte for byte. Trimmed pairs must hold the records fastp 1.4.0 keeps with the same settings, and every field of z-fastq's QC summary must equal fastp's (its own schema and tool names aside); fastp itself is that reference. QC tables must equal FastQC-Rust's \`fastqc_data.txt\` without its version and file name lines, tied overrepresented rows in sorted order."
word cells "A cell is the tool's verdict: \`accept\` (exit 0, outputs right), \`reject P001\` (with the rule it named), \`reject\`, \`wrong: accept\` (exit 0 on a broken pair or record), \`wrong: records differ\` or \`wrong: table R1\` (exit 0, other output), \`no table\` (no table written for a mate), \`refused\` (rejected a valid input), \`crash\`, or \`-\` (the tool has no such job)."

# Adapters: TruSeq for rnaFusion's settings, Nextera for FRAMEFusion's.
RF1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
RF2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
FF=CTGTCTCTTATACACATCT
# fastp's default filters, written out on both sides.
ZQ_FILTERS="--low-quality 15 --max-low-quality-percent 40 --max-n 5"
FP_FILTERS="-q 15 -u 40 -n 5"
FP_IO="{bin} -i {in1} -I {in2} -o {stem}_1.fq -O {stem}_2.fq"
FP_REPORTS="-w 1 -j {stem}.json -h {stem}.html"
ZQ_OUT="--out1 {stem}_1.fq --out2 {stem}_2.fq {in}"
ZQ_COPY="{zfastq} process $ZQ_OUT"
FP_COPY="$FP_IO -A -G -Q -L $FP_REPORTS"
ZQ_GZIP="{zfastq} process --gzip --out1 {stem}_1.fq.gz --out2 {stem}_2.fq.gz {in}"
FP_GZIP="{bin} -i {in1} -I {in2} -o {stem}_1.fq.gz -O {stem}_2.fq.gz -A -G -Q -L $FP_REPORTS"
ZQ_TRUSEQ="{zfastq} process --adapter1 $RF1 --adapter2 $RF2 $ZQ_FILTERS --min-length 15 --pair-duplication --qc-summary {stem}.json $ZQ_OUT"
FP_TRUSEQ="$FP_IO --adapter_sequence $RF1 --adapter_sequence_r2 $RF2 --disable_trim_poly_g $FP_FILTERS -l 15 $FP_REPORTS"
ZQ_NEXTERA="{zfastq} process --adapter1 $FF --adapter2 $FF --trim-poly-g 10 --trim-tail-window 4 --trim-tail-quality 20 $ZQ_FILTERS --min-length 30 --pair-duplication --qc-summary {stem}.json $ZQ_OUT"
FP_NEXTERA="$FP_IO --adapter_sequence $FF --adapter_sequence_r2 $FF --trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 $FP_FILTERS -l 30 $FP_REPORTS"
# z-fastq ships no contaminant list; its QC lane names the one FastQC 0.13.0
# ships, which both FastQC tools use by default.
CONTAMINANTS="$REPO/tools/.local/fastqc-current/runtime/Configuration/contaminant_list.txt"
ZQ_QC="{zfastq} process --qc-tables {stem}.qc --contaminants $CONTAMINANTS $ZQ_OUT"
FR_QC="{bin} --threads 1 -o {stem}.qc {in}"
FJ_QC="{bin} --threads 1 --quiet --outdir {stem}.qc {in}"

# Every lane writes new files under its stem, removed after each run: reads,
# fastp's JSON and HTML, and QC files inside the {stem}.qc directory.
OUTPUTS=outputs=_1.fq,_2.fq,_1.fq.gz,_2.fq.gz,.json,.html
lane z-fastq subject "$ZQ_COPY" "$OUTPUTS" outdir=.qc
lane fastp partial "$FP_COPY" "$OUTPUTS" note="writes its JSON and HTML report and measures duplication on every run; compression level 4"
lane fastqc-rust partial "$FR_QC" outdir=.qc note="QC tables only, no reads written and no pair checks; the QC tables' reference"
lane fastqc partial "$FJ_QC" outdir=.qc note="QC tables only, no reads written and no pair checks; a Java program whose JVM starts its own threads"
# Each peer again at 2, 4, 8, and 16 threads, on as many physical cores, in
# the same batches as z-fastq on one: at the 1 GiB inputs only.
THREADS="2 4 8 16"
for n in $THREADS; do
    lane "fastp-w$n" partial "${FP_COPY/-w 1/-w $n}" "$OUTPUTS" cores="$n" of=fastp headline-only
    lane "fastqc-rust-t$n" partial "${FR_QC/--threads 1/--threads $n}" outdir=.qc cores="$n" of=fastqc-rust headline-only
    lane "fastqc-t$n" partial "${FJ_QC/--threads 1/--threads $n}" outdir=.qc cores="$n" of=fastqc headline-only
done
FASTP_LANES="z-fastq,fastp$(printf ',fastp-w%s' $THREADS)"
QC_LANES="z-fastq,fastqc-rust$(printf ',fastqc-rust-t%s' $THREADS),fastqc$(printf ',fastqc-t%s' $THREADS)"
# fastp's command in the current family, at one worker and at each count.
fastp_command() {
    local n
    command fastp "$1"
    for n in $THREADS; do command "fastp-w$n" "${1/-w 1/-w $n}"; done
}

family trim-truseq "TruSeq adapters, quality and N filters, minimum length 15 (rnaFusion's settings)"
lanes "$FASTP_LANES"
look full
command z-fastq "$ZQ_TRUSEQ"
fastp_command "$FP_TRUSEQ"
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

family trim-nextera "Nextera adapters, poly-G and quality-tail trimming, filters, minimum length 30 (FRAMEFusion's settings)"
lanes "$FASTP_LANES"
look full
command z-fastq "$ZQ_NEXTERA"
fastp_command "$FP_NEXTERA"
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

family qc-tables "FastQC tables for both mates, no trimming"
lanes "$QC_LANES"
look full
command z-fastq "$ZQ_QC"
input pair "$PAIRED_SOURCES" medium
input pair paired tiny,xsmall,small,medium

# Copying measures what process costs before it does any work: reading,
# checking, and writing a pair. One pair at one thread answers that; the
# trimming and QC jobs carry the comparison across inputs and threads.
COPY_PAIR=illumina-drr1001190
COPY_LANES="z-fastq,fastp"

family copy-plain "two mate files into two plain files"
lanes "$COPY_LANES"
look full
input pair "$COPY_PAIR" medium

family copy-gzip "two mate files into two gzip files"
lanes "$COPY_LANES"
look full
command z-fastq "$ZQ_GZIP"
command fastp "$FP_GZIP"
input pair "$COPY_PAIR" medium

# The edge cases run each lane with its copy command, the QC tools with their
# tables; trimming and QC come first, so they name the sources and shapes.
edges process copy-plain qc-tables
# The report's sections: the cost of each step on one pair, then related jobs
# side by side in one figure per view.
SECTIONS="steps:What each step costs:steps=copy-plain/Copy to plain files,copy-gzip/Copy to gzip files,trim-truseq/TruSeq trim,trim-nextera/Nextera trim|trimming:Trimming:pooled=trim-truseq/TruSeq,trim-nextera/Nextera|qc:QC tables=qc-tables/FastQC tables"

# FastQC's name for an input: the file name without .gz, .bz2, .txt, .fastq,
# .fq, .csfastq, .sam, .bam, and .ubam, stripped in that order.
fastqc_stem() {
    local n="${1##*/}" s
    for s in .gz .bz2 .txt .fastq .fq .csfastq .sam .bam .ubam; do n="${n%"$s"}"; done
    printf '%s\n' "$n"
}
# What each lane writes into its {stem}.qc directory: z-fastq one folder per
# mate holding fastqc_data.txt, the FastQC tools an HTML page and a zip.
outdir_entries() { # lane family paths
    local p s
    for p in $3; do
        s="$(fastqc_stem "$p")"
        case "$1" in
            z-fastq) printf '%s_fastqc/fastqc_data.txt\n' "$s" ;;
            *) printf '%s_fastqc.html\n%s_fastqc.zip\n' "$s" "$s" ;;
        esac
    done
    [[ "$1" != z-fastq ]] || for p in $3; do printf '%s_fastqc/\n' "$(fastqc_stem "$p")"; done
    return 0
}

# --- edge cases ---------------------------------------------------------------
# Small pairs as R1 and R2. Rows run the command of one family: copy (p_),
# TruSeq trimming (t_), Nextera trimming (n_), or QC tables (q_). Trimming rows
# are the fastp-verified probe pairs of tests/data/synthetic/trim_*.fastq,
# grouped by what they probe.
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
t_overlap	accept	mates overlap past the insert (TruSeq)
t_adapter	accept	adapter by sequence, partial or with errors (TruSeq)
t_gap	accept	adapter after an insertion or deletion (TruSeq)
t_dimer	accept	adapter dimers and very short pairs (TruSeq)
n_overlap	accept	mates overlap past the insert (Nextera)
n_adapter	accept	adapter by sequence, gaps included (Nextera)
n_polyg	accept	poly-G tails (Nextera)
n_tail	accept	low-quality tails and N runs (Nextera)
q_ok	accept	valid pair (QC)
q_empty	accept	empty mate files (QC)
q_short	accept	reads of 12 bases or fewer (QC)
q_long	policy	a read longer than 1024 bases (QC)
q_names	reject P001	mate names differ (QC)
q_trunc	reject S004	truncated R1 record (QC)
EOF
}
edge_family() {
    case "$1" in
        p_*) echo copy-plain ;;
        t_*) echo trim-truseq ;;
        n_*) echo trim-nextera ;;
        *) echo qc-tables ;;
    esac
}
edge_input() {
    local base="$PROCESS_EDGES/$1"
    [[ "$1" != p_* ]] || base="$PROCESS_EDGES/${1#p_}"
    printf '%s_R1.fastq %s_R2.fastq\n' "$base" "$base"
}
# Writes NAME_R1.fastq and NAME_R2.fastq from two record lists (one record per
# argument, fields separated by |).
pair_files() { # name r1-records r2-records
    local name="$PROCESS_EDGES/$1" r
    : >"${name}_R1.fastq"; for r in $2; do tr '|' '\n' <<<"$r" >>"${name}_R1.fastq"; done
    : >"${name}_R2.fastq"; for r in $3; do tr '|' '\n' <<<"$r" >>"${name}_R2.fastq"; done
}
# The probe pairs whose names match PATTERN, from one fixture set (rf or ff).
probe_pairs() { # name set pattern
    local m
    for m in 1 2; do
        awk -v re="$3" 'NR % 4 == 1 { keep = ($0 ~ re) } keep' "$REPO/tests/data/synthetic/trim_$2_$m.fastq" >"$PROCESS_EDGES/$1_R$m.fastq"
    done
}
make_feature_edges() {
    [[ -f "$PROCESS_EDGES/q_trunc_R2.fastq" ]] && return 0
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
    probe_pairs t_overlap rf '^@ov_'
    probe_pairs t_adapter rf '^@ad_'
    probe_pairs t_gap rf '^@gap_'
    probe_pairs t_dimer rf '^@(dimer|short)_'
    probe_pairs n_overlap ff '^@ov_'
    probe_pairs n_adapter ff '^@(ad|gap)_'
    probe_pairs n_polyg ff '^@pg_'
    probe_pairs n_tail ff '^@qt_'
    local c1='@c1/1|ACGTACGTACGTACGTACGTACGT|+|IIIIIIIIIIIIIIIIIIIIIIII' c2='@c1/2|TTGATTGATTGATTGATTGATTGA|+|IIIIIIIIIIIIIIIIIIIIIIII'
    pair_files q_ok "$c1 $a2" "$c2 $b2"
    pair_files q_empty "" ""
    pair_files q_short '@s1/1|ACGTACGTACGT|+|IIIIIIIIIIII @s2/1|ACGTAC|+|IIIIII' '@s1/2|TTGATTGATTGA|+|IIIIIIIIIIII @s2/2|TTGATT|+|IIIIII'
    local long qual
    long="$(head -c 1500 /dev/zero | tr '\0' A)" qual="$(head -c 1500 /dev/zero | tr '\0' I)"
    pair_files q_long "@l1/1|$long|+|$qual $a2" "@l1/2|TTGA|+|IIII $b2"
    pair_files q_names "$c1 $a2" "@q1/2|TTGA|+|IIII $b2"
    pair_files q_trunc "$c1 @p2/1|GGCA" "$c2 $b2"
}

# --- truth --------------------------------------------------------------------
# The plain corpus file a path stands for: the copy beside a .gz.
plain_path() { printf '%s\n' "${1%.gz}"; }
decoded() { if [[ "$1" == *.gz ]]; then gzip -dc -- "$1"; else cat -- "$1"; fi; }
# A FASTQ file's decoded content as cksum prints it (CRC, bytes) and its
# record count; with SORTED, also the CRC of its records in sorted order (the
# same records in any order give the same value).
read_sums() { # file [sorted] -> crc bytes records [sorted-crc]
    local crc bytes lines sum sorted=""
    sum="$(decoded "$1" | cksum)" || return
    read -r crc bytes <<<"$sum"
    lines="$(decoded "$1" | wc -l)" || return
    if [[ -n "${2:-}" ]]; then
        sum="$(decoded "$1" | paste - - - - | sort -S 256M --parallel=1 -T /var/tmp | cksum)" || return
        read -r sorted _ <<<"$sum"
    fi
    printf '%s %s %s%s\n' "$crc" "$bytes" "$((lines / 4))" "${sorted:+ $sorted}"
}
# A table as compared: no version or file name line, tied overrepresented
# rows (same count) in sorted order.
table_norm() {
    awk -F'\t' 'NR == 1 || /^Filename\t/ { next }
        /^>>Overrepresented sequences/ { print; ov = 1; next }
        ov && /^#/ { print; next }
        ov && /^>>END_MODULE/ { ov = 0; close("sort -t \"\t\" -k2,2nr -k1,1"); print; next }
        ov { print | "sort -t \"\t\" -k2,2nr -k1,1"; next }
        { print }'
}
# A lane's normalized table for one input from its QC directory; fails when
# it wrote none.
table_of() { # lane dir input
    local s
    s="$(fastqc_stem "$3")"
    case "$1" in
        z-fastq) [[ -f "$2/${s}_fastqc/fastqc_data.txt" ]] && table_norm <"$2/${s}_fastqc/fastqc_data.txt" ;;
        *) [[ -f "$2/${s}_fastqc.zip" ]] && unzip -p "$2/${s}_fastqc.zip" "${s}_fastqc/fastqc_data.txt" | table_norm ;;
    esac
}
# A fastp-layout JSON file as sorted "field<TAB>value" lines, without the
# command line (it names the files, not the reads) and z-fastq's own schema,
# tool, and byte_strings fields. A fraction is written as fastp prints it:
# six significant digits, the C++ stream default; integers stay exact.
json_leaves() { # file
    jq -r 'del(.command, .schema, .tool, .byte_strings) | paths(scalars) as $p
        | [($p | map(tostring) | join(".")), (getpath($p) | tostring), (getpath($p) | type)] | @tsv' "$1" |
        awk -F'\t' -v OFS='\t' '$3 == "number" && $2 ~ /[.eE]/ { $2 = sprintf("%.6g", $2) } { print $1, $2 }' | sort
}
# The first field GOT writes that fastp's JSON (WANT) lacks or holds another
# value in; nothing when every one matches. fastp's fields that z-fastq does
# not write (the wiki's report reference lists them) are not compared.
json_diff() { # got want
    local got want diff
    got="$(json_leaves "$1" 2>/dev/null)" && [[ -n "$got" ]] || { printf 'unreadable\n'; return; }
    # No pipe to head: head would close early and comm would die with SIGPIPE
    # (141), which fails the checker under set -e with pipefail. Consume comm
    # fully, then take the first line's first field in the shell.
    want="$(json_leaves "$2")" || return
    diff="$(comm -23 <(printf '%s\n' "$got") <(printf '%s\n' "$want"))" || return
    if [[ -z "$diff" ]]; then return 0; fi
    diff="${diff%%$'\n'*}"
    printf '%s\n' "${diff%%$'\t'*}"
}

# The reference: fastp 1.4.0 at one worker on every trimming input, and
# FastQC-Rust on every QC input, run once per plain input pair and cached in
# results/reference.tsv while both binaries, the commands, and the inputs are
# unchanged. A gzip input shares its plain copy's reference.
REFERENCE="$RESULTS/reference.tsv"
REF_DIR="$RESULTS/reference"
ref_name() { local a="${1#"$BENCH"/data/}"; printf '%s\n' "${a//\//_}"; } # path -> file-name-safe id
reference_jobs() {
    cases | awk -F'\t' '$1 == "trim-truseq" || $1 == "trim-nextera" || $1 == "qc-tables" { print $1 "\t" $4 }' |
        while IFS=$'\t' read -r fam paths; do
            local -a p
            read -r -a p <<<"$paths"
            printf '%s\t%s %s\n' "$fam" "$(plain_path "${p[0]}")" "$(plain_path "${p[1]}")"
        done | sort -u
}
# One reference line: family, the two plain inputs, then for trimming the
# fastp outputs' read_sums per mate, mate 2 reads holding S, W, or D, and
# the kept JSON; for QC the two kept tables.
reference_one() {
    local fam paths stem dir
    IFS=$'\t' read -r fam paths <<<"$1"
    local -a p
    read -r -a p <<<"$paths"
    dir="$(scratch_dir "/var/tmp/bench-$FEATURE-ref.XXXXXX")" stem="$dir/o"
    bench_cleanup_path f "$dir/err" "$stem"_1.fq "$stem"_2.fq "$stem.json" "$stem.html"
    if [[ "$fam" == qc-tables ]]; then
        mkdir -- "$stem.qc"
        bench_cleanup_outdir "$stem.qc" fastqc-rust "$fam" "$paths"
        $(STEM="$stem" lane_cmd fastqc-rust "$paths" qc-tables) >/dev/null 2>"$dir/err" || { log reference FAIL "fastqc-rust on ${p[0]}: $(head -c 300 "$dir/err")"; rm -f -- "$dir/err"; clear_outdir "$stem.qc" fastqc-rust qc-tables "$paths"; rmdir -- "$dir"; return 1; }
        local m out=""
        for m in 0 1; do
            table_of fastqc-rust "$stem.qc" "${p[$m]}" >"$REF_DIR/$(ref_name "${p[$m]}").txt"
            out+=$'\t'"$REF_DIR/$(ref_name "${p[$m]}").txt"
        done
        clear_outdir "$stem.qc" fastqc-rust qc-tables "$paths"
        printf '%s\t%s\t%s%s\n' "$fam" "${p[0]}" "${p[1]}" "$out"
    else
        $(STEM="$stem" lane_cmd fastp "$paths" "$fam") >/dev/null 2>"$dir/err" || { log reference FAIL "fastp on ${p[0]}: $(head -c 300 "$dir/err")"; rm -f -- "$dir/err" "$stem"_1.fq "$stem"_2.fq "$stem.json" "$stem.html"; rmdir -- "$dir"; return 1; }
        local s1 s2 swd json="$REF_DIR/$fam-$(ref_name "${p[0]}").json"
        s1="$(read_sums "${stem}_1.fq" sorted)" s2="$(read_sums "${stem}_2.fq" sorted)"
        swd="$(awk 'NR % 4 == 2 && /[SWDswd]/ { n++ } END { print n + 0 }' "${p[1]}")"
        cp -- "$stem.json" "$json"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$fam" "${p[0]}" "${p[1]}" "$s1" "$s2" "$swd" "$json"
    fi
    rm -f -- "$dir/err" "$stem"_1.fq "$stem"_2.fq "$stem.json" "$stem.html"
    rmdir -- "$dir"
}
feature_prepare() {
    local needed available key
    needed="$(cases | awk -F'\t' '$6 > n { n = $6 } END { printf "%.0f", n * 6 + 1073741824 }')"
    available="$(df -PB1 /var/tmp | awk 'NR == 2 { print $4 }')"
    ((available >= needed)) || die "process verification needs at least $needed free bytes under /var/tmp; $available available"
    key="$(lane_id fastp trim-truseq; lane_id fastp trim-nextera; lane_id fastqc-rust qc-tables; reference_jobs | cut -f2 | tr ' ' '\n' | sort -u | TZ=UTC xargs stat -L -c "%n $FILE_ID_FORMAT"
        declare -f reference_one read_sums table_of table_norm fastqc_stem | cksum)"
    # New inputs or commands write separate sidecars. A failed refresh must
    # leave the JSON and tables named by the saved reference.tsv unchanged.
    REF_DIR="$RESULTS/reference/$(printf '%s' "$key" | sha256sum | cut -d' ' -f1)"
    mkdir -p -- "$REF_DIR"
    cached_step reference "$key" reference_jobs reference_one "$REFERENCE"
    awk -F'\t' '($1 ~ /^trim-/) && $6 > 0 { printf "%s: %d mate 2 reads with S, W, or D; fastp 1.4.0 is not a reference there\n", $2, $6 }' "$REFERENCE" |
        while IFS= read -r line; do log reference WARN "$line"; done
    return 0
}
# The reference line of FAMILY for the plain copies of PATHS; fails when none.
reference_of() { # family paths
    local -a p
    read -r -a p <<<"$2"
    awk -F'\t' -v f="$1" -v a="$(plain_path "${p[0]}")" -v b="$(plain_path "${p[1]}")" '$1 == f && $2 == a && $3 == b { print; found = 1; exit } END { exit !found }' "$REFERENCE"
}

# Facts: the verdict, then out_bytes (the bytes the lane wrote).
fact() {
    if [[ "${RC:-0}" != 0 ]]; then
        reject_code
        return
    fi
    case "${FAMILY%-threads}" in
        copy-*) copy_fact ;;
        trim-*) trim_fact ;;
        qc-tables) qc_fact ;;
        *) printf 'wrong: no fact for %s\n' "$FAMILY" ;;
    esac
}
# The lane's two read outputs, decoded or not, and their size on disk.
read_outputs() { O1="$OUT1" O2="$OUT2"; [[ -f "$O1" ]] || O1="$OUT1.gz" O2="$OUT2.gz"; }
out_bytes() { printf 'out_bytes %s\n' "$(stat -c %s -- "$@" 2>/dev/null | awk '{ s += $1 } END { print s + 0 }')"; }
copy_fact() {
    local -a in
    read -r -a in <<<"$IN"
    read_outputs
    if cmp -s <(decoded "$O1") "$(plain_path "${in[0]}")" && cmp -s <(decoded "$O2") "$(plain_path "${in[1]}")"; then printf 'accept\n'
    elif cmp -s <(decoded "$O1") <(tr -d '\r' <"$(plain_path "${in[0]}")") && cmp -s <(decoded "$O2") <(tr -d '\r' <"$(plain_path "${in[1]}")"); then printf 'CRLF written as LF\n'
    else printf 'wrong: bytes differ\n'; fi
    out_bytes "$O1" "$O2"
}
# Trimmed reads against fastp's records, then the lane's JSON against fastp's.
trim_fact() {
    local fam="${FAMILY%-threads}" ref w1 w2 swd json g1 g2 dir="" json_got
    read_outputs
    json_got="${OUT1%_1.fq}.json"
    if [[ -n "${EDGE:-}" ]]; then
        # Small edge inputs: fastp runs here, with the family's command.
        dir="$(scratch_dir "/var/tmp/bench-$FEATURE-edge.XXXXXX")"
        bench_cleanup_path f "$dir/o_1.fq" "$dir/o_2.fq" "$dir/o.json" "$dir/o.html"
        $(STEM="$dir/o" lane_cmd fastp "$IN" "$fam") >/dev/null 2>&1 || { printf 'wrong: no reference\n'; rm -f -- "$dir"/*; rmdir -- "$dir"; return; }
        w1="$(read_sums "$dir/o_1.fq" sorted)" w2="$(read_sums "$dir/o_2.fq" sorted)" swd=0 json="$dir/o.json"
    else
        ref="$(reference_of "$fam" "$IN")" || { printf 'wrong: no reference\n'; return; }
        IFS=$'\t' read -r _ _ _ w1 w2 swd json <<<"$ref"
    fi
    g1="$(read_sums "$O1")" g2="$(read_sums "$O2")"
    local -a w=($w1 $w2)
    if ((swd > 0)); then
        printf 'wrong: no reference (S, W, D in R2)\n'
    elif [[ "$g1 $g2" == "${w[0]} ${w[1]} ${w[2]} ${w[4]} ${w[5]} ${w[6]}" ]]; then
        local d="file missing"
        [[ ! -f "$json_got" ]] || d="$(json_diff "$json_got" "$json")"
        if [[ -z "$d" ]]; then
            printf 'accept\n'
        elif [[ "$LANE" == fastp* && "$d" != "file missing" ]]; then
            # fastp is the reference at one worker. With more, its records
            # stay the same but its own report does not (adapter sequence
            # counts, the insert size histogram), so that is a fact, not a
            # wrong output.
            printf 'accept\nsummary_differs %s\n' "$d"
        else
            printf 'wrong: summary %s\n' "$d"
        fi
    else
        # Sorted only now: the same records in another order are not wrong records.
        local -a g=($g1 $g2)
        if [[ "${g[1]} ${g[2]} ${g[4]} ${g[5]}" == "${w[1]} ${w[2]} ${w[5]} ${w[6]}" &&
            "$(read_sums "$O1" sorted | cut -d' ' -f4) $(read_sums "$O2" sorted | cut -d' ' -f4)" == "${w[3]} ${w[7]}" ]]; then
            printf 'wrong: same records, other order\n'
        else
            printf 'wrong: records differ\n'
        fi
    fi
    out_bytes "$O1" "$O2"
    if [[ -n "$dir" ]]; then rm -f -- "$dir"/*; rmdir -- "$dir"; fi
}
# Each mate's table against FastQC-Rust's; z-fastq's read copies against the input.
qc_fact() {
    local -a in
    read -r -a in <<<"$IN"
    local m want=() dir="" verdict=accept
    if [[ -n "${EDGE:-}" ]]; then
        dir="$(scratch_dir "/var/tmp/bench-$FEATURE-edge.XXXXXX")"
        mkdir -- "$dir/o.qc"
        bench_cleanup_path f "$dir/want0" "$dir/want1"
        bench_cleanup_outdir "$dir/o.qc" fastqc-rust qc-tables "$IN"
        $(STEM="$dir/o" lane_cmd fastqc-rust "$IN" qc-tables) >/dev/null 2>&1 || verdict="wrong: no reference"
        for m in 0 1; do table_of fastqc-rust "$dir/o.qc" "${in[$m]}" >"$dir/want$m" || true; want+=("$dir/want$m"); done
        clear_outdir "$dir/o.qc" fastqc-rust qc-tables "$IN"
    else
        local ref
        ref="$(reference_of qc-tables "$IN")" || { printf 'wrong: no reference\n'; return; }
        IFS=$'\t' read -r _ _ _ want[0] want[1] <<<"$ref"
    fi
    if [[ "$verdict" == accept ]]; then
        for m in 0 1; do
            local got
            if ! got="$(table_of "$LANE" "$OUTDIR" "${in[$m]}")"; then verdict="no table"; break; fi
            cmp -s <(printf '%s\n' "$got") "${want[$m]}" || { verdict="wrong: table R$((m + 1))"; break; }
        done
    fi
    if [[ "$verdict" == accept && "$LANE" == z-fastq ]]; then
        read_outputs
        cmp -s <(decoded "${in[0]}") "$O1" && cmp -s <(decoded "${in[1]}") "$O2" || verdict="wrong: copies differ"
    fi
    # A corpus input with no table is wrong; an edge cell shows it as it is.
    [[ "$verdict" != "no table" || -n "${EDGE:-}" ]] || verdict="wrong: no table"
    printf '%s\n' "$verdict"
    local copies=()
    [[ "$LANE" != z-fastq ]] || { read_outputs; copies=("$O1" "$O2"); }
    out_bytes $(find "$OUTDIR" -type f 2>/dev/null) "${copies[@]}"
    if [[ -n "$dir" ]]; then rm -f -- "$dir"/want*; rmdir -- "$dir"; fi
}

# --- report -------------------------------------------------------------------
declare -A FAMILY_TITLE=([copy-plain]="Copying into plain files" [copy-gzip]="Copying into gzip files" [trim-truseq]="Trimming, TruSeq settings"
    [trim-nextera]="Trimming, Nextera settings" [qc-tables]="QC tables")
tools_table() {
    local name fam line
    {
        printf 'Tool\tClass\tCopy to plain\tCopy to gzip\tTruSeq trim\tNextera trim\tQC tables\tThreads\tMalformed rejected\n'
        printf 'l\tl\tl\tl\tl\tl\tl\tl\tr\n'
        for name in "${LANES[@]}"; do
            [[ ! -v "LANE_OF[$name]" ]] || continue
            line="${REG_LABEL[$name]}"$'\t'"${LANE_CLASS[$name]}"
            for fam in copy-plain copy-gzip trim-truseq trim-nextera qc-tables; do
                if [[ " ${FAM_LANES[$fam]} " == *" $name "* ]]; then line+=$'\t'"$(mode_cell "$name" "$fam")"; else line+=$'\t-'; fi
            done
            case "$name" in z-fastq) line+=$'\t1' ;; fastp) line+=$'\t1 to 16 (-w)' ;; *) line+=$'\t1 to 16 (-t)' ;; esac
            printf '%s\t%s\n' "$line" "$(rejects_all "$name")"
        done
    } | mdtable
}
# One bullet per job: each other tool at one thread, its time and memory over
# z-fastq's on plain and gzip input; then its fastest thread count on plain
# input, with its time and CPU time there over z-fastq's on one core, and the
# fewest threads at which it is faster, if any. Geometric means over the
# headline inputs, each ratio within its batch.
job_bullet() { # family label
    awk -F'\t' -v fam="$1" -v lbl="$2" -v tier="$HEADLINE" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") fmt[$2] = $3; next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0; cores[$1] = $9 + 0; of[$1] = $10 } next }
        FNR == 1 || $1 != fam || !($2 in fmt) { next }
        { w[$2, $3] = $5; c[$2, $3] = $11; r[$2, $3] = $8; lanes[$3] = 1 }
        END {
            for (l in lanes) if (l != "z-fastq") for (x in fmt) {
                if (!((x, l) in w) || !((x, "z-fastq") in w) || c[x, l] <= 0 || c[x, "z-fastq"] <= 0) continue
                key = of[l] SUBSEP fmt[x] SUBSEP cores[l]; tools[of[l]] = 1
                k[key]++; sw[key] += log(w[x, l] / w[x, "z-fastq"]); sc[key] += log(c[x, l] / c[x, "z-fastq"]); sr[key] += log(r[x, l] / r[x, "z-fastq"])
            }
            n = 0; for (t in tools) list[++n] = t
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            out = ""
            for (i = 1; i <= n; i++) { t = list[i]; s = label[t]
                p1 = t SUBSEP "plain" SUBSEP 1; g1 = t SUBSEP "gzip" SUBSEP 1
                if (p1 in k) s = s sprintf(" takes %.1fx z-fastq%ss time on plain input", exp(sw[p1] / k[p1]), "\047")
                if (g1 in k) s = s sprintf(" and %.1fx on gzip", exp(sw[g1] / k[g1]))
                if (p1 in k) s = s sprintf(", with %.0fx its memory", exp(sr[p1] / k[p1]))
                best = 0; from = 0
                for (m = 2; m <= 1024; m *= 2) { q = t SUBSEP "plain" SUBSEP m; if (!(q in k)) continue
                    v = exp(sw[q] / k[q]); if (!best || v < bv) { best = m; bv = v; bc = exp(sc[q] / k[q]) }
                    if (!from && v < 1) from = m }
                if (best) s = s sprintf("; at its fastest, %d threads, %.2fx z-fastq%ss time for %.1fx its CPU time", best, bv, "\047", bc)
                if (from) s = s sprintf(" (faster than z-fastq from %d threads)", from)
                out = out (out == "" ? "" : "; ") s }
            if (out != "") printf "- **%s:** %s.\n", lbl, out
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
}
process_findings() {
    local others=0 name fam
    for name in "${LANES[@]}"; do [[ "$name" == z-fastq || -v "LANE_OF[$name]" ]] || others=$((others + 1)); done
    awk -F'\t' -v others="$others" 'FILENAME ~ /checks.tsv$/ { if ($1 == "z-fastq") { if ($3 == "pass") zp++; else zf++ } next }
        FNR > 1 && $6 != "policy" { if ($4 == "z-fastq") { if ($6 != "right") zbad++ } else if ($3 ~ /^reject/ && ($6 == "wrong" || $6 == "crash")) bad[$4] = 1 }
        END { line = (zp + zf == 0) ? "- **Right answers:** z-fastq has no current validation results" : zf ? sprintf("- **Right answers:** z-fastq gave wrong output on %d of %d checked inputs", zf, zp + zf) : sprintf("- **Right answers:** z-fastq wrote the right copies, fastp'"'"'s records, and FastQC'"'"'s tables on all %d checked inputs", zp)
              line = line (zbad ? "." : " and rejected every broken pair and record.")
              k = 0; for (b in bad) k++
              if (k) line = line sprintf(" %d of the %d other tools write output from a broken pair or record with exit 0.", k, others)
              print line }' "$REPORT_DIR/checks.tsv" "$REPORT_DIR/edges.tsv"
    valid_wrong_bullet
    steps_text bullet
    gzip_size_bullet
    for fam in trim-truseq trim-nextera qc-tables; do job_bullet "$fam" "${FAMILY_TITLE[$fam]}"; done
}
# The gzip files z-fastq and fastp wrote in the gzip copy, per plain headline
# input: one bullet, since a faster preset that writes more bytes is a cost.
gzip_size_bullet() {
    awk -F'\t' -v tier="$HEADLINE" -v size="${TIER_SIZE[$HEADLINE]}" 'FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == "copy-gzip" && $5 == tier && $3 == "plain") want[$2] = 1; next }
        $1 == "copy-gzip" && ($5 in want) && $3 == "pass" && ($2 == "z-fastq" || $2 == "fastp") { n = split($4, f, ";"); for (i = 1; i <= n; i++) if (f[i] ~ /^out_bytes /) { split(f[i], q, " "); s[$2] += q[2]; k[$2]++ } }
        END { if (k["z-fastq"] && k["fastp"]) { z = s["z-fastq"] / k["z-fastq"]; p = s["fastp"] / k["fastp"]
              printf "- **Gzip output size:** on that pair, z-fastq'"'"'s default `fast` preset writes %.1f MB, %.0f%% %s than fastp'"'"'s %.1f MB at its default level 4; `--gzip-preset even` or `dense` writes less, more slowly.\n", z / 1e6, (z > p ? z / p - 1 : 1 - z / p) * 100, (z > p ? "more" : "less"), p / 1e6 } }' \
        "$REPORT_DIR/cases.tsv" "$FAMILY_CHECKS"
}
# One report section: the jobs SLUG draws together.
readme_section() { # slug title family:job...
    local slug="$1" title="$2" spec fam jobs=""
    shift 2
    for spec; do fam="${spec%%:*}"; jobs+="${spec#*:}: ${FAM_DESC[$fam]}. "; done
    printf '## %s\n> <sub>%sThe real runs and the synthetic pair at %s, plain and gzip input, in the page cache, every tool on one CPU. The other tools at 2 to 16 threads are in the threads part below.</sub>\n\n' "$title" "$jobs" "${TIER_SIZE[$HEADLINE]}"
    picture "$slug-overview" "$title: time and peak RSS relative to z-fastq" && printf '\n'
    printf 'Geometric means over the %s inputs: z-fastq'"'"'s own time and memory, and every other tool'"'"'s as a multiple of them.\n\n' "${TIER_SIZE[$HEADLINE]}"
    MODE_HEADER=Job summary_table "$@"
    printf '\n'
    picture "$slug-speed" "$title: wall time by kind of input" && printf '\n'
    picture "$slug-memory" "$title: peak memory by kind of input" && printf '\n'
    for spec; do output_sentence "${spec%%:*}" "${spec#*:}"; done
    threads_part "$slug" "$title" "$@"
    details "Scaling"
    picture "$slug-scaling" "$title: throughput and peak memory by input size" && printf '\n'
    # Jobs whose scaling sentences match say them once, under all their names.
    local names="" text first="" same=1
    for spec; do
        text="$(scaling_sentences "${spec%%:*}")"
        names+="${names:+ and }${spec#*:}"
        if [[ -z "$first" ]]; then first="$text"; elif [[ "$text" != "$first" ]]; then same=0; fi
    done
    if ((same)); then printf '**%s.** %s\n' "$names" "$first"
    else for spec; do printf '**%s.** %s\n' "${spec#*:}" "$(scaling_sentences "${spec%%:*}")"; done; fi
    printf '</details>\n\n'
}
# The cost of each step on the one pair copying runs on: the copies, then the
# trims on the same input, z-fastq beside fastp at one thread.
STEP_FAMILIES="copy-plain copy-gzip trim-truseq trim-nextera"
declare -A STEP_NAME=([copy-plain]="Copy to plain files" [copy-gzip]="Copy to gzip files" [trim-truseq]="TruSeq trim" [trim-nextera]="Nextera trim")
steps_part() {
    printf '## What each step costs\n> <sub>Copying reads both mates, checks every pair, and writes them, with nothing else: what every job of `process` pays before its own work. One pair at %s, every tool on one CPU; each trim on the same pair, from its own job below.</sub>\n\n' "${TIER_SIZE[$HEADLINE]}"
    picture steps "What each step costs on one pair: wall time and peak memory" && printf '\n'
    steps_text table
    printf '\n'
    steps_text sentence
    printf '\n'
}
# MODE bullet (Findings), table, or sentence, from the medians on the pair.
steps_text() { # mode
    local names=""
    for fam in $STEP_FAMILIES; do names+="$fam=${STEP_NAME[$fam]}|"; done
    awk -F'\t' -v mode="$1" -v fams="$STEP_FAMILIES" -v names="${names%|}" -v tier="$HEADLINE" -v size="${TIER_SIZE[$HEADLINE]}" '
        function secs(ns) { return ns >= 1e9 ? sprintf("%.2f s", ns / 1e9) : sprintf("%.0f ms", ns / 1e6) }
        function mem(kib) { return kib >= 1048576 ? sprintf("%.1f GiB", kib / 1048576) : kib >= 10240 ? sprintf("%.1f MiB", kib / 1024) : sprintf("%.2f MiB", kib / 1024) }
        BEGIN { nf = split(fams, fl, " "); for (i = 1; i <= nf; i++) want[fl[i]] = 1
                n = split(names, nl, "|"); for (i = 1; i <= n; i++) { split(nl[i], q, "="); title[q[1]] = q[2] } }
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && ($1 in want) && $5 == tier && $6 != "many") { cs[$1, $4, $3] = $2; srcs[$4] = 1; lbl[$4] = $8 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
        FNR == 1 { next }
        { w[$1, $2, $3] = $5; m[$1, $2, $3] = $8 }
        END {
            # The one source every step ran.
            for (s in srcs) { ok = 1; for (i = 1; i <= nf; i++) if (!((fl[i], s, "plain") in cs)) ok = 0; if (ok) src = s }
            if (src == "") exit
            P = "plain"
            zc = w["copy-plain", cs["copy-plain", src, P], "z-fastq"]; fc = w["copy-plain", cs["copy-plain", src, P], "fastp"]
            zg = w["copy-gzip", cs["copy-gzip", src, P], "z-fastq"]; fg = w["copy-gzip", cs["copy-gzip", src, P], "fastp"]
            zt = w["trim-truseq", cs["trim-truseq", src, P], "z-fastq"]; ft = w["trim-truseq", cs["trim-truseq", src, P], "fastp"]
            zm = m["copy-plain", cs["copy-plain", src, P], "z-fastq"]; fm = m["copy-plain", cs["copy-plain", src, P], "fastp"]
            if (mode == "bullet") {
                printf "- **Reading and writing a pair:** z-fastq copies one %s pair to plain files in %s and to gzip files in %s; fastp takes %s and %s at one thread, holding %s.\n", size, secs(zc), secs(zg), secs(fc), secs(fg), mem(fm)
                exit
            }
            if (mode == "sentence") {
                split(lbl[src], q, ": "); name = q[length(q)]
                printf "On this pair (%s), copying to plain files takes z-fastq %s with %s: reading, checking, and writing are %.0f%% of its %s TruSeq trim. fastp needs %s with %s for the same copy, %.0f%% of its %s trim. Copying to gzip files instead takes z-fastq %s: compressing both mates adds %s, against %s that TruSeq trimming adds.\n", \
                    name, secs(zc), mem(zm), zc / zt * 100, secs(zt), secs(fc), mem(fm), fc / ft * 100, secs(ft), secs(zg), secs(zg - zc), secs(zt - zc)
                exit
            }
            print "Step\tTool\tPlain wall\tPlain peak memory\tGzip wall\tGzip peak memory"
            print "l\tl\tr\tr\tr\tr"
            split("z-fastq fastp", lanes, " ")
            for (i = 1; i <= nf; i++) for (j = 1; j <= 2; j++) { F = fl[i]; l = lanes[j]; line = title[F] "\t" label[l]
                for (f = 1; f <= 2; f++) { fmt = (f == 1 ? "plain" : "gzip"); c = cs[F, src, fmt]
                    line = line (((F, c, l) in w) ? "\t" secs(w[F, c, l]) "\t" mem(m[F, c, l]) : "\t-\t-") }
                print line }
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | { if [[ "$1" == table ]]; then mdtable; else cat; fi; }
}
# The other tools at 1 to 16 threads against z-fastq on one core: the figure,
# a table in absolute units, and where each tool first matches z-fastq.
threads_part() { # slug title family:job...
    local slug="$1" title="$2"
    shift 2
    printf '### %s: threads\n> <sub>Each other tool at 1, 2, 4, 8, and 16 threads on as many physical cores, z-fastq on one, in the same batches, at the %s inputs. Geometric means over those inputs. Cores used is CPU time over wall time.</sub>\n\n' "$title" "${TIER_SIZE[$HEADLINE]}"
    picture "$slug-threads" "$title: wall time, CPU time, and peak memory by thread count" && printf '\n'
    threads_table "$@"
    printf '\n'
    threads_sentences "$@"
}
threads_table() { # family:job...
    local fams
    fams="$(printf '%s|' "$@")"
    awk -F'\t' -v fams="${fams%|}" -v tier="$HEADLINE" '
        BEGIN { nf = split(fams, fl, "|"); for (i = 1; i <= nf; i++) { split(fl[i], q, ":"); fam[i] = q[1]; job[q[1]] = q[2]; at[q[1]] = 1 } }
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && ($1 in at) && $5 == tier && $6 != "many") { fmt[$1, $2] = $3; if ($3 == "gzip") hasgz = 1 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0; cores[$1] = ($9 > 0 ? $9 : 1); of[$1] = ($10 == "" ? $1 : $10) } next }
        FNR == 1 || !(($1, $2) in fmt) || $11 <= 0 { next }
        { k = $1 SUBSEP $3 SUBSEP fmt[$1, $2]; n[k]++; w[k] += log($5 / 1e9); c[k] += log($11 / 1e6); m[k] += log($8 / 1024); lanes[$1, $3] = 1 }
        END {
            print "Job\tTool\tThreads" "\tPlain wall s\tPlain CPU s\tPlain cores used\tPlain peak MiB" (hasgz ? "\tGzip wall s\tGzip CPU s\tGzip cores used\tGzip peak MiB" : "")
            print "l\tl\tr\tr\tr\tr\tr" (hasgz ? "\tr\tr\tr\tr" : "")
            split("plain gzip", fs, " ")
            for (x = 1; x <= nf; x++) { F = fam[x]; cnt = 0; delete list
                for (kk in lanes) { split(kk, q, SUBSEP); if (q[1] == F) list[++cnt] = q[2] }
                # By tool (its own registry order), then thread count.
                for (i = 1; i <= cnt; i++) for (j = i + 1; j <= cnt; j++) {
                    a = list[i]; b = list[j]
                    if (ord[of[b]] < ord[of[a]] || (of[b] == of[a] && cores[b] < cores[a])) { list[i] = b; list[j] = a } }
                for (i = 1; i <= cnt; i++) { l = list[i]; name = label[of[l]]; gsub("_", " ", name)
                    line = job[F] "\t" name "\t" cores[l]
                    for (f = 1; f <= (hasgz ? 2 : 1); f++) { k = F SUBSEP l SUBSEP fs[f]
                        if (!n[k]) { line = line "\t-\t-\t-\t-"; continue }
                        ww = exp(w[k] / n[k]); cc = exp(c[k] / n[k])
                        line = line sprintf("\t%.2f\t%.2f\t%.1f\t%.1f", ww, cc, cc / ww, exp(m[k] / n[k])) }
                    print line } }
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# Per job and tool: the fewest threads at which it matches z-fastq's wall
# time on one core, per input format, from the ratios within each batch.
threads_sentences() { # family:job...
    local spec
    for spec; do
        awk -F'\t' -v fam="${spec%%:*}" -v job="${spec#*:}" -v tier="$HEADLINE" '
            FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") fmt[$2] = $3; next }
            FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0; cores[$1] = ($9 > 0 ? $9 : 1); of[$1] = ($10 == "" ? $1 : $10) } next }
            FNR == 1 || $1 != fam || !($2 in fmt) { next }
            { w[$2, $3] = $5; c[$2, $3] = $11; lanes[$3] = 1 }
            END {
                for (l in lanes) if (l != "z-fastq") for (x in fmt) if ((x, l) in w && (x, "z-fastq") in w && c[x, l] > 0 && c[x, "z-fastq"] > 0) {
                    k = l SUBSEP fmt[x]; n[k]++; sw[k] += log(w[x, l] / w[x, "z-fastq"]); sc[k] += log(c[x, l] / c[x, "z-fastq"]); tools[of[l]] = 1 }
                nt = 0; for (t in tools) tl[++nt] = t
                for (i = 1; i <= nt; i++) for (j = i + 1; j <= nt; j++) if (ord[tl[j]] < ord[tl[i]]) { t = tl[i]; tl[i] = tl[j]; tl[j] = t }
                out = ""
                for (i = 1; i <= nt; i++) { t = tl[i]; name = label[t]; gsub("_", " ", name); hit = miss = ""; top = 0
                    split("plain gzip", fs, " ")
                    for (f = 1; f <= 2; f++) { best = 0; seen = 0
                        for (l in lanes) if (of[l] == t) { k = l SUBSEP fs[f]; if (!n[k]) continue
                            seen = 1; if (cores[l] > top) top = cores[l]
                            r = exp(sw[k] / n[k]); if (r < 1 && (!best || cores[l] < cores[best])) { best = l; br = r; bc = exp(sc[k] / n[k]) } }
                        if (!seen) continue
                        if (best) hit = hit (hit == "" ? "" : ", and ") sprintf("from %d threads on %s input (%.2fx z-fastq%ss time for %.1fx its CPU time)", cores[best], fs[f], br, "\047", bc)
                        else miss = miss (miss == "" ? "" : " or ") fs[f] }
                    if (hit != "" && miss != "") s = sprintf("%s matches z-fastq on one core %s; on %s input no thread count up to %d does.", name, hit, miss, top)
                    else if (hit != "") s = sprintf("%s matches z-fastq on one core %s.", name, hit)
                    else if (miss != "") s = sprintf("%s does not match z-fastq on one core at any thread count up to %d, on %s input.", name, top, miss)
                    else continue
                    out = out (out == "" ? "" : " ") s }
                if (out != "") printf "**%s.** %s\n\n", job, out
            }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
    done
}
feature_readme() {
    family_checks
    readme_head
    printf '## Findings\n\n'
    process_findings
    printf '\n## Compared tools\n> <sub>*partial* tools do less than z-fastq: fastp checks records less strictly, and the FastQC tools write no reads and check no pairs. A job is `yes` when the tool was right on every input of it, `partly` on some, `-` when the tool has no such job. Threads: z-fastq runs on one; in the trimming and QC jobs every other tool runs at 1 thread and again at 2, 4, 8, and 16 on as many physical cores, at the %s inputs, in the same batches as z-fastq; copying runs every tool on one. Malformed rejected: of the broken pairs and records in the edge cases, how many the tool rejects.</sub>\n\n' "${TIER_SIZE[$HEADLINE]}"
    tools_table
    printf '\n'
    local name
    for name in "${LANES[@]}"; do [[ -z "${LANE_NOTE[$name]}" ]] || printf -- '- `%s %s`: *%s*\n' "${REG_LABEL[$name]}" "$(lane_version "$name")" "${LANE_NOTE[$name]}"; done
    printf '\n> [!IMPORTANT]\n> fastp computes its QC report and pair duplication on every run, so z-fastq'"'"'s trimming lanes write the QC summary and measure pair duplication too. In the QC job z-fastq also writes both mates, which it cannot skip, while the FastQC tools write tables only.\n'
    readme_right
    steps_part
    readme_section trimming Trimming "trim-truseq:TruSeq" "trim-nextera:Nextera"
    readme_section qc "QC tables" "qc-tables:FastQC tables"
    readme_repro
}
# The bytes each tool wrote on the plain headline inputs, as one sentence.
output_sentence() { # family [job]
    awk -F'\t' -v fam="$1" -v job="${2:-}" -v tier="$HEADLINE" -v size="${TIER_SIZE[$HEADLINE]}" 'FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $3 == "plain") want[$2] = 1; next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0; base[$1] = ($10 == $1) } next }
        $1 == fam && ($5 in want) && $3 == "pass" && base[$2] { n = split($4, f, ";"); for (i = 1; i <= n; i++) if (f[i] ~ /^out_bytes /) { split(f[i], q, " "); s[$2] += q[2]; k[$2]++ } }
        END { m = 0; for (l in s) list[++m] = l
              for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
              # In the QC job z-fastq writes both mates as well as its tables; the FastQC tools write tables only.
              qc = fam == "qc-tables"
              out = ""; for (i = 1; i <= m; i++) out = out (i == 1 ? "" : ", ") sprintf("%s %.1f MB%s", label[list[i]], s[list[i]] / k[list[i]] / 1e6,
                  !qc ? "" : list[i] == "z-fastq" ? " (both mates in full, which this job cannot skip, and its tables)" : " (tables only)")
              if (out != "") printf "%sOutput per plain %s input, mean over %d inputs: %s.\n\n", (job == "" ? "" : "**" job ".** "), size, k[list[1]], out }' \
        "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$FAMILY_CHECKS"
}

bench_main "$@"

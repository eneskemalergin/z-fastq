#!/usr/bin/env bash
# Benchmark core, sourced by bench/FEATURE/run.sh. A feature declares its job,
# lanes, families, inputs, and a fact function; this file checks the setup,
# checks every lane against the truth, times lanes in zebrac batches through
# zexec, caches results, and writes the feature's README and figures.
# Modes: see usage below, or run.sh help.
set -euo pipefail
export LC_ALL=C

BENCH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(dirname -- "$BENCH")"
FEATURE_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
FEATURE="$(basename -- "$FEATURE_DIR")"
RESULTS="$FEATURE_DIR/results"
CACHE="$BENCH/cache"
ZEXEC="$CACHE/zexec"
# The sample feature's output checker (bench/lib/subset.zig); built when a
# feature sets SUBSET_CHECKER=1.
SUBSET="$CACHE/subset"
# The subject binary; ZFASTQ=PATH times a build kept elsewhere (never rebuilt here).
ZFASTQ_DEFAULT="$REPO/zig-out/bin/z-fastq"
ZFASTQ="${ZFASTQ:-$ZFASTQ_DEFAULT}"
BENCH_CPU="${BENCH_CPU:-3}"
BENCH_MEM="${BENCH_MEM:-4G}"
BENCH_JOBS="${BENCH_JOBS:-4}"
REAL="$BENCH/data/real"
SHAPES="$BENCH/data/synthetic"

# Every tool any feature can use: name, path, receipt (- none, adapter local
# build), pin variable, language, light color, dark color, marker shape, label.
# Color follows the implementation language: Zig orange for z-fastq, rust reds
# for Rust, pinks for C++, grays for C, each language's own hue otherwise.
# Order here is the order of every table and legend, in every feature: z-fastq,
# the parser libraries, the Rust tools, SeqFu, SeqKit, seqtk, fastp, fqtools
# (by far the slowest counter, so count's bars read from fast to slow), then
# fastQValidator.
REGISTRY='z-fastq zig-out/bin/z-fastq - - Zig #f7a41d #f7a41d circle z-fastq
needletail bench/cache/adapters/native/needletail-adapter adapter NEEDLETAIL_VERSION Rust #b7410e #e0703a diamond Needletail
helicase-checked bench/cache/adapters/native/helicase-adapter adapter HELICASE_VERSION Rust #7f2a0a #f29b72 triangle Helicase
fq tools/bin/fq fq FQ_VERSION Rust #8b2e0b #f5b08e star fq
fqkit tools/bin/fqkit fqkit FQKIT_VERSION Rust #d4561f #c84a14 triangle-down fqkit
rasusa tools/bin/rasusa rasusa RASUSA_VERSION Rust #a33a10 #e8865a plus Rasusa
irma-core tools/bin/irma-core irma-core IRMA_CORE_VERSION Rust #c2501a #f3a37c cross IRMA_Core
fasten tools/bin/fasten_sample fasten FASTEN_VERSION Rust #6e2408 #f7c1a6 square Fasten
seqfu tools/bin/seqfu seqfu SEQFU_VERSION Nim #8f7d1c #d8c55a star SeqFu
seqkit tools/bin/seqkit seqkit SEQKIT_VERSION Go #007d9c #00add8 diamond SeqKit
seqtk tools/bin/seqtk seqtk SEQTK_VERSION C #57606a #9198a1 cross seqtk
fastp tools/bin/fastp fastp FASTP_VERSION C++ #a61e4d #f9a3c0 triangle fastp
fqtools tools/bin/fqtools fqtools FQTOOLS_VERSION C #24292f #d1d7e0 square fqtools
fastqvalidator tools/bin/fastQValidator fastqvalidator FASTQ_VALIDATOR_VERSION C++ #d6336c #f478a2 square fastQValidator'
declare -A REG_PATH REG_RECEIPT REG_PIN REG_LANG REG_LIGHT REG_DARK REG_SHAPE REG_LABEL REG_ORDER
_n=0
while read -r name path receipt pin lang light dark shape label; do
    [[ "$path" == /* ]] || path="$REPO/$path"
    REG_PATH[$name]="$path" REG_RECEIPT[$name]="$receipt" REG_PIN[$name]="$pin" REG_LANG[$name]="$lang"
    REG_LIGHT[$name]="$light" REG_DARK[$name]="$dark" REG_SHAPE[$name]="$shape" REG_LABEL[$name]="${label//_/ }" REG_ORDER[$name]=$((_n++))
done <<<"$REGISTRY"
# z-fastq's identity, and so its checks and batches, follow ZFASTQ.
REG_PATH[z-fastq]="$ZFASTQ"

# Runs and warmups per tier; odd counts, since zebrac's even-count median is
# the upper middle value. One warmup: the input is read into the page cache
# before every batch. The headline 1 GiB tier takes 3 runs while the suites
# are in development; a publication run sets `runs medium 5`. The large tier
# (10 GiB) is for stress runs only: it adds nothing past 1 GiB and costs most
# of a run.
declare -A TIER_SIZE=([tiny]="1 MiB" [xsmall]="10 MiB" [small]="100 MiB" [medium]="1 GiB" [large]="10 GiB")
declare -A RUNS=([tiny]=31 [xsmall]=21 [small]=11 [medium]=3 [many]=11)
declare -A WARM=([tiny]=1 [xsmall]=1 [small]=1 [medium]=1 [many]=1)

# --- declarations ---------------------------------------------------------
JOB="" TRUTH="" HEADLINE=medium SCALING="" FAMILIES=() LANES=() CUR_FAMILY=""
declare -A LANE_CLASS LANE_TMPL LANE_NOTE LANE_HEADLINE LANE_OUTS DESCRIBE
declare -A FAM_DESC FAM_LANES FAM_STDIN FAM_INPUTS FAM_KIND FAM_PARAMS FAM_VARY FAM_FORMATS FAM_LOOK FAM_MEM
job() { JOB="$1"; }
# order LANES: this feature's table and legend order, ahead of the registry's.
order() { local n i=0; for n in ${1//,/ }; do REG_ORDER[$n]=$((i++ - 1000)); done; }
truth() { TRUTH="$1"; }
headline() { HEADLINE="$1"; }
# runs TIER N: measured runs for a tier, when a feature needs fewer than the default.
runs() { RUNS[$1]="$2"; }
scaling() { SCALING="${1//,/ }"; }
# Quick-check cases for peer and ab: picked from the last release (auto), a
# declared list on top (add), or the declared list alone (only).
PICK_MODE=auto PICK_CASES=""
pick() { PICK_MODE="$1"; shift; PICK_CASES="${*//,/ }"; }
describe() { DESCRIBE[$1]="$2"; }
# remark SECTION "TEXT": a feature's own note, shown as a NOTE callout in that
# README section (scaling).
declare -A REMARK=()
remark() { REMARK[$1]="$2"; }
# group "LABEL" SOURCES: the kinds of input the figures pool shapes into, in order.
declare -A GROUP_OF GROUP_RANK GROUP_MEMBERS
GROUP_LIST=()
group() {
    local src
    GROUP_LIST+=("$1") GROUP_MEMBERS[$1]="${2//,/ }"
    for src in ${2//,/ }; do GROUP_OF[$src]="$1" GROUP_RANK[$src]="${#GROUP_LIST[@]}"; done
}
# A source's name for text: the description for a real run, NAME: description
# for a synthetic shape.
# The lanes of FAMILY that run on a TIER: z-fastq alone for the verify set; a
# headline-only lane (too slow to sweep sizes) only at the headline tier.
lanes_at() { # family tier
    local name
    [[ "$1" != verify ]] || { printf 'z-fastq\n'; return; }
    for name in ${FAM_LANES[$1]}; do
        [[ "${LANE_HEADLINE[$name]}" == 0 || "$2" == "$HEADLINE" ]] && printf '%s\n' "$name"
    done
    return 0
}
label_of() { if [[ -d "$SHAPES/$1" ]]; then printf '%s: %s\n' "$1" "${DESCRIBE[$1]:-}"; else printf '%s\n' "${DESCRIBE[$1]:-$1}"; fi; }
lane() {
    local name="$1" class="$2" tmpl="$3" arg
    shift 3
    [[ -v "REG_PATH[$name]" ]] || die "lane $name is not in the registry"
    LANES+=("$name") LANE_CLASS[$name]="$class" LANE_TMPL[$name]="$tmpl" LANE_NOTE[$name]="" LANE_HEADLINE[$name]=0
    for arg; do
        case "$arg" in
            note=*) LANE_NOTE[$name]="${arg#note=}" ;;
            headline-only) LANE_HEADLINE[$name]=1 ;;
            outputs=*) LANE_OUTS[$name]="${arg#outputs=}" ;;
            *) die "lane $name: unknown option $arg" ;;
        esac
    done
}
family() {
    local name="$1" arg
    shift
    FAMILIES+=("$name") CUR_FAMILY="$name"
    FAM_STDIN[$name]=0 FAM_INPUTS[$name]="" FAM_KIND[$name]=file FAM_PARAMS[$name]="" FAM_VARY[$name]="" FAM_FORMATS[$name]="plain gzip" FAM_LOOK[$name]=""
    for arg; do
        case "$arg" in
            *) FAM_DESC[$name]="$arg" ;;
        esac
    done
}
# command LANE "TEMPLATE": the lane's command in the current family.
declare -A FAM_TMPL=()
command() { FAM_TMPL[$CUR_FAMILY|$1]="$2"; }
# param NAME VALUE: {NAME} in the current family's commands. VALUE is one value,
# or short:A,long:B for a value by the read class of the case's source.
param() { FAM_PARAMS[$CUR_FAMILY]+="$1=$2 "; }
# vary NAME V1,V2,...: one case per value on every input, the value in the case id.
vary() { FAM_VARY[$CUR_FAMILY]="$1 ${2//,/ }"; }
# formats plain: the current family's cases in these formats only.
formats() { FAM_FORMATS[$CUR_FAMILY]="${*//,/ }"; }
# look full|table|lines: how a feature's own README order shows the family.
look() { FAM_LOOK[$CUR_FAMILY]="$1"; }
# memory SIZE: the current family's memory cap per run, in place of BENCH_MEM.
memory() { FAM_MEM[$CUR_FAMILY]="$1"; }
# rival any: the quick-check rival is the fastest other tool of any class,
# not only the fastest complete one.
RIVAL_ANY=0
rival() { [[ "$1" == any ]] || die "rival: unknown value $1"; RIVAL_ANY=1; }
# facts piped: each lane's stdout streams into fact while the lane runs,
# instead of through a file; for outputs too large to keep.
FACT_PIPE=0
facts() { [[ "$1" == piped ]] || die "facts: unknown mode $1"; FACT_PIPE=1; }
# expect VERDICT: what every corpus case must print, when no manifest number is the truth.
EXPECT=""
expect() { EXPECT="$1"; }
lanes() {
    local list="$1"
    [[ "$list" != all ]] || list="$(printf '%s,' "${LANES[@]}")"
    FAM_LANES[$CUR_FAMILY]="${list//,/ }"
    [[ "${2:-}" != stdin ]] || FAM_STDIN[$CUR_FAMILY]=1
}
# input real SOURCES TIERS [interleaved] | input shape SHAPES TIERS [interleaved] | input pair SOURCES TIERS | input many TIER [interleaved]
input() { FAM_INPUTS[$CUR_FAMILY]+="$*"$'\n'; [[ "$1" != many ]] || FAM_KIND[$CUR_FAMILY]=many; }
# edges SET [FAMILY...] [PARAMS]: the edge cases of SET; each lane runs with
# the command of the first listed family it is in (else its first family),
# with PARAMS (k=1;p=1) for that family's placeholders.
EDGE_SET="" EDGE_FAMS="" EDGE_PARAMS=-
edges() {
    EDGE_SET="$1"
    shift
    local arg
    for arg; do if [[ "$arg" == *=* ]]; then EDGE_PARAMS="$arg"; else EDGE_FAMS+="$arg "; fi; done
}
# agree LANE FIELD: where no independent truth holds every field, LANE's
# fields on the same input are the truth, and LANE itself is checked on FIELD
# against the manifest count. Facts are then FIELD VALUE [TOLERANCE] lines.
AGREE_LANE="" AGREE_FIELD=""
agree() { AGREE_LANE="$1" AGREE_FIELD="$2"; }
# word KEY "TEXT": a README phrase; the defaults are count's.
declare -A WORD=(
    [title]="Counting records" [subject]="z-fastq count" [others]="counters" [timed]="its count was right"
    [right]="the right count" [wrong]="a wrong count" [missed]="missed the count" [verb]="counts" [past]="counted"
    [payer]="a counter" [same]="counting" [many]="Counting" [unseen]="the count is wrong" [legend]="wrong count, exit 0" [forms]="Other input forms" [forms_sub]=""
    [truth]="Each count is compared with the count recorded when the corpus was cut, without z-fastq."
    [cells]="A cell is the count printed, \`reject\` (exited with an error), \`wrong: N\` (printed another count, exit 0), \`refused\` (rejected a valid input), or \`crash\`."
    [probes]="The bold count is how many of three malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, and \`-\` instead of \`+\`."
)
word() { [[ -v "WORD[$1]" ]] || { printf 'error: word: unknown key %s\n' "$1" >&2; exit 2; }; WORD[$1]="$2"; }
# shellcheck source=sources.sh
source "$BENCH/lib/sources.sh"
SCALING="${SHAPE_SOURCES//,/ }"

# --- logging ------------------------------------------------------------------
# One line per event (time, stage, tag, message) to stderr, the tag colored on
# a terminal. During a run, stderr also feeds the run's log file (see
# run_start). Lines come from shell builtins alone: no process per line.
LOG_FILE="" LOG_COLOR="" STAGE=run RUN_MODE="" RUN_T0=0 RUN_LOADS="" N_TIMED=0 N_CACHED=0
[[ ! -t 2 ]] || LOG_COLOR=1
log() { # stage tag message
    local stamp color=""
    printf -v stamp '%(%H:%M:%S)T' -1
    if [[ -n "$LOG_COLOR" ]]; then
        case "$2" in FAIL | STOP) color=31 ;; WARN) color=33 ;; PASS | DONE | END) color=32 ;; CACHE | SKIP) color=2 ;; esac
    fi
    if [[ -n "$color" ]]; then
        printf '%s %-7s \e[%sm%-5s\e[0m %s\n' "$stamp" "$1" "$color" "$2" "$3" >&2
    else
        printf '%s %-7s %-5s %s\n' "$stamp" "$1" "$2" "$3" >&2
    fi
}
die() {
    log "$STAGE" STOP "$*"
    exit 1
}
note() { log "$STAGE" NOTE "$*"; }
duration() { if (($1 >= 60)); then printf '%d min %d s' $(($1 / 60)) $(($1 % 60)); else printf '%d s' "$1"; fi; }
# A run's own log file in results/logs/, with results/latest.log pointing at
# it, and its line in results/runs.tsv: start, end, mode, revision, z-fastq
# binary, batches timed and reused, highest load, and outcome. All of the
# run's stderr (log lines, build output, tool failures) goes through one awk
# that shows it and appends it, without color, to the log file; the exit
# trap writes the run's line and waits for that awk, also when a run stops.
run_start() { # mode
    local load
    RUN_MODE="$1" RUN_T0=$SECONDS
    mkdir -p -- "$RESULTS/logs"
    printf -v LOG_FILE '%s/logs/%(%Y%m%d-%H%M%S)T-%s.log' "$RESULTS" -1 "$1"
    ln -sfn -- "logs/${LOG_FILE##*/}" "$RESULTS/latest.log"
    exec 2> >(awk -v file="$LOG_FILE" '{ print > "/dev/stderr"; line = $0; gsub(/\033\[[0-9;]*m/, "", line); print line >> file; fflush() }')
    LOG_PID=$!
    printf -v RUN_STARTED '%(%Y-%m-%d %H:%M:%S)T' -1
    RUN_REV="$(git -C "$REPO" rev-parse --short HEAD)" SRC_REV="$RUN_REV"
    [[ -z "$(git -C "$REPO" status --porcelain -- src build.zig build.zig.zon | head -1)" ]] || SRC_REV+=", source modified"
    [[ -z "$(git -C "$REPO" status --porcelain -- src build.zig build.zig.zon bench tools | head -1)" ]] || RUN_REV+="+modified"
    read -r load _ </proc/loadavg
    trap 'run_end $?' EXIT
    local size mtime built=""
    if read -r size mtime < <(stat -L -c '%s %Y' -- "$ZFASTQ" 2>/dev/null); then printf -v built ', z-fastq %s bytes built %(%Y-%m-%d %H:%M)T' "$size" "$mtime"; fi
    log run START "$FEATURE $1 at $RUN_REV$built, load $load"
}
run_end() { # exit status
    local outcome=ok end loadmax
    ((${1:-0} == 0)) || outcome=stopped
    printf -v end '%(%Y-%m-%d %H:%M:%S)T' -1
    loadmax="$(printf '%s\n' $RUN_LOADS | sort -n | tail -1)"
    local batches=""
    ((N_TIMED + N_CACHED == 0)) || batches=", $N_TIMED batches timed and $N_CACHED reused"
    log run END "$RUN_MODE $outcome after $(duration $((SECONDS - RUN_T0)))$batches; log ${LOG_FILE#"$REPO"/}"
    [[ -f "$RESULTS/runs.tsv" ]] || printf 'start\tend\tmode\trevision\tz_fastq\ttimed\treused\tload_max\toutcome\tlog\n' >"$RESULTS/runs.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$RUN_STARTED" "$end" "$RUN_MODE" "$RUN_REV" "$(fileid "$ZFASTQ" 2>/dev/null)" \
        "$N_TIMED" "$N_CACHED" "${loadmax:--}" "$outcome" "${LOG_FILE##*/}" >>"$RESULTS/runs.tsv"
    exec 2>&-
    wait "$LOG_PID" 2>/dev/null || true
}

# --- identities -----------------------------------------------------------
# File facts identify binaries and inputs (size, mtime, inode); no hashing.
fileid() { stat -L -c '%s:%Y:%i' -- "$@" | paste -sd, -; }
lane_bin() { printf '%s\n' "${REG_PATH[$1]}"; }
# Lane identities are read once per run by cache_lane_ids; until then, stat.
declare -A LANE_IDS=()
# A lane's identity is its binary (size, mtime, inode) and the cksum of its
# declared command, so a rebuilt tool or a changed command line re-checks and
# re-times, and nothing else does.
lane_tmpl() { printf '%s\n' "${FAM_TMPL[${2:-}|$1]:-${LANE_TMPL[$1]}}"; } # name [family]
lane_id_now() { printf '%s#%s\n' "$(fileid "$(lane_bin "$1")")" "$(printf '%s' "$(lane_tmpl "$1" "${2:-}")" | cksum | cut -d' ' -f1)"; }
lane_id() { if [[ -v "LANE_IDS[$1|${2:-}]" ]]; then printf '%s\n' "${LANE_IDS[$1|${2:-}]}"; else lane_id_now "$1" "${2:-}"; fi; }
cache_lane_ids() {
    local name fam
    for name in "${LANES[@]}"; do
        for fam in "" verify "${FAMILIES[@]}"; do LANE_IDS[$name|$fam]="$(lane_id_now "$name" "$fam")"; done
    done
}

# Expands a command template: {zfastq}, {bin}, {in} (space-separated paths
# for many-files), each case parameter in PARAMS (NAME=VALUE;...; - for none),
# {out1} and {out2}, the paths in OUT1 and OUT2 (/dev/null when unset), and
# for a lane that names its outputs from a stem, {stem}, {dir}, and {name}:
# STEM, its directory, and its last part. {in1} and {in2} are the first and
# last input path.
expand() { # template lane input params
    local t="$1" kv
    if [[ -n "${STEM:-}" ]]; then t="${t//\{stem\}/$STEM}" t="${t//\{dir\}/${STEM%/*}}" t="${t//\{name\}/${STEM##*/}}"; fi
    t="${t//\{zfastq\}/$ZFASTQ}"
    t="${t//\{bin\}/$(lane_bin "$2")}"
    t="${t//\{in1\}/${3%% *}}" t="${t//\{in2\}/${3##* }}"
    t="${t//\{in\}/$3}"
    if [[ "$4" != - ]]; then for kv in ${4//;/ }; do t="${t//\{${kv%%=*}\}/${kv#*=}}"; done; fi
    t="${t//\{out1\}/${OUT1:-/dev/null}}"
    printf '%s\n' "${t//\{out2\}/${OUT2:-/dev/null}}"
}
lane_cmd() { expand "$(lane_tmpl "$1" "${3:-}")" "$1" "$2" "${4:--}"; } # name input [family] [params]
# z-fastq's command with another build in its place, for ab.
build_cmd() { local t; t="$(lane_tmpl z-fastq "${3:-}")"; expand "${t//\{zfastq\}/$1}" z-fastq "$2" "${4:--}"; } # bin paths [family] [params]
# The key of one check: lane, lane identity (with the case parameters when
# the case has any), form, and input identity.
status_key() { # lane family form input-id params
    local id="${LANE_IDS[$1|$2]:-}"
    [[ -n "$id" ]] || id="$(lane_id_now "$1" "$2")"
    [[ "${5:--}" == - ]] || id+="@$5"
    printf '%s|%s|%s|%s\n' "$1" "$id" "$3" "$4"
}

# --- cases and truth --------------------------------------------------------
# Every case: family, case id, format, paths, truth, plain bytes, source, tier,
# kind, and parameters (NAME=VALUE;..., - for none).
manifest_row() { # SOURCE FILE -> records bytes
    awk -F'\t' -v f="$2" '$3 == f { print $4, $5; exit }' "$REAL/$1/manifest.tsv"
}
real_file() { # SOURCE TIER VARIANT -> file name
    if [[ "$3" == interleaved && -n "$(manifest_row "$1" "$2_interleaved.fastq")" ]]; then
        printf '%s_interleaved.fastq\n' "$2"
    else
        printf '%s.fastq\n' "$2"
    fi
}
emit_case() { # family id format paths truth bytes source tier kind params
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@"
}
# The family's parameters for SOURCE into CASE_PARAMS, a short:A,long:B value
# resolved by the source's read class (LONG_SOURCES, sources.sh).
case_params() { # family source
    local kv name value
    CASE_PARAMS=""
    for kv in ${FAM_PARAMS[$1]}; do
        name="${kv%%=*}" value="${kv#*=}"
        if [[ "$value" == short:*,long:* ]]; then
            if [[ ",$LONG_SOURCES," == *",$2,"* ]]; then value="${value##*long:}"; else value="${value#short:}" value="${value%%,*}"; fi
        fi
        CASE_PARAMS+="${CASE_PARAMS:+;}$name=$value"
    done
}
# PARAMS in name order, so one parameter set always has one spelling, and a
# case of one family can reuse another family's check or batch.
sorted_params() { [[ "$1" == - || -z "$1" ]] && printf '%s\n' "${1:--}" || tr ';' '\n' <<<"$1" | sort | paste -sd';' -; }
# One case, or one per value of the family's varied parameter (named in the id).
emit_cases() { # family id-stem format paths truth bytes source tier kind
    local vname vals v
    case_params "$1" "$7"
    read -r vname vals <<<"${FAM_VARY[$1]}"
    if [[ -z "$vname" ]]; then
        emit_case "$1" "$2:$3" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$(sorted_params "${CASE_PARAMS:--}")"
        return 0
    fi
    for v in $vals; do emit_case "$1" "$2:$vname=$v:$3" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$(sorted_params "${CASE_PARAMS:+$CASE_PARAMS;}$vname=$v")"; done
}
cases() {
    local fam line kind sources tiers variant src tier f rec bytes fmt shape paths truths
    for fam in "${FAMILIES[@]}"; do
        while read -r kind sources tiers variant; do
            [[ -n "$kind" ]] || continue
            for fmt in ${FAM_FORMATS[$fam]}; do
                case "$kind" in
                    real)
                        for src in ${sources//,/ }; do
                            for tier in ${tiers//,/ }; do
                                f="$(real_file "$src" "$tier" "${variant:-}")"
                                read -r rec bytes < <(manifest_row "$src" "$f") || die "no manifest row: $src $f"
                                local p="$REAL/$src/$f"
                                [[ "$fmt" == plain ]] || p+=.gz
                                emit_cases "$fam" "$src:$tier" "$fmt" "$p" "$rec" "$bytes" "$src" "$tier" real
                            done
                        done
                        ;;
                    shape)
                        for shape in ${sources//,/ }; do
                            for tier in ${tiers//,/ }; do
                                read -r rec bytes < <(jq -r --arg t "$tier" '.entries[$t] | "\(.count) \(.decompressed_bytes)"' "$SHAPES/$shape/state.json")
                                local p="$SHAPES/$shape/$tier.fastq"
                                # A paired shape counts pairs; its interleaved file holds two records each.
                                [[ "${variant:-}" != interleaved ]] || { p="$SHAPES/$shape/${tier}_interleaved.fastq"; rec=$((rec * 2)); }
                                [[ "$fmt" == plain ]] || p+=.gz
                                emit_cases "$fam" "$shape:$tier" "$fmt" "$p" "$rec" "$bytes" "$shape" "$tier" shape
                            done
                        done
                        ;;
                    pair)
                        for src in ${sources//,/ }; do
                            for tier in ${tiers//,/ }; do
                                # The truth is the pair count: R1's records.
                                local base="$REAL/$src/$tier" r1 r2
                                if [[ -d "$SHAPES/$src" ]]; then
                                    base="$SHAPES/$src/$tier"
                                    rec="$(jq -r --arg t "$tier" '.entries[$t].count' "$SHAPES/$src/state.json")"
                                else
                                    read -r rec _ < <(manifest_row "$src" "${tier}_R1.fastq") || die "no manifest row: $src ${tier}_R1.fastq"
                                fi
                                r1="${base}_R1.fastq" r2="${base}_R2.fastq"
                                bytes=$(($(stat -c %s "$r1") + $(stat -c %s "$r2")))
                                [[ "$fmt" == plain ]] || { r1+=.gz; r2+=.gz; }
                                emit_cases "$fam" "$src:$tier" "$fmt" "$r1 $r2" "$rec" "$bytes" "$src" "$tier" pair
                            done
                        done
                        ;;
                    many)
                        paths="" truths="" bytes=0
                        tier="$sources" variant="$tiers"
                        for src in "$REAL"/*/; do
                            src="${src%/}" src="${src##*/}"
                            f="$(real_file "$src" "$tier" "$variant")"
                            read -r rec b < <(manifest_row "$src" "$f") || continue
                            local p="$REAL/$src/$f"
                            [[ "$fmt" == plain ]] || p+=.gz
                            paths+="$p " truths+="$rec;" bytes=$((bytes + b))
                        done
                        emit_cases "$fam" "all:$tier" "$fmt" "${paths% }" "${truths%;}" "$bytes" all "$tier" many
                        ;;
                esac
            done
        done <<<"${FAM_INPUTS[$fam]}"
    done
}
# z-fastq alone, on every corpus file and shape tier, plain and gzip; none
# with `verify off`, for a feature whose families already hold every source.
VERIFY=1
verify() { [[ "$1" == off ]] || die "verify: unknown value $1"; VERIFY=0; }
verify_cases() {
    ((VERIFY)) || return 0
    local m src f rec bytes tier shape
    # z-fastq runs the first family's command, with that family's parameters.
    for m in "$REAL"/*/manifest.tsv; do
        src="$(basename -- "$(dirname -- "$m")")"
        while IFS=$'\t' read -r _ tier f rec bytes _ _; do
            # Every file on disk; a tier removed to save space is skipped.
            [[ "$f" == *.fastq && -f "$REAL/$src/$f" ]] || continue
            case_params "${FAMILIES[0]}" "$src"
            emit_case verify "$src:$f:plain" plain "$REAL/$src/$f" "$rec" "$bytes" "$src" "$tier" real "${CASE_PARAMS:--}"
            emit_case verify "$src:$f:gzip" gzip "$REAL/$src/$f.gz" "$rec" "$bytes" "$src" "$tier" real "${CASE_PARAMS:--}"
        done < <(tail -n +2 "$m")
    done
    for shape in "$SHAPES"/*/; do
        shape="$(basename -- "$shape")"
        jq -r '.entries | to_entries[] | "\(.key) \(.value.count) \(.value.decompressed_bytes)"' "$SHAPES/$shape/state.json" |
            while read -r tier rec bytes; do
                [[ -f "$SHAPES/$shape/$tier.fastq" ]] || continue
                case_params "${FAMILIES[0]}" "$shape"
                emit_case verify "$shape:$tier:plain" plain "$SHAPES/$shape/$tier.fastq" "$rec" "$bytes" "$shape" "$tier" shape "${CASE_PARAMS:--}"
                [[ ! -f "$SHAPES/$shape/$tier.fastq.gz" ]] ||
                    emit_case verify "$shape:$tier:gzip" gzip "$SHAPES/$shape/$tier.fastq.gz" "$rec" "$bytes" "$shape" "$tier" shape "${CASE_PARAMS:--}"
            done
    done
}

# --- init -----------------------------------------------------------------
ISSUES=() STOPS=0
issue() { ISSUES+=("$1: $2"); [[ "$1" != stop ]] || STOPS=$((STOPS + 1)); }
run_init() {
    STAGE=init
    local c name bin receipt version recipe pin v
    for c in jq gawk taskset systemd-run cksum zebrac zig; do
        command -v "$c" >/dev/null || issue stop "$c not found"
    done
    [[ "$(zebrac --version 2>/dev/null)" == "zebrac 0.6."* ]] || issue stop "zebrac 0.6.x expected"
    local paranoid
    read -r paranoid </proc/sys/kernel/perf_event_paranoid
    ((paranoid <= 2)) || issue stop "perf_event_paranoid is $paranoid; zebrac needs 2 or lower"
    local uid ctl
    uid="$(id -u)"
    ctl="/sys/fs/cgroup/user.slice/user-$uid.slice/user@$uid.service/cgroup.subtree_control"
    [[ -S "${XDG_RUNTIME_DIR:-/run/user/$uid}/bus" && -r "$ctl" && "$(<"$ctl")" == *cpu*memory* ]] ||
        issue stop "user systemd scope cannot apply cpu and memory limits"
    [[ -x "$ZEXEC" && ! "$BENCH/lib/zexec.zig" -nt "$ZEXEC" ]] || issue note "zexec missing or stale; the mode rebuilds it"
    [[ -x "$ZFASTQ" ]] || issue stop "z-fastq not built: zig build -Dstatic=true -Doptimize=ReleaseFast"
    [[ -z "$(find "$REPO/src" "$REPO/build.zig" "$REPO/build.zig.zon" -newer "$ZFASTQ" -print -quit 2>/dev/null)" ]] ||
        issue note "z-fastq is older than its source; the mode rebuilds it"
    # shellcheck disable=SC1091
    source "$REPO/tools/versions.sh"
    for name in "${LANES[@]}"; do
        bin="$(lane_bin "$name")" receipt="${REG_RECEIPT[$name]}" pin="${REG_PIN[$name]}"
        if [[ "$receipt" == adapter ]]; then
            [[ -x "$bin" && -z "$(find "$REPO/tools/wrappers/src" -newer "$bin" -print -quit)" ]] ||
                issue note "$name adapter missing or stale; the mode rebuilds it"
            continue
        fi
        [[ -x "$bin" ]] || { issue lane "$name: $bin not executable (tools/install.sh)"; continue; }
        [[ "$receipt" != - ]] || continue
        version="" recipe=""
        if [[ -r "$REPO/tools/.local/$receipt-current/receipt.tsv" ]]; then
            while IFS=$'\t' read -r k v; do
                case "$k" in version) version="$v" ;; recipe) recipe="$v" ;; esac
            done <"$REPO/tools/.local/$receipt-current/receipt.tsv"
        fi
        [[ "$version" == "${!pin}" ]] || issue lane "$name: installed ${version:-none}, pinned ${!pin} (tools/install.sh $receipt)"
        local rpin="${pin%_VERSION}_RECIPE"
        [[ -z "${!rpin:-}" || "$recipe" == "${!rpin}" ]] || issue lane "$name: installed recipe ${recipe:-none}, pinned ${!rpin} (tools/install.sh $receipt)"
    done
    local missing
    missing="$({ cases; verify_cases; } | cut -f4 | tr ' ' '\n' | sort -u | while read -r p; do [[ -f "$p" ]] || printf '%s\n' "$p"; done | head -3)"
    [[ -z "$missing" ]] || issue stop "declared input missing: ${missing//$'\n'/, } (bench/lib/corpus/make.sh, bench/lib/shapes/make.sh)"
    local load
    read -r load _ </proc/loadavg
    awk -v l="$load" 'BEGIN { exit !(l > 8) }' && issue note "load average $load is above the 8-CPU budget"
    local i tag
    for i in "${ISSUES[@]}"; do
        case "$i" in stop:*) tag=STOP ;; lane:* | input:*) tag=WARN ;; *) tag=NOTE ;; esac
        log init "$tag" "${i#*: }"
    done
    log init "$( ((STOPS)) && echo STOP || echo PASS)" "${#LANES[@]} lanes, ${#FAMILIES[@]} families, ${#ISSUES[@]} issues"
}

# --- builds -----------------------------------------------------------------
ensure_builds() {
    STAGE=build
    local built=""
    mkdir -p -- "$CACHE" "$RESULTS"
    if [[ ! -x "$ZEXEC" || "$BENCH/lib/zexec.zig" -nt "$ZEXEC" ]]; then
        log build RUN zexec
        built+=" zexec"
        zig build-exe "$BENCH/lib/zexec.zig" -OReleaseSmall -fstrip -fsingle-threaded -target x86_64-linux-none -femit-bin="$ZEXEC"
    fi
    if [[ -n "${SUBSET_CHECKER:-}" && (! -x "$SUBSET" || "$BENCH/lib/subset.zig" -nt "$SUBSET") ]]; then
        log build RUN "subset checker"
        built+=" subset"
        zig build-exe "$BENCH/lib/subset.zig" -OReleaseFast -fstrip -fsingle-threaded -target x86_64-linux-none -femit-bin="$SUBSET"
    fi
    if [[ "$ZFASTQ" == "$ZFASTQ_DEFAULT" && -n "$(find "$REPO/src" "$REPO/build.zig" "$REPO/build.zig.zon" -newer "$ZFASTQ" -print -quit 2>/dev/null)" ]]; then
        log build RUN "z-fastq (static ReleaseFast)"
        built+=" z-fastq"
        (cd "$REPO" && zig build -Dstatic=true -Doptimize=ReleaseFast)
    fi
    # The source revision of each default binary, one line per build: the
    # check above leaves no source newer than the binary.
    local zid
    zid="$(fileid "$ZFASTQ")"
    if [[ "$ZFASTQ" == "$ZFASTQ_DEFAULT" ]] && ! grep -q "^$zid"$'\t' "$CACHE/z-fastq.revisions" 2>/dev/null; then
        printf '%s\t%s\n' "$zid" "$SRC_REV" >>"$CACHE/z-fastq.revisions"
    fi
    local out="$CACHE/adapters/native"
    if [[ ! -x "$out/needletail-adapter" || -n "$(find "$REPO/tools/wrappers/src" -newer "$out/needletail-adapter" -print -quit)" ]]; then
        log build RUN "adapters (native CPU)"
        built+=" adapters"
        mkdir -p -- "$out"
        env CARGO_HOME="$CACHE/cargo-home" CARGO_TARGET_DIR="$CACHE/cargo-target-native" RUSTFLAGS="-C target-cpu=native" \
            cargo build --quiet --locked --release --bins --jobs "$BENCH_JOBS" --manifest-path "$REPO/tools/wrappers/Cargo.toml" --features native-cpu
        install -m 755 "$CACHE/cargo-target-native/release/needletail-adapter" "$CACHE/cargo-target-native/release/helicase-adapter" "$out/"
    fi
    if [[ -n "$built" ]]; then log build DONE "built:$built"; else log build SKIP "z-fastq, zexec, and adapters are current"; fi
}

# --- checks -----------------------------------------------------------------
# One check: run LANE on a case (stdin or file form), map output to the fact,
# compare with the truth. fact sees LANE, FAMILY, IN (the case paths), TRUTH
# (the case's manifest value), PARAMS, and, once the lane has ended, RC and
# ERR. With piped facts the lane's stdout streams into fact, and each named
# output ({out1}, {out2}) is a FIFO that fact reads; this shell holds the
# FIFOs' write ends open while the lane runs, so fact always sees them end.
# A failed run's verdict for a fact: reject, with the rule (S001 to S006,
# P001, P002) the tool named in ERR, if any.
reject_code() {
    local code
    code="$(sed 's/\x1b\[[0-9;]*m//g' "$ERR" | grep -o -m1 -E ': \[?[SP]00[1-6]\]?:?' | grep -o -E '[SP]00[1-6]' | head -1)" || code=""
    printf 'reject%s\n' "${code:+ $code}"
}
# fact with the case's variables, read from the calling check_one's locals.
fact_of() { LANE="$name" FAMILY="$fam" IN="$paths" TRUTH="$truth" PARAMS="$params" fact; }
check_one() { # lane form fam id fmt paths want bytes params truth
    local name="$1" form="$2" fam="$3" paths="$6" want="$7" params="${9:--}" truth="${10:-$7}" out rc=0 value cmd input=- status=fail diff=""
    out="$(mktemp "/var/tmp/bench-$FEATURE.XXXXXX")"
    [[ "$form" == stdin ]] || input="$paths"
    if ((FACT_PIPE)) && uses_stem "$name" "$fam"; then
        # A lane that names its outputs from a stem writes new regular files
        # (it may refuse an existing path), read by fact once it has ended.
        local -a sfx
        read -r -a sfx <<<"${LANE_OUTS[$name]//,/ }"
        local stem="$out.o" o1="" o2=""
        o1="$stem${sfx[0]}"
        ((${#sfx[@]} < 2)) || o2="$stem${sfx[1]}"
        cmd="$(STEM="$stem" lane_cmd "$name" "$input" "$fam" "$params")"
        if [[ "$form" == stdin ]]; then $cmd <"$paths" >/dev/null 2>"$out.err" || rc=$?; else $cmd >/dev/null 2>"$out.err" || rc=$?; fi
        if ((rc == 0)); then
            value="$(OUT1="$o1" OUT2="$o2" fact_of </dev/null | paste -sd';' -)"
        else
            value="$(RC="$rc" ERR="$out.err" fact_of </dev/null | paste -sd';' -)"
        fi
        # Every file under the stem, also one a tool named its own way.
        rm -f -- "$o1" "$o2" "${stem:?}"*
    elif ((FACT_PIPE)); then
        local tmpl o1="" o2="" w1 w2 fpid
        tmpl="$(lane_tmpl "$name" "$fam")"
        [[ "$tmpl" != *"{out1}"* ]] || { o1="$out.1"; mkfifo -- "$o1"; }
        [[ "$tmpl" != *"{out2}"* ]] || { o2="$out.2"; mkfifo -- "$o2"; }
        cmd="$(OUT1="$o1" OUT2="$o2" lane_cmd "$name" "$input" "$fam" "$params")"
        if [[ -n "$o1" ]]; then
            OUT1="$o1" OUT2="$o2" fact_of </dev/null >"$out.fact" 2>>"$out.err" &
            fpid=$!
            exec {w1}>"$o1"
            [[ -z "$o2" ]] || exec {w2}>"$o2"
            if [[ "$form" == stdin ]]; then $cmd <"$paths" >/dev/null 2>"$out.err" || rc=$?; else $cmd >/dev/null 2>"$out.err" || rc=$?; fi
            exec {w1}>&-
            [[ -z "$o2" ]] || exec {w2}>&-
            wait "$fpid" || true
        else
            local -a st
            set +e
            if [[ "$form" == stdin ]]; then $cmd <"$paths" 2>"$out.err" | fact_of >"$out.fact"; else $cmd 2>"$out.err" | fact_of >"$out.fact"; fi
            st=("${PIPESTATUS[@]}")
            set -e
            rc="${st[0]}"
        fi
        value="$(paste -sd';' - <"$out.fact")"
        # A lane that failed is judged by its exit and stderr alone.
        ((rc == 0)) || value="$(RC="$rc" ERR="$out.err" fact_of </dev/null | paste -sd';' -)"
        rm -f -- "$out.fact" "$o1" "$o2"
    else
        cmd="$(lane_cmd "$name" "$input" "$fam" "$params")"
        if [[ "$form" == stdin ]]; then $cmd <"$paths" >"$out" 2>"$out.err" || rc=$?; else $cmd >"$out" 2>"$out.err" || rc=$?; fi
        value="$(RC="$rc" ERR="$out.err" fact_of <"$out" | paste -sd';' -)"
    fi
    rm -f -- "$out" "$out.err"
    if [[ -n "$EXPECT" ]]; then
        # A verdict is the first fact; the rest are facts the report reads.
        [[ "$rc" != 0 || "${value%%;*}" != "$want" ]] || status=pass
    elif [[ -z "$AGREE_LANE" ]]; then
        [[ "$rc" != 0 || "$value" != "$want" ]] || status=pass
    else
        # The agreed lane on a file is checked on its anchor field only; every
        # other check compares fields with the agreed lane's (TRUTH holds them).
        [[ "$name" != "$AGREE_LANE" || "$form" != file ]] || want="$AGREE_FIELD $want"
        diff="$(field_diff "$value" "$want")"
        [[ "$rc" != 0 || -n "$diff" ]] || status=pass
        # A failure records the first differing field instead of every field.
        [[ -z "$diff" ]] || value="$diff"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$(status_key "$name" "$fam" "$form" "$(fileid $paths)" "$params")" "$status" "$rc" "$value" "${paths//$REPO\//}"
}
# FIELD VALUE [TOLERANCE] facts joined by ";", one block per file with "--"
# between files. Prints the first field of GOT that differs from WANT by more
# than GOT's tolerance (none: exact), "empty" or "files N" when the blocks do
# not line up, "no field" when the two share none; nothing when they agree.
field_diff() { # got want
    awk -v a="$1" -v b="$2" 'BEGIN {
        na = blocks(a, A); nb = blocks(b, B)
        if (!na) { print "empty"; exit }
        if (na != nb) { print "files " na; exit }
        for (k = 1; k <= na; k++) {
            n = split(A[k], f, ";"); shared = 0
            for (i = 1; i <= n; i++) {
                m = split(f[i], p, " ")
                if (m < 2 || !((k, p[1]) in B)) continue
                shared++; d = p[2] - B[k, p[1]]
                if ((d < 0 ? -d : d) > (m > 2 ? p[3] : 0) + 1e-9) { print p[1] " " p[2] " (want " B[k, p[1]] ")"; exit }
            }
            if (!shared) { print "no field"; exit }
        }
    }
    function blocks(s, out,   e, n, i, k, q, open) {
        n = split(s, e, ";"); k = open = 0
        for (i = 1; i <= n; i++) {
            if (e[i] == "--") { open = 0; continue }
            if (e[i] == "") continue
            if (!open) { k++; open = 1 }
            if (split(e[i], q, " ") >= 2) out[k, q[1]] = q[2]
            out[k] = out[k] (out[k] == "" ? "" : ";") e[i]
        }
        return k
    }'
}
# What a lane must print for PATHS: the manifest TRUTH, or with agree, the
# agreed lane's facts on each file; fails when the agreed lane has none.
want_of() { # paths truth
    [[ -z "$EXPECT" ]] || { printf '%s\n' "$EXPECT"; return 0; }
    [[ -n "$AGREE_LANE" ]] || { printf '%s\n' "$2"; return 0; }
    local p key out=""
    for p in $1; do
        key="$(status_key "$AGREE_LANE" "" file "$(fileid "$p")" -)"
        [[ "${CHECK_STATUS[$key]:-}" == pass ]] || return 1
        out+="${CHECK_VALUE[$key]};"
    done
    printf '%s\n' "${out%;}"
}
same_answer() { # got want
    if [[ -n "$AGREE_LANE" ]]; then [[ -z "$(field_diff "$1" "$2")" ]]; else [[ "$1" == "$2" ]]; fi
}
CHECKS="$RESULTS/checks.tsv"
declare -A CHECK_STATUS CHECK_VALUE
load_checks() {
    [[ -f "$CHECKS" ]] || return 0
    local k s r v
    while IFS=$'\t' read -r k s r v _; do CHECK_STATUS[$k]="$s" CHECK_VALUE[$k]="$v"; done <"$CHECKS"
}
run_checks() {
    STAGE=check
    local t0=$SECONDS before=0
    [[ ! -f "$CHECKS" ]] || before="$(wc -l <"$CHECKS")"
    load_checks
    cache_lane_ids
    if [[ -n "$AGREE_LANE" ]]; then
        check_jobs agree | dispatch_checks
        load_checks
    fi
    check_jobs lanes | dispatch_checks
    load_checks
    # A FAIL line for each check this run made; the DONE line counts every
    # current check, cached ones included, so known refusals stay as numbers.
    local current line total fails
    while IFS= read -r line; do
        log check FAIL "$line"
    done < <(tail -n +"$((before + 1))" "$CHECKS" | awk -F'\t' -v repo="$REPO/" '$2 == "fail" { split($1, k, "|"); sub(repo, "", $5); printf "%s  %s  %s (exit %s)\n", k[1], $5, ($4 == "" ? "nothing printed" : $4), $3 }')
    current="$(current_checks)"
    total="$(grep -c . <<<"$current" || true)" fails="$(awk -F'\t' '$3 == "fail"' <<<"$current" | wc -l)"
    log check DONE "$total current: $((total - fails)) pass, $fails fail$( ((fails == 0)) || awk -F'\t' '$3 == "fail" { n[$1]++ } END { s = ""; for (l in n) s = s (s ? ", " : " (") l " " n[l]; print s ")" }' <<<"$current")  $(duration $((SECONDS - t0)))"
}
# One job line per unchecked lane, form, and input: the agreed lane on every
# single-file input (agree), or every declared lane against what it must
# print (lanes). An input checked under two families runs once.
check_jobs() { # agree | lanes
    local wave="$1" fam id fmt paths truth bytes src tier kind params inid name form want
    { cases; verify_cases; } | with_input_ids | while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params inid; do
        if [[ "$wave" == agree ]]; then
            [[ "$kind" == many || -v "CHECK_STATUS[$(status_key "$AGREE_LANE" "" file "$inid" -)]" ]] ||
                printf '%s\tfile\t%s\t%s\t%s\t%s\t%s\t%s\t-\t%s\n' "$AGREE_LANE" "$fam" "$id" "$fmt" "$paths" "$truth" "$bytes" "$truth"
            continue
        fi
        form=file
        [[ "$fam" == verify || "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
        want=""
        lanes_of "$fam" "$tier"
        for name in ${LANES_AT[$fam|$tier]}; do
            [[ ! -v "CHECK_STATUS[$(status_key "$name" "$fam" "$form" "$inid" "$params")]" ]] || continue
            [[ -n "$want" ]] || want="$(want_of "$paths" "$truth")" ||
                { log check WARN "$id not checked: ${REG_LABEL[$AGREE_LANE]} failed on it"; continue 2; }
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$form" "$fam" "$id" "$fmt" "$paths" "$want" "$bytes" "$params" "$truth"
        done
    done | awk -F'\t' '!seen[$1 FS $2 FS $6 FS $9]++'
}
# Calls FUNCTION LINE for every line on stdin, BENCH_JOBS at a time, each in
# a background subshell of this shell, so a job costs no new script start.
pool() { # function
    local line running=0
    while IFS= read -r line; do
        if ((running >= BENCH_JOBS)); then
            wait -n || true
            running=$((running - 1))
        fi
        "$1" "$line" &
        running=$((running + 1))
    done
    wait || true
}
check_line() { local f; IFS=$'\t' read -r -a f <<<"$1"; check_one "${f[@]}"; }
# A feature's own cached step: JOBS (a function printing one job line each)
# runs through the pool, each line through RUN (a function of one line), into
# OUT, unless KEY equals the key the last run of NAME stored.
cached_step() { # name key jobs-function run-function out
    local name="$1" key="$2" out="$5" t0=$SECONDS jobs n
    STAGE="$name"
    if [[ -f "$out" && "$(cat "$RESULTS/$name.key" 2>/dev/null)" == "$key" ]]; then
        log "$name" CACHE "$(wc -l <"$out") results, every lane unchanged"
        return 0
    fi
    jobs="$("$3")"
    n="$(grep -c . <<<"$jobs" || true)"
    log "$name" RUN "$n jobs, $BENCH_JOBS at a time"
    pool "$4" <<<"$jobs" >"$out.part"
    mv -- "$out.part" "$out"
    printf '%s\n' "$key" >"$RESULTS/$name.key"
    log "$name" DONE "$n jobs, $(wc -l <"$out") results  $(duration $((SECONDS - t0)))"
}
# Every lane's identity in every family it is in, one line each: the lane
# part of a feature step's key.
lane_keys() {
    local name fam
    for name in "${LANES[@]}"; do
        for fam in "${FAMILIES[@]}"; do [[ " ${FAM_LANES[$fam]} " != *" $name "* ]] || printf '%s %s %s\n' "$name" "$fam" "$(lane_id "$name" "$fam")"; done
    done
}
dispatch_checks() { # job lines on stdin
    local jobs n
    jobs="$(cat)"
    [[ -n "$jobs" ]] || return 0
    n="$(grep -c . <<<"$jobs")"
    log check RUN "$n lane and input pairs, $BENCH_JOBS at a time"
    pool check_line <<<"$jobs" >>"$CHECKS"
}

# --- edge cases -------------------------------------------------------------
# The edge inputs, made once and shared by every feature.
EDGES="$CACHE/edges"
# Every edge row: name, then what each edge set (count, stats, check, sample) expects ("-": not in the set;
# a count, FIELD VALUE facts, reject, or policy), "probe" for the malformed
# inputs the tool list reports, and the row label.
edge_table() {
    cat <<'EOF'
lf	2	reads 2;bases 6	accept	accept	-	two LF records
nofinal	2	reads 2;bases 6	accept	accept	-	no final newline
crlf	2	reads 2;bases 6	accept	accept	-	CRLF line endings, final CRLF
empty	0	reads 0;bases 0	accept	accept	-	empty file
plusname	2	reads 2;bases 6	accept	accept	-	annotated plus line (`+r1`)
atqual	2	reads 2;bases 6	accept	accept	-	quality line starting with `@`
emptyseq	1	reads 1;bases 0	policy	policy	-	empty sequence and quality
lower	2	reads 2;bases 6	accept	accept	-	lowercase bases and `N`
tabhead	2	reads 2;bases 6	accept	accept	-	tab inside the header
longline	1	reads 1;bases 1048576	accept	accept	-	one 1 MiB sequence line
members	4	reads 4;bases 12	accept	accept	-	two gzip members
qrange	-	reads 1;bases 4;quality_sum 143;mean_quality 35.75;q20_bases 3;q30_bases 2;q20_pct 75;q30_pct 50	accept	accept	-	Q0, Q20, Q30, and Q93 in one record
iupac	-	reads 1;bases 7;a 1;c 1;g 1;t 1;n 1;other_bases 2	policy	policy	-	lowercase bases and IUPAC codes
blank	reject	reject	reject S003	reject S003	-	blank line between records
lenmis	reject	reject	reject S005	reject S005	probe	sequence and quality lengths differ
trunc	reject	reject	reject S004	reject S004	probe	truncated last record
dash	reject	reject	reject S001	reject S001	probe	`-` instead of `+`
qlow	-	reject	reject S006	reject S006	probe	quality byte below 33
qhigh	-	reject	reject S006	reject S006	-	quality byte above 126
noat	reject	reject	reject S003	reject S003	-	header without `@`
gztrunc	reject	reject	reject	reject	-	truncated gzip member
gzgarbage	reject	reject	reject	reject	-	non-zero bytes after the last gzip member
crlfnofinal	policy	policy	policy	policy	-	CRLF with no final LF
gzzeros	policy	policy	policy	policy	-	zero bytes after the last gzip member
dup	-	-	-	accept	-	the same record twice
EOF
}
edge_col() { case "$EDGE_SET" in count) echo 2 ;; stats) echo 3 ;; check) echo 4 ;; sample) echo 5 ;; *) declare -F feature_edge_list >/dev/null || die "unknown edge set $EDGE_SET" ;; esac; }
# name<TAB>expected<TAB>row label, for the declared edge set. A feature with
# its own rows defines feature_edge_list, make_feature_edges (its inputs),
# edge_input EDGE (the paths a row runs on), and edge_family EDGE (the family
# whose commands run the row; a lane outside it shows "-").
edge_list() {
    if declare -F feature_edge_list >/dev/null; then feature_edge_list; return; fi
    edge_table | awk -F'\t' -v c="$(edge_col)" -v OFS='\t' '$c != "-" { print $1, $c, $7 }'
}
probe_rows() { edge_table | awk -F'\t' -v c="$(edge_col)" '$c != "-" && $6 == "probe" { print $1 }'; }
make_edges() {
    mkdir -p -- "$EDGES"
    if declare -F make_feature_edges >/dev/null; then make_feature_edges; fi
    local r='@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n' e="$EDGES"
    [[ -f "$e/qhigh.fq" && -f "$e/gzzeros.fq.gz" && -f "$e/dup.fq" ]] && return 0
    printf '@r1\nACGT\n+\n!5?~\n' >"$e/qrange.fq"
    printf '@r1\nacgtRYn\n+\nIIIIIII\n' >"$e/iupac.fq"
    printf '@r1\nACGT\n+\nII I\n' >"$e/qlow.fq"
    printf '@r1\nACGT\n+\nII\177I\n' >"$e/qhigh.fq"
    printf "$r" >"$e/lf.fq"
    printf '@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII' >"$e/nofinal.fq"
    printf '@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nAC\r\n+\r\nII\r\n' >"$e/crlf.fq"
    : >"$e/empty.fq"
    printf '@r1\nACGT\n+r1\nIIII\n@r2\nAC\n+r2\nII\n' >"$e/plusname.fq"
    printf '@r1\nACGT\n+\n@III\n@r2\nAC\n+\nII\n' >"$e/atqual.fq"
    printf '@r1\n\n+\n\n' >"$e/emptyseq.fq"
    printf '@r1\nacgN\n+\nIIII\n@r2\nnn\n+\nII\n' >"$e/lower.fq"
    printf '@r1\tx=1\nACGT\n+\nIIII\n@r2\ty=2\nAC\n+\nII\n' >"$e/tabhead.fq"
    { printf '@r1\n'; head -c 1048576 /dev/zero | tr '\0' A; printf '\n+\n'; head -c 1048576 /dev/zero | tr '\0' I; printf '\n'; } >"$e/longline.fq"
    { printf "$r" | gzip -n; printf "$r" | gzip -n; } >"$e/members.fq.gz"
    printf '@r1\nACGT\n+\nIIII\n\n@r2\nAC\n+\nII\n' >"$e/blank.fq"
    printf '@r1\nACGT\n+\nIII\n@r2\nAC\n+\nII\n' >"$e/lenmis.fq"
    printf '@r1\nACGT\n+\nIIII\n@r2\nAC\n' >"$e/trunc.fq"
    printf '@r1\nACGT\n-\nIIII\n' >"$e/dash.fq"
    printf 'r1\nACGT\n+\nIIII\n' >"$e/noat.fq"
    printf "$r" | gzip -n | head -c 30 >"$e/gztrunc.fq.gz"
    { printf "$r" | gzip -n; printf 'garbage'; } >"$e/gzgarbage.fq.gz"
    printf '@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nAC\r\n+\r\nII\r' >"$e/crlfnofinal.fq"
    { printf "$r" | gzip -n; head -c 8 /dev/zero; } >"$e/gzzeros.fq.gz"
    printf '@r1\nACGT\n+\nIIII\n@r1\nACGT\n+\nIIII\n' >"$e/dup.fq"
}
edge_path() {
    if declare -F edge_input >/dev/null; then edge_input "$1"; return; fi
    if [[ -f "$EDGES/$1.fq" ]]; then printf '%s\n' "$EDGES/$1.fq"; else printf '%s\n' "$EDGES/$1.fq.gz"; fi
}
# Cell text: the count when right, "reject" when refused, "wrong: V" when a
# count differs, "refused" when a valid input was rejected, "crash" when the
# lane died by a signal.
run_edges() {
    [[ -n "$EDGE_SET" ]] || return 0
    make_edges
    local out="$RESULTS/edges.tsv" name edge expected label path value rc cell key fam efam
    declare -A edge_fam=()
    # Each lane runs with the command of the first edge family it is in, or
    # with its own command when edge families are named and it is in none;
    # with none named, the first family it is in.
    local search
    for name in "${LANES[@]}"; do
        efam="" search="$EDGE_FAMS"
        [[ -n "$search" ]] || search="${FAMILIES[*]}"
        for fam in $search; do [[ " ${FAM_LANES[$fam]} " != *" $name "* ]] || { efam="$fam"; break; }; done
        edge_fam[$name]="$efam"
    done
    # Cached while every lane's command and every edge file is unchanged.
    key="$(for name in "${LANES[@]}"; do printf '%s %s\n' "$name" "$(lane_id "$name" "${edge_fam[$name]}")"; done
        [[ "$EDGE_PARAMS" == - ]] || printf 'params %s\n' "$EDGE_PARAMS"
        if declare -F edge_family >/dev/null; then lane_keys; declare -f edge_family edge_input | cksum; fi
        edge_list | while IFS=$'\t' read -r edge _; do fileid $(edge_path "$edge"); done)"
    STAGE=edges
    if [[ -f "$out" && "$(cat "$RESULTS/edges.key" 2>/dev/null)" == "$key" ]]; then
        log edges CACHE "$(edge_list | wc -l) rows, every lane unchanged"
        return 0
    fi
    local t0=$SECONDS
    : >"$out"
    local o="$RESULTS/edges.out" tmpl cmd params="$EDGE_PARAMS" paths truth=""
    local stem=""
    for name in "${LANES[@]}"; do
        while IFS=$'\t' read -r edge expected label; do
            fam="${edge_fam[$name]}"
            if declare -F edge_family >/dev/null; then
                fam="$(edge_family "$edge")"
                if [[ " ${FAM_LANES[$fam]} " != *" $name "* ]]; then printf '%s\t%s\t-\n' "$name" "$edge" >>"$out"; continue; fi
            fi
            tmpl="$(lane_tmpl "$name" "$fam")"
            local o1="" o2="" stem=""
            [[ "$tmpl" != *"{out1}"* ]] || o1="$o.1"
            [[ "$tmpl" != *"{out2}"* ]] || o2="$o.2"
            if [[ -n "$fam" ]] && uses_stem "$name" "$fam"; then
                local -a sfx
                read -r -a sfx <<<"${LANE_OUTS[$name]//,/ }"
                stem="$o.s" o1="$stem${sfx[0]}" o2=""
                ((${#sfx[@]} < 2)) || o2="$stem${sfx[1]}"
                rm -f -- "$o1" "$o2"
            fi
            path="$(edge_path "$edge")" paths="$path"
            rc=0
            # A lane whose command names no input ({in}, {in1}, {in2}) reads it on stdin.
            if [[ (-n "$fam" && "${FAM_STDIN[$fam]}" == 1) || "$tmpl" != *"{in"* ]]; then
                cmd="$(STEM="$stem" OUT1="$o1" OUT2="$o2" lane_cmd "$name" - "$fam" "$params")"
                $cmd <"$path" >"$o" 2>"$o.err" || rc=$?
            else
                cmd="$(STEM="$stem" OUT1="$o1" OUT2="$o2" lane_cmd "$name" "$path" "$fam" "$params")"
                $cmd >"$o" 2>"$o.err" || rc=$?
            fi
            value="$(RC="$rc" ERR="$o.err" OUT1="$o1" OUT2="$o2" EDGE=1 fact_of <"$o" | paste -sd';' -)"
            # A verdict is the first fact.
            [[ -z "$EXPECT" ]] || value="${value%%;*}"
            # A signal death is 128 + the signal; exit 255 is a tool's own error.
            if ((rc > 128 && rc <= 192)); then
                cell="crash (signal $((rc - 128)))"
            elif [[ -n "$EXPECT" ]]; then
                # Verdicts: the lane's own verdict, marked wrong when it passes
                # a malformed input and refused when it fails a valid one.
                if [[ "$expected" == reject* && "$rc" == 0 ]]; then cell="wrong: accept"
                elif [[ "$expected" == policy || ("$expected" == reject* && "$value" != accept) ]]; then cell="$value"
                elif [[ "$expected" == reject* ]]; then cell="wrong: accept"
                elif [[ "$value" == accept ]]; then cell=accept
                elif ((rc == 0)); then cell="$value"
                else cell=refused; fi
            elif ((rc != 0)); then
                [[ "$expected" == reject || "$expected" == policy ]] && cell=reject || cell=refused
            elif [[ -n "$AGREE_LANE" ]]; then
                # Facts: ok, or the first wrong field; a malformed input that
                # exits 0 is wrong whatever it printed.
                local diff=""
                [[ "$expected" == reject || "$expected" == policy ]] || diff="$(field_diff "$value" "$expected")"
                if [[ "$expected" == reject ]]; then cell="wrong: exit 0"
                elif [[ "$expected" == policy || -z "$diff" ]]; then cell=ok
                elif [[ "$diff" == "no field" ]]; then cell="no field"
                else cell="wrong: ${diff%% *}"; fi
            elif [[ "$expected" == reject ]]; then
                cell="wrong: ${value:-empty}"
            elif [[ "$expected" == policy || "$value" == "$expected" ]]; then
                cell="${value:-empty}"
            else
                cell="wrong: ${value:-empty}"
            fi
            printf '%s\t%s\t%s\n' "$name" "$edge" "$cell" >>"$out"
            # Every file under the stem, also one a tool named its own way.
            [[ -z "$stem" ]] || rm -f -- "$o1" "$o2" "${stem:?}"*
        done < <(edge_list)
    done
    rm -f -- "$o" "$o.err" "$o.1" "$o.2"
    printf '%s\n' "$key" >"$RESULTS/edges.key"
    log edges DONE "$(edge_list | wc -l) rows  $(duration $((SECONDS - t0)))"
}

# --- timing -----------------------------------------------------------------
batch_settings() { # tier kind -> runs warmups
    local key="$1"
    [[ "$2" != many ]] || key=many
    [[ -v "RUNS[$key]" ]] || die "tier $key has no run settings"
    printf '%s %s\n' "${RUNS[$key]}" "${WARM[$key]}"
}
# A batch's identity, one "key value" line each: form, input files, runs,
# warmups, CPU, zexec, then one "lane NAME ID" line per lane timed.
batch_meta() { # form input-id runs warm lanes... -> sets META
    local form="$1" inid="$2" runs="$3" warm="$4" name
    shift 4
    printf -v META 'form %s\ninputs %s\nruns %s\nwarm %s\ncpu %s\nzexec %s\n' "$form" "$inid" "$runs" "$warm" "$BENCH_CPU" "$ZEXEC_ID"
    for name; do META+="lane $name ${LANE_IDS[$name|${BATCH_FAM:-}]}"$'\n'; done
}
# The newest cached batch that can stand for META: same form, inputs, runs,
# CPU, and zexec; at least as many warmups; every wanted lane timed with the
# same identity. Extra lanes in the batch do not matter: dropping a lane never
# re-times a case, adding one does (lanes compare only inside one batch).
find_batch() { # glob-prefix; reads META, sets BATCH
    local dir newest="" line have hw want=()
    BATCH=""
    mapfile -t want <<<"${META%$'\n'}"
    for dir in "$1"*/; do
        dir="${dir%/}"
        [[ -f "$dir/summary.tsv" && -f "$dir/meta.tsv" ]] || continue
        have=$'\n'"$(<"$dir/meta.tsv")"$'\n'
        for line in "${want[@]}"; do
            if [[ "$line" == "warm "* ]]; then
                hw="${have#*$'\n'warm }" hw="${hw%%$'\n'*}"
                ((hw >= ${line#warm })) || continue 2
            else
                [[ "$have" == *$'\n'"$line"$'\n'* ]] || continue 2
            fi
        done
        [[ -z "$newest" || "$dir/summary.tsv" -nt "$newest" ]] || continue
        newest="$dir/summary.tsv" BATCH="$dir"
    done
    [[ -n "$BATCH" ]]
}
# Whether LANE's command in FAMILY names its output files from a stem.
uses_stem() { local t; t="$(lane_tmpl "$1" "$2")"; [[ -n "${LANE_OUTS[$1]:-}" && ("$t" == *"{stem}"* || "$t" == *"{name}"*) ]]; }
# Timed runs of a lane that writes named files write them under OUT_DIR with
# zexec's @RUN@ token in the stem, so every run writes new files, and zexec
# removes them when the run ends.
OUT_DIR="/var/tmp/bench-outputs"
declare -A UNLINKS=()
# Appends one timed command (NAME=command) to SPECS, with the lane's output
# files set up and their -u options in UNLINKS[NAME]; NAME is the batch's name
# for the run (the lane, or a1, a2, b for ab's builds), BIN a z-fastq build.
SPECS=()
timed_spec() { # name lane input family params [bin]
    local name="$1" lane="$2" stem="" s cmd
    UNLINKS[$name]=""
    if uses_stem "$lane" "$4"; then
        mkdir -p -- "$OUT_DIR"
        stem="$OUT_DIR/$FEATURE-$name-@RUN@"
        for s in ${LANE_OUTS[$lane]//,/ }; do UNLINKS[$name]+=" -u $stem$s"; done
    fi
    if [[ -n "${6:-}" ]]; then cmd="$(STEM="$stem" build_cmd "$6" "$3" "$4" "$5")"; else cmd="$(STEM="$stem" lane_cmd "$lane" "$3" "$4" "$5")"; fi
    SPECS+=("$name=$cmd")
}
# Times every passing lane of one case in one zebrac batch, unless a cached
# batch can stand for it.
run_batch() { # fam id fmt paths truth bytes src tier kind params
    local fam="$1" id="$2" fmt="$3" paths="$4" tier="$8" kind="$9" params="${10}" name form key runs warm cmds=()
    form=file
    [[ "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
    local timed=() inid
    inid="$(fileid $paths)"
    for name in $(lanes_at "$fam" "$tier"); do
        [[ "${CHECK_STATUS[$(status_key "$name" "$fam" "$form" "$inid" "$params")]:-}" == pass ]] || continue
        timed+=("$name")
    done
    if ((${#timed[@]} == 0)); then
        log time SKIP "${PROGRESS:-} $fam $id: no lane passed its check"
        return 0
    fi
    read -r runs warm <<<"$(batch_settings "$tier" "$kind")"
    local bkey dir prefix="$RESULTS/batches/$fam-${id//:/-}-"
    BATCH_FAM="$fam"
    batch_meta "$form" "$inid" "$runs" "$warm" "${timed[@]}"
    [[ "$params" == - ]] || META+="params $params"$'\n'
    # The family's own batch, else any batch with the same identity: the same
    # commands on the same input time the same work, whatever family asked.
    if [[ -z "${RETIME:-}" ]] && { find_batch "$prefix" || find_batch "$RESULTS/batches/"; }; then
        printf '%s\t%s\t%s\n' "$fam" "$id" "$BATCH" >>"$RESULTS/current.tsv.part"
        log time CACHE "${PROGRESS:-} $fam $id"
        N_CACHED=$((N_CACHED + 1))
        return 0
    fi
    bkey="$(printf '%s' "$META" | cksum | cut -d' ' -f1)"
    dir="$prefix$bkey"
    printf '%s\t%s\t%s\n' "$fam" "$id" "$dir" >>"$RESULTS/current.tsv.part"
    SPECS=()
    for name in "${timed[@]}"; do timed_spec "$name" "$name" "$([[ "$form" == stdin ]] && echo - || echo "$paths")" "$fam" "$params"; done
    time_lanes "$dir" "$form" "$paths" "$runs" "$warm" "$fam" "$id" "${SPECS[@]}"
    N_TIMED=$((N_TIMED + 1))
}
# Times one zebrac batch into DIR, every NAME=COMMAND through zexec -r (stdin
# from PATHS when FORM is stdin), checks every run, and writes DIR/summary.tsv.
# META, when set, is written as DIR/meta.tsv first.
time_lanes() { # dir form paths runs warm fam id name=command...
    local dir="$1" form="$2" paths="$3" runs="$4" warm="$5" fam="$6" id="$7" spec name cmds=() names=() p load0 load1 lines
    shift 7
    mkdir -p -- "$dir"
    rm -f -- "$dir"/*.tsv "$dir/zebrac.json"
    [[ -z "${META:-}" ]] || printf '%s' "$META" >"$dir/meta.tsv"
    for spec; do
        name="${spec%%=*}"
        names+=("$name")
        if [[ "$form" == stdin ]]; then
            cmds+=("$ZEXEC -r $dir/$name.tsv -i $paths${UNLINKS[$name]:-} -- ${spec#*=}")
        else
            cmds+=("$ZEXEC -r $dir/$name.tsv${UNLINKS[$name]:-} -- ${spec#*=}")
        fi
    done
    for p in $paths; do cat -- "$p" >/dev/null; done
    read -r load0 _ </proc/loadavg
    log time RUN "${PROGRESS:-} $fam $id: ${#names[@]} lanes, $runs runs + $warm warmup${ETA_TEXT:+, $ETA_TEXT}"
    local t0=$SECONDS
    systemd-run --user --scope -q -p CPUQuota=100% -p "MemoryMax=${FAM_MEM[$fam]:-$BENCH_MEM}" \
        taskset -c "$BENCH_CPU" zebrac -q -w "$warm" -i "$runs" -a "$runs" -d 1 --json="$dir/zebrac.json" "${cmds[@]}" 2>"$dir/zebrac.err" ||
        die "zebrac failed for $fam $id: $(head -c 400 "$dir/zebrac.err")"
    read -r load1 _ </proc/loadavg
    # Output files a run left behind (zexec removes each run's own).
    rm -f -- "${OUT_DIR:?}/$FEATURE"-*
    # Each lane's zexec file must hold warmups + runs lines, every status 0.
    for name in "${names[@]}"; do
        lines="$(wc -l <"$dir/$name.tsv")"
        ((lines == runs + warm)) || die "$fam $id $name: $lines usage lines, expected $((runs + warm))"
        awk -F'\t' '$5 != 0 { exit 1 }' "$dir/$name.tsv" || die "$fam $id $name: a run exited nonzero"
    done
    summarize "$dir" "$fam" "$id" "$warm" "$load0" "$load1" "${names[@]}" >"$dir/summary.tsv.part"
    mv -- "$dir/summary.tsv.part" "$dir/summary.tsv"
    # The source the timed z-fastq came from; outside meta.tsv, so not part of the batch's identity.
    BATCH_SECS=$((SECONDS - t0)) RUN_LOADS+=" $load0 $load1"
    log time DONE "${PROGRESS:-} $(duration "$BATCH_SECS")  $(awk -F'\t' '{ ms = $5 / 1e6; printf "%s%s %s", (NR > 1 ? " | " : ""), $3, (ms >= 1000 ? sprintf("%.2f s", ms / 1000) : sprintf("%.1f ms", ms)) }' "$dir/summary.tsv")  load $load0-$load1"
}
median_col() { # column file skip -> median of the column after skipping the warmup lines
    tail -n +"$(($3 + 1))" "$2" | cut -f"$1" | sort -n | awk '{ v[NR] = $1 } END { print v[int((NR + 1) / 2)] }'
}
summarize() {
    local dir="$1" fam="$2" id="$3" warm="$4" l0="$5" l1="$6" i=0 name
    shift 6
    for name; do
        local t="$dir/$name.tsv" rss_med rss_min rss_max cpu_med usr_med sys_med flt
        rss_med="$(median_col 1 "$t" "$warm")"
        rss_min="$(tail -n +"$((warm + 1))" "$t" | cut -f1 | sort -n | head -1)"
        rss_max="$(tail -n +"$((warm + 1))" "$t" | cut -f1 | sort -n | tail -1)"
        cpu_med="$(tail -n +"$((warm + 1))" "$t" | awk -F'\t' '{ print $2 + $3 }' | sort -n | awk '{ v[NR] = $1 } END { print v[int((NR + 1) / 2)] }')"
        usr_med="$(median_col 2 "$t" "$warm")" sys_med="$(median_col 3 "$t" "$warm")" flt="$(median_col 4 "$t" "$warm")"
        jq -r --argjson i "$i" --arg fam "$fam" --arg id "$id" --arg lane "$name" --arg class "${LANE_CLASS[$name]:-subject}" \
            --arg rss "$rss_med" --arg rmin "$rss_min" --arg rmax "$rss_max" --arg cpu "$cpu_med" --arg usr "$usr_med" \
            --arg sys "$sys_med" --arg flt "$flt" --arg l0 "$l0" --arg l1 "$l1" '
            .results[$i] | [$fam, $id, $lane, $class, .wall_time.median, .wall_time.q1, .wall_time.q3,
            $rss, $rmin, $rmax, $cpu, $usr, $sys, $flt, .instructions.median, .cpu_cycles.median,
            .cache_references.median, .cache_misses.median, .branch_misses.median, .sample_count, $l0, $l1] | @tsv' "$dir/zebrac.json"
        i=$((i + 1))
    done
}
SUMMARY_HEAD=$'family\tcase\tlane\tclass\twall_ns\twall_q1_ns\twall_q3_ns\trss_kib\trss_min_kib\trss_max_kib\tcpu_us\tuser_us\tsys_us\tminor_faults\tinstructions\tcycles\tcache_refs\tcache_misses\tbranch_misses\tsamples\tload_start\tload_end'
run_timing() { # [case id pattern] [family]
    local filter="${1:-}" only="${2:-}" fam id fmt paths truth bytes src tier kind params
    cache_lane_ids
    ZEXEC_ID="$(fileid "$ZEXEC")"
    mkdir -p -- "$RESULTS/batches"
    : >"$RESULTS/current.tsv.part"
    STAGE=time
    local list=() tiers=() stored=() line k i before left
    declare -A tier_secs=() tier_runs=()
    while IFS= read -r line; do
        IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params <<<"$line"
        [[ -z "$filter" ]] || [[ "$id" =~ $filter ]] || continue
        [[ -z "$only" || "$fam" == "$only" ]] || continue
        list+=("$line") tiers+=("$tier")
        # Any stored batch of the case: a glob, so no process per case.
        local found=("$RESULTS/batches/$fam-${id//:/-}-"*/summary.tsv)
        [[ -e "${found[0]}" ]] && stored+=(1) || stored+=("")
    done < <(cases)
    for ((k = 0; k < ${#list[@]}; k++)); do
        IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params <<<"${list[$k]}"
        # Time still to go: each remaining case with no batch on disk yet, at
        # the mean of its tier so far (of any tier when its own has none).
        ETA_TEXT="" left=0
        if ((${#tier_runs[@]})); then
            local all=0 alln=0 t
            for t in "${!tier_runs[@]}"; do all=$((all + tier_secs[$t])) alln=$((alln + tier_runs[$t])); done
            for ((i = k; i < ${#list[@]}; i++)); do
                [[ -z "${stored[$i]}" ]] || continue
                t="${tiers[$i]}"
                if [[ -v "tier_runs[$t]" ]]; then left=$((left + tier_secs[$t] / tier_runs[$t])); else left=$((left + all / alln)); fi
            done
            ((left == 0)) || ETA_TEXT="ETA about $(duration "$left")"
        fi
        PROGRESS="[$((k + 1))/${#list[@]}]" before=$N_TIMED
        run_batch "$fam" "$id" "$fmt" "$paths" "$truth" "$bytes" "$src" "$tier" "$kind" "$params"
        if ((N_TIMED > before)); then tier_secs[$tier]=$((${tier_secs[$tier]:-0} + BATCH_SECS)) tier_runs[$tier]=$((${tier_runs[$tier]:-0} + 1)); fi
    done
    PROGRESS="" ETA_TEXT=""
    mv -- "$RESULTS/current.tsv.part" "$RESULTS/current.tsv"
}
# The batch of every declared case, lanes declared for its tier only, so a
# report never needs a timing run: the batch the last timing run chose for
# it (current.tsv), else the newest of its own. A batch timed for another
# family's case shows under this case's family and id.
collect() {
    printf '%s\n' "$SUMMARY_HEAD"
    local declared fam id fmt paths truth bytes src tier kind params newest
    declare -A chosen=()
    if [[ -f "$RESULTS/current.tsv" ]]; then
        while IFS=$'\t' read -r fam id newest; do chosen[$fam|$id]="$newest/summary.tsv"; done <"$RESULTS/current.tsv"
    fi
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params; do
        declared=" $(lanes_at "$fam" "$tier" | tr '\n' ' ')"
        newest="${chosen[$fam|$id]:-}"
        if [[ ! -f "$newest" ]]; then
            newest=""
            local f
            for f in "$RESULTS/batches/$fam-${id//:/-}-"*/summary.tsv; do
                [[ -f "$f" && (-z "$newest" || "$f" -nt "$newest") ]] && newest="$f"
            done
        fi
        [[ -n "$newest" ]] || continue
        awk -F'\t' -v d="$declared" -v fam="$fam" -v id="$id" -v OFS='\t' 'index(d, " " $3 " ") { $1 = fam; $2 = id; print }' "$newest"
    done < <(cases)
}

# --- report -----------------------------------------------------------------
# Aligns a TSV (header, alignment row of l/r, rows) into a Markdown table.
mdtable() {
    awk -F'\t' '
        NR == 1 { n = NF; for (i = 1; i <= NF; i++) { h[i] = $i; w[i] = length($i) } next }
        NR == 2 { for (i = 1; i <= NF; i++) a[i] = $i; next }
        { r++; for (i = 1; i <= n; i++) { c[r, i] = $i; if (length($i) > w[i]) w[i] = length($i) } }
        END {
            for (i = 1; i <= n; i++) if (w[i] < 3) w[i] = 3
            line = "|"; for (i = 1; i <= n; i++) line = line " " pad(h[i], w[i], a[i]) " |"; print line
            line = "|"; for (i = 1; i <= n; i++) line = line " " (a[i] == "r" ? dashes(w[i] - 1) ":" : dashes(w[i])) " |"; print line
            for (j = 1; j <= r; j++) { line = "|"; for (i = 1; i <= n; i++) line = line " " pad(c[j, i], w[i], a[i]) " |"; print line }
        }
        function dashes(k,   s) { s = ""; while (k-- > 0) s = s "-"; return s }
        function pad(s, k, al,   g) { g = k - length(s); while (g-- > 0) s = (al == "r") ? " " s : s " "; return s }'
}
picture() { # stem alt; nothing when the figure was not drawn
    [[ -e "$RESULTS/figures/$1-light.png" ]] || return 0
    cat <<EOF
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/$1-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/$1-light.png">
  <img src="results/figures/$1-light.png" alt="$2" width="100%">
</picture>
EOF
}
style_tsv() {
    printf 'lane\tlabel\tclass\tlight\tdark\tshape\torder\tlanguage\n'
    local name
    for name in "${LANES[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "${REG_LABEL[$name]}" "${LANE_CLASS[$name]}" "${REG_LIGHT[$name]}" "${REG_DARK[$name]}" "${REG_SHAPE[$name]}" "${REG_ORDER[$name]}" "${REG_LANG[$name]}"
    done
}
cases_tsv() {
    printf 'family\tcase\tformat\tsource\ttier\tkind\tbytes\tlabel\torder\tgroup\tgroup_order\tparams\n'
    local fam id fmt paths truth bytes src tier kind params n=0
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params; do
        # A pair cut from a synthetic shape is a shape in the figures.
        [[ "$kind" != pair || ! -d "$SHAPES/$src" ]] || kind=shape
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$fam" "$id" "$fmt" "$src" "$tier" "$kind" "$bytes" \
            "$(label_of "$src")" $((n++)) "${GROUP_OF[$src]:-Other}" "${GROUP_RANK[$src]:-99}" "$params"
    done < <(cases)
}

lane_version() {
    local name="$1" bin receipt pin v
    bin="$(lane_bin "$name")" receipt="${REG_RECEIPT[$name]}" pin="${REG_PIN[$name]}"
    case "$name" in
        z-fastq) "$bin" --version | awk '{ print $2 }' ;;
        *)
            if [[ "$receipt" == adapter ]]; then
                "$bin" --version | sed -n 's/.*engine=[a-z]* \([^;]*\);.*/\1/p'
            else
                awk -F'\t' '$1 == "version" { print $2 }' "$REPO/tools/.local/$receipt-current/receipt.tsv"
            fi
            ;;
    esac
}
# Per lane, the geometric mean over the headline shapes of FAMILY of its time
# and peak RSS as multiples of z-fastq's, on the shapes both ran. With several
# FAMILY:MODE arguments, one table with a Mode column, a row per mode and lane.
# With THROUGH=1, a last column counts the malformed inputs each tool lets through (lets_through).
summary_table() { # family | family:mode...
    local fams
    fams="$(printf '%s|' "$@")"
    awk -F'\t' -v fams="${fams%|}" -v tier="$HEADLINE" -v through="${THROUGH:-0}" '
        BEGIN { nf = split(fams, fl, "|"); for (i = 1; i <= nf; i++) { split(fl[i], q, ":"); fam[i] = q[1]; mode[i] = q[2]; at[q[1]] = i } }
        FILENAME ~ /through.tsv$/ { lets[$2, $1] = $4; next }
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && ($1 in at) && $5 == tier && $6 != "many") { fmt[$1, $2] = $3; if ($3 == "gzip") hasgz = 1 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        FNR == 1 || !(($1, $2) in fmt) { next }
        { w[$1, $2, $3] = $5; m[$1, $2, $3] = $8; lanes[$1, $3] = 1; cs[$1, $2] = 1 }
        END {
            moded = mode[1] != ""; nfmt = hasgz ? 2 : 1
            print (moded ? "Mode\t" : "") "Tool\tPlain time\tPlain memory" (hasgz ? "\tGzip time\tGzip memory" : "") (through ? "\tPasses malformed" : "")
            print (moded ? "l\t" : "") "l\tr\tr" (hasgz ? "\tr\tr" : "") (through ? "\tl" : "")
            for (x = 1; x <= nf; x++) {
                F = fam[x]; n = 0; delete list
                for (k in lanes) { split(k, q, SUBSEP); if (q[1] == F) list[++n] = q[2] }
                for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
                for (i = 1; i <= n; i++) {
                    l = list[i]; line = (moded ? mode[x] "\t" : "") label[l]; gsub("_", " ", line)
                    split("plain gzip", fs, " ")
                    for (f = 1; f <= nfmt; f++) {
                        k = sw = sm = 0
                        for (c in cs) { split(c, q, SUBSEP); if (q[1] != F) continue; c2 = q[2]
                            if (fmt[F, c2] == fs[f] && (F, c2, l) in w && (F, c2, "z-fastq") in w) {
                                k++; sw += log(w[F, c2, l] / w[F, c2, "z-fastq"]); sm += log(m[F, c2, l] / m[F, c2, "z-fastq"]) } }
                        line = line "\t" (k ? sprintf("%.2fx", exp(sw / k)) : "-") "\t" (k ? sprintf("%.2fx", exp(sm / k)) : "-")
                    }
                    print line (through ? "\t" (((F, l) in lets) ? lets[F, l] : (("", l) in lets) ? lets["", l] : "-") : "")
                }
            }
        }' "$REPORT_DIR/through.tsv" "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# z-fastq's own geometric means over the headline inputs, as one sentence.
subject_sentence() { # family
    awk -F'\t' -v fam="$1" -v tier="$HEADLINE" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") fmt[$2] = $3; next }
        FNR == 1 || $1 != fam || !($2 in fmt) || $3 != "z-fastq" { next }
        { k[fmt[$2]]++; w[fmt[$2]] += log($5); m[fmt[$2]] += log($8) }
        END {
            if (!k["gzip"]) {
                printf "z-fastq itself takes %.1f ms and %.2f MiB per plain input (geometric mean over %d inputs); the table gives every tool as a multiple of that.\n",
                    exp(w["plain"] / k["plain"]) / 1e6, exp(m["plain"] / k["plain"]) / 1024, k["plain"]
                exit
            }
            printf "z-fastq itself takes %.1f ms and %.2f MiB per plain input, %.1f ms and %.2f MiB per gzip input (geometric means over %s); the table gives every tool as a multiple of that.\n",
                exp(w["plain"] / k["plain"]) / 1e6, exp(m["plain"] / k["plain"]) / 1024, exp(w["gzip"] / k["gzip"]) / 1e6, exp(m["gzip"] / k["gzip"]) / 1024,
                (k["plain"] == k["gzip"] ? k["plain"] " inputs each" : k["plain"] " and " k["gzip"] " inputs")
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/summary.tsv"
}
# Per lane, user-space counters per input byte, geometric mean over the plain headline inputs.
counters_table() { # family
    awk -F'\t' -v fam="$1" -v tier="$HEADLINE" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many" && $3 == "plain") bytes[$2] = $7; next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        FNR == 1 || $1 != fam || !($2 in bytes) || $15 <= 0 || $16 <= 0 { next }
        { mib = bytes[$2] / 1048576; k[$3]++; ins[$3] += log($15 / bytes[$2]); ipc[$3] += log($15 / $16)
          cm[$3] += log(($18 + 1) / mib); bm[$3] += log(($19 + 1) / mib) }
        END {
            print "Tool\tInstructions per byte\tInstructions per cycle\tCache misses per MiB\tBranch misses per MiB"
            print "l\tr\tr\tr\tr"
            n = 0; for (l in k) list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= n; i++) { l = list[i]; name = label[l]; gsub("_", " ", name)
                printf "%s\t%.2f\t%.2f\t%.0f\t%.0f\n", name, exp(ins[l] / k[l]), exp(ipc[l] / k[l]), exp(cm[l] / k[l]), exp(bm[l] / k[l]) }
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# The headline findings, computed from the report tables so they stay true
# after every run: correctness, speed against the fastest full-check tool,
# memory against every tool, and the kinds of input where z-fastq is not
# fastest. PARTS picks the bullets (right speed memory where); LABEL prefixes
# each bullet's name; ANY=1 names the fastest tool of any class before the
# fastest full-check one.
findings() { # family [parts] [label]
    awk -F'\t' -v fam="$1" -v tier="$HEADLINE" -v right="${WORD[right]}" -v missed="${WORD[missed]}" -v wrong="${WORD[wrong]}" \
        -v parts=" ${2:-right speed memory where} " -v lbl="${3:-}" -v any="${ANY:-0}" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") { fmt[$2] = $3; grp[$2] = $10 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; cls[$1] = $3; ord[$1] = $7 } next }
        FILENAME ~ /checks.tsv$/ { if ($1 == "z-fastq") { if ($3 == "pass") zp++; else zf++ } next }
        FILENAME ~ /edges.tsv$/ {
            if (FNR == 1 || $6 == "policy") next
            if ($4 == "z-fastq") { if ($6 != "right") zbad++ } else if ($6 == "wrong" || $6 == "crash") bad[$4] = 1
            next
        }
        FNR == 1 || $1 != fam || !($2 in fmt) { next }
        { w[$2, $3] = $5; m[$2, $3] = $8; lanes[$3] = 1 }
        END {
            n = 0; for (l in lanes) if (l != "z-fastq") list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            # Correctness.
            line = sprintf("- **Right answers:** z-fastq returned %s on all %d checked inputs", right, zp)
            if (zf) line = sprintf("- **Right answers:** z-fastq %s on %d of %d checked inputs", missed, zf, zp + zf)
            line = line (zbad ? "." : " and rejected every malformed edge case.")
            k = 0; for (i = 1; i <= n; i++) if (label[list[i]] in bad) k++
            if (k) line = line sprintf(" %d of the %d other tools print %s or crash on at least one.", k, n, wrong)
            if (index(parts, " right ")) print line
            # Speed against the fastest complete lane, per format.
            split("plain gzip", fs, " ")
            for (f = 1; f <= 2; f++) {
                best[f] = ""; br[f] = 0
                for (i = 1; i <= n; i++) if (cls[list[i]] == "complete") { r = ratio(list[i], fs[f], w); if (r && (best[f] == "" || r < br[f])) { best[f] = list[i]; br[f] = r } }
            }
            if (index(parts, " speed ")) {
                head = "- **" (lbl == "" ? "Speed" : lbl " speed") ":** "
                if (any) {
                    # The fastest tool of any class first, per format.
                    for (f = 1; f <= 2; f++) { fa[f] = ""; fr[f] = 0
                        for (i = 1; i <= n; i++) { r = ratio(list[i], fs[f], w); if (r && (fa[f] == "" || r < fr[f])) { fa[f] = list[i]; fr[f] = r } } }
                    head = head sprintf("the fastest other tool takes %.2fx z-fastq%ss time on plain input (%s)", fr[1], "\047", label[fa[1]]) (fr[2] ? sprintf(" and %.2fx on gzip (%s)", fr[2], label[fa[2]]) : "") "; "
                }
                if (best[1] == "") print head "no other full-check tool runs it."
                else if (br[2] == 0)
                    printf "%s%s, the fastest other full-check tool, takes %.2fx z-fastq%ss time on plain input.\n", head, label[best[1]], br[1], "\047"
                else if (best[1] == best[2])
                    printf "%s%s, the fastest other full-check tool, takes %.2fx z-fastq%ss time on plain input and %.2fx on gzip.\n", head, label[best[1]], br[1], "\047", br[2]
                else
                    printf "%sthe fastest other full-check tool takes %.2fx z-fastq%ss time on plain input (%s) and %.2fx on gzip (%s).\n", head, br[1], "\047", label[best[1]], br[2], label[best[2]]
            }
            # Memory against every other lane, both formats.
            lo = hi = 0
            for (i = 1; i <= n; i++) for (f = 1; f <= 2; f++) { r = ratio(list[i], fs[f], m); if (r) { if (!lo || r < lo) lo = r; if (r > hi) hi = r } }
            for (c in fmt) if ((c, "z-fastq") in m) { zm += log(m[c, "z-fastq"]); zk++; if (fmt[c] == "gzip") zg = 1 }
            if (index(parts, " memory ")) printf "- **Memory%s:** every other tool peaks at %.1fx to %.0fx z-fastq%ss %.2f MiB.\n", (lbl == "" ? "" : ", " tolower(lbl)), lo, hi, "\047", exp(zm / zk) / 1024
            # Kinds of input where some tool is faster than z-fastq.
            out = ""
            for (c in fmt) gs[grp[c]] = 1
            for (g in gs) for (f = 1; f <= 2; f++) for (i = 1; i <= n; i++) {
                k = s = 0
                for (c in fmt) if (grp[c] == g && fmt[c] == fs[f] && (c, list[i]) in w && (c, "z-fastq") in w) { k++; s += log(w[c, list[i]] / w[c, "z-fastq"]) }
                if (k && exp(s / k) < 0.99) out = out sprintf("%s%s on %s %s (%.2fx)", (out == "" ? "" : "; "), label[list[i]], fs[f], tolower(g), exp(s / k))
            }
            where = "- **Where z-fastq is not fastest" (lbl == "" ? "" : ", " tolower(lbl)) ":** "
            if (!index(parts, " where ")) { }
            else if (out == "") print where "on no kind of input; it is the fastest tool on every kind of input (geometric mean over the group), " (zg ? "plain and gzip." : "plain.")
            else print where out "."
        }
        function ratio(l, f, a,   c, k, s) { for (c in fmt) if (fmt[c] == f && (c, l) in a && (c, "z-fastq") in a) { k++; s += log(a[c, l] / a[c, "z-fastq"]) } return k ? exp(s / k) : 0 }
    ' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/checks.tsv" "$REPORT_DIR/edges.tsv" "$REPORT_DIR/summary.tsv"
}
# Two sentences under the scaling figures, computed from the scaling cases:
# z-fastq's throughput at the largest size, who passes it there, and whose
# peak memory grows with the input.
scaling_sentences() { # family
    local range
    range="$(scaling_range)"
    awk -F'\t' -v fam="$1" -v srcs=" $SCALING " -v toplabel="${range##* to }" -v verb="${WORD[verb]}" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && index(srcs, " " $4 " ")) { fmt[$2] = $3; src[$2] = $4; size[$2] = $7 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        FNR == 1 || $1 != fam || !($2 in fmt) { next }
        { tp = size[$2] / ($5 / 1e9); v[$2, $3] = tp; r[$2, $3] = $8; lanes[$3] = 1
          if (size[$2] > top[src[$2], fmt[$2]]) { top[src[$2], fmt[$2]] = size[$2]; big[src[$2], fmt[$2]] = $2 } }
        END {
            n = 0; for (l in lanes) if (l != "z-fastq") list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            split("plain gzip", fs, " ")
            line = "At " toplabel ", z-fastq " verb
            for (f = 1; f <= 2; f++) {
                lo = hi = 0
                for (k in big) { split(k, p, SUBSEP); if (p[2] != fs[f]) continue; c = big[k]; t = v[c, "z-fastq"] / 2^30
                    if (!lo || t < lo) lo = t; if (t > hi) hi = t }
                line = line sprintf("%s %.1f to %.1f GiB/s on %s input", (f > 1 ? " and" : ""), lo, hi, fs[f])
            }
            out = ""
            for (i = 1; i <= n; i++) for (k in big) { c = big[k]; split(k, p, SUBSEP); l = list[i]
                if (v[c, l] > v[c, "z-fastq"] * 1.01)
                    out = out (out == "" ? "" : "; ") sprintf("%s on %s %s (%.2fx)", label[l], p[2], p[1], v[c, l] / v[c, "z-fastq"]) }
            print line (out == "" ? "; no other tool is faster there." : ". Faster on a single shape there, as throughput over z-fastq'"'"'s: " out ".")
            zlo = zhi = 0; grow = ""
            for (c in fmt) if ((c, "z-fastq") in r) { x = r[c, "z-fastq"]; if (!zlo || x < zlo) zlo = x; if (x > zhi) zhi = x }
            for (i = 1; i <= n; i++) { l = list[i]; g = 0
                for (k in big) { split(k, p, SUBSEP); small = 0
                    for (c in fmt) if (src[c] == p[1] && fmt[c] == p[2] && (c, l) in r && (!small || size[c] < size[sm])) { small = 1; sm = c }
                    if (small && (big[k], l) in r) { q = r[big[k], l] / r[sm, l]; if (q > g) { g = q; gk = p[2] " " p[1] } } }
                if (g >= 1.25) grow = grow (grow == "" ? "" : "; ") sprintf("%s by %.1fx on %s", label[l], g, gk) }
            if (zhi - zlo < 256) printf "z-fastq peaks at %.2f MiB at every size and shape.", zhi / 1024
            else printf "z-fastq peaks at %.2f to %.2f MiB across sizes and shapes.", zlo / 1024, zhi / 1024
            print (grow == "" ? " No other tool grows by more than a quarter either." : " Memory that grows with the input: " grow ".")
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
}
# A CAUTION callout naming the tools that print a wrong answer with exit 0 or
# crash on malformed edge cases; nothing when none does.
caution() {
    awk -F'\t' -v wrong="${WORD[wrong]}" -v unseen="${WORD[unseen]}" 'FILENAME ~ /style.tsv$/ { if (FNR > 1) ord[$2] = $7 + 0; next }
        FNR > 1 && $6 == "wrong" { w[$4] = 1 } FNR > 1 && $6 == "crash" { c[$4] = 1 }
        END {
            for (t in w) k1++
            for (t in c) k2++
            wl = names(w); cl = names(c)
            if (wl == "" && cl == "") exit
            printf "> [!CAUTION]\n>"
            if (wl != "") printf " %s %s %s and %s 0 on some malformed inputs, so a script cannot tell %s.", wl, (k1 > 1 ? "print" : "prints"), wrong, (k1 > 1 ? "exit" : "exits"), unseen
            if (cl != "") printf " %s %s on some malformed gzip input.", cl, (k2 > 1 ? "crash" : "crashes")
            printf "\n\n"
        }
        # The tools in a set, in registry order, joined with commas and "and".
        function names(set,   t, list, i, j, n, out) {
            n = 0; for (t in set) list[++n] = t
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= n; i++) out = out (i == 1 ? "" : i == n ? " and " : ", ") list[i]
            return out
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/edges.tsv"
}
# Per lane, its time and memory reading stdin as a multiple of its own on the
# same file by path, geometric mean over the shared shapes.
stdin_table() { # family
    awk -F'\t' -v fam="$1" -v base="${FAMILIES[0]}" '
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        FNR == 1 { next }
        $1 == fam { sw[$2, $3] = $5; sm[$2, $3] = $8; lanes[$3] = 1; cs[$2] = 1; next }
        $1 == base { fw[$2, $3] = $5; fm[$2, $3] = $8 }
        END {
            print "Tool\tPlain time\tPlain memory\tGzip time\tGzip memory"
            print "l\tr\tr\tr\tr"
            n = 0; for (l in lanes) list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= n; i++) {
                l = list[i]; line = label[l]
                split("plain gzip", fs, " ")
                for (f = 1; f <= 2; f++) {
                    k = a = b = 0
                    for (c in cs) { p = split(c, q, ":"); if (q[p] == fs[f] && (c, l) in sw && (c, l) in fw) { k++; a += log(sw[c, l] / fw[c, l]); b += log(sm[c, l] / fm[c, l]) } }
                    line = line "\t" (k ? sprintf("%.2fx", exp(a / k)) : "-") "\t" (k ? sprintf("%.2fx", exp(b / k)) : "-")
                }
                print line
            }
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# One sentence for a family timed as one call over many files.
many_sentence() { # family
    local files bytes
    read -r files bytes < <(cases | awk -F'\t' -v f="$1" '$1 == f { print split($4, p, " "), $6; exit }')
    awk -F'\t' -v fam="$1" -v files="$files" -v bytes="$bytes" -v gerund="${WORD[many]}" -v verb="${WORD[verb]}" '
        FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
        FNR == 1 || $1 != fam { next }
        { n = split($2, p, ":"); f = p[n]; w[$3, f] = $5; m[$3, f] = $8; if ($3 != "z-fastq") other[$3] = 1 }
        END {
            if (!(("z-fastq", "plain") in w) || !(("z-fastq", "gzip") in w)) { print "Not timed yet."; exit }
            printf "%s all %d xsmall corpus files (%.0f MiB) in one call, z-fastq takes %.1f ms on plain input and %.1f ms on gzip, at %.2f MiB peak.", gerund, files, bytes / 1048576, w["z-fastq", "plain"] / 1e6, w["z-fastq", "gzip"] / 1e6, m["z-fastq", "plain"] / 1024
            k = 0; for (o in other) k++
            for (o in other) printf " %s%s takes %.1fx and %.1fx as long with %.0fx and %.0fx the memory.", label[o], (k == 1 ? ", the only other tool that " verb " many files in one call," : ""), w[o, "plain"] / w["z-fastq", "plain"], w[o, "gzip"] / w["z-fastq", "gzip"], m[o, "plain"] / m["z-fastq", "plain"], m[o, "gzip"] / m["z-fastq", "gzip"]
            printf "\n"
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
}
# Per lane (and per family, for a feature whose edge rows each run in one
# family), the malformed edge inputs it lets through with exit 0, as short
# names ("nothing" when it rejects them all), and their count: "N of M" or
# "none". One line per lane: lane, family (empty for all), names, count.
lets_through() {
    local name fam fams=""
    declare -F edge_family >/dev/null && fams="${FAMILIES[*]}"
    for name in "${LANES[@]}"; do
        for fam in "" $fams; do
            awk -F'\t' -v l="$name" -v f="$fam" -v rowfam="$(edge_list | while IFS=$'\t' read -r e _; do printf '%s=%s ' "$e" "$(declare -F edge_family >/dev/null && edge_family "$e")"; done)" 'BEGIN {
                    short["blank"] = "blank line"; short["lenmis"] = "length mismatch"; short["trunc"] = "truncated record"
                    short["dash"] = "bad plus line"; short["qlow"] = "quality range"; short["qhigh"] = "quality range"
                    short["noat"] = "bad header"; short["gztrunc"] = "broken gzip"; short["gzgarbage"] = "broken gzip"
                    n = split(rowfam, rf, " "); for (i = 1; i <= n; i++) { split(rf[i], kv, "="); of[kv[1]] = kv[2] } }
                FNR == NR { if ($2 ~ /^reject/) { want[$1] = 1; label[$1] = $3 } next }
                $1 == l && ($2 in want) && $3 != "-" && (f == "" || of[$2] == f) { total++
                    if ($3 ~ /^wrong/) { passed++; k = ($2 in short) ? short[$2] : label[$2]; if (!(k in seen)) { seen[k] = 1; out = out (out == "" ? "" : ", ") k } } }
                END { if (f == "" || total) printf "%s\t%s\t%s\t%s\n", l, f, (out == "" ? "nothing" : out), (passed ? passed " of " total : "none") }' <(edge_list) "$RESULTS/edges.tsv"
        done
    done
}
# One findings bullet: each valid edge input on which some tool exits 0 but
# writes the wrong output, with those tools (in table order) and how; nothing
# when none does.
valid_wrong_bullet() {
    awk -F'\t' 'FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0 } next }
        FNR == NR { if ($2 !~ /^reject/ && $2 != "policy") { ok[$1] = 1; what[$1] = $3; row[$1] = ++rows; name[rows] = $1 } next }
        ($2 in ok) && $3 ~ /^wrong/ { s = $3; sub(/^wrong: /, "", s); if (s ~ /^[0-9]+$/) s = "printed " s; else if (s ~ /^[a-z_]+$/) s = s " wrong"
            n = ++cnt[$2]; tool[$2, n] = $1; how[$2, n] = s }
        END {
            for (r = 1; r <= rows; r++) { e = name[r]; if (!(e in cnt)) continue
                m = cnt[e]
                for (a = 1; a <= m; a++) for (b = a + 1; b <= m; b++) if (ord[tool[e, b]] < ord[tool[e, a]]) {
                    t = tool[e, a]; tool[e, a] = tool[e, b]; tool[e, b] = t; t = how[e, a]; how[e, a] = how[e, b]; how[e, b] = t }
                part = what[e] ": "
                for (a = 1; a <= m; a++) part = part (a > 1 ? ", " : "") label[tool[e, a]] ((a == m || how[e, a] != how[e, a + 1]) ? " (" how[e, a] ")" : "")
                out = out (out == "" ? "" : "; ") part }
            if (out != "") { gsub(/`/, "", out); print "- **Wrong on valid input, exit 0:** " out "." } }' \
        <(edge_list) "$REPORT_DIR/style.tsv" "$RESULTS/edges.tsv"
}
# Per lane, "N of M": the malformed edge inputs it rejects, over every one it ran.
rejects_all() { # lane
    awk -F'\t' -v l="$1" 'FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next } $1 == l && ($2 in want) && $3 != "-" { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' <(edge_list) "$RESULTS/edges.tsv"
}
# yes when every check of LANE in FAMILY passed, partly when some did, no when none.
mode_cell() { # lane family
    awk -F'\t' -v l="$1" -v f="$2" '$1 == f && $2 == l { n++; if ($3 == "pass") p++ } END { print (n == 0) ? "-" : (p == n) ? "yes" : (p > 0) ? "partly" : "no" }' "$FAMILY_CHECKS"
}
# Every current check of every family into FAMILY_CHECKS: family, lane,
# status, value, case id. Written once per report.
FAMILY_CHECKS=""
family_checks() {
    local fam id fmt paths truth bytes src tier kind params inid name form key
    FAMILY_CHECKS="$REPORT_DIR/family-checks.tsv"
    load_checks
    cache_lane_ids
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params inid; do
        form=file
        [[ "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
        for name in ${FAM_LANES[$fam]}; do
            key="$(status_key "$name" "$fam" "$form" "$inid" "$params")"
            [[ ! -v "CHECK_STATUS[$key]" ]] || printf '%s\t%s\t%s\t%s\t%s\n' "$fam" "$name" "${CHECK_STATUS[$key]}" "${CHECK_VALUE[$key]}" "$id"
        done
    done < <(cases | with_input_ids) >"$FAMILY_CHECKS"
}
# One row per lane: each kind of malformed edge input (SPEC is
# Label=edge,edge|...) as rejects, or **passes** when the lane exits 0 on any
# of its edge inputs; then the columns EXTRA (tab-separated headers) that
# FUNCTION prints for the lane, or none.
kinds_table() { # spec [extra-headers function]
    local name spec="$1" extra="${2:-}" fn="${3:-}" head="" align="" tail
    head="$(tr '|' '\n' <<<"$spec" | cut -d= -f1 | paste -sd$'\t' -)"
    align="$(tr '|' '\n' <<<"$spec" | sed 's/.*/l/' | paste -sd$'\t' -)"
    {
        printf 'Tool\t%s%s\n' "$head" "${extra:+$'\t'$extra}"
        printf 'l\t%s%s\n' "$align" "${extra:+$'\t'$(tr '\t' '\n' <<<"$extra" | sed 's/.*/l/' | paste -sd$'\t' -)}"
        for name in "${LANES[@]}"; do
            tail=""
            [[ -z "$fn" ]] || tail="$("$fn" "$name")"
            awk -F'\t' -v l="$name" -v label="${REG_LABEL[$name]}" -v spec="$spec" -v tail="$tail" '
                BEGIN { nk = split(spec, kinds, "|"); for (i = 1; i <= nk; i++) sub(/^[^=]*=/, "", kinds[i]) }
                FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next }
                $1 == l && ($2 in want) { cell[$2] = $3 }
                END {
                    line = label
                    for (i = 1; i <= nk; i++) { n = split(kinds[i], e, ","); passed = 0; seen = 0
                        for (j = 1; j <= n; j++) { if (e[j] in cell) seen = 1; if (cell[e[j]] ~ /^wrong/) passed = 1 }
                        line = line "\t" (!seen || cell[e[1]] == "-" ? "-" : passed ? "**passes**" : "rejects") }
                    print line (tail == "" ? "" : "\t" tail)
                }' <(edge_list) "$RESULTS/edges.tsv"
        done
    } | mdtable
}
# Per lane, "N of M": the malformed edge inputs it rejects, over every one in the set.
rejects_all() { # lane
    awk -F'\t' -v l="$1" 'FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next } $1 == l && ($2 in want) { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' <(edge_list) "$RESULTS/edges.tsv"
}
probes() { # lane -> "N of M", over the edge set's probe rows
    awk -F'\t' -v l="$1" -v rows=" $(probe_rows | tr '\n' ' ')" '$1 == l && index(rows, " " $2 " ") { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' "$RESULTS/edges.tsv"
}
# Appends each case's input identity (every path's fileid, joined by ",") as
# a tenth column, from one stat call for all the files instead of one each.
with_input_ids() {
    local rows
    rows="$(cat)"
    [[ -n "$rows" ]] || return 0
    awk -F'\t' 'FNR == NR { id[$1] = $2; next }
        { n = split($4, p, " "); s = ""; for (i = 1; i <= n; i++) s = s (i > 1 ? "," : "") id[p[i]]; print $0 "\t" s }' \
        <(cut -f4 <<<"$rows" | tr ' ' '\n' | sort -u | xargs -d '\n' stat -L -c $'%n\t%s:%Y:%i') - <<<"$rows"
}
# The lanes of FAMILY at TIER, as a space-separated list cached per family and tier.
declare -A LANES_AT=()
lanes_of() { [[ -v "LANES_AT[$1|$2]" ]] || LANES_AT[$1|$2]="$(lanes_at "$1" "$2" | tr '\n' ' ')"; }
# Every check of a declared input with the current binaries, one line per
# lane and input: lane, source, status, exit, printed, path. Older check lines
# (rebuilt binaries, inputs no longer declared) stay cached but never count.
current_checks() {
    local fam id fmt paths truth bytes src tier kind params inid name form
    cache_lane_ids
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params inid; do
        form=file
        [[ "$fam" == verify || "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
        lanes_of "$fam" "$tier"
        for name in ${LANES_AT[$fam|$tier]}; do
            printf '%s\t%s\t%s\t%s\n' "$(status_key "$name" "$fam" "$form" "$inid" "$params")" "$name" "$src" "${paths%% *}"
        done
    done < <({ cases; verify_cases; } | with_input_ids) | awk -F'\t' 'FNR == NR { if (!($1 in want)) { want[$1] = $2 "\t" $3 "\t" $4 } next }
        ($1 in want) && !($1 in done) { done[$1] = 1; split(want[$1], w, "\t")
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", w[1], w[2], $2, $3, substr($4, 1, 40), w[3] }' - "$CHECKS"
}
# One sentence on z-fastq, then the other tools wrong on inputs they were
# given, as a table (tool, wrong, of, inputs), or a sentence when none was.
correctness_sentence() {
    local name pass fail rows="" zf="" where
    for name in "${LANES[@]}"; do
        read -r pass fail < <(awk -F'\t' -v l="$name" '$1 == l { if ($3 == "pass") a++; else b++ } END { print a + 0, b + 0 }' "$REPORT_DIR/checks.tsv")
        if [[ "$name" == z-fastq ]]; then
            local t top=""
            for t in tiny xsmall small medium large; do verify_cases | cut -f8 | grep -qx "$t" && top="${TIER_SIZE[$t]}"; done
            if ((VERIFY)); then zf="z-fastq returned ${WORD[right]} on all $pass inputs: every corpus file and synthetic shape on disk, ${TIER_SIZE[tiny]} to $top, plain and gzip."
            else zf="z-fastq returned ${WORD[right]} on all $pass inputs it was given."; fi
            ((fail == 0)) || zf="z-fastq was wrong on $fail of $((pass + fail)) inputs."
            continue
        fi
        ((fail > 0)) || continue
        where="$(awk -F'\t' -v l="$name" '$1 == l && $3 == "fail" { print $2 }' "$REPORT_DIR/checks.tsv" | sort -u |
            while read -r src; do if [[ -d "$SHAPES/$src" ]]; then printf '%s, ' "$src"; else printf '%s, ' "${DESCRIBE[$src]:-$src}"; fi; done)"
        rows+="${REG_LABEL[$name]}"$'\t'"$fail"$'\t'"$((pass + fail))"$'\t'"${where%, }"$'\n'
    done
    printf '%s ' "$zf"
    if [[ -z "$rows" ]]; then
        printf 'Every other tool was right on every input it was given.\n'
    else
        printf 'Every other tool was right on every input it was given except these, which are not timed where they %s:\n\n' "${WORD[missed]}"
        { printf 'Tool\tWrong\tOf\tInputs\n'; printf 'l\tr\tr\tl\n'; printf '%s' "$rows"; } | mdtable
    fi
}
# The edge-case matrix for plot.py: row, label, expected, tool, cell, verdict.
edges_tsv() {
    printf 'row\tlabel\texpected\ttool\tcell\tverdict\n'
    local row=0 edge expected label name cell verdict
    while IFS=$'\t' read -r edge expected label; do
        for name in "${LANES[@]}"; do
            cell="$(awk -F'\t' -v l="$name" -v e="$edge" '$1 == l && $2 == e { print $3 }' "$RESULTS/edges.tsv")"
            case "$cell" in
                wrong*) verdict=wrong ;;
                crash*) verdict=crash ;;
                refused) verdict=refused ;;
                -) verdict=policy ;;
                *) verdict=right ;;
            esac
            [[ "$expected" != policy ]] || verdict=policy
            # Expected facts show as ok, as their right cells do.
            [[ -z "$AGREE_LANE" || "$expected" != *" "* ]] || expected=ok
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$row" "${label//\`/}" "${expected/policy/tool policy}" "${REG_LABEL[$name]}" "$cell" "$verdict"
        done
        row=$((row + 1))
    done < <(edge_list)
}
details() { printf '<details><summary>%s</summary>\n\n' "$1"; }
# The size span of the scaling inputs, as "1 MiB to 10 GiB".
scaling_range() {
    local tiers
    tiers="$(cases | awk -F'\t' -v s=" $SCALING " 'index(s, " " $7 " ") { print $8 }' | sort -u)"
    local t first="" last="" name
    for t in tiny xsmall small medium large; do
        grep -qx "$t" <<<"$tiers" || continue
        name="${TIER_SIZE[$t]}"
        [[ -n "$first" ]] || first="$name"
        last="$name"
    done
    printf '%s to %s\n' "$first" "$last"
}
# Lanes timed at the headline tier only, as "fqtools".
headline_only() { local name out=""; for name in "${LANES[@]}"; do [[ "${LANE_HEADLINE[$name]}" == 0 ]] || out+="${REG_LABEL[$name]}, "; done; printf '%s\n' "${out%, }"; }
# The inputs behind every figure: one line per kind, then the scaling sources.
data_list() {
    local kind src reals=0 shapes=0 lines="" synth="" used
    used=" $(cases | awk -F'\t' -v f="${FAMILIES[0]}" -v t="$HEADLINE" '$1 == f && $8 == t { print $7 }' | sort -u | tr '\n' ' ')"
    for kind in "${GROUP_LIST[@]}"; do
        local real=""
        for src in ${GROUP_MEMBERS[$kind]}; do
            [[ "$used" == *" $src "* ]] || continue
            if [[ -d "$SHAPES/$src" ]]; then
                synth+="    - \`$src\`: *${DESCRIBE[$src]}*"$'\n' shapes=$((shapes + 1))
            else
                real+="${DESCRIBE[$src]}; " reals=$((reals + 1))
            fi
        done
        [[ -z "$real" ]] || lines+="    - \`$kind\`: *${real%; }*"$'\n'
    done
    printf -- '- `data`:\n'
    printf -- '  - **Real**: *%d public sequencing runs, %s each, pooled in the bars by the kind of reads they hold*\n%s' "$reals" "${TIER_SIZE[$HEADLINE]}" "$lines"
    printf -- '  - **Synthetic**: *%d %s generated from real reads, each changing one thing; one group in the bars at %s, and the scaling lines from %s*\n%s' "$shapes" "$( ((shapes == 1)) && echo shape || echo shapes)" "${TIER_SIZE[$HEADLINE]}" "$(scaling_range)" "$synth"
}
# The source revisions of the z-fastq binaries the current batches timed, as
# recorded when each was built, with batch counts when there are several.
zfastq_revision() {
    local dir
    awk -F'\t' '{ print $3 }' "$RESULTS/current.tsv" | sort -u | while IFS= read -r dir; do
        sed -n 's/^lane z-fastq \([^#]*\)#.*/\1/p' "$dir/meta.tsv"
    done | awk -F'\t' -v revs="$CACHE/z-fastq.revisions" 'BEGIN { while ((getline line < revs) > 0) { split(line, f, "\t"); rev[f[1]] = f[2] } }
        { n[($1 in rev) ? rev[$1] : "revision not recorded"]++ }
        END { for (r in n) { k++; list[k] = r }
            if (k == 1) { print (list[1] == "revision not recorded" ? "" : "at ") list[1]; exit }
            for (i = 1; i <= k; i++) printf "%s%s (%d batches)", (i > 1 ? "; " : "at "), list[i], n[list[i]]; print "" }'
}
# The README's sections, each printing its Markdown; readme() puts them in
# the default order, and a feature with its own feature_readme orders them.
readme_head() {
    local first="${FAMILIES[0]}" name
    local date rev
    date="$(date -u +%Y-%m-%d)"
    rev="$(zfastq_revision)"
    local counters reals shapes
    counters=0
    for name in "${LANES[@]}"; do
        case "${LANE_CLASS[$name]}" in complete | partial) counters=$((counters + 1)) ;; esac
    done
    # Sources of the first family at the headline tier: a shape has a synthetic directory.
    local used src fmts
    used="$(cases | awk -F'\t' -v f="$first" -v t="$HEADLINE" '$1 == f && $8 == t && $9 != "many" { print $7 }' | sort -u)"
    reals=0 shapes=0
    for src in $used; do if [[ -d "$SHAPES/$src" ]]; then shapes=$((shapes + 1)); else reals=$((reals + 1)); fi; done
    fmts="$(cases | awk -F'\t' '$3 == "gzip" { g = 1 } END { print g ? "plain and gzip" : "plain" }')"
    local cpu
    cpu="$(awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo | sed 's/^AMD //; s/^Intel(R) Core(TM) //; s/ [0-9]*-Core Processor$//')"
    cat <<EOF
<!-- markdownlint-disable MD013 MD033 MD041 -->

# ${WORD[title]}
> <sub>\`${WORD[subject]}\` against $counters other ${WORD[others]}, on $reals real sequencing run$( ((reals == 1)) || echo s) and $shapes synthetic shape$( ((shapes == 1)) || echo s), $fmts. A tool is timed only where ${WORD[timed]}; time and memory come from the same run, on one CPU.</sub>

- \`z-fastq\`: *$("$ZFASTQ" --version | awk '{ print $2 }') $rev*
- \`machine\`: *$cpu, Linux $(uname -r | cut -d- -f1); one CPU per run*
- \`runner\`: *$(zebrac --version) with zexec*
$(data_list)
- \`date\`: *$date*

EOF
}
readme_tools() {
    local name
    cat <<EOF
## Compared tools
> <sub>*complete* tools check as much as z-fastq; *partial* tools check less, so they do less work. ${WORD[probes]}</sub>

EOF
    for name in "${LANES[@]}"; do
        local what="${LANE_NOTE[$name]}" probe=""
        case "${LANE_CLASS[$name]}" in
            subject | complete) what="checks every record's structure${what:+; $what}" ;;
        esac
        probe="; **$(probes "$name") rejected**"
        printf -- '- `%s %s` **(%s)**; *%s*%s\n' "${REG_LABEL[$name]}" "$(lane_version "$name")" "${LANE_CLASS[$name]}" "$what" "$probe"
    done
    if [[ " ${LANE_CLASS[*]} " == *" partial "* ]]; then
        printf '\n> [!IMPORTANT]\n> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.\n'
    fi
}
readme_right() {
    printf '\n## Right answers\n> <sub>%s The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>\n\n' "${WORD[truth]}"
    correctness_sentence
    printf '\n'
    caution
    picture "edges" "What each tool does on small edge-case inputs"
    printf '\n'
    printf '> [!NOTE]\n> %s "tool policy" rows have no single right answer.\n\n' "${WORD[cells]}"
    if awk -F'\t' '$3 == "fail" { found = 1 } END { exit !found }' "$REPORT_DIR/checks.tsv"; then
        details "Inputs where a tool ${WORD[missed]}"
        {
            printf 'Tool\tInput\tExit\tPrinted\n'
            printf 'l\tl\tr\tl\n'
            # A verdict feature prints its verdict, the first fact; a row seen in
            # two families shows once.
            awk -F'\t' -v repo="$REPO/" -v verdict="${EXPECT:+1}" 'FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
                $3 == "fail" { sub(repo, "", $6); v = $5; if (verdict) sub(/;.*/, "", v); printf "%s\t%s\t%s\t%s\n", label[$1], $6, $4, (v == "" ? "nothing" : v) }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/checks.tsv" | sort -u
        } | mdtable
        printf '\n</details>\n\n'
    fi
}
readme_speed() { # family [title]
    local fam="$1"
    printf '## %s\n' "${2:-Speed and memory}"
    printf '> <sub>The real runs and synthetic shapes at the %s tier, about %s of FASTQ each, pooled into kinds of input. Every tool runs on one CPU with the file in the page cache.</sub>\n\n' "$HEADLINE" "${TIER_SIZE[$HEADLINE]}"
    picture "$fam-overview" "Time and peak RSS relative to z-fastq"
    printf '\n'
    subject_sentence "$fam"
    printf '\n'
    summary_table "$fam"
    printf '\n'
    picture "$fam-speed" "Wall time by kind of input" && printf '\n'
    picture "$fam-memory" "Peak memory by kind of input" && printf '\n'
    if [[ -e "$RESULTS/figures/$fam-faults-light.png" ]]; then
        details "Minor page faults by kind of input"
        picture "$fam-faults" "Minor page faults by kind of input"
        printf '\n</details>\n\n'
    fi
}
readme_scaling() { # family [heading level]
    local fam="$1" slow
    slow="$(headline_only)"
    if [[ "$SCALING" == *" "* ]]; then
        printf '%s Scaling\n> <sub>Each synthetic shape from %s of FASTQ: throughput and peak memory against input size. A shape changes one thing %s pays for, so its line shows what that one thing costs each tool.</sub>\n\n' "${2:-##}" "$(scaling_range)" "${WORD[payer]}"
    else
        printf '%s Scaling\n> <sub>The `%s` shape (%s) from %s of FASTQ: throughput and peak memory against input size.</sub>\n\n' "${2:-##}" "$SCALING" "${DESCRIBE[$SCALING]}" "$(scaling_range)"
    fi
    picture "$fam-scaling" "Throughput and peak memory by input size" && printf '\n'
    picture "$fam-scaling-speed" "Throughput by input size" && printf '\n'
    picture "$fam-scaling-memory" "Peak memory by input size" && printf '\n'
    scaling_sentences "$fam"
    printf '\n'
    if [[ -n "$slow" || -n "${REMARK[scaling]:-}" ]]; then
        printf '> [!NOTE]\n> %s\n\n' "${slow:+$slow is timed at the $HEADLINE tier only: it is slow enough that sweeping sizes with it would take most of the run. }${REMARK[scaling]:-}"
    fi
}
readme_forms() {
    local fam
    printf '## %s\n> <sub>%s</sub>\n\n' "${WORD[forms]}" "${WORD[forms_sub]:-The same ${WORD[same]}, with the input on stdin or as many files in one call.}"
    for fam in "${FAMILIES[@]:1}"; do
        printf '### %s\n\n%s.\n\n' "$(tr '-' ' ' <<<"${fam^}")" "${FAM_DESC[$fam]^}"
        if [[ "${FAM_KIND[$fam]}" == many ]]; then
            many_sentence "$fam"
            printf '\n'
        elif [[ "${FAM_STDIN[$fam]}" == 1 ]]; then
            printf "Each tool's time and peak memory reading the file from a pipe, as a multiple of its own with the path (geometric mean over the shared shapes):\n\n"
            stdin_table "$fam"
            printf '\n'
        else
            subject_sentence "$fam"
            printf '\n'
            summary_table "$fam"
            printf '\n'
        fi
    done
}
# The families with their own memory cap, as ", 8G for deinterleave".
mem_caps() { local fam out=""; for fam in "${FAMILIES[@]}"; do [[ -z "${FAM_MEM[$fam]:-}" ]] || out+=", ${FAM_MEM[$fam]} for $fam"; done; printf '%s' "$out"; }
readme_repro() {
    local name fam date rev
    date="$(date -u +%Y-%m-%d)"
    rev="$(zfastq_revision)"
    printf '## Reproducibility\n\n'
    details "How it was measured"
    cat <<EOF
Every tool runs through \`zexec -r\` inside one zebrac batch per case, held to one CPU core by a systemd scope with a CPU quota and a $BENCH_MEM memory cap$(mem_caps), after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. The kernel reports peak memory in 256 KiB steps, so values closer than that are equal. Runs per tier, each after one warmup: tiny ${RUNS[tiny]}, xsmall ${RUNS[xsmall]}, small ${RUNS[small]}, medium ${RUNS[medium]}; sizes 1 MiB to 1 GiB. Generated by \`bench/$FEATURE/run.sh release\` on $date.

</details>

EOF
    details "Identity"
    printf -- '- `z-fastq`: *%s %s*\n' "$("$ZFASTQ" --version | awk '{ print $2 }')" "$rev"
    printf -- '- `zebrac`: *%s*\n' "$(zebrac --version | awk '{ print $2 }')"
    printf -- '- `zexec`: *%s bytes, static, no libc*\n' "$(stat -c %s "$ZEXEC")"
    printf -- '- `kernel`: *%s*\n' "$(uname -r)"
    printf -- '- `cpu`: *%s*\n' "$(awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo)"
    printf -- '- `load at batch start`: *%s to %s*\n' "$(awk -F'\t' 'NR > 1 { print $21 }' "$REPORT_DIR/summary.tsv" | sort -n | head -1)" "$(awk -F'\t' 'NR > 1 { print $21 }' "$REPORT_DIR/summary.tsv" | sort -n | tail -1)"
    printf '\n</details>\n\n'
    details "Commands"
    local any=""
    for fam in "${FAMILIES[@]}"; do any+="${FAM_PARAMS[$fam]}${FAM_VARY[$fam]}"; done
    if [[ -z "$any" ]]; then
        for name in "${LANES[@]}"; do printf -- '- `%s`: *%s*\n' "${REG_LABEL[$name]}" "$(lane_cmd "$name" IN | sed "s#$REPO/##g")"; done
        for fam in "${FAMILIES[@]}"; do
            for name in ${FAM_LANES[$fam]}; do
                [[ -v "FAM_TMPL[$fam|$name]" ]] || continue
                printf -- '- `%s`, %s: *%s*\n' "${REG_LABEL[$name]}" "$fam" "$(lane_cmd "$name" IN "$fam" | sed "s#$REPO/##g")"
            done
        done
    else
        # Every family's command, with the parameters of its first case; a
        # named output is the null device when timed.
        local params
        for fam in "${FAMILIES[@]}"; do
            params="$(cases | awk -F'\t' -v f="$fam" '$1 == f && !p { p = $10 } END { print p }')"
            for name in "${LANES[@]}"; do
                [[ " ${FAM_LANES[$fam]} " == *" $name "* ]] || continue
                printf -- '- `%s`, %s: *%s*\n' "${REG_LABEL[$name]}" "$fam" "$(lane_cmd "$name" IN "$fam" "$params" | sed "s#$REPO/##g")"
            done
        done
    fi
    printf '\n</details>\n\n'
    if [[ -f "$RESULTS/picks.tsv" ]]; then
        details "Quick-check cases for peer and ab"
        if ((RIVAL_ANY)); then
            printf 'Picked from these results: per family, group, and format, the case where z-fastq is closest to, or furthest behind, its fastest rival of any class. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.\n\n'
        else
            printf 'Picked from these results: per group and format, the case where z-fastq leads its fastest full-check rival by the least. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.\n\n'
        fi
        awk -F'\t' -v many="$(($(wc -w <<<"$(full_families)") > 1))" 'FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
            FILENAME ~ /cases.tsv$/ { if (FNR > 1) name[$1 "\t" $2] = $8; next }
            FNR > 1 { k = $1 "\t" $2
                printf "- `%s%s, %s`: *%s*", (many ? $1 ": " : ""), $3, $4, (k in name ? name[k] : $2)
                if ($5 != "-") printf "; rival %s at %.2fx z-fastq%ss time", label[$5], $7, "\047"
                if ($6 != "-") printf "; memory rival %s", label[$6]
                if ($8 == "declared") printf "; declared"
                printf "\n" }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/cases.tsv" "$RESULTS/picks.tsv"
        printf '\n</details>\n\n'
    fi
    details "Hardware counters (user space, plain $HEADLINE-tier inputs, geometric means)"
    counters_table "${FAMILIES[0]}"
    printf '\n</details>\n'
}
readme() {
    if declare -F feature_readme >/dev/null; then feature_readme; return; fi
    readme_head
    printf '## Findings\n\n%s\n' "$(findings "${FAMILIES[0]}")"
    valid_wrong_bullet
    printf '\n'
    readme_tools
    readme_right
    if declare -F feature_report >/dev/null; then feature_report; fi
    readme_speed "${FAMILIES[0]}"
    readme_scaling "${FAMILIES[0]}"
    readme_forms
    readme_repro
}
report() {
    STAGE=report
    REPORT_DIR="$RESULTS/report"
    local figures="$RESULTS/figures" stamp key f
    mkdir -p -- "$REPORT_DIR" "$figures"
    [[ -d "$RESULTS/batches" ]] || die "no timing run yet: run.sh release"
    collect >"$REPORT_DIR/summary.tsv"
    style_tsv >"$REPORT_DIR/style.tsv"
    cases_tsv >"$REPORT_DIR/cases.tsv"
    edges_tsv >"$REPORT_DIR/edges.tsv"
    current_checks >"$REPORT_DIR/checks.tsv"
    if [[ -f "$RESULTS/edges.tsv" ]]; then lets_through >"$REPORT_DIR/through.tsv"; else : >"$REPORT_DIR/through.tsv"; fi
    printf 'wrong\t%s\ndone\t%s\nsize\t%s\n' "${WORD[legend]}" "${WORD[past]}" "${TIER_SIZE[$HEADLINE]}" >"$REPORT_DIR/words.tsv"
    # Figures depend only on these files, plot.py, and the arguments; redraw
    # only when one of them changed or a figure is missing.
    local mains lines="" fam
    mains="$(full_families)" mains="${mains% }" mains="${mains// /,}"
    for fam in "${FAMILIES[@]}"; do [[ "${FAM_LOOK[$fam]}" != lines ]] || lines+="$fam:${FAM_VARY[$fam]%% *},"; done
    key="$(cat "$REPORT_DIR"/{summary,style,cases,edges,words}.tsv "$BENCH/lib/plot.py" | cksum) $HEADLINE $SCALING $mains${lines:+ ${lines%,}}"
    if [[ "$(cat "$REPORT_DIR/figures.key" 2>/dev/null)" != "$key" ]] || ! ls "$figures"/*.png >/dev/null 2>&1; then
        stamp="$REPORT_DIR/stamp"
        : >"$stamp"
        "$REPO/tools/venv/bin/python" -I "$BENCH/lib/plot.py" "$REPORT_DIR" "$figures" "$HEADLINE" "$SCALING" "$mains" ${lines:+"${lines%,}"}
        # A figure the plot step did not just write belongs to nothing declared now.
        for f in "$figures"/*.svg "$figures"/*.png; do [[ "$f" -nt "$stamp" || ! -e "$f" ]] || rm -- "$f"; done
        printf '%s\n' "$key" >"$REPORT_DIR/figures.key"
    fi
    write_picks
    readme >"$FEATURE_DIR/README.md.part"
    mv -- "$FEATURE_DIR/README.md.part" "$FEATURE_DIR/README.md"
    log report DONE "${FEATURE_DIR#"$REPO"/}/README.md, $(ls "$figures"/*.png | wc -l) figures"
}

# --- quick checks: picks, peer, ab ------------------------------------------
# The machine a cached reference was timed on: a reference from another host,
# CPU, or kernel is never reused.
host_id() { printf '%s|%s|%s\n' "$(hostname)" "$(awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo)" "$(uname -r)"; }

# The families with the full report look: every family declared look full,
# else the first family.
full_families() {
    local fam out=""
    for fam in "${FAMILIES[@]}"; do [[ "${FAM_LOOK[$fam]}" != full ]] || out+="$fam "; done
    printf '%s\n' "${out:-${FAMILIES[0]}}"
}
# picks.tsv: the quick-check cases, written by every report from the release
# results. Per bar group and format at the headline tier of each full family,
# the case where z-fastq's lead over its speed rival (the fastest complete
# lane, or the fastest lane when the family has no complete one) is smallest; the memory rival is the complete lane with the lowest
# peak RSS, kept only where it is below z-fastq's. Columns: family, case,
# group, format, speed rival, memory rival, rival time over z-fastq's, reason.
write_picks() {
    local out="$RESULTS/picks.tsv" fam="${FAMILIES[0]}" pf
    {
        printf 'family\tcase\tgroup\tformat\trival\tmemory_rival\tlead\treason\n'
        if [[ "$PICK_MODE" != only ]]; then
            for pf in $(full_families); do
            awk -F'\t' -v fam="$pf" -v tier="$HEADLINE" -v any="$RIVAL_ANY" '
                FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") { fmt[$2] = $3; grp[$2] = $10; ord[$2] = $9 + 0 } next }
                FILENAME ~ /style.tsv$/ { if (FNR > 1) cls[$1] = $3; next }
                FNR == 1 || $1 != fam || !($2 in fmt) { next }
                { w[$2, $3] = $5; m[$2, $3] = $8; lanes[$3] = 1 }
                END {
                    for (c in fmt) {
                        if (!((c, "z-fastq") in w)) continue
                        r = mr = ""
                        for (l in lanes) if (cls[l] == "complete" && (c, l) in w) {
                            if (r == "" || w[c, l] < w[c, r]) r = l
                            if (m[c, l] < m[c, "z-fastq"] && (mr == "" || m[c, l] < m[c, mr])) mr = l
                        }
                        # With no complete lane in the family, or with rival any, the rival is the fastest other lane.
                        if (any) r = mr = ""
                        if (r == "") for (l in lanes) if (l != "z-fastq" && (c, l) in w) {
                            if (r == "" || w[c, l] < w[c, r]) r = l
                            if (m[c, l] < m[c, "z-fastq"] && (mr == "" || m[c, l] < m[c, mr])) mr = l
                        }
                        if (r == "") continue
                        lead = w[c, r] / w[c, "z-fastq"]; k = grp[c] SUBSEP fmt[c]
                        if (!(k in best) || lead < bl[k] || (lead == bl[k] && ord[c] < ord[best[k]])) { best[k] = c; bl[k] = lead; br[k] = r; bm[k] = mr }
                    }
                    for (k in best) { split(k, p, SUBSEP); printf "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\tsmallest lead\n", fam, best[k], p[1], p[2], br[k], (bm[k] == "" ? "-" : bm[k]), bl[k] }
                }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | sort -t$'\t' -k4,4r -k3,3
            done
        fi
        local c
        for c in $PICK_CASES; do printf '%s\t%s\t-\t%s\t-\t-\t-\tdeclared\n' "$fam" "$c" "${c##*:}"; done
    } >"$out.part"
    mv -- "$out.part" "$out"
}
# The picks as case rows: family, case, and the case's fields from cases.
pick_rows() {
    [[ -f "$RESULTS/picks.tsv" ]] || die "no quick-check cases yet: run.sh release first (or declare pick only CASES)"
    awk -F'\t' 'FNR == NR { if (FNR > 1) { want[$1 "\t" $2] = 1; rival[$1 "\t" $2] = $5; mem[$1 "\t" $2] = $6 } next }
        ($1 "\t" $2) in want { print $0 "\t" rival[$1 "\t" $2] "\t" mem[$1 "\t" $2] }' "$RESULTS/picks.tsv" <(cases)
}
# Checks z-fastq (or BIN) on one case before it is timed; dies on a wrong count.
quick_check() { # name fam id fmt paths truth bytes params [bin]
    local name="$1" fam="$2" id="$3" paths="$5" params="$8" bin="${9:-}" form=file out truth
    [[ "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
    truth="$(want_of "$paths" "$6")" || die "no truth for $id: run.sh check first"
    if [[ -n "$bin" ]]; then
        out="$($(build_cmd "$bin" "$paths" "$fam" "$params") 2>/dev/null | LANE=z-fastq FAMILY="$fam" IN="$paths" TRUTH="$6" PARAMS="$params" fact | paste -sd';' -)" || die "$name failed on $id"
        same_answer "$out" "$truth" || die "$name printed a wrong answer on $id"
        return 0
    fi
    local key
    key="$(status_key "$name" "$fam" "$form" "$(fileid $paths)" "$params")"
    if [[ ! -v "CHECK_STATUS[$key]" ]]; then
        check_one "$name" "$form" "$fam" "$id" "$4" "$paths" "$truth" "$7" "$params" "$6" >>"$CHECKS"
        load_checks
    fi
    [[ "${CHECK_STATUS[$key]}" == pass ]] || die "$name gives a wrong answer on $id"
}

# peer: z-fastq against each pick's rivals. A rival's reference is timed once
# with long settings and kept while its binary, the input, zexec, the CPU, and
# the host are unchanged; z-fastq is timed with the headline settings and kept
# while its binary is unchanged. A rival's reference and z-fastq are timed at
# different moments, so the verdict compares interquartile ranges, and the
# load at both timings is shown.
PEER_RUNS=31 PEER_WARM=3
run_peer() {
    cache_lane_ids
    load_checks
    ZEXEC_ID="$(fileid "$ZEXEC")"
    local host fam id fmt paths truth bytes src tier kind params rival mem form inid runs warm dir lanes
    host="$(host_id)"
    local out="$RESULTS/peer.txt"
    printf '%-40s %-11s %9s %9s %7s %-12s %8s %-8s %s\n' case rival z-fastq rival "time" verdict memory verdict load | tee "$out"
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params rival mem; do
        [[ "$rival" != - ]] || continue
        form=file; [[ "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
        inid="$(fileid $paths)"
        quick_check z-fastq "$fam" "$id" "$fmt" "$paths" "$truth" "$bytes" "$params"
        lanes="$rival"; [[ "$mem" == - || "$mem" == "$rival" ]] || lanes+=" $mem"
        BATCH_FAM="$fam"
        # The rivals' reference: long settings, kept per host.
        batch_meta "$form" "$inid" "$PEER_RUNS" "$PEER_WARM" $lanes
        [[ "$params" == - ]] || META+="params $params"$'\n'
        META+="host $host"$'\n'
        if ! find_batch "$RESULTS/peer/ref/$fam-${id//:/-}-"; then
            local l
            SPECS=()
            for l in $lanes; do timed_spec "$l" "$l" "$([[ "$form" == stdin ]] && echo - || echo "$paths")" "$fam" "$params"; done
            BATCH="$RESULTS/peer/ref/$fam-${id//:/-}-$(printf '%s' "$META" | cksum | cut -d' ' -f1)"
            time_lanes "$BATCH" "$form" "$paths" "$PEER_RUNS" "$PEER_WARM" "$fam" "$id" "${SPECS[@]}"
        fi
        local ref="$BATCH/summary.tsv"
        # z-fastq now: headline settings, kept while its binary is unchanged.
        read -r runs warm <<<"$(batch_settings "$tier" "$kind")"
        batch_meta "$form" "$inid" "$runs" "$warm" z-fastq
        [[ "$params" == - ]] || META+="params $params"$'\n'
        META+="host $host"$'\n'
        if ! find_batch "$RESULTS/peer/now/$fam-${id//:/-}-"; then
            BATCH="$RESULTS/peer/now/$fam-${id//:/-}-$(printf '%s' "$META" | cksum | cut -d' ' -f1)"
            SPECS=()
            timed_spec z-fastq z-fastq "$([[ "$form" == stdin ]] && echo - || echo "$paths")" "$fam" "$params"
            time_lanes "$BATCH" "$form" "$paths" "$runs" "$warm" "$fam" "$id" "${SPECS[@]}"
        fi
        META=""
        awk -F'\t' -v id="$id" -v rival="$rival" -v mem="$mem" '
            { w[$3] = $5; q1[$3] = $6; q3[$3] = $7; r[$3] = $8; l0[$3] = $21 }
            END {
                v = (q3["z-fastq"] < q1[rival]) ? "faster" : (q1["z-fastq"] > q3[rival]) ? "slower" : "within noise"
                mr = (mem == "-") ? rival : mem
                rv = (r["z-fastq"] + 256 <= r[mr]) ? "smaller" : (r["z-fastq"] >= r[mr] + 256) ? "larger" : "same"
                printf "%-40s %-11s %7.2fms %7.2fms %6.2fx %-12s %7.2fx %-8s %s/%s\n", id, rival, w["z-fastq"] / 1e6, w[rival] / 1e6,
                    w[rival] / w["z-fastq"], v, r[mr] / r["z-fastq"], rv, l0["z-fastq"], l0[rival]
            }' "$ref" "$BATCH/summary.tsv" | tee -a "$out"
    done < <(pick_rows)
    note "peer: time is the rival's over z-fastq's (above 1x: z-fastq faster); faster or slower only when the interquartile ranges do not overlap; memory is the memory rival's (or the speed rival's) peak RSS over z-fastq's; load is at z-fastq's and the reference's timing"
}

# ab: z-fastq at REV_A against REV_B (default: the working tree) on the picks,
# with A timed twice as the A/A control. Always timed fresh: A and B are only
# comparable inside one batch. The last run is kept in results/ab/.
build_rev() { # rev -> binary path
    local sha bin wt="$CACHE/worktree"
    sha="$(git -C "$REPO" rev-parse --verify "$1^{commit}")" || die "unknown revision $1"
    bin="$CACHE/bin/$sha/z-fastq"
    if [[ ! -x "$bin" ]]; then
        if [[ ! -d "$wt" ]]; then git -C "$REPO" worktree add --detach "$wt" "$sha" >&2; else git -C "$wt" checkout -q --detach "$sha"; fi
        log build RUN "z-fastq at ${sha:0:9}"
        (cd "$wt" && systemd-run --user --scope -q -p CPUQuota=400% -p MemoryMax=4G zig build -Dstatic=true -Doptimize=ReleaseFast -j4) >&2
        mkdir -p -- "$(dirname -- "$bin")"
        install -m 755 "$wt/zig-out/bin/z-fastq" "$bin"
    fi
    printf '%s\n' "$bin"
}
run_ab() {
    local rev_a="$1" rev_b="${2:-}" a b fam id fmt paths truth bytes src tier kind params rival mem runs warm dir
    cache_lane_ids
    load_checks
    a="$(build_rev "$rev_a")"
    if [[ -n "$rev_b" ]]; then b="$(build_rev "$rev_b")"; else b="$ZFASTQ"; fi
    local out="$RESULTS/ab.txt"
    printf '%-40s %10s %10s %8s %8s %-14s %10s %10s %-8s\n' case "A" "B" change "A/A" verdict "A RSS" "B RSS" "RSS" | tee "$out"
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params rival mem; do
        quick_check A "$fam" "$id" "$fmt" "$paths" "$truth" "$bytes" "$params" "$a"
        quick_check B "$fam" "$id" "$fmt" "$paths" "$truth" "$bytes" "$params" "$b"
        read -r runs warm <<<"$(batch_settings "$tier" "$kind")"
        dir="$RESULTS/ab/$fam-${id//:/-}"
        META=""
        SPECS=()
        timed_spec a1 z-fastq "$paths" "$fam" "$params" "$a"
        timed_spec a2 z-fastq "$paths" "$fam" "$params" "$a"
        timed_spec b z-fastq "$paths" "$fam" "$params" "$b"
        time_lanes "$dir" file "$paths" "$runs" "$warm" "$fam" "$id" "${SPECS[@]}"
        awk -F'\t' -v id="$id" '
            { w[$3] = $5; r[$3] = $8 }
            END {
                a = (w["a1"] + w["a2"]) / 2; ch = (w["b"] - a) / a * 100; sp = (w["a1"] > w["a2"] ? w["a1"] - w["a2"] : w["a2"] - w["a1"]) / a * 100
                # One A/A pair underestimates the noise on a busy machine: a change
                # counts only beyond twice that spread and at least 3%.
                v = (ch < -2 * sp && ch < -3) ? "faster" : (ch > 2 * sp && ch > 3) ? "slower" : "within noise"
                ra = (r["a1"] + r["a2"]) / 2; rs = (r["a1"] > r["a2"] ? r["a1"] - r["a2"] : r["a2"] - r["a1"])
                rv = (r["b"] < ra - rs - 256) ? "smaller" : (r["b"] > ra + rs + 256) ? "larger" : "same"
                printf "%-40s %8.2fms %8.2fms %+7.1f%% %7.1f%% %-14s %7dKiB %7dKiB %-8s\n", id, a / 1e6, w["b"] / 1e6, ch, sp, v, ra, r["b"], rv
            }' "$dir/summary.tsv" | tee -a "$out"
    done < <(pick_rows)
    note "ab: A=${rev_a} B=${rev_b:-working tree}; a change is real only beyond twice the A/A spread and 3% (time), or one 256 KiB step beyond the A/A spread (RSS)"
}

# --- guards -------------------------------------------------------------------
# Declaration mistakes stop every mode at once, before any work: an unknown
# name, a source without a description or group, a scaling source with one
# size, a pick that names no case.
validate() {
    local errs=() fam name src t list used
    [[ -n "$JOB" && -n "$TRUTH" ]] || errs+=("job and truth must be declared")
    [[ -v "RUNS[$HEADLINE]" ]] || errs+=("headline tier $HEADLINE has no run settings")
    ((${#FAMILIES[@]})) || errs+=("no family declared")
    for fam in "${FAMILIES[@]}"; do
        [[ -n "${FAM_LANES[$fam]:-}" ]] || errs+=("family $fam declares no lanes")
        [[ -n "${FAM_INPUTS[$fam]:-}" ]] || errs+=("family $fam declares no inputs")
        for name in ${FAM_LANES[$fam]:-}; do [[ -v "LANE_TMPL[$name]" ]] || errs+=("family $fam: lane $name is not declared"); done
        [[ "${FAM_LOOK[$fam]}" =~ ^(|full|table|lines)$ ]] || errs+=("family $fam: unknown look ${FAM_LOOK[$fam]}")
        [[ "${FAM_LOOK[$fam]}" != lines || -n "${FAM_VARY[$fam]}" ]] || errs+=("family $fam: look lines needs a vary line")
        # Every placeholder in a family's commands is a core one or a family parameter.
        local known=" zfastq bin in in1 in2 out1 out2 stem dir name ${FAM_VARY[$fam]%% *} " kv ph tmpl
        for kv in ${FAM_PARAMS[$fam]}; do known+="${kv%%=*} "; done
        for name in ${FAM_LANES[$fam]:-}; do
            [[ -v "LANE_TMPL[$name]" ]] || continue
            tmpl="$(lane_tmpl "$name" "$fam")"
            for ph in $(grep -o '{[a-z0-9_]*}' <<<"$tmpl" | tr -d '{}' || true); do
                [[ "$known" == *" $ph "* ]] || errs+=("family $fam: lane $name uses {$ph}, which no param line declares")
            done
        done
    done
    list="$(cases)" || errs+=("the declared inputs do not resolve (see the error above)")
    used="$(awk -F'\t' -v f="${FAMILIES[0]}" -v t="$HEADLINE" '$1 == f && $8 == t && $9 != "many" { print $7 }' <<<"$list" | sort -u)"
    for src in $used; do
        [[ -v "DESCRIBE[$src]" ]] || errs+=("source $src has no describe line")
        [[ -v "GROUP_OF[$src]" ]] || errs+=("source $src is in no group")
    done
    for src in $SCALING; do
        t="$(awk -F'\t' -v f="${FAMILIES[0]}" -v s="$src" '$1 == f && $7 == s { print $8 }' <<<"$list" | sort -u | wc -l)"
        ((t >= 2)) || errs+=("scaling source $src is declared at $t size(s); scaling needs two or more")
    done
    [[ "$PICK_MODE" =~ ^(auto|add|only)$ ]] || errs+=("pick must be auto, add CASES, or only CASES")
    [[ "$PICK_MODE" != only || -n "$PICK_CASES" ]] || errs+=("pick only needs cases")
    for t in $PICK_CASES; do grep -q $'\t'"$t"$'\t' <<<"$list" || errs+=("pick: no case $t"); done
    for t in "${!REMARK[@]}"; do [[ "$t" == scaling ]] || errs+=("remark: unknown section $t"); done
    [[ -z "$AGREE_LANE" || -v "LANE_TMPL[$AGREE_LANE]" ]] || errs+=("agree: lane $AGREE_LANE is not declared")
    [[ -z "$EDGE_SET" || "$EDGE_SET" =~ ^(count|stats|check|sample)$ ]] || declare -F feature_edge_list >/dev/null || errs+=("edges: unknown set $EDGE_SET")
    for fam in $EDGE_FAMS; do [[ -v "FAM_LANES[$fam]" ]] || errs+=("edges: no family $fam"); done
    ((${#errs[@]} == 0)) && return 0
    printf 'run.sh: %s\n' "${errs[@]}" >&2
    exit 2
}
# One timing job at a time on this machine: every feature times on the same
# CPU, so two runs at once would slow each other. The lock is released when
# the process exits.
take_lock() {
    mkdir -p -- "$CACHE"
    # One timing run per CPU: runs pinned to different CPUs may overlap.
    exec 9>"$CACHE/timing-cpu$BENCH_CPU.lock"
    flock -n 9 || die "another benchmark run holds $CACHE/timing-cpu$BENCH_CPU.lock; wait for it or use another BENCH_CPU"
}
usage() {
    cat <<USAGE
bench/$FEATURE/run.sh MODE: $JOB

  init           check tools, versions, data, and the machine; no timing (seconds)
  check          check every lane's output against the truth, and the edge cases;
                 only new lane and input pairs run (cached in results/checks.tsv)
  release        check, time every case not yet cached, then report; a repeat with
                 nothing new takes seconds, a rebuilt z-fastq re-times every case
  report         README.md and figures from cached results only; never times
  peer           z-fastq against each quick-check case's rival (results/picks.tsv,
                 written by report): rival references are timed once and cached,
                 z-fastq is re-timed only when rebuilt (results/peer.txt);
                 a pointer, not proof: the two are timed at different moments
  ab REV [REV2]  two z-fastq builds and an A/A control on the quick-check cases,
                 always timed fresh in one batch (results/ab.txt); REV2 defaults
                 to the working tree; the acceptance check for a z-fastq change
  list           every timed case: family, id, format, paths, truth, bytes, ...
  list-verify    the extra inputs z-fastq alone is checked on
  help           this text

Timing modes (check, release, peer, ab) hold one machine-wide lock. Results
and caches live in bench/$FEATURE/results/ and bench/cache/; both are ignored
by git and reused, never deleted by these modes. RETIME=1 forces new batches.
USAGE
}
bench_main() {
    local mode="${1:-}"
    mapfile -t LANES < <(for n in "${LANES[@]}"; do printf '%s %s\n' "${REG_ORDER[$n]}" "$n"; done | sort -n | cut -d' ' -f2)
    case "$mode" in
        "" | help | -h | --help) usage; return 0 ;;
        init | check | release | report | list | list-verify | peer | ab) validate ;;
        *) usage >&2; exit 2 ;;
    esac
    case "$mode" in check | release | peer | ab) take_lock ;; esac
    case "$mode" in check | release | report | peer | ab) run_start "$mode" ;; esac
    case "$mode" in
        init)
            run_init
            ((STOPS == 0)) || exit 2
            # Notes inform; a lane issue (a peer missing or off its pin) fails.
            [[ "${ISSUES[*]}" != *"lane: "* ]] || exit 1
            ;;
        check)
            run_init
            ((STOPS == 0)) || exit 2
            ensure_builds
            run_checks
            run_edges
            if declare -F feature_check >/dev/null; then feature_check; fi
            ;;
        release)
            run_init
            ((STOPS == 0)) || exit 2
            ensure_builds
            run_checks
            run_edges
            if declare -F feature_check >/dev/null; then feature_check; fi
            run_timing
            report
            ;;
        report) report ;;
        list) cases ;;
        list-verify) verify_cases ;;
        peer)
            run_init
            ((STOPS == 0)) || exit 2
            ensure_builds
            run_peer
            ;;
        ab)
            [[ -n "${2:-}" ]] || die "usage: run.sh ab REV_A [REV_B]"
            run_init
            ((STOPS == 0)) || exit 2
            ensure_builds
            run_ab "$2" "${3:-}"
            ;;
        *) usage >&2; exit 2 ;;
    esac
}

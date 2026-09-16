#!/usr/bin/env bash
# Sample benchmark runner: contract checks, capability probes, then zebrac.
#
# Workloads come from features.tsv:
#   sample/<publication|small>/<se|paired|interleaved>
# Each catalog role is timed twice: --fraction and --count.
#
# seqtk is the screened SE exact-reference peer. Other samplers are descriptive.
# File-writing peers sink to /dev/null so zebrac still measures the tool, not a shell.
# IRMA Core is timed on SE and paired plain only. One interleaved file is SE
# records, not pair units. IRMA gzip is probed and not timed: gzip and plain
# select different records.
#
# Usage:
#   bash bench/sample/run.sh
#   bash bench/sample/run.sh --small-real --runs 5 --warmup 3
#   bash bench/sample/run.sh --skip-tests --skip-report
#   SAMPLE_RUN_TIMESTAMP=<ts> bash bench/sample/run.sh --skip-tests

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
PLAIN_DIR="$BENCH_ROOT/shared/cache/plain"
FIXTURE_DIR="$PROJECT_ROOT/tests/data/synthetic"

# shellcheck disable=SC1091
source "$BENCH_ROOT/shared/tools.sh"
catalog_load

RUNS=25
WARMUP=5
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_FULL=true
DO_REPORT=true
SMALL_REAL=false
ALLOW_INCOMPLETE=false

SAMPLE_SEED=11
SAMPLE_FRACTION=0.1
SAMPLE_COUNT=1000
# Fixture probes use these, not the timed SAMPLE_* values.
FIXTURE_FRACTION=0.5
FIXTURE_COUNT=3
FIXTURE_PERCENT=50

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-full) DO_FULL=false; shift ;;
        --skip-report) DO_REPORT=false; shift ;;
        --small-real) SMALL_REAL=true; shift ;;
        --allow-incomplete) ALLOW_INCOMPLETE=true; shift ;;
        -h|--help)
            sed -n '2,/^set -euo pipefail$/p' "$0" | head -n -1
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

SAMPLE_PERCENT="$(awk -v p="$SAMPLE_FRACTION" 'BEGIN {
    pct = int(p * 100 + 0.5)
    if (pct < 1 || pct > 100) exit 1
    printf "%d", pct
}')" || {
    echo "error: SAMPLE_FRACTION=$SAMPLE_FRACTION is not an IRMA -p percent in 1..100" >&2
    exit 1
}

ZEBRAC_MIN_SAMPLES="$RUNS"
ZEBRAC_WARMUP="$WARMUP"
export ZEBRAC_DURATION_MS ZEBRAC_MIN_SAMPLES ZEBRAC_WARMUP

if [[ -n "${SAMPLE_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$SAMPLE_RUN_TIMESTAMP"
else
    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
fi

mkdir -p "$RESULTS_DIR"
METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
VERIFY_LOG="$RESULTS_DIR/verify_${TIMESTAMP}.log"
MANIFEST_PATH="$RESULTS_DIR/run_${TIMESTAMP}.json"
: >"$METADATA_JSONL"

PYTHON="$(report_python)" || exit 1
CHECK_DIR="$(mktemp -d "$RESULTS_DIR/.check.XXXXXX")"
trap 'rm -rf -- "$CHECK_DIR"' EXIT

SAMPLE_SET="publication"
if $SMALL_REAL; then
    SAMPLE_SET="small"
fi

declare -a SE_IDS=()
declare -a PAIRED_IDS=()
declare -a INTERLEAVED_IDS=()
declare -a WORKLOAD_KEYS=()
declare -a WORKLOAD_ROWS=()
declare -A DATA_GZ=()
declare -A DATA_PLAIN=()
declare -A DATA_EXPECTED=()
declare -A DATA_DECODED=()
declare -A LANE_ENABLED=()
declare -A LANE_REASON=()
declare -a PEER_FIXTURE_ROWS=()

SAMPLE_PEER_ORDER=(seqtk seqkit rasusa fq fqkit irma bbtools)

json_string_array() {
    local sep="" value
    printf '['
    for value in "$@"; do
        printf '%s%s' "$sep" "$(zebrac_json_string "$value")"
        sep=','
    done
    printf ']'
}

json_number_array() {
    local sep="" value
    printf '['
    for value in "$@"; do
        printf '%s%s' "$sep" "$value"
        sep=','
    done
    printf ']'
}

workload_key() {
    local role="$1"
    shift
    local joined
    joined="$(IFS=+; printf '%s' "$*")"
    printf '%s__%s' "$role" "$joined"
}

register_workload() {
    local role="$1"
    shift
    local -a ids=("$@")
    local key id
    key="$(workload_key "$role" "${ids[@]}")"
    WORKLOAD_KEYS+=("$key")
    local -a accessions=()
    local -a records=()
    for id in "${ids[@]}"; do
        accessions+=("${CATALOG_ACCESSION[$id]}")
        records+=("${CATALOG_EXPECTED[$id]}")
    done
    WORKLOAD_ROWS+=("$(printf \
        '{"key":%s,"role":%s,"ids":%s,"accessions":%s,"records":%s,"plain":true,"gzip":true}' \
        "$(zebrac_json_string "$key")" \
        "$(zebrac_json_string "$role")" \
        "$(json_string_array "${ids[@]}")" \
        "$(json_string_array "${accessions[@]}")" \
        "$(json_number_array "${records[@]}")")")
}

workloads_json() {
    local sep="" row
    printf '['
    for row in "${WORKLOAD_ROWS[@]}"; do
        printf '%s%s' "$sep" "$row"
        sep=','
    done
    printf ']'
}

peer_fixtures_json() {
    local sep="" row
    printf '['
    for row in "${PEER_FIXTURE_ROWS[@]}"; do
        printf '%s%s' "$sep" "$row"
        sep=','
    done
    printf ']'
}

bind_dataset() {
    local id="$1"
    local filename
    filename="$(bench_catalog field "$id" filename)"
    DATA_GZ["$id"]="$DATA_DIR/$filename"
    DATA_PLAIN["$id"]="$PLAIN_DIR/${filename%.gz}"
    DATA_EXPECTED["$id"]="$(bench_catalog expected "$id")"
    DATA_DECODED["$id"]="$(file_size_bytes "${DATA_PLAIN[$id]}")"
}

load_role_ids() {
    local role="$1"
    local -n destination="$2"
    mapfile -t destination < <(
        bench_catalog ids --suite sample --set "$SAMPLE_SET" --role "$role" --no-expand
    )
    ((${#destination[@]} > 0)) || {
        echo "error: sample/$SAMPLE_SET/$role has no datasets" >&2
        exit 1
    }
}

ensure_real_data() {
    local download_args=(--suite sample)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (sample/$SAMPLE_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    load_role_ids se SE_IDS
    load_role_ids paired PAIRED_IDS
    load_role_ids interleaved INTERLEAVED_IDS

    local -a all_ids=(
        "${SE_IDS[@]}"
        "${PAIRED_IDS[@]}"
        "${INTERLEAVED_IDS[@]}"
    )
    local id
    for id in "${all_ids[@]}"; do
        bind_dataset "$id"
    done

    local se_id
    for se_id in "${SE_IDS[@]}"; do
        register_workload se_fraction "$se_id"
        register_workload se_count "$se_id"
    done
    register_workload paired_fraction "${PAIRED_IDS[@]}"
    register_workload paired_count "${PAIRED_IDS[@]}"
    register_workload interleaved_fraction "${INTERLEAVED_IDS[@]}"
    register_workload interleaved_count "${INTERLEAVED_IDS[@]}"

    local -a unique_ids=()
    declare -A seen_ids=()
    for id in "${all_ids[@]}"; do
        [[ -n "${seen_ids[$id]+present}" ]] && continue
        seen_ids["$id"]=1
        unique_ids+=("$id")
    done
    for id in "${unique_ids[@]}"; do
        [[ -f "${DATA_GZ[$id]}" ]] || {
            echo "error: missing gzip dataset ${DATA_GZ[$id]}" >&2
            exit 1
        }
        [[ -f "${DATA_PLAIN[$id]}" ]] || {
            echo "error: missing plain cache ${DATA_PLAIN[$id]}" >&2
            exit 1
        }
        echo "  $id: ${DATA_GZ[$id]}  (expected ${DATA_EXPECTED[$id]})"
    done
}

log_verify() {
    printf '%s\n' "$*" | tee -a "$VERIFY_LOG"
}

summarize_output() {
    local file="$1"
    if [[ ! -s "$file" ]]; then
        printf '(empty)'
        return
    fi
    local bytes text
    bytes="$(file_size_bytes "$file")"
    text="$(tr -d '\0' <"$file" | head -c 400)"
    if ((bytes > 400)); then
        printf '%s... (%s bytes)' "$text" "$bytes"
    else
        printf '%s' "$text"
    fi
}

sample_fail() {
    local label="$1" expected="$2" actual="$3"
    {
        echo "SAMPLE CONTRACT FAIL"
        echo "  case:     $label"
        echo "  expected: $expected"
        echo "  actual:   $actual"
    } | tee -a "$VERIFY_LOG" >&2
    exit 1
}

capture_command() {
    local stdout_file="$1"
    local stderr_file="$2"
    shift 2
    local status=0
    "$@" >"$stdout_file" 2>"$stderr_file" || status=$?
    printf '%s' "$status"
}

dataset_path() {
    local compression="$1" id="$2"
    if [[ "$compression" == gzip ]]; then
        printf '%s' "${DATA_GZ[$id]}"
    else
        printf '%s' "${DATA_PLAIN[$id]}"
    fi
}

resolve_dataset_paths() {
    local -n result="$1"
    local compression="$2"
    shift 2
    local id
    result=()
    for id in "$@"; do
        result+=("$(dataset_path "$compression" "$id")")
    done
}

min_u64() {
    local a="$1" b="$2"
    if ((a < b)); then
        printf '%s' "$a"
    else
        printf '%s' "$b"
    fi
}

exact_output_records() {
    local mode="$1"
    shift
    local -a ids=("$@")
    local n units
    case "$mode" in
        se_count)
            n="${DATA_EXPECTED[${ids[0]}]}"
            min_u64 "$SAMPLE_COUNT" "$n"
            ;;
        paired_count)
            n="${DATA_EXPECTED[${ids[0]}]}"
            units="$(min_u64 "$SAMPLE_COUNT" "$n")"
            printf '%s' "$((units * 2))"
            ;;
        interleaved_count)
            n="${DATA_EXPECTED[${ids[0]}]}"
            units="$(min_u64 "$SAMPLE_COUNT" "$((n / 2))")"
            printf '%s' "$((units * 2))"
            ;;
        *)
            echo "error: exact_output_records for $mode" >&2
            return 1
            ;;
    esac
}

seqtk_byte_exact_id() {
    local id="$1"
    [[ "$id" != Dense ]]
}

tool_ready() {
    local tool="$1"
    bench_has_tool "$tool" || return 1
    if [[ "$tool" == bbtools ]] && ! command -v java >/dev/null 2>&1; then
        return 1
    fi
    return 0
}

tool_supports_mode() {
    local tool="$1" mode="$2"
    case "$tool:$mode" in
        seqtk:se_fraction|seqtk:se_count) return 0 ;;
        seqkit:se_fraction|seqkit:se_count) return 0 ;;
        rasusa:se_fraction|rasusa:se_count|rasusa:paired_fraction|rasusa:paired_count) return 0 ;;
        fq:se_fraction|fq:se_count|fq:paired_fraction|fq:paired_count) return 0 ;;
        fqkit:se_count) return 0 ;;
        irma:se_fraction|irma:se_count|irma:paired_fraction|irma:paired_count) return 0 ;;
        bbtools:se_fraction|bbtools:se_count|bbtools:paired_fraction|bbtools:paired_count|bbtools:interleaved_fraction|bbtools:interleaved_count) return 0 ;;
        *) return 1 ;;
    esac
}

build_sample_args() {
    local -n result="$1"
    local tool="$2"
    local mode="$3"
    shift 3
    local -a paths=("$@")
    result=()

    case "$tool" in
        z-fastq|z-fastq-native)
            local binary="$ZFASTQ"
            [[ "$tool" == z-fastq-native ]] && binary="$ZFASTQ_NATIVE"
            result=("$binary" sample --seed "$SAMPLE_SEED")
            case "$mode" in
                se_fraction) result+=(--fraction "$SAMPLE_FRACTION" "${paths[0]}") ;;
                se_count) result+=(--count "$SAMPLE_COUNT" "${paths[0]}") ;;
                paired_fraction) result+=(--paired --fraction "$SAMPLE_FRACTION" "${paths[0]}" "${paths[1]}") ;;
                paired_count) result+=(--paired --count "$SAMPLE_COUNT" "${paths[0]}" "${paths[1]}") ;;
                interleaved_fraction) result+=(--interleaved --fraction "$SAMPLE_FRACTION" "${paths[0]}") ;;
                interleaved_count) result+=(--interleaved --count "$SAMPLE_COUNT" "${paths[0]}") ;;
                *) echo "error: unknown sample mode $mode" >&2; return 1 ;;
            esac
            ;;
        seqtk)
            case "$mode" in
                se_fraction) result=("$SEQTK" sample -s "$SAMPLE_SEED" "${paths[0]}" "$SAMPLE_FRACTION") ;;
                se_count) result=("$SEQTK" sample -2 -s "$SAMPLE_SEED" "${paths[0]}" "$SAMPLE_COUNT") ;;
                *) echo "error: seqtk does not support $mode" >&2; return 1 ;;
            esac
            ;;
        seqkit)
            case "$mode" in
                se_fraction)
                    result=("$SEQKIT" sample --quiet -j 1 -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" -o - "${paths[0]}")
                    ;;
                se_count)
                    result=("$SEQKIT" sample2 --quiet -j 1 -2 -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" -o - "${paths[0]}")
                    ;;
                *) echo "error: seqkit does not support $mode" >&2; return 1 ;;
            esac
            ;;
        rasusa)
            case "$mode" in
                se_fraction) result=("$RASUSA" reads -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" "${paths[0]}") ;;
                se_count) result=("$RASUSA" reads -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" "${paths[0]}") ;;
                paired_fraction)
                    result=("$RASUSA" reads -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" -o /dev/null -o /dev/null "${paths[0]}" "${paths[1]}")
                    ;;
                paired_count)
                    result=("$RASUSA" reads -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" -o /dev/null -o /dev/null "${paths[0]}" "${paths[1]}")
                    ;;
                *) echo "error: rasusa does not support $mode" >&2; return 1 ;;
            esac
            ;;
        fq)
            case "$mode" in
                se_fraction)
                    result=("$FQ" subsample -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" --r1-dst /dev/null "${paths[0]}")
                    ;;
                se_count)
                    result=("$FQ" subsample -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" --r1-dst /dev/null "${paths[0]}")
                    ;;
                paired_fraction)
                    result=("$FQ" subsample -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" --r1-dst /dev/null --r2-dst /dev/null "${paths[0]}" "${paths[1]}")
                    ;;
                paired_count)
                    result=("$FQ" subsample -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" --r1-dst /dev/null --r2-dst /dev/null "${paths[0]}" "${paths[1]}")
                    ;;
                *) echo "error: fq does not support $mode" >&2; return 1 ;;
            esac
            ;;
        fqkit)
            result=("$FQKIT" subfq -q -@ 1 -2 -n "$SAMPLE_COUNT" -s "$SAMPLE_SEED" "${paths[0]}")
            ;;
        irma)
            result=("$IRMA_CORE" sampler -s "$SAMPLE_SEED")
            case "$mode" in
                se_fraction|paired_fraction) result+=(-p "$SAMPLE_PERCENT") ;;
                se_count|paired_count) result+=(-t "$SAMPLE_COUNT") ;;
                *) echo "error: irma does not support $mode" >&2; return 1 ;;
            esac
            result+=("${paths[@]}")
            ;;
        bbtools)
            result=("$REFORMAT" -Xmx200m threads=1 sampleseed="$SAMPLE_SEED" qin=33 qout=33 changequality=f overwrite=t)
            case "$mode" in
                se_fraction)
                    result+=(in="${paths[0]}" out=/dev/null samplerate="$SAMPLE_FRACTION")
                    ;;
                se_count)
                    result+=(in="${paths[0]}" out=/dev/null samplereadstarget="$SAMPLE_COUNT")
                    ;;
                paired_fraction)
                    result+=(in="${paths[0]}" in2="${paths[1]}" out=/dev/null samplerate="$SAMPLE_FRACTION")
                    ;;
                paired_count)
                    result+=(in="${paths[0]}" in2="${paths[1]}" out=/dev/null samplereadstarget="$SAMPLE_COUNT")
                    ;;
                interleaved_fraction)
                    result+=(in="${paths[0]}" int=t out=/dev/null samplerate="$SAMPLE_FRACTION")
                    ;;
                interleaved_count)
                    result+=(in="${paths[0]}" int=t out=/dev/null samplereadstarget="$SAMPLE_COUNT")
                    ;;
                *) echo "error: bbtools does not support $mode" >&2; return 1 ;;
            esac
            ;;
        *)
            echo "error: unknown sample tool $tool" >&2
            return 1
            ;;
    esac
}

workload_input_bytes() {
    local compression="$1"
    shift
    local total=0 id path
    for id in "$@"; do
        path="$(dataset_path "$compression" "$id")"
        total=$((total + $(file_size_bytes "$path")))
    done
    printf '%s' "$total"
}

workload_decoded_bytes() {
    local total=0 id
    for id in "$@"; do
        total=$((total + DATA_DECODED[$id]))
    done
    printf '%s' "$total"
}

parse_single_int() {
    local text="$1"
    local line
    line="$(printf '%s' "$text" | tr -d '\r' | awk 'NF{print; exit}')"
    if [[ ! "$line" =~ ^[0-9]+$ ]]; then
        return 1
    fi
    printf '%s' "$line"
}

count_fastq_records() {
    local file="$1"
    local out="$CHECK_DIR/count.out" err="$CHECK_DIR/count.err"
    local status got
    status="$(capture_command "$out" "$err" "$ZFASTQ" count "$file")"
    got="$(parse_single_int "$(cat "$out")" || true)"
    [[ "$status" == 0 && -n "$got" ]] || return 1
    printf '%s' "$got"
}

write_screened_fastq() {
    local path="$1"
    local n="${2:-16}"
    local i
    : >"$path"
    for ((i = 1; i <= n; i++)); do
        printf '@r%d\nACGT\n+\nIIII\n' "$i" >>"$path"
    done
}

write_pair_fastq() {
    local r1="$1" r2="$2" interleaved="$3"
    local n="${4:-16}"
    local i
    : >"$r1"
    : >"$r2"
    : >"$interleaved"
    for ((i = 1; i <= n; i++)); do
        printf '@cluster%d/1\nAC\n+\n!!\n' "$i" >>"$r1"
        printf '@cluster%d/2\nGT\n+\n##\n' "$i" >>"$r2"
        printf '@cluster%d/1\nAC\n+\n!!\n' "$i" >>"$interleaved"
        printf '@cluster%d/2\nGT\n+\n##\n' "$i" >>"$interleaved"
    done
}

record_peer_fixture() {
    local fixture="$1" tool="$2" result="$3" detail="${4:-}"
    PEER_FIXTURE_ROWS+=("$(printf \
        '{"fixture":%s,"tool":%s,"result":%s,"detail":%s}' \
        "$(zebrac_json_string "$fixture")" \
        "$(zebrac_json_string "$tool")" \
        "$(zebrac_json_string "$result")" \
        "$(zebrac_json_string "$detail")")")
}

peer_writes_stdout() {
    local tool="$1" mode="$2"
    case "$tool" in
        fq|bbtools) return 1 ;;
        rasusa)
            [[ "$mode" == paired_* ]] && return 1
            return 0
            ;;
        *) return 0 ;;
    esac
}

probe_peer_fixture() {
    local fixture="$1" tool="$2" mode="$3"
    shift 3
    local -a paths=("$@")
    if ! tool_supports_mode "$tool" "$mode"; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "not supported for mode"
        log_verify "  fixture $fixture $tool unsupported"
        return 0
    fi
    if ! tool_ready "$tool"; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "tool unavailable"
        log_verify "  fixture $fixture $tool unavailable"
        return 0
    fi
    local -a args=()
    build_sample_args args "$tool" "$mode" "${paths[@]}"
    local stem="fixture_${fixture}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local out="$CHECK_DIR/${stem}.out" err="$CHECK_DIR/${stem}.err"
    local status n expected
    status="$(capture_command "$out" "$err" "${args[@]}")"
    if [[ "$status" != 0 ]]; then
        record_peer_fixture "$fixture" "$tool" "fail" "exit $status"
        log_verify "  fixture FAIL $fixture $tool (exit $status)"
        return 0
    fi
    if peer_writes_stdout "$tool" "$mode"; then
        n="$(count_fastq_records "$out" || true)"
        if [[ -z "$n" ]]; then
            record_peer_fixture "$fixture" "$tool" "fail" "unreadable stdout"
            log_verify "  fixture FAIL $fixture $tool (unreadable stdout)"
            return 0
        fi
        case "$mode" in
            se_count)
                expected="$FIXTURE_COUNT"
                if [[ "$n" != "$expected" ]]; then
                    record_peer_fixture "$fixture" "$tool" "fail" "got $n records want $expected"
                    log_verify "  fixture FAIL $fixture $tool (count $n != $expected)"
                    return 0
                fi
                ;;
            paired_count|interleaved_count)
                expected="$((FIXTURE_COUNT * 2))"
                if [[ "$n" != "$expected" ]]; then
                    record_peer_fixture "$fixture" "$tool" "fail" "got $n records want $expected"
                    log_verify "  fixture FAIL $fixture $tool (count $n != $expected)"
                    return 0
                fi
                ;;
            paired_fraction|interleaved_fraction)
                if ((n % 2 != 0)); then
                    record_peer_fixture "$fixture" "$tool" "fail" "odd record count $n"
                    log_verify "  fixture FAIL $fixture $tool (odd pair output $n)"
                    return 0
                fi
                ;;
        esac
    fi
    record_peer_fixture "$fixture" "$tool" "pass" "exit 0"
    log_verify "  fixture PASS $fixture $tool"
}

run_peer_fixture_probes() {
    local tool
    local screened="$CHECK_DIR/screened.fastq"
    local pair_r1="$CHECK_DIR/pair_r1.fastq"
    local pair_r2="$CHECK_DIR/pair_r2.fastq"
    local interleaved="$CHECK_DIR/interleaved.fastq"
    local saved_fraction="$SAMPLE_FRACTION"
    local saved_count="$SAMPLE_COUNT"
    local saved_percent="$SAMPLE_PERCENT"
    SAMPLE_FRACTION="$FIXTURE_FRACTION"
    SAMPLE_COUNT="$FIXTURE_COUNT"
    SAMPLE_PERCENT="$FIXTURE_PERCENT"
    log_verify "--- peer fixture probes (descriptive; do not fail the run) ---"
    log_verify "  fixture parameters: fraction=$SAMPLE_FRACTION count=$SAMPLE_COUNT seed=$SAMPLE_SEED"
    for tool in "${SAMPLE_PEER_ORDER[@]}"; do
        probe_peer_fixture screened_se_fraction "$tool" se_fraction "$screened"
        probe_peer_fixture screened_se_count "$tool" se_count "$screened"
        probe_peer_fixture slash_pair_fraction "$tool" paired_fraction "$pair_r1" "$pair_r2"
        probe_peer_fixture slash_pair_count "$tool" paired_count "$pair_r1" "$pair_r2"
        probe_peer_fixture interleaved_fraction "$tool" interleaved_fraction "$interleaved"
        probe_peer_fixture interleaved_count "$tool" interleaved_count "$interleaved"
    done
    SAMPLE_FRACTION="$saved_fraction"
    SAMPLE_COUNT="$saved_count"
    SAMPLE_PERCENT="$saved_percent"
}

expect_ok_binary() {
    local binary="$1" label="$2"
    shift 2
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == 0 ]] ||
        sample_fail "$label" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_ok_stdin() {
    local binary="$1" label="$2" infile="$3"
    shift 3
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status=0
    "$binary" "$@" <"$infile" >"$out" 2>"$err" || status=$?
    [[ "$status" == 0 ]] ||
        sample_fail "$label" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

real_sample_out() {
    local compression="$1" mode="$2" kind="${3:-isa-l}"
    if [[ "$kind" == native ]]; then
        printf '%s' "$CHECK_DIR/real_${compression}_${mode}_native.out"
    else
        printf '%s' "$CHECK_DIR/real_${compression}_${mode}.out"
    fi
}

expect_empty_stdout() {
    local binary="$1" label="$2"
    shift 2
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == 0 && ! -s "$out" ]] ||
        sample_fail "$label" "exit 0 and empty stdout" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_record_count() {
    local label="$1" expected="$2" file="$3"
    local got
    got="$(count_fastq_records "$file" || true)"
    [[ "$got" == "$expected" ]] ||
        sample_fail "$label" "$expected records" "count=${got:-?}"
}

cmp_or_fail() {
    local label="$1" left="$2" right="$3"
    cmp -s -- "$left" "$right" ||
        sample_fail "$label" "byte-identical files" \
            "left=$(file_size_bytes "$left") bytes right=$(file_size_bytes "$right") bytes"
}

run_zfastq_to() {
    local out="$1"
    shift
    local err="$CHECK_DIR/zfastq.err"
    local status=0
    "$@" >"$out" 2>"$err" || status=$?
    [[ "$status" == 0 ]] ||
        sample_fail "$*" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

run_contract_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== sample contract ${TIMESTAMP} ==="
    log_verify "parameters: fraction=$SAMPLE_FRACTION count=$SAMPLE_COUNT seed=$SAMPLE_SEED"

    local screened="$CHECK_DIR/screened.fastq"
    local screened_gz="$CHECK_DIR/screened.fastq.gz"
    local pair_r1="$CHECK_DIR/pair_r1.fastq"
    local pair_r2="$CHECK_DIR/pair_r2.fastq"
    local interleaved="$CHECK_DIR/interleaved.fastq"
    write_screened_fastq "$screened" 16
    gzip -c -- "$screened" >"$screened_gz"
    write_pair_fastq "$pair_r1" "$pair_r2" "$interleaved" 16

    local binary out err status
    log_verify "--- valid fixtures ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_empty_stdout "$binary" "$(basename "$binary") fraction 0" \
            sample --fraction 0 --seed "$SAMPLE_SEED" "$FIXTURE_DIR/acgtn_valid.fastq"
        expect_ok_binary "$binary" "$(basename "$binary") fraction 1" \
            sample --fraction 1 --seed "$SAMPLE_SEED" "$FIXTURE_DIR/acgtn_valid.fastq"
        expect_empty_stdout "$binary" "$(basename "$binary") count 0" \
            sample --count 0 --seed "$SAMPLE_SEED" "$FIXTURE_DIR/acgtn_valid.fastq"
        expect_ok_binary "$binary" "$(basename "$binary") count 1" \
            sample --count 1 --seed "$SAMPLE_SEED" "$FIXTURE_DIR/acgtn_valid.fastq"
        expect_ok_binary "$binary" "$(basename "$binary") gzip fraction" \
            sample --fraction "$SAMPLE_FRACTION" --seed "$SAMPLE_SEED" "$screened_gz"
        expect_ok_stdin "$binary" "$(basename "$binary") stdin fraction 1" \
            "$FIXTURE_DIR/acgtn_valid.fastq" \
            sample --fraction 1 --seed "$SAMPLE_SEED" -
    done

    log_verify "--- screened seqtk exact reference ---"
    bench_require_tool seqtk
    out="$CHECK_DIR/zf.frac.out"
    run_zfastq_to "$out" "$ZFASTQ" sample --fraction 0.5 --seed 11 "$screened"
    status="$(capture_command "$CHECK_DIR/seqtk.frac.out" "$CHECK_DIR/seqtk.frac.err" \
        "$SEQTK" sample -s 11 "$screened" 0.5)"
    [[ "$status" == 0 ]] ||
        sample_fail "seqtk fraction 0.5" "exit 0" "exit=$status"
    cmp_or_fail "z-fastq vs seqtk screened fraction 0.5 seed 11" \
        "$out" "$CHECK_DIR/seqtk.frac.out"

    out="$CHECK_DIR/zf.count.out"
    run_zfastq_to "$out" "$ZFASTQ" sample --count 3 --seed 11 "$screened"
    status="$(capture_command "$CHECK_DIR/seqtk.count.out" "$CHECK_DIR/seqtk.count.err" \
        "$SEQTK" sample -2 -s 11 "$screened" 3)"
    [[ "$status" == 0 ]] ||
        sample_fail "seqtk count 3" "exit 0" "exit=$status"
    cmp_or_fail "z-fastq vs seqtk screened count 3 seed 11" \
        "$out" "$CHECK_DIR/seqtk.count.out"
    expect_record_count "screened exact count 3" 3 "$out"

    log_verify "--- paired and interleaved select the same pairs ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        run_zfastq_to "$CHECK_DIR/paired.frac.out" "$binary" sample --paired --fraction 0.5 --seed 11 \
            "$pair_r1" "$pair_r2"
        run_zfastq_to "$CHECK_DIR/inter.frac.out" "$binary" sample --interleaved --fraction 0.5 --seed 11 \
            "$interleaved"
        cmp_or_fail "$(basename "$binary") paired vs interleaved fraction 0.5" \
            "$CHECK_DIR/paired.frac.out" "$CHECK_DIR/inter.frac.out"
        local n
        n="$(count_fastq_records "$CHECK_DIR/paired.frac.out")"
        ((n % 2 == 0)) ||
            sample_fail "$(basename "$binary") paired fraction pair completeness" "even record count" "$n"

        run_zfastq_to "$CHECK_DIR/paired.count.out" "$binary" sample --paired --count 3 --seed 11 \
            "$pair_r1" "$pair_r2"
        run_zfastq_to "$CHECK_DIR/inter.count.out" "$binary" sample --interleaved --count 3 --seed 11 \
            "$interleaved"
        cmp_or_fail "$(basename "$binary") paired vs interleaved count 3" \
            "$CHECK_DIR/paired.count.out" "$CHECK_DIR/inter.count.out"
        expect_record_count "$(basename "$binary") paired exact 3 pairs" 6 "$CHECK_DIR/paired.count.out"
    done

    log_verify "--- invocation rejects ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local reject_out="$CHECK_DIR/reject.out" reject_err="$CHECK_DIR/reject.err"
        status="$(capture_command "$reject_out" "$reject_err" "$binary" sample --count 1 -)"
        [[ "$status" == 2 ]] ||
            sample_fail "$(basename "$binary") stdin exact-count" "exit 2" "exit=$status"
        status="$(capture_command "$reject_out" "$reject_err" "$binary" sample --json "$screened")"
        [[ "$status" == 2 ]] ||
            sample_fail "$(basename "$binary") --json" "exit 2" "exit=$status"
    done

    run_peer_fixture_probes

    log_verify "--- real workloads ---"
    prepare_real_lanes
    log_verify "ALL PASSED"
}

preflight_peer() {
    local section="$1" compression="$2" mode="$3" key="$4" tool="$5"
    shift 5
    if ! tool_supports_mode "$tool" "$mode"; then
        LANE_REASON["$section|$key|$tool"]="not supported for mode"
        log_verify "  skip $section $key $tool (mode not supported)"
        return 0
    fi
    if ! tool_ready "$tool"; then
        LANE_REASON["$section|$key|$tool"]="tool unavailable"
        log_verify "  skip $section $key $tool (tool unavailable)"
        return 0
    fi
    if [[ "$tool" == irma && "$compression" == gzip ]]; then
        LANE_REASON["$section|$key|$tool"]="probed, not timed: gzip and plain select different records"
        log_verify "  peer PROBE $section $key $tool (gzip not timed)"
        local -a probe_args=() probe_paths=()
        resolve_dataset_paths probe_paths "$compression" "$@"
        build_sample_args probe_args "$tool" "$mode" "${probe_paths[@]}"
        local stem="${section}_${key}_${tool}_probe"
        stem="${stem//[^A-Za-z0-9_.+-]/_}"
        capture_command /dev/null "$CHECK_DIR/${stem}.err" "${probe_args[@]}" >/dev/null || true
        return 0
    fi
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "$@"
    build_sample_args args "$tool" "$mode" "${paths[@]}"
    local stem="${section}_${key}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local out="$CHECK_DIR/${stem}.out" err="$CHECK_DIR/${stem}.err"
    local status
    status="$(capture_command /dev/null "$err" "${args[@]}")"
    if [[ "$status" == 0 ]]; then
        LANE_ENABLED["$section|$key|$tool"]=1
        LANE_REASON["$section|$key|$tool"]="accepted positive workload"
        log_verify "  peer PASS $section $key $tool"
    else
        LANE_REASON["$section|$key|$tool"]="exit $status on positive workload"
        log_verify "  peer SKIP $section $key $tool (exit $status)"
    fi
}

qualify_zfastq_exact() {
    local compression="$1" mode="$2"
    shift 2
    local -a ids=("$@")
    local expected out err status got native_out native_err
    expected="$(exact_output_records "$mode" "${ids[@]}")"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"
    build_sample_args args z-fastq "$mode" "${paths[@]}"
    out="$(real_sample_out "$compression" "$mode")"
    err="$CHECK_DIR/zfastq_exact.err"
    status="$(capture_command "$out" "$err" "${args[@]}")"
    [[ "$status" == 0 ]] ||
        sample_fail "$compression $mode z-fastq" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
    got="$(count_fastq_records "$out" || true)"
    [[ "$got" == "$expected" ]] ||
        sample_fail "$compression $mode z-fastq exact size" "$expected records" "count=${got:-?}"

    if [[ "$compression" == gzip ]]; then
        build_sample_args args z-fastq-native "$mode" "${paths[@]}"
        native_out="$(real_sample_out "$compression" "$mode" native)"
        native_err="$CHECK_DIR/native_exact.err"
        status="$(capture_command "$native_out" "$native_err" "${args[@]}")"
        [[ "$status" == 0 ]] ||
            sample_fail "$compression $mode z-fastq-native" "exit 0" \
                "exit=$status stderr=$(summarize_output "$native_err")"
        got="$(count_fastq_records "$native_out" || true)"
        [[ "$got" == "$expected" ]] ||
            sample_fail "$compression $mode z-fastq-native exact size" "$expected records" "count=${got:-?}"
        cmp_or_fail "$compression $mode ISA-L vs native exact" "$out" "$native_out"
    fi
}

fraction_output_is_large() {
    local n=0 id
    for id in "$@"; do
        n=$((n + DATA_EXPECTED[$id]))
    done
    ((n >= 1000000))
}

qualify_zfastq_fraction() {
    local compression="$1" mode="$2"
    shift 2
    local -a ids=("$@")
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"
    build_sample_args args z-fastq "$mode" "${paths[@]}"
    local out err status n native_out
    out="$(real_sample_out "$compression" "$mode")"
    err="$CHECK_DIR/zfastq_frac.err"
    if fraction_output_is_large "${ids[@]}"; then
        status="$(capture_command /dev/null "$err" "${args[@]}")"
        [[ "$status" == 0 ]] ||
            sample_fail "$compression $mode z-fastq" "exit 0" \
                "exit=$status stderr=$(summarize_output "$err")"
        rm -f -- "$out"
        log_verify "  fraction output not retained ($mode; input records >= 1e6)"
    else
        status="$(capture_command "$out" "$err" "${args[@]}")"
        [[ "$status" == 0 ]] ||
            sample_fail "$compression $mode z-fastq" "exit 0" \
                "exit=$status stderr=$(summarize_output "$err")"
        n="$(count_fastq_records "$out" || true)"
        [[ -n "$n" ]] || sample_fail "$compression $mode z-fastq fraction count" "integer" "unreadable"
        if [[ "$mode" == paired_fraction || "$mode" == interleaved_fraction ]]; then
            ((n % 2 == 0)) ||
                sample_fail "$compression $mode pair completeness" "even record count" "$n"
        fi
    fi

    if [[ "$compression" == gzip ]]; then
        build_sample_args args z-fastq-native "$mode" "${paths[@]}"
        native_out="$(real_sample_out "$compression" "$mode" native)"
        if fraction_output_is_large "${ids[@]}"; then
            status="$(capture_command /dev/null "$CHECK_DIR/native_frac.err" "${args[@]}")"
            rm -f -- "$native_out"
        else
            status="$(capture_command "$native_out" "$CHECK_DIR/native_frac.err" "${args[@]}")"
        fi
        [[ "$status" == 0 ]] ||
            sample_fail "$compression $mode z-fastq-native" "exit 0" \
                "exit=$status stderr=$(summarize_output "$CHECK_DIR/native_frac.err")"
        if [[ -f "$out" && -f "$native_out" ]]; then
            cmp_or_fail "$compression $mode ISA-L vs native fraction" "$out" "$native_out"
        fi
    fi
}

assert_pair_layout_identity() {
    local compression="$1" kind="$2"
    local paired interleaved
    paired="$(real_sample_out "$compression" "paired_${kind}")"
    interleaved="$(real_sample_out "$compression" "interleaved_${kind}")"
    if [[ -f "$paired" && -f "$interleaved" ]]; then
        cmp_or_fail "$compression $kind paired vs interleaved real files" "$paired" "$interleaved"
        log_verify "  pair layout identical $compression $kind"
    else
        log_verify "  pair layout identity skipped ($compression $kind; output not retained)"
    fi
}

qualify_seqtk_bytes() {
    local compression="$1" mode="$2" id="$3"
    seqtk_byte_exact_id "$id" || {
        log_verify "  seqtk byte-cmp skipped on $id (named plus / unscreened)"
        return 0
    }
    if [[ "$mode" == se_fraction ]] && fraction_output_is_large "$id"; then
        log_verify "  seqtk byte-cmp skipped on $id (fraction output would be large)"
        return 0
    fi
    tool_ready seqtk || return 0
    local -a z_args=() tk_args=() paths=()
    resolve_dataset_paths paths "$compression" "$id"
    build_sample_args z_args z-fastq "$mode" "${paths[@]}"
    build_sample_args tk_args seqtk "$mode" "${paths[@]}"
    local z_out="$CHECK_DIR/seqtk_cmp_z.out" tk_out="$CHECK_DIR/seqtk_cmp_tk.out"
    local z_err="$CHECK_DIR/seqtk_cmp_z.err" tk_err="$CHECK_DIR/seqtk_cmp_tk.err"
    local status
    status="$(capture_command "$z_out" "$z_err" "${z_args[@]}")"
    [[ "$status" == 0 ]] ||
        sample_fail "$compression $mode $id z-fastq for seqtk cmp" "exit 0" "exit=$status"
    status="$(capture_command "$tk_out" "$tk_err" "${tk_args[@]}")"
    [[ "$status" == 0 ]] ||
        sample_fail "$compression $mode $id seqtk for cmp" "exit 0" "exit=$status"
    cmp_or_fail "$compression $mode $id z-fastq vs seqtk" "$z_out" "$tk_out"
    log_verify "  seqtk byte-identical $compression $mode $id"
}

prepare_workload_lanes() {
    local section="$1" compression="$2" mode="$3"
    shift 3
    local -a ids=("$@")
    local key
    key="$(workload_key "$mode" "${ids[@]}")"

    case "$mode" in
        *_count) qualify_zfastq_exact "$compression" "$mode" "${ids[@]}" ;;
        *) qualify_zfastq_fraction "$compression" "$mode" "${ids[@]}" ;;
    esac

    if [[ "$mode" == se_fraction || "$mode" == se_count ]]; then
        if (( ${#ids[@]} == 1 )); then
            qualify_seqtk_bytes "$compression" "$mode" "${ids[0]}"
        fi
    fi

    local tool
    for tool in "${SAMPLE_PEER_ORDER[@]}"; do
        preflight_peer "$section" "$compression" "$mode" "$key" "$tool" "${ids[@]}"
    done
}

prepare_real_lanes() {
    local id
    for id in "${SE_IDS[@]}"; do
        prepare_workload_lanes perf_plain plain se_fraction "$id"
        prepare_workload_lanes perf_gzip gzip se_fraction "$id"
        prepare_workload_lanes perf_plain plain se_count "$id"
        prepare_workload_lanes perf_gzip gzip se_count "$id"
    done
    prepare_workload_lanes perf_plain plain paired_fraction "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip paired_fraction "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_plain plain paired_count "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip paired_count "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_plain plain interleaved_fraction "${INTERLEAVED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip interleaved_fraction "${INTERLEAVED_IDS[@]}"
    prepare_workload_lanes perf_plain plain interleaved_count "${INTERLEAVED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip interleaved_count "${INTERLEAVED_IDS[@]}"

    assert_pair_layout_identity plain fraction
    assert_pair_layout_identity gzip fraction
    assert_pair_layout_identity plain count
    assert_pair_layout_identity gzip count
}

sample_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "sample" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$command"
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    sample_add_command "$section" "$workload" "$tool" "$family" \
        "$json_out" "$command" "$input_bytes" "$decoded_bytes"
    if ! bench_group "$json_out"; then
        echo "error: zebrac failed for $section $workload $tool" >&2
        echo "  command: $command" >&2
        return 1
    fi
    t1="$(date +%s)"
    elapsed=$((t1 - t0))
    echo "  << $section $workload $tool  ${elapsed}s"
}

run_timed_workload() {
    local section="$1" compression="$2" mode="$3"
    shift 3
    local -a ids=("$@")
    local key
    key="$(workload_key "$mode" "${ids[@]}")"
    local input_bytes decoded_bytes json command tool
    input_bytes="$(workload_input_bytes "$compression" "${ids[@]}")"
    decoded_bytes="$(workload_decoded_bytes "${ids[@]}")"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"

    build_sample_args args z-fastq "$mode" "${paths[@]}"
    command="$(zebrac_command "${args[@]}")"
    json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq.json"
    run_zebrac_tool "$section" "$key" z-fastq z-fastq "$json" \
        "$command" "$input_bytes" "$decoded_bytes"

    if [[ "$compression" == gzip ]]; then
        build_sample_args args z-fastq-native "$mode" "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq-native.json"
        run_zebrac_tool "$section" "$key" z-fastq-native z-fastq "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    fi

    for tool in "${SAMPLE_PEER_ORDER[@]}"; do
        [[ "${LANE_ENABLED[$section|$key|$tool]:-0}" == 1 ]] || continue
        build_sample_args args "$tool" "$mode" "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__${tool}.json"
        run_zebrac_tool "$section" "$key" "$tool" "$tool" "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    done
}

run_perf() {
    mkdir -p "$RESULTS_DIR/perf_plain_${TIMESTAMP}" "$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
    local id
    echo "=== perf_plain ==="
    for id in "${SE_IDS[@]}"; do
        run_timed_workload perf_plain plain se_fraction "$id"
        run_timed_workload perf_plain plain se_count "$id"
    done
    run_timed_workload perf_plain plain paired_fraction "${PAIRED_IDS[@]}"
    run_timed_workload perf_plain plain paired_count "${PAIRED_IDS[@]}"
    run_timed_workload perf_plain plain interleaved_fraction "${INTERLEAVED_IDS[@]}"
    run_timed_workload perf_plain plain interleaved_count "${INTERLEAVED_IDS[@]}"

    echo "=== perf_gzip ==="
    for id in "${SE_IDS[@]}"; do
        run_timed_workload perf_gzip gzip se_fraction "$id"
        run_timed_workload perf_gzip gzip se_count "$id"
    done
    run_timed_workload perf_gzip gzip paired_fraction "${PAIRED_IDS[@]}"
    run_timed_workload perf_gzip gzip paired_count "${PAIRED_IDS[@]}"
    run_timed_workload perf_gzip gzip interleaved_fraction "${INTERLEAVED_IDS[@]}"
    run_timed_workload perf_gzip gzip interleaved_count "${INTERLEAVED_IDS[@]}"
}

write_manifest() {
    local verify_skipped_json=false
    [[ "$VERIFY_SKIPPED" == 1 ]] && verify_skipped_json=true
    local verify_pass_json=null
    [[ -n "$VERIFY_PASS" ]] && verify_pass_json="$(zebrac_json_string "$VERIFY_PASS")"
    if [[ -f "$VERIFY_LOG" ]] && rg -Fq 'ALL PASSED' "$VERIFY_LOG"; then
        verify_skipped_json=false
        verify_pass_json='"ALL PASSED"'
    fi

    local skip_full_json=false
    $DO_FULL || skip_full_json=true
    local temporary_manifest="${MANIFEST_PATH}.tmp.$$"
    {
        printf '{\n'
        printf '  "schema_version": "sample-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "sample",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$SAMPLE_SET")"
        printf '  "sample_fraction": %s,\n' "$(zebrac_json_string "$SAMPLE_FRACTION")"
        printf '  "sample_count": %s,\n' "$(zebrac_json_number_or_null "$SAMPLE_COUNT")"
        printf '  "sample_seed": %s,\n' "$(zebrac_json_number_or_null "$SAMPLE_SEED")"
        printf '  "workloads": %s,\n' "$(workloads_json)"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${SAMPLE_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${SAMPLE_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${SAMPLE_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${SAMPLE_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${SAMPLE_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$SAMPLE_TOOLS_JSON"
        printf '  "peer_fixtures": %s,\n' "$(peer_fixtures_json)"
        printf '  "lane_reasons": {'
        local lane_sep="" lane_key
        for lane_key in "${!LANE_REASON[@]}"; do
            printf '%s%s:%s' "$lane_sep" \
                "$(zebrac_json_string "$lane_key")" \
                "$(zebrac_json_string "${LANE_REASON[$lane_key]}")"
            lane_sep=','
        done
        printf '},\n'
        printf '  "sections": {'
        local section_key section_prefix section_sep=""
        for section_key in perf_plain perf_gzip; do
            case "$section_key" in
                perf_plain) section_prefix=perf_plain_ ;;
                perf_gzip) section_prefix=perf_gzip_ ;;
            esac
            if [[ -d "$RESULTS_DIR/${section_prefix}${TIMESTAMP}" ]]; then
                printf '%s\n    %s: %s' "$section_sep" \
                    "$(zebrac_json_string "$section_key")" \
                    "$(zebrac_json_string "${section_prefix}${TIMESTAMP}")"
                section_sep=,
            fi
        done
        printf '\n  },\n'
        printf '  "skip_full": %s\n' "$skip_full_json"
        printf '}\n'
    } >"$temporary_manifest"
    mv -- "$temporary_manifest" "$MANIFEST_PATH"

    local temporary_latest="$RESULTS_DIR/LATEST.tmp.$$"
    printf '%s\n' "$TIMESTAMP" >"$temporary_latest"
    mv -- "$temporary_latest" "$RESULTS_DIR/LATEST"
}

echo "z-fastq sample bench  $TIMESTAMP"
echo

build_subjects() {
    echo "Building z-fastq native inflate, then ISA-L product binary..."
    (
        cd "$PROJECT_ROOT"
        zig build -Doptimize=ReleaseFast -Disa-l=false
    )
    cp -f -- "$PROJECT_ROOT/zig-out/bin/z-fastq" "$ZFASTQ_NATIVE"
    chmod +x "$ZFASTQ_NATIVE"
    (
        cd "$PROJECT_ROOT"
        zig build -Doptimize=ReleaseFast
    )
    bench_require_tool z-fastq
    bench_require_tool z-fastq-native
    echo "  ISA-L:  $ZFASTQ"
    echo "  native: $ZFASTQ_NATIVE"
}

build_subjects
ensure_real_data

VERIFY_PASS=""
VERIFY_SKIPPED=1
if $DO_TESTS; then
    VERIFY_SKIPPED=0
    bench_require_tool seqtk
    run_contract_tests
else
    : >"$VERIFY_LOG"
    log_verify "WARNING: skipping the sample contract preflight (--skip-tests); not for a published report"
    if $DO_BENCHMARKS; then
        prepare_real_lanes
    fi
fi

if $DO_BENCHMARKS; then
    bench_require_tool zebrac
    bench_require_tool seqtk
    if $DO_FULL; then
        run_perf
    fi
fi

SAMPLE_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
SAMPLE_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
SAMPLE_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
SAMPLE_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
SAMPLE_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
SAMPLE_TOOLS_JSON="$(
    printf '{'
    printf '"seqtk":%s,' "$(zebrac_json_string "$(bench_tool_version seqtk || true)")"
    printf '"seqkit":%s,' "$(zebrac_json_string "$(bench_tool_version seqkit || true)")"
    printf '"rasusa":%s,' "$(zebrac_json_string "$(bench_tool_version rasusa || true)")"
    printf '"fq":%s,' "$(zebrac_json_string "$(bench_tool_version fq || true)")"
    printf '"fqkit":%s,' "$(zebrac_json_string "$(bench_tool_version fqkit || true)")"
    printf '"irma":%s,' "$(zebrac_json_string "$(bench_tool_version irma || true)")"
    printf '"bbtools":%s' "$(zebrac_json_string "$(bench_tool_version bbtools || true)")"
    printf '}'
)"
write_manifest

if $DO_REPORT; then
    echo "Writing REPORT.md ..."
    allow="$ALLOW_INCOMPLETE"
    if ! $DO_BENCHMARKS || ! $DO_FULL; then
        allow=true
    fi
    "$PYTHON" "$SCRIPT_DIR/generate_report.py" --allow-incomplete "$allow"
fi

echo "done. results: $RESULTS_DIR (LATEST=$TIMESTAMP)"

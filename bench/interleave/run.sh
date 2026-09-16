#!/usr/bin/env bash
# Interleave benchmark runner: contract checks, capability probes, then zebrac.
#
# Workloads come from features.tsv:
#   interleave/<publication|small>/time
# The catalog role is two mate files. The timed job is interleaved stdout.
#
# seqtk mergepe is the screened exact layout reference (LF, nonempty, bare-plus).
# Other interleavers are descriptive. BBTools writes a file; timed runs use
# out=/dev/null so zebrac still measures the tool, not a shell.
#
# Usage:
#   bash bench/interleave/run.sh
#   bash bench/interleave/run.sh --small-real --runs 5 --warmup 3
#   bash bench/interleave/run.sh --skip-tests --skip-report
#   INTERLEAVE_RUN_TIMESTAMP=<ts> bash bench/interleave/run.sh --skip-tests

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
PLAIN_DIR="$BENCH_ROOT/shared/cache/plain"

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

ZEBRAC_MIN_SAMPLES="$RUNS"
ZEBRAC_WARMUP="$WARMUP"
export ZEBRAC_DURATION_MS ZEBRAC_MIN_SAMPLES ZEBRAC_WARMUP

if [[ -n "${INTERLEAVE_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$INTERLEAVE_RUN_TIMESTAMP"
else
    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
fi

mkdir -p "$RESULTS_DIR"
METADATA_JSONL="$RESULTS_DIR/metadata_${TIMESTAMP}.jsonl"
VERIFY_LOG="$RESULTS_DIR/verify_${TIMESTAMP}.log"
MANIFEST_PATH="$RESULTS_DIR/run_${TIMESTAMP}.json"
if $DO_BENCHMARKS && $DO_FULL; then
    : >"$METADATA_JSONL"
fi

PYTHON="$(report_python)" || exit 1
CHECK_DIR="$(mktemp -d "$RESULTS_DIR/.check.XXXXXX")"
trap 'rm -rf -- "$CHECK_DIR"' EXIT

INTERLEAVE_SET="publication"
if $SMALL_REAL; then
    INTERLEAVE_SET="small"
fi

declare -a PAIR_IDS=()
declare -a WORKLOAD_KEYS=()
declare -a WORKLOAD_ROWS=()
declare -A DATA_GZ=()
declare -A DATA_PLAIN=()
declare -A DATA_EXPECTED=()
declare -A DATA_DECODED=()
declare -A LANE_ENABLED=()
declare -A LANE_REASON=()
declare -a PEER_FIXTURE_ROWS=()

INTERLEAVE_PEER_ORDER=(seqtk seqfu fqkit irma bbtools)

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
        bench_catalog ids --suite interleave --set "$INTERLEAVE_SET" --role "$role" --no-expand
    )
    ((${#destination[@]} > 0)) || {
        echo "error: interleave/$INTERLEAVE_SET/$role has no datasets" >&2
        exit 1
    }
}

ensure_real_data() {
    local download_args=(--suite interleave)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (interleave/$INTERLEAVE_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    load_role_ids time PAIR_IDS
    ((${#PAIR_IDS[@]} == 2)) || {
        echo "error: interleave/$INTERLEAVE_SET/time must contain exactly 2 datasets" >&2
        exit 1
    }

    local id
    for id in "${PAIR_IDS[@]}"; do
        bind_dataset "$id"
    done
    register_workload time "${PAIR_IDS[@]}"

    for id in "${PAIR_IDS[@]}"; do
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

interleave_fail() {
    local label="$1" expected="$2" actual="$3"
    {
        echo "INTERLEAVE CONTRACT FAIL"
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

seqtk_byte_exact_ids() {
    local id
    for id in "$@"; do
        [[ "$id" != Dense ]] || return 1
    done
    return 0
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
        seqtk:time|seqfu:time|fqkit:time|irma:time|bbtools:time) return 0 ;;
        *) return 1 ;;
    esac
}

peer_writes_stdout() {
    local tool="$1"
    case "$tool" in
        bbtools) return 1 ;;
        *) return 0 ;;
    esac
}

build_interleave_args() {
    local -n result="$1"
    local tool="$2"
    shift 2
    local -a paths=("$@")
    result=()
    case "$tool" in
        z-fastq|z-fastq-native)
            local binary="$ZFASTQ"
            [[ "$tool" == z-fastq-native ]] && binary="$ZFASTQ_NATIVE"
            result=("$binary" interleave "${paths[0]}" "${paths[1]}")
            ;;
        seqtk)
            result=("$SEQTK" mergepe "${paths[0]}" "${paths[1]}")
            ;;
        seqfu)
            result=("$SEQFU" interleave -c -1 "${paths[0]}" -2 "${paths[1]}")
            ;;
        fqkit)
            result=("$FQKIT" merge -q -@ 1 --read1 "${paths[0]}" --read2 "${paths[1]}")
            ;;
        irma)
            result=("$IRMA_CORE" xleave "${paths[0]}" "${paths[1]}")
            ;;
        bbtools)
            result=("$REFORMAT" -Xmx200m threads=1 qin=33 qout=33 changequality=f overwrite=t
                in="${paths[0]}" in2="${paths[1]}" out="${BBTOOLS_OUT:-/dev/null}")
            ;;
        *)
            echo "error: unknown interleave tool $tool" >&2
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

write_slash_pair() {
    local r1="$1" r2="$2"
    local n="${3:-4}"
    local nl="${4:-$'\n'}"
    local i
    : >"$r1"
    : >"$r2"
    for ((i = 1; i <= n; i++)); do
        printf '@cluster%d/1%sAC%s+%s!!%s' "$i" "$nl" "$nl" "$nl" "$nl" >>"$r1"
        printf '@cluster%d/2%sGT%s+%s##%s' "$i" "$nl" "$nl" "$nl" "$nl" >>"$r2"
    done
}

write_casava_pair() {
    local r1="$1" r2="$2"
    printf '@cluster 1:N:0:index\nAC\n+\n!!\n' >"$r1"
    printf '@cluster 2:Y:0:index\nGT\n+\n##\n' >"$r2"
}

real_interleave_out() {
    local compression="$1" kind="${2:-isa-l}"
    if [[ "$kind" == native ]]; then
        printf '%s' "$CHECK_DIR/real_${compression}_native.out"
    else
        printf '%s' "$CHECK_DIR/real_${compression}.out"
    fi
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

probe_peer_fixture() {
    local fixture="$1" tool="$2"
    shift 2
    local -a paths=("$@")
    if ! tool_supports_mode "$tool" time; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "not supported for mode"
        log_verify "  fixture $fixture $tool unsupported"
        return 0
    fi
    if ! tool_ready "$tool"; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "tool unavailable"
        log_verify "  fixture $fixture $tool unavailable"
        return 0
    fi
    local stem="fixture_${fixture}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local out="$CHECK_DIR/${stem}.out" err="$CHECK_DIR/${stem}.err"
    local -a args=()
    if [[ "$tool" == bbtools ]]; then
        BBTOOLS_OUT="$out" build_interleave_args args "$tool" "${paths[@]}"
    else
        build_interleave_args args "$tool" "${paths[@]}"
    fi
    local status n
    if [[ "$tool" == bbtools ]]; then
        status="$(capture_command /dev/null "$err" "${args[@]}")"
    else
        status="$(capture_command "$out" "$err" "${args[@]}")"
    fi
    if [[ "$status" != 0 ]]; then
        local detail bytes
        bytes="$(file_size_bytes "$out")"
        n="$(count_fastq_records "$out" || true)"
        if [[ -s "$out" && -n "$n" ]]; then
            detail="exit $status wrote $n records"
        elif [[ -s "$out" ]]; then
            detail="exit $status non-FASTQ stdout ($bytes bytes)"
        else
            detail="exit $status empty stdout"
        fi
        record_peer_fixture "$fixture" "$tool" "fail" "$detail"
        log_verify "  fixture FAIL $fixture $tool ($detail) stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
        return 0
    fi
    if peer_writes_stdout "$tool" || [[ "$tool" == bbtools ]]; then
        n="$(count_fastq_records "$out" || true)"
        if [[ -z "$n" ]]; then
            record_peer_fixture "$fixture" "$tool" "fail" "unreadable output"
            log_verify "  fixture FAIL $fixture $tool (unreadable output)"
            return 0
        fi
        if ((n % 2 != 0)); then
            record_peer_fixture "$fixture" "$tool" "fail" "odd record count $n"
            log_verify "  fixture FAIL $fixture $tool (odd pair output $n)"
            return 0
        fi
    fi
    record_peer_fixture "$fixture" "$tool" "pass" "exit 0"
    log_verify "  fixture PASS $fixture $tool"
}

run_peer_fixture_probes() {
    local tool
    local slash_r1="$CHECK_DIR/slash_r1.fastq"
    local slash_r2="$CHECK_DIR/slash_r2.fastq"
    local casava_r1="$CHECK_DIR/casava_r1.fastq"
    local casava_r2="$CHECK_DIR/casava_r2.fastq"
    local empty_r1="$CHECK_DIR/empty_r1.fastq"
    local empty_r2="$CHECK_DIR/empty_r2.fastq"
    local crlf_r1="$CHECK_DIR/peer_crlf_r1.fastq"
    local crlf_r2="$CHECK_DIR/peer_crlf_r2.fastq"
    write_slash_pair "$slash_r1" "$slash_r2" 4
    write_casava_pair "$casava_r1" "$casava_r2"
    write_slash_pair "$crlf_r1" "$crlf_r2" 4 $'\r\n'
    : >"$empty_r1"
    : >"$empty_r2"
    log_verify "--- peer fixture probes (descriptive; do not fail the run) ---"
    for tool in "${INTERLEAVE_PEER_ORDER[@]}"; do
        probe_peer_fixture screened_slash "$tool" "$slash_r1" "$slash_r2"
        probe_peer_fixture casava_names "$tool" "$casava_r1" "$casava_r2"
        probe_peer_fixture crlf_slash "$tool" "$crlf_r1" "$crlf_r2"
        probe_peer_fixture empty_pair "$tool" "$empty_r1" "$empty_r2"
    done
}

expect_ok_binary() {
    local binary="$1" label="$2"
    shift 2
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == 0 ]] ||
        interleave_fail "$label" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_status() {
    local binary="$1" label="$2" expected="$3"
    shift 3
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == "$expected" ]] ||
        interleave_fail "$label" "exit $expected" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_empty_stdout() {
    local binary="$1" label="$2"
    shift 2
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == 0 && ! -s "$out" ]] ||
        interleave_fail "$label" "exit 0 and empty stdout" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_ok_stdin() {
    local binary="$1" label="$2" infile="$3"
    shift 3
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status=0
    "$binary" "$@" <"$infile" >"$out" 2>"$err" || status=$?
    [[ "$status" == 0 ]] ||
        interleave_fail "$label" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_record_count() {
    local label="$1" expected="$2" file="$3"
    local got
    got="$(count_fastq_records "$file" || true)"
    [[ "$got" == "$expected" ]] ||
        interleave_fail "$label" "$expected records" "count=${got:-?}"
}

cmp_or_fail() {
    local label="$1" left="$2" right="$3"
    cmp -s -- "$left" "$right" ||
        interleave_fail "$label" "byte-identical files" \
            "left=$(file_size_bytes "$left") bytes right=$(file_size_bytes "$right") bytes"
}

run_zfastq_to() {
    local out="$1"
    shift
    local err="$CHECK_DIR/zfastq.err"
    local status=0
    "$@" >"$out" 2>"$err" || status=$?
    [[ "$status" == 0 ]] ||
        interleave_fail "$*" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

run_contract_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== interleave contract ${TIMESTAMP} ==="

    local slash_r1="$CHECK_DIR/slash_r1.fastq"
    local slash_r2="$CHECK_DIR/slash_r2.fastq"
    local slash_gz1="$CHECK_DIR/slash_r1.fastq.gz"
    local slash_gz2="$CHECK_DIR/slash_r2.fastq.gz"
    local casava_r1="$CHECK_DIR/casava_r1.fastq"
    local casava_r2="$CHECK_DIR/casava_r2.fastq"
    local empty_r1="$CHECK_DIR/empty_r1.fastq"
    local empty_r2="$CHECK_DIR/empty_r2.fastq"
    local exact_r1="$CHECK_DIR/exact_r1.fastq"
    local exact_r2="$CHECK_DIR/exact_r2.fastq"
    write_slash_pair "$slash_r1" "$slash_r2" 4
    gzip -c -- "$slash_r1" >"$slash_gz1"
    gzip -c -- "$slash_r2" >"$slash_gz2"
    write_casava_pair "$casava_r1" "$casava_r2"
    : >"$empty_r1"
    : >"$empty_r2"
    printf '@same left\nA\n+\n!\n' >"$exact_r1"
    printf '@same right\nT\n+\n#\n' >"$exact_r2"

    local binary
    log_verify "--- valid fixtures ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_empty_stdout "$binary" "$(basename "$binary") empty pair" \
            interleave "$empty_r1" "$empty_r2"
        expect_ok_binary "$binary" "$(basename "$binary") slash pair" \
            interleave "$slash_r1" "$slash_r2"
        expect_ok_binary "$binary" "$(basename "$binary") casava pair" \
            interleave "$casava_r1" "$casava_r2"
        expect_ok_binary "$binary" "$(basename "$binary") gzip pair" \
            interleave "$slash_gz1" "$slash_gz2"
        expect_ok_binary "$binary" "$(basename "$binary") mixed gzip R2" \
            interleave "$slash_r1" "$slash_gz2"
        expect_ok_binary "$binary" "$(basename "$binary") mixed gzip R1" \
            interleave "$slash_gz1" "$slash_r2"
        expect_ok_binary "$binary" "$(basename "$binary") --pair-names exact" \
            interleave --pair-names exact "$exact_r1" "$exact_r2"
        expect_ok_stdin "$binary" "$(basename "$binary") stdin R1" \
            "$slash_r1" interleave - "$slash_r2"
        expect_ok_stdin "$binary" "$(basename "$binary") stdin R2" \
            "$slash_r2" interleave "$slash_r1" -
    done

    log_verify "--- screened seqtk exact layout ---"
    bench_require_tool seqtk
    local z_out="$CHECK_DIR/zf.slash.out"
    run_zfastq_to "$z_out" "$ZFASTQ" interleave "$slash_r1" "$slash_r2"
    local status
    status="$(capture_command "$CHECK_DIR/seqtk.slash.out" "$CHECK_DIR/seqtk.slash.err" \
        "$SEQTK" mergepe "$slash_r1" "$slash_r2")"
    [[ "$status" == 0 ]] ||
        interleave_fail "seqtk mergepe slash" "exit 0" "exit=$status"
    cmp_or_fail "z-fastq vs seqtk screened slash pair" \
        "$z_out" "$CHECK_DIR/seqtk.slash.out"
    expect_record_count "screened slash 4 pairs" 8 "$z_out"

    run_zfastq_to "$CHECK_DIR/zf.empty.out" "$ZFASTQ" interleave "$empty_r1" "$empty_r2"
    status="$(capture_command "$CHECK_DIR/seqtk.empty.out" "$CHECK_DIR/seqtk.empty.err" \
        "$SEQTK" mergepe "$empty_r1" "$empty_r2")"
    [[ "$status" == 0 ]] ||
        interleave_fail "seqtk mergepe empty pair" "exit 0" "exit=$status"
    cmp_or_fail "z-fastq vs seqtk empty pair" \
        "$CHECK_DIR/zf.empty.out" "$CHECK_DIR/seqtk.empty.out"

    status="$(capture_command "$CHECK_DIR/seqtk.gz.out" "$CHECK_DIR/seqtk.gz.err" \
        "$SEQTK" mergepe "$slash_gz1" "$slash_gz2")"
    [[ "$status" == 0 ]] ||
        interleave_fail "seqtk mergepe gzip slash" "exit 0" "exit=$status"
    run_zfastq_to "$CHECK_DIR/zf.gz.out" "$ZFASTQ" interleave "$slash_gz1" "$slash_gz2"
    cmp_or_fail "z-fastq vs seqtk gzip slash pair" \
        "$CHECK_DIR/zf.gz.out" "$CHECK_DIR/seqtk.gz.out"
    cmp_or_fail "gzip vs plain slash interleave" "$z_out" "$CHECK_DIR/zf.gz.out"
    run_zfastq_to "$CHECK_DIR/zf.native.gz.out" "$ZFASTQ_NATIVE" interleave "$slash_gz1" "$slash_gz2"
    cmp_or_fail "native gzip slash vs LF" "$z_out" "$CHECK_DIR/zf.native.gz.out"

    if tool_ready irma; then
        status="$(capture_command "$CHECK_DIR/irma.slash.out" "$CHECK_DIR/irma.slash.err" \
            "$IRMA_CORE" xleave "$slash_r1" "$slash_r2")"
        local irma_gz_status
        irma_gz_status="$(capture_command "$CHECK_DIR/irma.gz.out" "$CHECK_DIR/irma.gz.err" \
            "$IRMA_CORE" xleave "$slash_gz1" "$slash_gz2")"
        if [[ "$status" == 0 && "$irma_gz_status" == 0 ]] &&
            cmp -s -- "$CHECK_DIR/irma.slash.out" "$CHECK_DIR/irma.gz.out" &&
            cmp -s -- "$z_out" "$CHECK_DIR/irma.slash.out"; then
            log_verify "  IRMA xleave slash gzip == plain == z-fastq (descriptive)"
        else
            log_verify "  IRMA xleave slash gzip/plain/z-fastq differed (descriptive; not a z-fastq failure)"
        fi
    fi

    log_verify "--- mixed gzip, stdin, and native match LF ---"
    run_zfastq_to "$CHECK_DIR/zf.native.slash.out" "$ZFASTQ_NATIVE" interleave "$slash_r1" "$slash_r2"
    cmp_or_fail "native vs ISA-L slash pair" "$z_out" "$CHECK_DIR/zf.native.slash.out"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag
        tag="$(basename "$binary")"
        run_zfastq_to "$CHECK_DIR/${tag}.mixed_r2.out" "$binary" interleave "$slash_r1" "$slash_gz2"
        cmp_or_fail "$tag mixed gzip R2 vs LF" "$z_out" "$CHECK_DIR/${tag}.mixed_r2.out"
        run_zfastq_to "$CHECK_DIR/${tag}.mixed_r1.out" "$binary" interleave "$slash_gz1" "$slash_r2"
        cmp_or_fail "$tag mixed gzip R1 vs LF" "$z_out" "$CHECK_DIR/${tag}.mixed_r1.out"
        expect_ok_stdin "$binary" "$tag stdin R1 bytes" "$slash_r1" interleave - "$slash_r2"
        cmp_or_fail "$tag stdin R1 vs files" "$z_out" "$CHECK_DIR/contract.out"
        expect_ok_stdin "$binary" "$tag stdin R2 bytes" "$slash_r2" interleave "$slash_r1" -
        cmp_or_fail "$tag stdin R2 vs files" "$z_out" "$CHECK_DIR/contract.out"
    done

    log_verify "--- CRLF to LF ---"
    local crlf_r1="$CHECK_DIR/crlf_r1.fastq"
    local crlf_r2="$CHECK_DIR/crlf_r2.fastq"
    write_slash_pair "$crlf_r1" "$crlf_r2" 4 $'\r\n'
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        run_zfastq_to "$CHECK_DIR/$(basename "$binary").crlf.out" "$binary" interleave "$crlf_r1" "$crlf_r2"
        cmp_or_fail "$(basename "$binary") CRLF slash vs LF slash" \
            "$z_out" "$CHECK_DIR/$(basename "$binary").crlf.out"
    done
    status="$(capture_command "$CHECK_DIR/seqtk.crlf.out" "$CHECK_DIR/seqtk.crlf.err" \
        "$SEQTK" mergepe "$crlf_r1" "$crlf_r2")"
    [[ "$status" == 0 ]] ||
        interleave_fail "seqtk mergepe CRLF slash" "exit 0" "exit=$status"
    cmp_or_fail "z-fastq vs seqtk CRLF slash pair" \
        "$z_out" "$CHECK_DIR/seqtk.crlf.out"

    log_verify "--- name policy ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local exact_out="$CHECK_DIR/exact_slash.out" exact_err="$CHECK_DIR/exact_slash.err"
        status="$(capture_command "$exact_out" "$exact_err" \
            "$binary" interleave --pair-names exact "$slash_r1" "$slash_r2")"
        [[ "$status" == 1 && ! -s "$exact_out" ]] ||
            interleave_fail "$(basename "$binary") exact on slash names" \
                "exit 1 and empty stdout" \
                "exit=$status stdout=$(summarize_output "$exact_out")"
    done

    log_verify "--- invocation rejects ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_status "$binary" "$(basename "$binary") no inputs" 2 interleave
        expect_status "$binary" "$(basename "$binary") one input" 2 interleave "$slash_r1"
        expect_status "$binary" "$(basename "$binary") three inputs" 2 \
            interleave "$slash_r1" "$slash_r2" "$slash_r1"
        expect_status "$binary" "$(basename "$binary") double stdin" 2 interleave - -
        expect_status "$binary" "$(basename "$binary") --json" 2 \
            interleave --json "$slash_r1" "$slash_r2"
        expect_status "$binary" "$(basename "$binary") --paired" 2 \
            interleave --paired "$slash_r1" "$slash_r2"
        expect_status "$binary" "$(basename "$binary") --pair-names other" 2 \
            interleave --pair-names other "$slash_r1" "$slash_r2"
        expect_status "$binary" "$(basename "$binary") --alphabet dna" 2 \
            interleave --alphabet dna "$slash_r1" "$slash_r2"
    done

    log_verify "--- pair mismatch writes no records ---"
    local bad1="$CHECK_DIR/bad_r1.fastq" bad2="$CHECK_DIR/bad_r2.fastq"
    printf '@left/1\nA\n+\n!\n' >"$bad1"
    printf '@right/2\nT\n+\n#\n' >"$bad2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local out="$CHECK_DIR/mismatch.out" err="$CHECK_DIR/mismatch.err"
        status="$(capture_command "$out" "$err" "$binary" interleave "$bad1" "$bad2")"
        [[ "$status" == 1 && ! -s "$out" ]] ||
            interleave_fail "$(basename "$binary") P001 stdout" \
                "exit 1 and empty stdout" \
                "exit=$status stdout=$(summarize_output "$out")"
        grep -Fq 'P001' "$err" ||
            interleave_fail "$(basename "$binary") P001 diagnostic" \
                "stderr contains P001" "stderr=$(summarize_output "$err")"
    done

    log_verify "--- unequal counts keep complete earlier pairs ---"
    local extra1="$CHECK_DIR/p002_r1.fastq" extra2="$CHECK_DIR/p002_r2.fastq"
    local first_pair="$CHECK_DIR/p002_first.fastq"
    printf '@ok/1\nA\n+\n!\n@extra/1\nA\n+\n!\n' >"$extra1"
    printf '@ok/2\nT\n+\n#\n' >"$extra2"
    printf '@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n' >"$first_pair"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local out="$CHECK_DIR/p002.out" err="$CHECK_DIR/p002.err"
        status="$(capture_command "$out" "$err" "$binary" interleave "$extra1" "$extra2")"
        [[ "$status" == 1 ]] ||
            interleave_fail "$(basename "$binary") P002 extra R1 exit" \
                "exit 1" "exit=$status stdout=$(summarize_output "$out")"
        cmp_or_fail "$(basename "$binary") P002 kept first pair" "$first_pair" "$out"
        grep -Fq 'P002' "$err" ||
            interleave_fail "$(basename "$binary") P002 extra R1 diagnostic" \
                "stderr contains P002" "stderr=$(summarize_output "$err")"
    done
    printf '@ok/1\nA\n+\n!\n' >"$extra1"
    printf '@ok/2\nT\n+\n#\n@extra/2\nT\n+\n#\n' >"$extra2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local out="$CHECK_DIR/p002_extra_r2.out" err="$CHECK_DIR/p002_extra_r2.err"
        status="$(capture_command "$out" "$err" "$binary" interleave "$extra1" "$extra2")"
        [[ "$status" == 1 ]] ||
            interleave_fail "$(basename "$binary") P002 extra R2 exit" \
                "exit 1" "exit=$status stdout=$(summarize_output "$out")"
        cmp_or_fail "$(basename "$binary") P002 extra R2 kept first pair" "$first_pair" "$out"
        grep -Fq 'P002' "$err" ||
            interleave_fail "$(basename "$binary") P002 extra R2 diagnostic" \
                "stderr contains P002" "stderr=$(summarize_output "$err")"
    done
    : >"$extra1"
    printf '@ok/2\nT\n+\n#\n' >"$extra2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local out="$CHECK_DIR/p002_r2.out" err="$CHECK_DIR/p002_r2.err"
        status="$(capture_command "$out" "$err" "$binary" interleave "$extra1" "$extra2")"
        [[ "$status" == 1 && ! -s "$out" ]] ||
            interleave_fail "$(basename "$binary") P002 empty R1 stdout" \
                "exit 1 and empty stdout" \
                "exit=$status stdout=$(summarize_output "$out")"
        grep -Fq 'P002' "$err" ||
            interleave_fail "$(basename "$binary") P002 empty R1 diagnostic" \
                "stderr contains P002" "stderr=$(summarize_output "$err")"
    done

    run_peer_fixture_probes

    log_verify "--- real workloads ---"
    prepare_real_lanes
    log_verify "ALL PASSED"
}

preflight_peer() {
    local section="$1" compression="$2" key="$3" tool="$4"
    shift 4
    if ! tool_supports_mode "$tool" time; then
        LANE_REASON["$section|$key|$tool"]="not supported for mode"
        log_verify "  skip $section $key $tool (mode not supported)"
        return 0
    fi
    if ! tool_ready "$tool"; then
        LANE_REASON["$section|$key|$tool"]="tool unavailable"
        log_verify "  skip $section $key $tool (tool unavailable)"
        return 0
    fi
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "$@"
    build_interleave_args args "$tool" "${paths[@]}"
    local stem="${section}_${key}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local err="$CHECK_DIR/${stem}.err"
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

qualify_zfastq_real() {
    local compression="$1"
    shift
    local -a ids=("$@")
    local expected n out err status native_out
    n="${DATA_EXPECTED[${ids[0]}]}"
    [[ "${DATA_EXPECTED[${ids[1]}]}" == "$n" ]] ||
        interleave_fail "$compression pair counts" "equal mate counts" \
            "${DATA_EXPECTED[${ids[0]}]} vs ${DATA_EXPECTED[${ids[1]}]}"
    expected="$((n * 2))"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"
    build_interleave_args args z-fastq "${paths[@]}"
    out="$(real_interleave_out "$compression")"
    err="$CHECK_DIR/zfastq_real.err"
    status="$(capture_command "$out" "$err" "${args[@]}")"
    [[ "$status" == 0 ]] ||
        interleave_fail "$compression z-fastq" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
    expect_record_count "$compression z-fastq pair-complete" "$expected" "$out"
    local got
    got="$(count_fastq_records "$out")"
    ((got % 2 == 0)) ||
        interleave_fail "$compression pair completeness" "even record count" "$got"

    if [[ "$compression" == gzip ]]; then
        build_interleave_args args z-fastq-native "${paths[@]}"
        native_out="$(real_interleave_out "$compression" native)"
        status="$(capture_command "$native_out" "$CHECK_DIR/native_real.err" "${args[@]}")"
        [[ "$status" == 0 ]] ||
            interleave_fail "$compression z-fastq-native" "exit 0" \
                "exit=$status stderr=$(summarize_output "$CHECK_DIR/native_real.err")"
        expect_record_count "$compression z-fastq-native pair-complete" "$expected" "$native_out"
        cmp_or_fail "$compression ISA-L vs native" "$out" "$native_out"
    fi
}

qualify_seqtk_bytes() {
    local compression="$1"
    shift
    local -a ids=("$@")
    seqtk_byte_exact_ids "${ids[@]}" || {
        log_verify "  seqtk byte-cmp skipped (${ids[*]}; named plus / unscreened)"
        return 0
    }
    tool_ready seqtk || return 0
    local -a z_args=() tk_args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"
    build_interleave_args z_args z-fastq "${paths[@]}"
    build_interleave_args tk_args seqtk "${paths[@]}"
    local z_out="$CHECK_DIR/seqtk_cmp_z.out" tk_out="$CHECK_DIR/seqtk_cmp_tk.out"
    local status
    status="$(capture_command "$z_out" "$CHECK_DIR/seqtk_cmp_z.err" "${z_args[@]}")"
    [[ "$status" == 0 ]] ||
        interleave_fail "$compression seqtk cmp z-fastq" "exit 0" "exit=$status"
    status="$(capture_command "$tk_out" "$CHECK_DIR/seqtk_cmp_tk.err" "${tk_args[@]}")"
    [[ "$status" == 0 ]] ||
        interleave_fail "$compression seqtk cmp seqtk" "exit 0" "exit=$status"
    cmp_or_fail "$compression ${ids[*]} z-fastq vs seqtk" "$z_out" "$tk_out"
    log_verify "  seqtk byte-identical $compression ${ids[*]}"
}

prepare_workload_lanes() {
    local section="$1" compression="$2"
    shift 2
    local -a ids=("$@")
    local key
    key="$(workload_key time "${ids[@]}")"
    qualify_zfastq_real "$compression" "${ids[@]}"
    qualify_seqtk_bytes "$compression" "${ids[@]}"
    local tool
    for tool in "${INTERLEAVE_PEER_ORDER[@]}"; do
        preflight_peer "$section" "$compression" "$key" "$tool" "${ids[@]}"
    done
}

prepare_real_lanes() {
    prepare_workload_lanes perf_plain plain "${PAIR_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip "${PAIR_IDS[@]}"
    cmp_or_fail "catalog gzip vs plain interleave" \
        "$(real_interleave_out plain)" "$(real_interleave_out gzip)"
}

interleave_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "interleave" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$command"
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    interleave_add_command "$section" "$workload" "$tool" "$family" \
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
    local section="$1" compression="$2"
    shift 2
    local -a ids=("$@")
    local key
    key="$(workload_key time "${ids[@]}")"
    local input_bytes decoded_bytes json command tool
    input_bytes="$(workload_input_bytes "$compression" "${ids[@]}")"
    decoded_bytes="$(workload_decoded_bytes "${ids[@]}")"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"

    build_interleave_args args z-fastq "${paths[@]}"
    command="$(zebrac_command "${args[@]}")"
    json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq.json"
    run_zebrac_tool "$section" "$key" z-fastq z-fastq "$json" \
        "$command" "$input_bytes" "$decoded_bytes"

    if [[ "$compression" == gzip ]]; then
        build_interleave_args args z-fastq-native "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq-native.json"
        run_zebrac_tool "$section" "$key" z-fastq-native z-fastq "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    fi

    for tool in "${INTERLEAVE_PEER_ORDER[@]}"; do
        [[ "${LANE_ENABLED[$section|$key|$tool]:-0}" == 1 ]] || continue
        build_interleave_args args "$tool" "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__${tool}.json"
        run_zebrac_tool "$section" "$key" "$tool" "$tool" "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    done
}

run_perf() {
    mkdir -p "$RESULTS_DIR/perf_plain_${TIMESTAMP}" "$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
    echo "=== perf_plain ==="
    run_timed_workload perf_plain plain "${PAIR_IDS[@]}"
    echo "=== perf_gzip ==="
    run_timed_workload perf_gzip gzip "${PAIR_IDS[@]}"
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
        printf '  "schema_version": "interleave-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "interleave",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$INTERLEAVE_SET")"
        printf '  "workloads": %s,\n' "$(workloads_json)"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${INTERLEAVE_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${INTERLEAVE_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${INTERLEAVE_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${INTERLEAVE_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${INTERLEAVE_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$INTERLEAVE_TOOLS_JSON"
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

echo "z-fastq interleave bench  $TIMESTAMP"
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
    log_verify "WARNING: skipping the interleave contract preflight (--skip-tests); not for a published report"
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

INTERLEAVE_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
INTERLEAVE_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
INTERLEAVE_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
INTERLEAVE_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
INTERLEAVE_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
INTERLEAVE_TOOLS_JSON="$(
    printf '{'
    printf '"seqtk":%s,' "$(zebrac_json_string "$(bench_tool_version seqtk || true)")"
    printf '"seqfu":%s,' "$(zebrac_json_string "$(bench_tool_version seqfu || true)")"
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

#!/usr/bin/env bash
# Deinterleave benchmark runner: contract checks, capability probes, then zebrac.
#
# Workloads come from features.tsv:
#   deinterleave/<publication|small>/time
# The catalog role is one interleaved file. The timed job writes two mate files.
#
# seqtk has no deinterleave command. `seqtk seq -l 0 -1` / `seqtk seq -l 0 -2` is the
# screened positional layout reference (two processes; contract-only, not timed).
# Other splitters are descriptive. z-fastq creates outputs exclusively, so
# timed runs write to a tmpfs sink whose known output names a companion
# reaper unlinks so exclusive-create z-fastq can be sampled on the same argv.
# Open file descriptors still complete write(); this is not "leave two files
# on disk." Do not wrap the timed argv in a shell.
#
# Usage:
#   bash bench/deinterleave/run.sh
#   bash bench/deinterleave/run.sh --small-real --runs 5 --warmup 3
#   bash bench/deinterleave/run.sh --skip-tests --skip-report
#   DEINTERLEAVE_RUN_TIMESTAMP=<ts> bash bench/deinterleave/run.sh --skip-tests

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

if [[ -n "${DEINTERLEAVE_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$DEINTERLEAVE_RUN_TIMESTAMP"
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
SINK_DIR=""
REAPER_PID=""

stop_sink_reaper() {
    if [[ -n "${REAPER_PID}" ]]; then
        kill "$REAPER_PID" 2>/dev/null || true
        wait "$REAPER_PID" 2>/dev/null || true
        REAPER_PID=""
    fi
    if [[ -n "$SINK_DIR" ]]; then
        rm -rf -- "$SINK_DIR"
        SINK_DIR=""
    fi
}

cleanup() {
    stop_sink_reaper
    rm -rf -- "$CHECK_DIR"
}
trap cleanup EXIT

DEINTERLEAVE_SET="publication"
if $SMALL_REAL; then
    DEINTERLEAVE_SET="small"
fi

declare -a INPUT_IDS=()
declare -a MATE_IDS=()
declare -a WORKLOAD_KEYS=()
declare -a WORKLOAD_ROWS=()
declare -A DATA_GZ=()
declare -A DATA_PLAIN=()
declare -A DATA_EXPECTED=()
declare -A DATA_DECODED=()
declare -A LANE_ENABLED=()
declare -A LANE_REASON=()
declare -a PEER_FIXTURE_ROWS=()

DEINTERLEAVE_PEER_ORDER=(seqtk seqfu fqkit irma bbtools)

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
        bench_catalog ids --suite deinterleave --set "$DEINTERLEAVE_SET" --role "$role" --no-expand
    )
    ((${#destination[@]} > 0)) || {
        echo "error: deinterleave/$DEINTERLEAVE_SET/$role has no datasets" >&2
        exit 1
    }
}

ensure_real_data() {
    local download_args=(--suite deinterleave)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (deinterleave/$DEINTERLEAVE_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    load_role_ids time INPUT_IDS
    ((${#INPUT_IDS[@]} == 1)) || {
        echo "error: deinterleave/$DEINTERLEAVE_SET/time must contain exactly 1 dataset" >&2
        exit 1
    }

    local id="${INPUT_IDS[0]}"
    bind_dataset "$id"
    mapfile -t MATE_IDS < <(bench_catalog field "$id" derived_deps | tr ',' '\n')
    ((${#MATE_IDS[@]} == 2)) || {
        echo "error: $id must have exactly 2 derived mate dependencies" >&2
        exit 1
    }
    local mate
    for mate in "${MATE_IDS[@]}"; do
        bind_dataset "$mate"
    done
    register_workload time "$id"

    for id in "${INPUT_IDS[@]}" "${MATE_IDS[@]}"; do
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

deinterleave_fail() {
    local label="$1" expected="$2" actual="$3"
    {
        echo "DEINTERLEAVE CONTRACT FAIL"
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
        seqfu:time|fqkit:time|irma:time|bbtools:time) return 0 ;;
        *) return 1 ;;
    esac
}

# Timed peers write two files. seqtk seq -l 0 -1/-2 is contract-only (two processes).
tool_is_timed() {
    local tool="$1"
    case "$tool" in
        seqtk) return 1 ;;
        *) return 0 ;;
    esac
}

build_deinterleave_args() {
    local -n result="$1"
    local tool="$2"
    local input="$3"
    local out1="$4"
    local out2="$5"
    local extra="$6"
    result=()
    case "$tool" in
        z-fastq|z-fastq-native)
            local binary="$ZFASTQ"
            [[ "$tool" == z-fastq-native ]] && binary="$ZFASTQ_NATIVE"
            result=("$binary" deinterleave --out1 "$out1" --out2 "$out2" "$input")
            ;;
        seqfu)
            result=("$SEQFU" deinterleave -c -o "$extra" "$input")
            ;;
        fqkit)
            result=("$FQKIT" split -q -@ 1 -p demo -o "$extra" "$input")
            ;;
        irma)
            result=("$IRMA_CORE" xleave -1 "$out1" -2 "$out2" "$input")
            ;;
        bbtools)
            result=("$REFORMAT" -Xmx200m threads=1 qin=33 qout=33 changequality=f overwrite=t
                int=t in="$input" out="$out1" out2="$out2")
            ;;
        *)
            echo "error: unknown deinterleave tool $tool" >&2
            return 1
            ;;
    esac
}

peer_output_paths() {
    local -n r1_var="$1"
    local -n r2_var="$2"
    local tool="$3" out1="$4" out2="$5" extra="$6"
    case "$tool" in
        seqfu)
            r1_var="${extra}_R1.fq"
            r2_var="${extra}_R2.fq"
            ;;
        fqkit)
            r1_var="${extra}/demo_r1.fq"
            r2_var="${extra}/demo_r2.fq"
            ;;
        *)
            r1_var="$out1"
            r2_var="$out2"
            ;;
    esac
}

timed_sink_spec() {
    local tool="$1"
    local out1="$SINK_DIR/r1.fastq"
    local out2="$SINK_DIR/r2.fastq"
    local extra="$SINK_DIR"
    case "$tool" in
        seqfu) extra="$SINK_DIR/out" ;;
        fqkit) extra="$SINK_DIR" ;;
    esac
    printf '%s\t%s\t%s' "$out1" "$out2" "$extra"
}

workload_input_bytes() {
    local compression="$1" id="$2"
    file_size_bytes "$(dataset_path "$compression" "$id")"
}

workload_decoded_bytes() {
    local id="$1"
    printf '%s' "${DATA_DECODED[$id]}"
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

seqtk_split_to() {
    local input="$1" r1="$2" r2="$3" label="$4"
    local stem="${label//[^A-Za-z0-9_.+-]/_}"
    local err="$CHECK_DIR/seqtk.${stem}.err"
    local status
    status="$(capture_command "$r1" "$err" "$SEQTK" seq -l 0 -1 "$input")"
    [[ "$status" == 0 ]] ||
        deinterleave_fail "seqtk seq -1 $label" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
    status="$(capture_command "$r2" "$err" "$SEQTK" seq -l 0 -2 "$input")"
    [[ "$status" == 0 ]] ||
        deinterleave_fail "seqtk seq -2 $label" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
}

write_slash_interleaved() {
    local out="$1"
    local n="${2:-4}"
    local nl="${3:-$'\n'}"
    local i
    : >"$out"
    for ((i = 1; i <= n; i++)); do
        printf '@cluster%d/1%sAC%s+%s!!%s' "$i" "$nl" "$nl" "$nl" "$nl" >>"$out"
        printf '@cluster%d/2%sGT%s+%s##%s' "$i" "$nl" "$nl" "$nl" "$nl" >>"$out"
    done
}

write_casava_interleaved() {
    local out="$1"
    printf '@cluster 1:N:0:index\nAC\n+\n!!\n@cluster 2:Y:0:index\nGT\n+\n##\n' >"$out"
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

probe_seqtk_fixture() {
    local fixture="$1" input="$2"
    local stem="fixture_${fixture}_seqtk"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local r1="$CHECK_DIR/${stem}_r1.out" r2="$CHECK_DIR/${stem}_r2.out"
    local err="$CHECK_DIR/${stem}.err"
    local status n1 n2
    status="$(capture_command "$r1" "$err" "$SEQTK" seq -l 0 -1 "$input")"
    if [[ "$status" != 0 ]]; then
        record_peer_fixture "$fixture" seqtk "fail" "seq -1 exit $status"
        log_verify "  fixture FAIL $fixture seqtk (seq -1 exit $status) stderr=$(summarize_output "$err")"
        return 0
    fi
    status="$(capture_command "$r2" "$err" "$SEQTK" seq -l 0 -2 "$input")"
    if [[ "$status" != 0 ]]; then
        record_peer_fixture "$fixture" seqtk "fail" "seq -2 exit $status"
        log_verify "  fixture FAIL $fixture seqtk (seq -2 exit $status) stderr=$(summarize_output "$err")"
        return 0
    fi
    n1="$(count_fastq_records "$r1" || true)"
    n2="$(count_fastq_records "$r2" || true)"
    if [[ -z "$n1" || -z "$n2" ]]; then
        record_peer_fixture "$fixture" seqtk "fail" "unreadable output"
        log_verify "  fixture FAIL $fixture seqtk (unreadable output)"
        return 0
    fi
    if [[ "$n1" != "$n2" ]]; then
        record_peer_fixture "$fixture" seqtk "fail" "unequal mate counts $n1 $n2"
        log_verify "  fixture FAIL $fixture seqtk (unequal mate counts $n1 $n2)"
        return 0
    fi
    record_peer_fixture "$fixture" seqtk "pass" "exit 0"
    log_verify "  fixture PASS $fixture seqtk"
}

probe_peer_fixture() {
    local fixture="$1" tool="$2" input="$3"
    if [[ "$tool" == seqtk ]]; then
        if ! tool_ready seqtk; then
            record_peer_fixture "$fixture" seqtk "unsupported" "tool unavailable"
            log_verify "  fixture $fixture seqtk unavailable"
            return 0
        fi
        probe_seqtk_fixture "$fixture" "$input"
        return 0
    fi
    if ! tool_supports_mode "$tool" time || ! tool_is_timed "$tool"; then
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
    local outdir="$CHECK_DIR/${stem}.dir"
    mkdir -p "$outdir"
    local out1="$outdir/r1.fastq" out2="$outdir/r2.fastq" extra="$outdir"
    case "$tool" in
        seqfu) extra="$outdir/out" ;;
        fqkit) extra="$outdir" ;;
    esac
    local -a args=()
    build_deinterleave_args args "$tool" "$input" "$out1" "$out2" "$extra"
    local err="$CHECK_DIR/${stem}.err"
    local status r1 r2 n1 n2 bytes
    status="$(capture_command /dev/null "$err" "${args[@]}")"
    peer_output_paths r1 r2 "$tool" "$out1" "$out2" "$extra"
    if [[ "$status" != 0 ]]; then
        local detail
        bytes="$(file_size_bytes "$r1")"
        n1="$(count_fastq_records "$r1" || true)"
        if [[ -s "$r1" && -n "$n1" ]]; then
            detail="exit $status wrote $n1 R1 records"
        elif [[ -s "$r1" ]]; then
            detail="exit $status non-FASTQ R1 ($bytes bytes)"
        else
            detail="exit $status empty R1"
        fi
        record_peer_fixture "$fixture" "$tool" "fail" "$detail"
        log_verify "  fixture FAIL $fixture $tool ($detail) stderr=$(summarize_output "$err")"
        return 0
    fi
    n1="$(count_fastq_records "$r1" || true)"
    n2="$(count_fastq_records "$r2" || true)"
    if [[ -z "$n1" || -z "$n2" ]]; then
        record_peer_fixture "$fixture" "$tool" "fail" "unreadable output"
        log_verify "  fixture FAIL $fixture $tool (unreadable output)"
        return 0
    fi
    if [[ "$n1" != "$n2" ]]; then
        record_peer_fixture "$fixture" "$tool" "fail" "unequal mate counts $n1 $n2"
        log_verify "  fixture FAIL $fixture $tool (unequal mate counts $n1 $n2)"
        return 0
    fi
    record_peer_fixture "$fixture" "$tool" "pass" "exit 0"
    log_verify "  fixture PASS $fixture $tool"
}

run_peer_fixture_probes() {
    local tool
    local slash="$CHECK_DIR/slash.interleaved.fastq"
    local casava="$CHECK_DIR/casava.interleaved.fastq"
    local empty="$CHECK_DIR/empty.interleaved.fastq"
    local crlf="$CHECK_DIR/crlf.interleaved.fastq"
    local p001="$CHECK_DIR/fixture.p001.interleaved.fastq"
    local odd="$CHECK_DIR/fixture.odd.interleaved.fastq"
    write_slash_interleaved "$slash" 4
    write_casava_interleaved "$casava"
    write_slash_interleaved "$crlf" 4 $'\r\n'
    : >"$empty"
    printf '@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n' >"$p001"
    printf '@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@extra/1\nA\n+\n!\n' >"$odd"
    log_verify "--- peer fixture probes (descriptive; do not fail the run) ---"
    for tool in "${DEINTERLEAVE_PEER_ORDER[@]}"; do
        probe_peer_fixture screened_slash "$tool" "$slash"
        probe_peer_fixture casava_names "$tool" "$casava"
        probe_peer_fixture crlf_slash "$tool" "$crlf"
        probe_peer_fixture empty_pair "$tool" "$empty"
        probe_peer_fixture p001_mismatch "$tool" "$p001"
        probe_peer_fixture odd_record "$tool" "$odd"
    done
}

expect_status() {
    local binary="$1" label="$2" expected="$3"
    shift 3
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == "$expected" ]] ||
        deinterleave_fail "$label" "exit $expected" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

cmp_or_fail() {
    local label="$1" left="$2" right="$3"
    cmp -s -- "$left" "$right" ||
        deinterleave_fail "$label" "byte-identical files" \
            "left=$(file_size_bytes "$left") bytes right=$(file_size_bytes "$right") bytes"
}

expect_empty_files() {
    local label="$1"
    shift
    local file
    for file in "$@"; do
        [[ -f "$file" ]] ||
            deinterleave_fail "$label" "created empty files" "$file missing"
        [[ ! -s "$file" ]] ||
            deinterleave_fail "$label" "empty files" \
                "$file=$(summarize_output "$file")"
    done
}

run_deinterleave_to() {
    local binary="$1" input="$2" out1="$3" out2="$4"
    shift 4
    local err="$CHECK_DIR/zfastq.err"
    local status=0
    "$binary" deinterleave "$@" --out1 "$out1" --out2 "$out2" "$input" >"$CHECK_DIR/contract.out" 2>"$err" || status=$?
    [[ "$status" == 0 ]] ||
        deinterleave_fail "$binary deinterleave $input" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
}

start_sink_reaper() {
    if [[ -d /dev/shm && -w /dev/shm ]]; then
        SINK_DIR="/dev/shm/zfastq-deinterleave-${TIMESTAMP}"
    else
        SINK_DIR="$RESULTS_DIR/.sink.${TIMESTAMP}"
    fi
    rm -rf -- "$SINK_DIR"
    mkdir -p "$SINK_DIR"
    python3 - "$SINK_DIR" <<'PY' &
import os
import sys
import time

path = sys.argv[1]
known = {
    "r1.fastq",
    "r2.fastq",
    "out_R1.fq",
    "out_R2.fq",
    "demo_r1.fq",
    "demo_r2.fq",
}
while True:
    try:
        for name in os.listdir(path):
            if name not in known:
                continue
            target = os.path.join(path, name)
            try:
                if os.path.isfile(target) or os.path.islink(target):
                    os.unlink(target)
            except FileNotFoundError:
                pass
            except IsADirectoryError:
                pass
    except FileNotFoundError:
        break
    time.sleep(0.001)
PY
    REAPER_PID=$!
}

run_contract_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== deinterleave contract ${TIMESTAMP} ==="

    local slash="$CHECK_DIR/slash.interleaved.fastq"
    local slash_gz="$CHECK_DIR/slash.interleaved.fastq.gz"
    local casava="$CHECK_DIR/casava.interleaved.fastq"
    local empty="$CHECK_DIR/empty.interleaved.fastq"
    local exact="$CHECK_DIR/exact.interleaved.fastq"
    local slash_r1="$CHECK_DIR/slash.r1.ref.fastq"
    local slash_r2="$CHECK_DIR/slash.r2.ref.fastq"
    local casava_r1="$CHECK_DIR/casava.r1.ref.fastq"
    local casava_r2="$CHECK_DIR/casava.r2.ref.fastq"
    local exact_r1="$CHECK_DIR/exact.r1.ref.fastq"
    local exact_r2="$CHECK_DIR/exact.r2.ref.fastq"
    write_slash_interleaved "$slash" 4
    cp -- "$slash" "$CHECK_DIR/slash.interleaved.ref.fastq"
    gzip -c -- "$slash" >"$slash_gz"
    write_casava_interleaved "$casava"
    : >"$empty"
    printf '@same left\nA\n+\n!\n@same right\nT\n+\n#\n' >"$exact"
    printf '@cluster1/1\nAC\n+\n!!\n@cluster2/1\nAC\n+\n!!\n@cluster3/1\nAC\n+\n!!\n@cluster4/1\nAC\n+\n!!\n' >"$slash_r1"
    printf '@cluster1/2\nGT\n+\n##\n@cluster2/2\nGT\n+\n##\n@cluster3/2\nGT\n+\n##\n@cluster4/2\nGT\n+\n##\n' >"$slash_r2"
    printf '@cluster 1:N:0:index\nAC\n+\n!!\n' >"$casava_r1"
    printf '@cluster 2:Y:0:index\nGT\n+\n##\n' >"$casava_r2"
    printf '@same left\nA\n+\n!\n' >"$exact_r1"
    printf '@same right\nT\n+\n#\n' >"$exact_r2"

    local binary out1 out2 status
    log_verify "--- valid fixtures ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag
        tag="$(basename "$binary")"
        out1="$CHECK_DIR/${tag}.empty.r1" out2="$CHECK_DIR/${tag}.empty.r2"
        run_deinterleave_to "$binary" "$empty" "$out1" "$out2"
        expect_empty_files "$tag empty pair" "$out1" "$out2"
        out1="$CHECK_DIR/${tag}.slash.r1" out2="$CHECK_DIR/${tag}.slash.r2"
        run_deinterleave_to "$binary" "$slash" "$out1" "$out2"
        cmp_or_fail "$tag slash R1" "$slash_r1" "$out1"
        cmp_or_fail "$tag slash R2" "$slash_r2" "$out2"
        out1="$CHECK_DIR/${tag}.casava.r1" out2="$CHECK_DIR/${tag}.casava.r2"
        run_deinterleave_to "$binary" "$casava" "$out1" "$out2"
        cmp_or_fail "$tag Casava R1" "$casava_r1" "$out1"
        cmp_or_fail "$tag Casava R2" "$casava_r2" "$out2"
        out1="$CHECK_DIR/${tag}.gzip.r1" out2="$CHECK_DIR/${tag}.gzip.r2"
        run_deinterleave_to "$binary" "$slash_gz" "$out1" "$out2"
        cmp_or_fail "$tag gzip slash vs LF R1" "$slash_r1" "$out1"
        cmp_or_fail "$tag gzip slash vs LF R2" "$slash_r2" "$out2"
        out1="$CHECK_DIR/${tag}.exact.r1" out2="$CHECK_DIR/${tag}.exact.r2"
        run_deinterleave_to "$binary" "$exact" "$out1" "$out2" --pair-names exact
        cmp_or_fail "$tag exact R1" "$exact_r1" "$out1"
        cmp_or_fail "$tag exact R2" "$exact_r2" "$out2"
        out1="$CHECK_DIR/${tag}.stdin.r1" out2="$CHECK_DIR/${tag}.stdin.r2"
        local status=0
        "$binary" deinterleave --out1 "$out1" --out2 "$out2" - <"$slash" >"$CHECK_DIR/contract.out" 2>"$CHECK_DIR/contract.err" || status=$?
        [[ "$status" == 0 ]] ||
            deinterleave_fail "$tag stdin" "exit 0" "exit=$status stderr=$(summarize_output "$CHECK_DIR/contract.err")"
        cmp_or_fail "$tag stdin R1 vs files" "$slash_r1" "$out1"
        cmp_or_fail "$tag stdin R2 vs files" "$slash_r2" "$out2"
    done

    log_verify "--- screened seqtk positional layout ---"
    bench_require_tool seqtk
    local zf_slash_r1="$CHECK_DIR/zf.slash.r1" zf_slash_r2="$CHECK_DIR/zf.slash.r2"
    run_deinterleave_to "$ZFASTQ" "$slash" "$zf_slash_r1" "$zf_slash_r2"
    seqtk_split_to "$slash" "$CHECK_DIR/seqtk.slash.r1" "$CHECK_DIR/seqtk.slash.r2" "slash"
    cmp_or_fail "z-fastq vs seqtk screened slash R1" "$zf_slash_r1" "$CHECK_DIR/seqtk.slash.r1"
    cmp_or_fail "z-fastq vs seqtk screened slash R2" "$zf_slash_r2" "$CHECK_DIR/seqtk.slash.r2"

    run_deinterleave_to "$ZFASTQ" "$empty" "$CHECK_DIR/zf.empty.r1" "$CHECK_DIR/zf.empty.r2"
    seqtk_split_to "$empty" "$CHECK_DIR/seqtk.empty.r1" "$CHECK_DIR/seqtk.empty.r2" "empty"
    cmp_or_fail "z-fastq vs seqtk empty R1" "$CHECK_DIR/zf.empty.r1" "$CHECK_DIR/seqtk.empty.r1"
    cmp_or_fail "z-fastq vs seqtk empty R2" "$CHECK_DIR/zf.empty.r2" "$CHECK_DIR/seqtk.empty.r2"

    seqtk_split_to "$slash_gz" "$CHECK_DIR/seqtk.gz.r1" "$CHECK_DIR/seqtk.gz.r2" "gzip slash"
    cmp_or_fail "z-fastq vs seqtk gzip slash R1" "$zf_slash_r1" "$CHECK_DIR/seqtk.gz.r1"
    cmp_or_fail "z-fastq vs seqtk gzip slash R2" "$zf_slash_r2" "$CHECK_DIR/seqtk.gz.r2"

    if tool_ready irma; then
        status="$(capture_command /dev/null "$CHECK_DIR/irma.slash.err" \
            "$IRMA_CORE" xleave -1 "$CHECK_DIR/irma.slash.r1" -2 "$CHECK_DIR/irma.slash.r2" "$slash")"
        local irma_gz_status
        irma_gz_status="$(capture_command /dev/null "$CHECK_DIR/irma.gz.err" \
            "$IRMA_CORE" xleave -1 "$CHECK_DIR/irma.gz.r1" -2 "$CHECK_DIR/irma.gz.r2" "$slash_gz")"
        if [[ "$status" == 0 && "$irma_gz_status" == 0 ]] &&
            cmp -s -- "$CHECK_DIR/irma.slash.r1" "$CHECK_DIR/irma.gz.r1" &&
            cmp -s -- "$CHECK_DIR/irma.slash.r2" "$CHECK_DIR/irma.gz.r2" &&
            cmp -s -- "$zf_slash_r1" "$CHECK_DIR/irma.slash.r1" &&
            cmp -s -- "$zf_slash_r2" "$CHECK_DIR/irma.slash.r2"; then
            log_verify "  IRMA xleave slash gzip == plain == z-fastq (descriptive)"
        else
            log_verify "  IRMA xleave slash gzip/plain/z-fastq differed (descriptive; not a z-fastq failure)"
        fi
    fi

    if tool_ready bbtools; then
        local bb1="$CHECK_DIR/bb.slash.r1" bb2="$CHECK_DIR/bb.slash.r2"
        status="$(capture_command /dev/null "$CHECK_DIR/bb.slash.err" \
            "$REFORMAT" -Xmx50m threads=1 qin=33 qout=33 changequality=f overwrite=t \
            int=t in="$slash" out="$bb1" out2="$bb2")"
        if [[ "$status" == 0 ]] && grep -Fq '##' "$bb1" && ! grep -Fq '!!' "$bb1"; then
            log_verify "  BBTools changequality=f rewrote slash R1 !! to ## (descriptive)"
        else
            log_verify "  BBTools slash quality rewrite not observed (descriptive; exit=$status)"
        fi
    fi

    if tool_ready seqfu; then
        local sf_extra="$CHECK_DIR/seqfu_slash_layout"
        status="$(capture_command /dev/null "$CHECK_DIR/seqfu.slash.err" \
            "$SEQFU" deinterleave -c -o "$sf_extra" "$slash")"
        if [[ "$status" == 0 ]] && grep -Fq '@cluster1/1 ' "${sf_extra}_R1.fq"; then
            log_verify "  SeqFu deinterleave -c rewrote slash headers with a trailing space (descriptive)"
        else
            log_verify "  SeqFu slash header rewrite not observed (descriptive; exit=$status)"
        fi
    fi

    log_verify "--- CRLF to LF ---"
    local crlf="$CHECK_DIR/crlf.interleaved.fastq"
    write_slash_interleaved "$crlf" 4 $'\r\n'
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag
        tag="$(basename "$binary")"
        run_deinterleave_to "$binary" "$crlf" "$CHECK_DIR/${tag}.crlf.r1" "$CHECK_DIR/${tag}.crlf.r2"
        cmp_or_fail "$tag CRLF slash vs LF R1" "$slash_r1" "$CHECK_DIR/${tag}.crlf.r1"
        cmp_or_fail "$tag CRLF slash vs LF R2" "$slash_r2" "$CHECK_DIR/${tag}.crlf.r2"
    done
    seqtk_split_to "$crlf" "$CHECK_DIR/seqtk.crlf.r1" "$CHECK_DIR/seqtk.crlf.r2" "CRLF"
    cmp_or_fail "z-fastq vs seqtk CRLF R1" "$slash_r1" "$CHECK_DIR/seqtk.crlf.r1"
    cmp_or_fail "z-fastq vs seqtk CRLF R2" "$slash_r2" "$CHECK_DIR/seqtk.crlf.r2"

    log_verify "--- name policy ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag policy_r1 policy_r2 policy_err
        tag="$(basename "$binary")"
        policy_r1="$CHECK_DIR/${tag}.exact_slash.r1"
        policy_r2="$CHECK_DIR/${tag}.exact_slash.r2"
        policy_err="$CHECK_DIR/${tag}.exact_slash.err"
        status="$(capture_command /dev/null "$policy_err" \
            "$binary" deinterleave --pair-names exact --out1 "$policy_r1" --out2 "$policy_r2" "$slash")"
        [[ "$status" == 1 ]] ||
            deinterleave_fail "$tag exact on slash names" "exit 1" "exit=$status"
        expect_empty_files "$tag exact on slash names" "$policy_r1" "$policy_r2"
    done

    log_verify "--- invocation rejects ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag
        tag="$(basename "$binary")"
        expect_status "$binary" "$tag missing --out1" 2 deinterleave "$slash"
        expect_status "$binary" "$tag missing --out2" 2 deinterleave --out1 "$CHECK_DIR/x.r1" "$slash"
        expect_status "$binary" "$tag no input" 2 deinterleave --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2"
        expect_status "$binary" "$tag two inputs" 2 \
            deinterleave --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2" "$slash" "$slash"
        expect_status "$binary" "$tag stdout --out1" 2 \
            deinterleave --out1 - --out2 "$CHECK_DIR/x.r2" "$slash"
        expect_status "$binary" "$tag stdout --out2" 2 \
            deinterleave --out1 "$CHECK_DIR/x.r1" --out2 - "$slash"
        expect_status "$binary" "$tag identical outputs" 2 \
            deinterleave --out1 "$CHECK_DIR/same.fq" --out2 "$CHECK_DIR/same.fq" "$slash"
        expect_status "$binary" "$tag identical /dev/null" 2 \
            deinterleave --out1 /dev/null --out2 /dev/null "$slash"
        expect_status "$binary" "$tag --json" 2 \
            deinterleave --json --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2" "$slash"
        expect_status "$binary" "$tag --paired" 2 \
            deinterleave --paired --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2" "$slash"
        expect_status "$binary" "$tag --pair-names other" 2 \
            deinterleave --pair-names other --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2" "$slash"
        expect_status "$binary" "$tag --alphabet dna" 2 \
            deinterleave --alphabet dna --out1 "$CHECK_DIR/x.r1" --out2 "$CHECK_DIR/x.r2" "$slash"
    done

    log_verify "--- existing output is refused ---"
    printf 'keep\n' >"$CHECK_DIR/exists.r1"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag err
        tag="$(basename "$binary")"
        err="$CHECK_DIR/${tag}.exists.err"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$CHECK_DIR/exists.r1" --out2 "$CHECK_DIR/${tag}.exists.r2" "$slash")"
        [[ "$status" == 3 ]] ||
            deinterleave_fail "$tag existing --out1" "exit 3" "exit=$status"
        grep -Fq 'already exists' "$err" ||
            deinterleave_fail "$tag existing --out1 diagnostic" \
                "stderr contains already exists" "stderr=$(summarize_output "$err")"
        [[ -f "$CHECK_DIR/${tag}.exists.r2" ]] &&
            deinterleave_fail "$tag existing --out1" "--out2 not created" "created"
        [[ "$(cat "$CHECK_DIR/exists.r1")" == "keep" ]] ||
            deinterleave_fail "$tag existing --out1" "preserved contents" "$(cat "$CHECK_DIR/exists.r1")"
    done

    printf 'keep-two\n' >"$CHECK_DIR/exists.r2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag err
        tag="$(basename "$binary")"
        err="$CHECK_DIR/${tag}.exists2.err"
        rm -f -- "$CHECK_DIR/${tag}.exists2.r1"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$CHECK_DIR/${tag}.exists2.r1" --out2 "$CHECK_DIR/exists.r2" "$slash")"
        [[ "$status" == 3 ]] ||
            deinterleave_fail "$tag existing --out2" "exit 3" "exit=$status"
        grep -Fq 'already exists' "$err" ||
            deinterleave_fail "$tag existing --out2 diagnostic" \
                "stderr contains already exists" "stderr=$(summarize_output "$err")"
        expect_empty_files "$tag existing --out2 created --out1" "$CHECK_DIR/${tag}.exists2.r1"
        [[ "$(cat "$CHECK_DIR/exists.r2")" == "keep-two" ]] ||
            deinterleave_fail "$tag existing --out2" "preserved contents" "$(cat "$CHECK_DIR/exists.r2")"
    done

    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag err
        tag="$(basename "$binary")"
        err="$CHECK_DIR/${tag}.input_as_out.err"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$slash" --out2 "$CHECK_DIR/${tag}.input_as_out.r2" "$slash")"
        [[ "$status" == 3 ]] ||
            deinterleave_fail "$tag input as --out1" "exit 3" "exit=$status"
        grep -Fq 'already exists' "$err" ||
            deinterleave_fail "$tag input as --out1 diagnostic" \
                "stderr contains already exists" "stderr=$(summarize_output "$err")"
        [[ -f "$CHECK_DIR/${tag}.input_as_out.r2" ]] &&
            deinterleave_fail "$tag input as --out1" "--out2 not created" "created"
        cmp_or_fail "$tag input as --out1 preserved input" \
            "$CHECK_DIR/slash.interleaved.ref.fastq" "$slash"
    done

    log_verify "--- pair mismatch writes no records ---"
    local bad="$CHECK_DIR/p001.interleaved.fastq"
    printf '@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n' >"$bad"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag out1 out2 err
        tag="$(basename "$binary")"
        out1="$CHECK_DIR/${tag}.p001.r1" out2="$CHECK_DIR/${tag}.p001.r2"
        err="$CHECK_DIR/${tag}.p001.err"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$out1" --out2 "$out2" "$bad")"
        [[ "$status" == 1 ]] ||
            deinterleave_fail "$tag P001 exit" "exit 1" "exit=$status"
        expect_empty_files "$tag P001 outputs" "$out1" "$out2"
        grep -Fq 'P001' "$err" ||
            deinterleave_fail "$tag P001 diagnostic" \
                "stderr contains P001" "stderr=$(summarize_output "$err")"
    done

    log_verify "--- odd count preserves earlier complete pairs ---"
    local odd="$CHECK_DIR/p002.interleaved.fastq"
    printf '@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@extra/1\nA\n+\n!\n' >"$odd"
    printf '@ok/1\nA\n+\n!\n' >"$CHECK_DIR/p002.expected.r1"
    printf '@ok/2\nT\n+\n#\n' >"$CHECK_DIR/p002.expected.r2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag out1 out2 err
        tag="$(basename "$binary")"
        out1="$CHECK_DIR/${tag}.p002.r1" out2="$CHECK_DIR/${tag}.p002.r2"
        err="$CHECK_DIR/${tag}.p002.err"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$out1" --out2 "$out2" "$odd")"
        [[ "$status" == 1 ]] ||
            deinterleave_fail "$tag P002 extra R1 exit" "exit 1" "exit=$status"
        cmp_or_fail "$tag P002 earlier R1" "$CHECK_DIR/p002.expected.r1" "$out1"
        cmp_or_fail "$tag P002 earlier R2" "$CHECK_DIR/p002.expected.r2" "$out2"
        grep -Fq 'P002' "$err" ||
            deinterleave_fail "$tag P002 extra R1 diagnostic" \
                "stderr contains P002" "stderr=$(summarize_output "$err")"
    done
    printf '@ok/1\nA\n+\n!\n' >"$odd"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        local tag out1 out2 err
        tag="$(basename "$binary")"
        out1="$CHECK_DIR/${tag}.p002odd.r1" out2="$CHECK_DIR/${tag}.p002odd.r2"
        err="$CHECK_DIR/${tag}.p002odd.err"
        status="$(capture_command /dev/null "$err" \
            "$binary" deinterleave --out1 "$out1" --out2 "$out2" "$odd")"
        [[ "$status" == 1 ]] ||
            deinterleave_fail "$tag P002 lone R1 exit" "exit 1" "exit=$status"
        expect_empty_files "$tag P002 lone R1 outputs" "$out1" "$out2"
        grep -Fq 'P002' "$err" ||
            deinterleave_fail "$tag P002 lone R1 diagnostic" \
                "stderr contains P002" "stderr=$(summarize_output "$err")"
    done

    run_peer_fixture_probes

    log_verify "--- real workloads ---"
    prepare_real_lanes
    log_verify "ALL PASSED"
}

preflight_peer() {
    local section="$1" compression="$2" key="$3" tool="$4" id="$5"
    if ! tool_is_timed "$tool"; then
        LANE_REASON["$section|$key|$tool"]="not timed (two processes)"
        return 0
    fi
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
    local input outdir out1 out2 extra
    input="$(dataset_path "$compression" "$id")"
    outdir="$CHECK_DIR/preflight_${section}_${tool}"
    mkdir -p "$outdir"
    out1="$outdir/r1.fastq" out2="$outdir/r2.fastq" extra="$outdir"
    case "$tool" in
        seqfu) extra="$outdir/out" ;;
        fqkit) extra="$outdir" ;;
    esac
    local -a args=()
    build_deinterleave_args args "$tool" "$input" "$out1" "$out2" "$extra"
    local err="$CHECK_DIR/preflight_${section}_${tool}.err"
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
    local compression="$1" id="$2"
    local expected_pairs n1 n2 out1 out2 err status native1 native2
    expected_pairs="$((${DATA_EXPECTED[$id]} / 2))"
    [[ "$((expected_pairs * 2))" == "${DATA_EXPECTED[$id]}" ]] ||
        deinterleave_fail "$compression $id record count" "even interleaved count" \
            "${DATA_EXPECTED[$id]}"
    out1="$(real_mate_out "$compression" r1)"
    out2="$(real_mate_out "$compression" r2)"
    err="$CHECK_DIR/zfastq_real.err"
    rm -f -- "$out1" "$out2"
    status="$(capture_command /dev/null "$err" \
        "$ZFASTQ" deinterleave --out1 "$out1" --out2 "$out2" "$(dataset_path "$compression" "$id")")"
    [[ "$status" == 0 ]] ||
        deinterleave_fail "$compression z-fastq" "exit 0" \
            "exit=$status stderr=$(summarize_output "$err")"
    n1="$(count_fastq_records "$out1" || true)"
    n2="$(count_fastq_records "$out2" || true)"
    [[ "$n1" == "$expected_pairs" && "$n2" == "$expected_pairs" ]] ||
        deinterleave_fail "$compression z-fastq pair-complete" \
            "$expected_pairs + $expected_pairs" "r1=$n1 r2=$n2"
    cmp_or_fail "$compression R1 vs catalog mate" "${DATA_PLAIN[${MATE_IDS[0]}]}" "$out1"
    cmp_or_fail "$compression R2 vs catalog mate" "${DATA_PLAIN[${MATE_IDS[1]}]}" "$out2"

    if [[ "$compression" == gzip ]]; then
        native1="$(real_mate_out "$compression" r1 native)"
        native2="$(real_mate_out "$compression" r2 native)"
        rm -f -- "$native1" "$native2"
        status="$(capture_command /dev/null "$CHECK_DIR/native_real.err" \
            "$ZFASTQ_NATIVE" deinterleave --out1 "$native1" --out2 "$native2" "$(dataset_path gzip "$id")")"
        [[ "$status" == 0 ]] ||
            deinterleave_fail "$compression z-fastq-native" "exit 0" \
                "exit=$status stderr=$(summarize_output "$CHECK_DIR/native_real.err")"
        n1="$(count_fastq_records "$native1" || true)"
        n2="$(count_fastq_records "$native2" || true)"
        [[ "$n1" == "$expected_pairs" && "$n2" == "$expected_pairs" ]] ||
            deinterleave_fail "$compression z-fastq-native pair-complete" \
                "$expected_pairs + $expected_pairs" "r1=$n1 r2=$n2"
        cmp_or_fail "$compression ISA-L vs native R1" "$out1" "$native1"
        cmp_or_fail "$compression ISA-L vs native R2" "$out2" "$native2"
    fi
}

real_mate_out() {
    local compression="$1" mate="$2" kind="${3:-isa-l}"
    if [[ "$kind" == native ]]; then
        printf '%s' "$CHECK_DIR/real_${compression}_native_${mate}.out"
    else
        printf '%s' "$CHECK_DIR/real_${compression}_${mate}.out"
    fi
}

qualify_seqtk_bytes() {
    local compression="$1" id="$2"
    seqtk_byte_exact_ids "$id" "${MATE_IDS[@]}" || {
        log_verify "  seqtk byte-cmp skipped ($id; named plus / unscreened)"
        return 0
    }
    tool_ready seqtk || return 0
    local input tk1 tk2
    input="$(dataset_path "$compression" "$id")"
    tk1="$CHECK_DIR/seqtk_cmp_${compression}_r1.out"
    tk2="$CHECK_DIR/seqtk_cmp_${compression}_r2.out"
    seqtk_split_to "$input" "$tk1" "$tk2" "$compression $id"
    cmp_or_fail "$compression $id z-fastq vs seqtk R1" "$(real_mate_out "$compression" r1)" "$tk1"
    cmp_or_fail "$compression $id z-fastq vs seqtk R2" "$(real_mate_out "$compression" r2)" "$tk2"
    log_verify "  seqtk byte-identical $compression $id"
}

prepare_workload_lanes() {
    local section="$1" compression="$2" id="$3"
    local key
    key="$(workload_key time "$id")"
    qualify_zfastq_real "$compression" "$id"
    qualify_seqtk_bytes "$compression" "$id"
    local tool
    for tool in "${DEINTERLEAVE_PEER_ORDER[@]}"; do
        preflight_peer "$section" "$compression" "$key" "$tool" "$id"
    done
}

prepare_real_lanes() {
    local id="${INPUT_IDS[0]}"
    prepare_workload_lanes perf_plain plain "$id"
    prepare_workload_lanes perf_gzip gzip "$id"
    cmp_or_fail "catalog gzip vs plain R1" \
        "$(real_mate_out plain r1)" "$(real_mate_out gzip r1)"
    cmp_or_fail "catalog gzip vs plain R2" \
        "$(real_mate_out plain r2)" "$(real_mate_out gzip r2)"
}

deinterleave_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "deinterleave" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$command"
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    deinterleave_add_command "$section" "$workload" "$tool" "$family" \
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
    local section="$1" compression="$2" id="$3"
    local key
    key="$(workload_key time "$id")"
    local input_bytes decoded_bytes json command tool
    input_bytes="$(workload_input_bytes "$compression" "$id")"
    decoded_bytes="$(workload_decoded_bytes "$id")"
    local input out1 out2 extra spec
    input="$(dataset_path "$compression" "$id")"

    spec="$(timed_sink_spec z-fastq)"
    IFS=$'\t' read -r out1 out2 extra <<<"$spec"
    local -a args=()
    build_deinterleave_args args z-fastq "$input" "$out1" "$out2" "$extra"
    command="$(zebrac_command "${args[@]}")"
    json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq.json"
    run_zebrac_tool "$section" "$key" z-fastq z-fastq "$json" \
        "$command" "$input_bytes" "$decoded_bytes"

    if [[ "$compression" == gzip ]]; then
        spec="$(timed_sink_spec z-fastq-native)"
        IFS=$'\t' read -r out1 out2 extra <<<"$spec"
        build_deinterleave_args args z-fastq-native "$input" "$out1" "$out2" "$extra"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq-native.json"
        run_zebrac_tool "$section" "$key" z-fastq-native z-fastq "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    fi

    for tool in "${DEINTERLEAVE_PEER_ORDER[@]}"; do
        [[ "${LANE_ENABLED[$section|$key|$tool]:-0}" == 1 ]] || continue
        spec="$(timed_sink_spec "$tool")"
        IFS=$'\t' read -r out1 out2 extra <<<"$spec"
        build_deinterleave_args args "$tool" "$input" "$out1" "$out2" "$extra"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__${tool}.json"
        run_zebrac_tool "$section" "$key" "$tool" "$tool" "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    done
}

run_perf() {
    mkdir -p "$RESULTS_DIR/perf_plain_${TIMESTAMP}" "$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
    start_sink_reaper
    echo "=== perf_plain (sink=$SINK_DIR) ==="
    run_timed_workload perf_plain plain "${INPUT_IDS[0]}"
    echo "=== perf_gzip ==="
    run_timed_workload perf_gzip gzip "${INPUT_IDS[0]}"
    stop_sink_reaper
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
        printf '  "schema_version": "deinterleave-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "deinterleave",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$DEINTERLEAVE_SET")"
        printf '  "workloads": %s,\n' "$(workloads_json)"
        printf '  "mate_ids": %s,\n' "$(json_string_array "${MATE_IDS[@]}")"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${DEINTERLEAVE_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${DEINTERLEAVE_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${DEINTERLEAVE_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${DEINTERLEAVE_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${DEINTERLEAVE_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$DEINTERLEAVE_TOOLS_JSON"
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

echo "z-fastq deinterleave bench  $TIMESTAMP"
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
    log_verify "WARNING: skipping the deinterleave contract preflight (--skip-tests); not for a published report"
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

DEINTERLEAVE_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
DEINTERLEAVE_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
DEINTERLEAVE_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
DEINTERLEAVE_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
DEINTERLEAVE_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
DEINTERLEAVE_TOOLS_JSON="$(
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

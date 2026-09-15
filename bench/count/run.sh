#!/usr/bin/env bash
# Count benchmark runner: same-job count check, then zebrac.
#
# Usage:
#   bash bench/count/run.sh [options]
#
# Defaults: count-agreement check first, then perf (runs=25, warmup=5, duration=5000).
# A count mismatch or unexpected exit stops the run before any timed command.
#
#   bash bench/count/run.sh
#   bash bench/count/run.sh --skip-tests --skip-report
#   bash bench/count/run.sh --skip-scale
#   bash bench/count/run.sh --small-real --skip-scale
#   COUNT_RUN_TIMESTAMP=<ts> bash bench/count/run.sh --skip-tests
#
#   -h|--help  print this header

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
SCALING_DIR="$BENCH_ROOT/shared/cache/scaling"
DATA_DIR="$BENCH_ROOT/shared/data"
PLAIN_DIR="$BENCH_ROOT/shared/cache/plain"
ORACLE="$SCRIPT_DIR/oracle.py"
FIXTURE_DIR="$PROJECT_ROOT/tests/data/synthetic"

source "$BENCH_ROOT/shared/tools.sh"

SIZE_MBS=(1 5 10 25 50 100 250 500)
READS_FIXED_COUNTS=(100000 250000 500000 1000000)

RUNS=25
WARMUP=5
ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
DO_TESTS=true
DO_BENCHMARKS=true
DO_FULL=true
DO_SCALE=true
DO_SCALE_SIZE=true
DO_SCALE_READS=true
DO_SCALE_GZIP=true
DO_REPORT=true
SMALL_REAL=false
REGENERATE_FIXTURES=false
ALLOW_INCOMPLETE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --runs) RUNS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --duration) ZEBRAC_DURATION_MS="$2"; shift 2 ;;
        --skip-tests|--skip-verify) DO_TESTS=false; shift ;;
        --skip-benchmarks|--skip-perf) DO_BENCHMARKS=false; shift ;;
        --skip-full) DO_FULL=false; shift ;;
        --skip-scale) DO_SCALE=false; shift ;;
        --skip-size) DO_SCALE_SIZE=false; shift ;;
        --skip-seqs|--skip-reads) DO_SCALE_READS=false; shift ;;
        --skip-gzip-scale) DO_SCALE_GZIP=false; shift ;;
        --skip-report) DO_REPORT=false; shift ;;
        --small-real) SMALL_REAL=true; shift ;;
        --regenerate-fixtures) REGENERATE_FIXTURES=true; shift ;;
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

if [[ -n "${COUNT_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$COUNT_RUN_TIMESTAMP"
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

DENSE_ID="Dense"
DENSE_GZ_NAME="REAL_Dense.fastq.gz"
if $SMALL_REAL; then
    DENSE_ID="DenseSmall"
    DENSE_GZ_NAME="REAL_Dense_small.fastq.gz"
fi

REAL_ORDER=(Dense Variable Long)
declare -A REAL_GZ=()
declare -A REAL_PLAIN=()
declare -A REAL_EXPECTED=()
declare -A REAL_DECODED=()

count_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "count" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$script"
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    count_add_command "$section" "$workload" "$tool" "$family" \
        "$json_out" "$script" "$input_bytes" "$decoded_bytes"
    if ! bench_group "$json_out"; then
        echo "error: zebrac failed for $section $workload $tool" >&2
        echo "  command: $script" >&2
        return 1
    fi
    t1="$(date +%s)"
    elapsed=$((t1 - t0))
    echo "  << $section $workload $tool  ${elapsed}s"
}

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

ensure_real_data() {
    local download_args=()
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    REAL_GZ["Dense"]="$DATA_DIR/$DENSE_GZ_NAME"
    REAL_PLAIN["Dense"]="$PLAIN_DIR/${DENSE_GZ_NAME%.gz}"
    REAL_GZ["Variable"]="$DATA_DIR/REAL_Variable.fastq.gz"
    REAL_PLAIN["Variable"]="$PLAIN_DIR/REAL_Variable.fastq"
    REAL_GZ["Long"]="$DATA_DIR/REAL_Long.fastq.gz"
    REAL_PLAIN["Long"]="$PLAIN_DIR/REAL_Long.fastq"

    local id
    for id in Dense Variable Long; do
        [[ -f "${REAL_GZ[$id]}" ]] || {
            echo "error: missing gzip dataset ${REAL_GZ[$id]}" >&2
            exit 1
        }
        [[ -f "${REAL_PLAIN[$id]}" ]] || {
            echo "error: missing plain cache ${REAL_PLAIN[$id]}" >&2
            exit 1
        }
        REAL_DECODED["$id"]="$(file_size_bytes "${REAL_PLAIN[$id]}")"
    done

    REAL_EXPECTED["Dense"]="$(awk -F'\t' -v id="$DENSE_ID" '$1==id{print $3; exit}' "$BENCH_ROOT/shared/datasets.manifest")"
    REAL_EXPECTED["Variable"]="$(awk -F'\t' '$1=="Variable"{print $3; exit}' "$BENCH_ROOT/shared/datasets.manifest")"
    REAL_EXPECTED["Long"]="$(awk -F'\t' '$1=="Long"{print $3; exit}' "$BENCH_ROOT/shared/datasets.manifest")"
    echo "  Dense gzip: ${REAL_GZ[Dense]}  (manifest id $DENSE_ID, expected ${REAL_EXPECTED[Dense]})"
}

ensure_scaling_fixtures() {
    if ! $DO_SCALE; then
        return 0
    fi
    local args=()
    if $REGENERATE_FIXTURES; then
        args+=(--force)
    fi
    echo "  Ensuring scaling FASTQ under $SCALING_DIR ..."
    if $DO_SCALE_SIZE; then
        bench_ensure_scaling --mode size "${args[@]}"
    fi
    if $DO_SCALE_READS || $DO_SCALE_GZIP; then
        bench_ensure_scaling --mode reads "${args[@]}"
    fi
    if $DO_SCALE_GZIP; then
        bench_ensure_scaling --mode gzip "${args[@]}"
    fi
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
    if (( bytes > 400 )); then
        printf '%s... (%s bytes)' "$text" "$bytes"
    else
        printf '%s' "$text"
    fi
}

gold_fail() {
    local file="$1" command="$2" status="$3" expected="$4" actual="$5"
    {
        echo "COUNT CHECK FAIL"
        echo "  file:     $file"
        echo "  command:  $command"
        echo "  exit:     $status"
        echo "  expected: $expected"
        echo "  actual:   $actual"
    } | tee -a "$VERIFY_LOG" >&2
    exit 1
}

capture_count() {
    local out="$1"
    shift
    local status=0
    "$@" >"$out" 2>"$out.err" || status=$?
    printf '%s' "$status"
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

parse_seqtk_reads() {
    local text="$1"
    local line field
    line="$(printf '%s' "$text" | tr -d '\r' | awk 'NF{print; exit}')"
    field="${line%%$'\t'*}"
    field="${field%% *}"
    if [[ ! "$field" =~ ^[0-9]+$ ]]; then
        return 1
    fi
    printf '%s' "$field"
}

should_run_oracle() {
    local expected="$1"
    [[ "$expected" =~ ^[0-9]+$ ]] && [[ "$expected" -lt 1000000 ]]
}

check_same_count() {
    local file="$1" expected="$2" run_native="$3" run_oracle="$4"
    local out status got
    log_verify "  file $file  expected $expected"

    if [[ "$run_oracle" == "true" ]]; then
        out="$CHECK_DIR/oracle.out"
        status="$(capture_count "$out" "$PYTHON" "$ORACLE" "$file")"
        got="$(parse_single_int "$(cat "$out")" || true)"
        if [[ "$status" != "0" || "$got" != "$expected" ]]; then
            gold_fail "$file" "oracle.py $file" "$status" "$expected" \
                "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
        fi
    fi

    out="$CHECK_DIR/zfastq.out"
    status="$(capture_count "$out" "$ZFASTQ" count "$file")"
    got="$(parse_single_int "$(cat "$out")" || true)"
    if [[ "$status" != "0" || "$got" != "$expected" ]]; then
        gold_fail "$file" "z-fastq count $file" "$status" "$expected" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi

    if [[ "$run_native" == "true" ]]; then
        out="$CHECK_DIR/native.out"
        status="$(capture_count "$out" "$ZFASTQ_NATIVE" count "$file")"
        got="$(parse_single_int "$(cat "$out")" || true)"
        if [[ "$status" != "0" || "$got" != "$expected" ]]; then
            gold_fail "$file" "z-fastq-native count $file" "$status" "$expected" \
                "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
        fi
    fi

    out="$CHECK_DIR/seqtk.out"
    status="$(capture_count "$out" "$SEQTK" size "$file")"
    got="$(parse_seqtk_reads "$(cat "$out")" || true)"
    if [[ "$status" != "0" || "$got" != "$expected" ]]; then
        gold_fail "$file" "seqtk size $file" "$status" "$expected" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi

    out="$CHECK_DIR/needletail.out"
    status="$(capture_count "$out" "$NEEDLETAIL" count "$file")"
    got="$(parse_single_int "$(cat "$out")" || true)"
    if [[ "$status" != "0" || "$got" != "$expected" ]]; then
        gold_fail "$file" "needletail-adapter count $file" "$status" "$expected" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi

    out="$CHECK_DIR/helicase.out"
    status="$(capture_count "$out" "$HELICASE" count "$file")"
    got="$(parse_single_int "$(cat "$out")" || true)"
    if [[ "$status" != "0" || "$got" != "$expected" ]]; then
        gold_fail "$file" "helicase-adapter count $file" "$status" "$expected" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi

    out="$CHECK_DIR/fqtools.out"
    status="$(capture_count "$out" "$FQTOOLS" count "$file")"
    got="$(parse_single_int "$(cat "$out")" || true)"
    if [[ "$status" != "0" || "$got" != "$expected" ]]; then
        gold_fail "$file" "fqtools count $file" "$status" "$expected" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi
}

expect_count_fail() {
    local file="$1" binary="$2" label="$3"
    local out="$CHECK_DIR/reject.out" status=0
    "$binary" count "$file" >"$out" 2>"$out.err" || status=$?
    if [[ "$status" == "0" ]]; then
        gold_fail "$file" "$label count $file" "0" "non-zero exit" "$(summarize_output "$out")"
    fi
    log_verify "  PASS reject $label on $(basename "$file") (exit $status)"
}

run_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== count agreement check $TIMESTAMP ==="
    bench_require_tool seqtk
    bench_require_tool fqtools
    bench_require_tool needletail
    bench_require_tool helicase

    local fixture="$FIXTURE_DIR/basic_valid.fastq"
    local fixture_gz="$CHECK_DIR/basic_valid.fastq.gz"
    gzip -c -- "$fixture" >"$fixture_gz"
    local fixture_n
    fixture_n="$("$PYTHON" "$ORACLE" "$fixture")"
    fixture_n="$(parse_single_int "$fixture_n")"
    log_verify "--- fixtures ---"
    check_same_count "$fixture" "$fixture_n" true true
    check_same_count "$fixture_gz" "$fixture_n" true true

    log_verify "--- malformed (z-fastq must reject; seqtk is not compared) ---"
    local bad
    for bad in bad_header.fastq bad_plus.fastq truncated_record.fastq bad_qual_length.fastq; do
        expect_count_fail "$FIXTURE_DIR/$bad" "$ZFASTQ" "z-fastq"
        expect_count_fail "$FIXTURE_DIR/$bad" "$ZFASTQ_NATIVE" "z-fastq-native"
    done

    local name
    log_verify "--- REAL plain ---"
    for name in "${REAL_ORDER[@]}"; do
        check_same_count "${REAL_PLAIN[$name]}" "${REAL_EXPECTED[$name]}" false \
            "$(should_run_oracle "${REAL_EXPECTED[$name]}" && echo true || echo false)"
    done
    log_verify "--- REAL gzip ---"
    for name in "${REAL_ORDER[@]}"; do
        check_same_count "${REAL_GZ[$name]}" "${REAL_EXPECTED[$name]}" true \
            "$(should_run_oracle "${REAL_EXPECTED[$name]}" && echo true || echo false)"
    done

    if $DO_SCALE && $DO_BENCHMARKS; then
        local mb count path
        if $DO_SCALE_SIZE; then
            log_verify "--- scale size ---"
            for mb in "${SIZE_MBS[@]}"; do
                path="$SCALING_DIR/size_${mb}mb.fastq"
                check_same_count "$path" "100000" false true
            done
        fi
        if $DO_SCALE_READS; then
            log_verify "--- scale reads ---"
            for count in "${READS_FIXED_COUNTS[@]}"; do
                path="$SCALING_DIR/reads_fixed_${count}.fastq"
                check_same_count "$path" "$count" false \
                    "$(should_run_oracle "$count" && echo true || echo false)"
            done
        fi
        if $DO_SCALE_GZIP; then
            log_verify "--- scale gzip ---"
            for count in "${READS_FIXED_COUNTS[@]}"; do
                path="$SCALING_DIR/reads_fixed_${count}.fastq.gz"
                check_same_count "$path" "$count" true \
                    "$(should_run_oracle "$count" && echo true || echo false)"
            done
        fi
    fi

    log_verify "ALL PASSED"
}

run_count_tools() {
    local section="$1" workload="$2" file="$3" out_dir="$4"
    local include_native="$5" include_hatched="$6"
    local decoded_bytes="$7"

    local nbytes json
    nbytes="$(file_size_bytes "$file")"

    json="$out_dir/${workload}__z-fastq.json"
    run_zebrac_tool "$section" "$workload" z-fastq z-fastq "$json" \
        "$(zebrac_command "$ZFASTQ" count "$file")" "$nbytes" "$decoded_bytes"

    if [[ "$include_native" == "true" ]]; then
        json="$out_dir/${workload}__z-fastq-native.json"
        run_zebrac_tool "$section" "$workload" z-fastq-native z-fastq "$json" \
            "$(zebrac_command "$ZFASTQ_NATIVE" count "$file")" "$nbytes" "$decoded_bytes"
    fi

    json="$out_dir/${workload}__needletail.json"
    run_zebrac_tool "$section" "$workload" needletail needletail "$json" \
        "$(zebrac_command "$NEEDLETAIL" count "$file")" "$nbytes" "$decoded_bytes"

    json="$out_dir/${workload}__helicase.json"
    run_zebrac_tool "$section" "$workload" helicase helicase "$json" \
        "$(zebrac_command "$HELICASE" count "$file")" "$nbytes" "$decoded_bytes"

    json="$out_dir/${workload}__seqtk.json"
    run_zebrac_tool "$section" "$workload" seqtk seqtk "$json" \
        "$(zebrac_command "$SEQTK" size "$file")" "$nbytes" "$decoded_bytes"

    if [[ "$include_hatched" == "true" ]] && bench_has_tool seqfu; then
        json="$out_dir/${workload}__seqfu.json"
        run_zebrac_tool "$section" "$workload" seqfu seqfu "$json" \
            "$(zebrac_command "$SEQFU" count "$file")" "$nbytes" "$decoded_bytes"
    fi

    json="$out_dir/${workload}__fqtools.json"
    run_zebrac_tool "$section" "$workload" fqtools fqtools "$json" \
        "$(zebrac_command "$FQTOOLS" count "$file")" "$nbytes" "$decoded_bytes"
}

run_perf() {
    local name mb count path decoded
    if $DO_FULL; then
        local plain_dir gzip_dir
        plain_dir="$RESULTS_DIR/perf_plain_${TIMESTAMP}"
        gzip_dir="$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
        mkdir -p "$plain_dir" "$gzip_dir"
        echo "=== perf_plain ==="
        for name in "${REAL_ORDER[@]}"; do
            decoded="${REAL_DECODED[$name]}"
            run_count_tools perf_plain "$name" "${REAL_PLAIN[$name]}" "$plain_dir" false true "$decoded"
        done
        echo "=== perf_gzip ==="
        for name in "${REAL_ORDER[@]}"; do
            decoded="${REAL_DECODED[$name]}"
            run_count_tools perf_gzip "$name" "${REAL_GZ[$name]}" "$gzip_dir" true true "$decoded"
        done
    fi

    if $DO_SCALE && $DO_SCALE_SIZE; then
        local size_dir="$RESULTS_DIR/scale_size_${TIMESTAMP}"
        mkdir -p "$size_dir"
        echo "=== scale_size ==="
        for mb in "${SIZE_MBS[@]}"; do
            path="$SCALING_DIR/size_${mb}mb.fastq"
            decoded="$(file_size_bytes "$path")"
            run_count_tools scale_size "${mb}mb" "$path" "$size_dir" false false "$decoded"
        done
    fi

    if $DO_SCALE && $DO_SCALE_READS; then
        local reads_dir="$RESULTS_DIR/scale_reads_${TIMESTAMP}"
        mkdir -p "$reads_dir"
        echo "=== scale_reads ==="
        for count in "${READS_FIXED_COUNTS[@]}"; do
            path="$SCALING_DIR/reads_fixed_${count}.fastq"
            decoded="$(file_size_bytes "$path")"
            run_count_tools scale_reads "$count" "$path" "$reads_dir" false false "$decoded"
        done
    fi

    if $DO_SCALE && $DO_SCALE_GZIP; then
        local gz_dir="$RESULTS_DIR/scale_gzip_${TIMESTAMP}"
        mkdir -p "$gz_dir"
        echo "=== scale_gzip ==="
        for count in "${READS_FIXED_COUNTS[@]}"; do
            path="$SCALING_DIR/reads_fixed_${count}.fastq.gz"
            decoded="$(file_size_bytes "$SCALING_DIR/reads_fixed_${count}.fastq")"
            run_count_tools scale_gzip "$count" "$path" "$gz_dir" true false "$decoded"
        done
    fi
}

write_manifest() {
    local python_bin
    python_bin="$(report_python)"
    "$python_bin" - "$MANIFEST_PATH" <<'PY'
import json, os, sys
from pathlib import Path

manifest = Path(sys.argv[1])
ts = os.environ["COUNT_TS"]
results = Path(os.environ["COUNT_RESULTS"])
out = {
    "schema_version": "count-run.v1",
    "timestamp": ts,
    "runner": "zebrac",
    "mode": "warm",
    "suite": "count",
    "real_set": os.environ.get("COUNT_REAL_SET", "full"),
    "zebrac": os.environ.get("COUNT_ZEBRAC_VER", ""),
    "z_fastq": os.environ.get("COUNT_ZFASTQ_VER", ""),
    "z_fastq_native": os.environ.get("COUNT_ZFASTQ_NATIVE_VER", ""),
    "z_fastq_bytes": int(os.environ.get("COUNT_ZFASTQ_BYTES", "0") or 0),
    "z_fastq_native_bytes": int(os.environ.get("COUNT_ZFASTQ_NATIVE_BYTES", "0") or 0),
    "runs": int(os.environ["COUNT_RUNS"]),
    "warmup": int(os.environ["COUNT_WARMUP"]),
    "duration_ms": int(os.environ["COUNT_DURATION"]),
    "metadata": f"metadata_{ts}.jsonl",
    "verify_log": f"verify_{ts}.log",
    "verify_skipped": os.environ.get("COUNT_VERIFY_SKIPPED") == "1",
    "verify_pass": os.environ.get("COUNT_VERIFY_PASS") or None,
    "tools": json.loads(os.environ.get("COUNT_TOOLS_JSON", "{}")),
    "sections": {},
    "skip_full": os.environ.get("COUNT_SKIP_FULL") == "1",
    "skip_scale": os.environ.get("COUNT_SKIP_SCALE") == "1",
}
for key, prefix in (
    ("perf_plain", "perf_plain_"),
    ("perf_gzip", "perf_gzip_"),
    ("scale_size", "scale_size_"),
    ("scale_reads", "scale_reads_"),
    ("scale_gzip", "scale_gzip_"),
):
    path = results / f"{prefix}{ts}"
    if path.is_dir():
        out["sections"][key] = path.name
log = results / f"verify_{ts}.log"
if log.is_file() and "ALL PASSED" in log.read_text(encoding="utf-8", errors="replace"):
    out["verify_skipped"] = False
    out["verify_pass"] = "ALL PASSED"
manifest.write_text(json.dumps(out, indent=2) + "\n")
PY
    printf '%s\n' "$TIMESTAMP" >"$RESULTS_DIR/LATEST"
}

# --- run ---
echo "z-fastq count bench  $TIMESTAMP"
echo

build_subjects
ensure_real_data
if $DO_SCALE && $DO_BENCHMARKS; then
    ensure_scaling_fixtures
fi

VERIFY_PASS=""
VERIFY_SKIPPED=1
if $DO_TESTS; then
    VERIFY_SKIPPED=0
    run_tests
    VERIFY_PASS="ALL PASSED"
elif $DO_BENCHMARKS; then
    echo "warning: skipping the count-agreement check (--skip-tests); not for a published report" >&2
fi

if $DO_BENCHMARKS; then
    bench_require_tool zebrac
    bench_require_tool seqtk
    bench_require_tool fqtools
    bench_require_tool needletail
    bench_require_tool helicase
    run_perf
fi

export COUNT_TS="$TIMESTAMP"
export COUNT_RESULTS="$RESULTS_DIR"
export COUNT_REAL_SET="$($SMALL_REAL && echo small || echo full)"
export COUNT_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
export COUNT_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
export COUNT_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
export COUNT_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
export COUNT_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
export COUNT_RUNS="$RUNS"
export COUNT_WARMUP="$WARMUP"
export COUNT_DURATION="$ZEBRAC_DURATION_MS"
export COUNT_VERIFY_SKIPPED="$VERIFY_SKIPPED"
export COUNT_VERIFY_PASS="$VERIFY_PASS"
export COUNT_SKIP_FULL="$($DO_FULL && echo 0 || echo 1)"
export COUNT_SKIP_SCALE="$($DO_SCALE && echo 0 || echo 1)"
export COUNT_SEQTK_VER="$(bench_tool_version seqtk || true)"
export COUNT_FQTOOLS_VER="$(bench_tool_version fqtools || true)"
export COUNT_NEEDLETAIL_VER="$(bench_tool_version needletail || true)"
export COUNT_HELICASE_VER="$(bench_tool_version helicase || true)"
export COUNT_SEQFU_VER="$(bench_tool_version seqfu || true)"
COUNT_TOOLS_JSON="$("$PYTHON" - <<'PY'
import json, os
print(json.dumps({
    "needletail": os.environ.get("COUNT_NEEDLETAIL_VER", ""),
    "helicase": os.environ.get("COUNT_HELICASE_VER", ""),
    "seqtk": os.environ.get("COUNT_SEQTK_VER", ""),
    "seqfu": os.environ.get("COUNT_SEQFU_VER", ""),
    "fqtools": os.environ.get("COUNT_FQTOOLS_VER", ""),
}))
PY
)"
export COUNT_TOOLS_JSON
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

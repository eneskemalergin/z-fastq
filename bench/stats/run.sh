#!/usr/bin/env bash
# Stats benchmark runner: same-job stats check, then zebrac.
#
# Timed files come from features.tsv (stats/<publication|small>/time).
#
# Usage:
#   bash bench/stats/run.sh [options]
#
# Defaults: stats-agreement check first, then perf (runs=25, warmup=5, duration=5000).
# A field mismatch or unexpected exit stops the run before any timed command.
#
#   bash bench/stats/run.sh
#   bash bench/stats/run.sh --skip-tests --skip-report
#   bash bench/stats/run.sh --small-real --runs 5 --warmup 3
#   STATS_RUN_TIMESTAMP=<ts> bash bench/stats/run.sh --skip-tests
#
#   -h|--help  print this header

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$BENCH_ROOT")"
RESULTS_DIR="$SCRIPT_DIR/results"
DATA_DIR="$BENCH_ROOT/shared/data"
PLAIN_DIR="$BENCH_ROOT/shared/cache/plain"
ORACLE="$SCRIPT_DIR/oracle.py"
FIXTURE_DIR="$PROJECT_ROOT/tests/data/synthetic"

# shellcheck disable=SC1091
source "$BENCH_ROOT/shared/tools.sh"

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
    case $1 in
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

if [[ -n "${STATS_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$STATS_RUN_TIMESTAMP"
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

STATS_SET="publication"
if $SMALL_REAL; then
    STATS_SET="small"
fi

TIME_IDS=()
declare -A REAL_GZ=()
declare -A REAL_PLAIN=()
declare -A REAL_EXPECTED=()
declare -A REAL_DECODED=()

stats_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "stats" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$script"
}

stats_dataset_json() {
    local id="$1"
    catalog_has_id "$id" || return 1
    printf '{'
    printf '"manifest_id":%s,' "$(zebrac_json_string "$id")"
    printf '"filename":%s,' "$(zebrac_json_string "${CATALOG_FILENAME[$id]}")"
    printf '"expected_records":%s,' "${CATALOG_EXPECTED[$id]}"
    printf '"accession":%s' "$(zebrac_json_string "${CATALOG_ACCESSION[$id]}")"
    printf '}'
}

stats_datasets_json() {
    local id sep=""
    printf '{'
    for id in "${TIME_IDS[@]}"; do
        printf '%s%s:' "${sep}" "$(zebrac_json_string "${id}")"
        stats_dataset_json "${id}" || return 1
        sep=','
    done
    printf '}'
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    stats_add_command "$section" "$workload" "$tool" "$family" \
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

bind_dataset() {
    local id="$1"
    local filename
    filename="$(bench_catalog field "$id" filename)"
    REAL_GZ["$id"]="$DATA_DIR/$filename"
    REAL_PLAIN["$id"]="$PLAIN_DIR/${filename%.gz}"
    REAL_EXPECTED["$id"]="$(bench_catalog expected "$id")"
}

ensure_real_data() {
    local download_args=(--suite stats)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (stats/$STATS_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    local id
    mapfile -t TIME_IDS < <(bench_catalog ids --suite stats --set "$STATS_SET" --role time --no-expand)
    if [[ ${#TIME_IDS[@]} -lt 1 ]]; then
        echo "error: stats/$STATS_SET/time must list at least one id" >&2
        exit 1
    fi
    for id in "${TIME_IDS[@]}"; do
        bind_dataset "$id"
        [[ -f "${REAL_GZ[$id]}" ]] || {
            echo "error: missing gzip dataset ${REAL_GZ[$id]}" >&2
            exit 1
        }
        [[ -f "${REAL_PLAIN[$id]}" ]] || {
            echo "error: missing plain cache ${REAL_PLAIN[$id]}" >&2
            exit 1
        }
        REAL_DECODED["$id"]="$(file_size_bytes "${REAL_PLAIN[$id]}")"
        echo "  $id gzip: ${REAL_GZ[$id]}  (expected ${REAL_EXPECTED[$id]})"
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
    if (( bytes > 400 )); then
        printf '%s... (%s bytes)' "$text" "$bytes"
    else
        printf '%s' "$text"
    fi
}

gold_fail() {
    local file="$1" command="$2" status="$3" expected="$4" actual="$5"
    {
        echo "STATS CHECK FAIL"
        echo "  file:     $file"
        echo "  command:  $command"
        echo "  exit:     $status"
        echo "  expected: $expected"
        echo "  actual:   $actual"
    } | tee -a "$VERIFY_LOG" >&2
    exit 1
}

capture_out() {
    local out="$1"
    shift
    local status=0
    "$@" >"$out" 2>"$out.err" || status=$?
    printf '%s' "$status"
}

should_run_oracle() {
    local expected="$1"
    [[ "$expected" =~ ^[0-9]+$ ]] && [[ "$expected" -lt 1000000 ]]
}

agree_or_fail() {
    local file="$1" command="$2" status="$3" gold="$4" peer="$5"
    if [[ "$status" != "0" ]]; then
        gold_fail "$file" "$command" "$status" "exit 0 and matching fields" \
            "$(summarize_output "$peer"; printf '\n'; summarize_output "$peer.err")"
    fi
    if ! "$PYTHON" "$ORACLE" --agree-human "$gold" "$peer" 2>"$CHECK_DIR/agree.err"; then
        gold_fail "$file" "$command" "$status" "human fields matching $gold" \
            "$(summarize_output "$peer"; printf '\n'; cat "$CHECK_DIR/agree.err")"
    fi
}

overlap_or_fail() {
    local flag="$1" file="$2" command="$3" status="$4" gold="$5" peer="$6"
    if [[ "$status" != "0" ]]; then
        gold_fail "$file" "$command" "$status" "exit 0 and overlapping fields" \
            "$(summarize_output "$peer"; printf '\n'; summarize_output "$peer.err")"
    fi
    if ! "$PYTHON" "$ORACLE" "$flag" "$gold" "$peer" 2>"$CHECK_DIR/agree.err"; then
        gold_fail "$file" "$command" "$status" "overlapping fields vs z-fastq stats" \
            "$(summarize_output "$peer"; printf '\n'; cat "$CHECK_DIR/agree.err")"
    fi
}

check_same_stats() {
    local file="$1" expected="$2" run_native="$3" run_oracle="$4"
    local out status gold
    log_verify "  file $file  expected_records $expected"

    out="$CHECK_DIR/zfastq.out"
    status="$(capture_out "$out" "$ZFASTQ" stats "$file")"
    if [[ "$status" != "0" ]]; then
        gold_fail "$file" "z-fastq stats $file" "$status" "exit 0" \
            "$(summarize_output "$out"; printf '\n'; summarize_output "$out.err")"
    fi
    gold="$out"

    if [[ "$run_oracle" == "true" ]]; then
        out="$CHECK_DIR/oracle.out"
        status="$(capture_out "$out" "$PYTHON" "$ORACLE" "$file")"
        agree_or_fail "$file" "oracle.py $file" "$status" "$gold" "$out"
    fi

    if [[ "$run_native" == "true" ]]; then
        out="$CHECK_DIR/native.out"
        status="$(capture_out "$out" "$ZFASTQ_NATIVE" stats "$file")"
        agree_or_fail "$file" "z-fastq-native stats $file" "$status" "$gold" "$out"
    fi

    out="$CHECK_DIR/needletail.out"
    status="$(capture_out "$out" "$NEEDLETAIL" stats "$file")"
    agree_or_fail "$file" "needletail-adapter stats $file" "$status" "$gold" "$out"

    out="$CHECK_DIR/helicase.out"
    status="$(capture_out "$out" "$HELICASE" stats "$file")"
    agree_or_fail "$file" "helicase-adapter stats $file" "$status" "$gold" "$out"

    out="$CHECK_DIR/seqkit.out"
    status="$(capture_out "$out" "$SEQKIT" stats -a -T -j 1 "$file")"
    overlap_or_fail --overlap-seqkit "$file" "seqkit stats -a -T -j 1 $file" \
        "$status" "$gold" "$out"

    if bench_has_tool seqfu; then
        out="$CHECK_DIR/seqfu.out"
        status="$(capture_out "$out" "$SEQFU" stats --threads 1 --gc "$file")"
        overlap_or_fail --overlap-seqfu "$file" "seqfu stats --threads 1 --gc $file" \
            "$status" "$gold" "$out"
    fi
}

expect_stats_fail() {
    local file="$1" binary="$2" label="$3"
    local out="$CHECK_DIR/reject.out" status=0
    "$binary" stats "$file" >"$out" 2>"$out.err" || status=$?
    if [[ "$status" == "0" ]]; then
        gold_fail "$file" "$label stats $file" "0" "non-zero exit" "$(summarize_output "$out")"
    fi
    log_verify "  PASS reject $label on $(basename "$file") (exit $status)"
}

run_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== stats agreement check $TIMESTAMP ==="
    bench_require_tool needletail
    bench_require_tool helicase
    bench_require_tool seqkit

    local fixture="$FIXTURE_DIR/basic_valid.fastq"
    local fixture_gz="$CHECK_DIR/basic_valid.fastq.gz"
    gzip -c -- "$fixture" >"$fixture_gz"
    log_verify "--- fixtures ---"
    check_same_stats "$fixture" 5 true true
    check_same_stats "$fixture_gz" 5 true true

    log_verify "--- malformed (z-fastq must reject; SeqKit/SeqFu are not compared) ---"
    local bad
    for bad in bad_header.fastq bad_plus.fastq truncated_record.fastq bad_qual_length.fastq; do
        expect_stats_fail "$FIXTURE_DIR/$bad" "$ZFASTQ" "z-fastq"
        expect_stats_fail "$FIXTURE_DIR/$bad" "$ZFASTQ_NATIVE" "z-fastq-native"
    done

    local id
    log_verify "--- REAL plain ---"
    for id in "${TIME_IDS[@]}"; do
        check_same_stats "${REAL_PLAIN[$id]}" "${REAL_EXPECTED[$id]}" false \
            "$(should_run_oracle "${REAL_EXPECTED[$id]}" && echo true || echo false)"
    done
    log_verify "--- REAL gzip ---"
    for id in "${TIME_IDS[@]}"; do
        check_same_stats "${REAL_GZ[$id]}" "${REAL_EXPECTED[$id]}" true \
            "$(should_run_oracle "${REAL_EXPECTED[$id]}" && echo true || echo false)"
    done

    log_verify "ALL PASSED"
}

run_stats_tools() {
    local section="$1" workload="$2" file="$3" out_dir="$4"
    local include_native="$5" include_hatched="$6"
    local decoded_bytes="$7"

    local nbytes json
    nbytes="$(file_size_bytes "$file")"

    json="$out_dir/${workload}__z-fastq.json"
    run_zebrac_tool "$section" "$workload" z-fastq z-fastq "$json" \
        "$(zebrac_command "$ZFASTQ" stats "$file")" "$nbytes" "$decoded_bytes"

    if [[ "$include_native" == "true" ]]; then
        json="$out_dir/${workload}__z-fastq-native.json"
        run_zebrac_tool "$section" "$workload" z-fastq-native z-fastq "$json" \
            "$(zebrac_command "$ZFASTQ_NATIVE" stats "$file")" "$nbytes" "$decoded_bytes"
    fi

    json="$out_dir/${workload}__needletail.json"
    run_zebrac_tool "$section" "$workload" needletail needletail "$json" \
        "$(zebrac_command "$NEEDLETAIL" stats "$file")" "$nbytes" "$decoded_bytes"

    json="$out_dir/${workload}__helicase.json"
    run_zebrac_tool "$section" "$workload" helicase helicase "$json" \
        "$(zebrac_command "$HELICASE" stats "$file")" "$nbytes" "$decoded_bytes"

    json="$out_dir/${workload}__seqkit.json"
    run_zebrac_tool "$section" "$workload" seqkit seqkit "$json" \
        "$(zebrac_command "$SEQKIT" stats -a -T -j 1 "$file")" "$nbytes" "$decoded_bytes"

    if [[ "$include_hatched" == "true" ]] && bench_has_tool seqfu; then
        json="$out_dir/${workload}__seqfu.json"
        run_zebrac_tool "$section" "$workload" seqfu seqfu "$json" \
            "$(zebrac_command "$SEQFU" stats --threads 1 --gc "$file")" "$nbytes" "$decoded_bytes"
    fi
}

run_perf() {
    local id decoded
    if $DO_FULL; then
        local plain_dir gzip_dir
        plain_dir="$RESULTS_DIR/perf_plain_${TIMESTAMP}"
        gzip_dir="$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
        mkdir -p "$plain_dir" "$gzip_dir"
        echo "=== perf_plain ==="
        for id in "${TIME_IDS[@]}"; do
            decoded="${REAL_DECODED[$id]}"
            run_stats_tools perf_plain "$id" "${REAL_PLAIN[$id]}" "$plain_dir" false true "$decoded"
        done
        echo "=== perf_gzip ==="
        for id in "${TIME_IDS[@]}"; do
            decoded="${REAL_DECODED[$id]}"
            run_stats_tools perf_gzip "$id" "${REAL_GZ[$id]}" "$gzip_dir" true true "$decoded"
        done
    fi
}

write_manifest() {
    local real_set="full"
    $SMALL_REAL && real_set="small"

    local verify_skipped_json=false
    [[ "$VERIFY_SKIPPED" == 1 ]] && verify_skipped_json=true
    local verify_pass_json=null
    [[ -n "$VERIFY_PASS" ]] && verify_pass_json="$(zebrac_json_string "$VERIFY_PASS")"
    if [[ -f "$VERIFY_LOG" ]] && grep -Fq 'ALL PASSED' "$VERIFY_LOG"; then
        verify_skipped_json=false
        verify_pass_json='"ALL PASSED"'
    fi

    local skip_full_json=false
    $DO_FULL || skip_full_json=true

    local temporary_manifest="${MANIFEST_PATH}.tmp.$$"
    {
        printf '{\n'
        printf '  "schema_version": "stats-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "stats",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$real_set")"
        printf '  "datasets": %s,\n' "$STATS_DATASETS_JSON"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${STATS_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${STATS_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${STATS_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${STATS_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${STATS_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$STATS_TOOLS_JSON"
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

# --- run ---
echo "z-fastq stats bench  $TIMESTAMP"
echo

build_subjects
ensure_real_data

VERIFY_PASS=""
VERIFY_SKIPPED=1
if $DO_TESTS; then
    VERIFY_SKIPPED=0
    run_tests
    VERIFY_PASS="ALL PASSED"
elif $DO_BENCHMARKS; then
    echo "warning: skipping the stats-agreement check (--skip-tests); not for a published report" >&2
fi

if $DO_BENCHMARKS; then
    bench_require_tool zebrac
    bench_require_tool needletail
    bench_require_tool helicase
    bench_require_tool seqkit
    run_perf
fi

STATS_DATASETS_JSON="$(stats_datasets_json)"
STATS_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
STATS_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
STATS_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
STATS_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
STATS_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
STATS_NEEDLETAIL_VER="$(bench_tool_version needletail || true)"
STATS_HELICASE_VER="$(bench_tool_version helicase || true)"
STATS_SEQKIT_VER="$(bench_tool_version seqkit || true)"
STATS_SEQFU_VER="$(bench_tool_version seqfu || true)"
STATS_TOOLS_JSON="$(
    printf '{'
    printf '"needletail":%s,' "$(zebrac_json_string "$STATS_NEEDLETAIL_VER")"
    printf '"helicase":%s,' "$(zebrac_json_string "$STATS_HELICASE_VER")"
    printf '"seqkit":%s,' "$(zebrac_json_string "$STATS_SEQKIT_VER")"
    printf '"seqfu":%s' "$(zebrac_json_string "$STATS_SEQFU_VER")"
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

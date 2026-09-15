#!/usr/bin/env bash
# Count benchmark runner: same-job count check, then zebrac.
#
# Timed files and extra agreement files come from features.tsv
# (count/<publication|small>/time and .../check).
#
# Usage:
#   bash bench/count/run.sh [options]
#
# Defaults: count-agreement check first, then perf (runs=25, warmup=5, duration=5000).
# A count mismatch or unexpected exit stops the run before any timed command.
#
#   bash bench/count/run.sh
#   bash bench/count/run.sh --skip-tests --skip-report
#   bash bench/count/run.sh --small-real --runs 5 --warmup 3
#   COUNT_RUN_TIMESTAMP=<ts> bash bench/count/run.sh --skip-tests
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

COUNT_SET="publication"
if $SMALL_REAL; then
    COUNT_SET="small"
fi

DENSE_ID=""
VARIABLE_ID=""
LONG_ID=""
REAL_ORDER=(Dense Variable Long)
declare -A REAL_GZ=()
declare -A REAL_PLAIN=()
declare -A REAL_EXPECTED=()
declare -A REAL_DECODED=()
declare -A EXTRA_GZ=()
declare -A EXTRA_PLAIN=()
declare -A EXTRA_EXPECTED=()
EXTRA_ORDER=()

count_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" script="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "count" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$script"
}

count_dataset_json() {
    local id="$1"
    catalog_has_id "$id" || return 1
    printf '{'
    printf '"manifest_id":%s,' "$(zebrac_json_string "$id")"
    printf '"filename":%s,' "$(zebrac_json_string "${CATALOG_FILENAME[$id]}")"
    printf '"expected_records":%s,' "${CATALOG_EXPECTED[$id]}"
    printf '"accession":%s' "$(zebrac_json_string "${CATALOG_ACCESSION[$id]}")"
    printf '}'
}

count_datasets_json() {
    local dense_id="$1"
    local variable_id="$2"
    local long_id="$3"
    local extras_csv="${4:-}"
    local -a extra_ids=()
    local id key sep=""
    if [[ -n "${extras_csv}" ]]; then
        IFS=',' read -r -a extra_ids <<< "${extras_csv}"
    fi

    printf '{'
    for key in Dense Variable Long; do
        case "${key}" in
            Dense) id="${dense_id}" ;;
            Variable) id="${variable_id}" ;;
            Long) id="${long_id}" ;;
        esac
        printf '%s%s:' "${sep}" "$(zebrac_json_string "${key}")"
        count_dataset_json "${id}" || return 1
        sep=','
    done
    for id in "${extra_ids[@]}"; do
        [[ -n "${id}" ]] || continue
        printf '%s%s:' "${sep}" "$(zebrac_json_string "${id}")"
        count_dataset_json "${id}" || return 1
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

bind_dataset() {
    local slot="$1" id="$2"
    local filename
    filename="$(bench_catalog field "$id" filename)"
    REAL_GZ["$slot"]="$DATA_DIR/$filename"
    REAL_PLAIN["$slot"]="$PLAIN_DIR/${filename%.gz}"
    REAL_EXPECTED["$slot"]="$(bench_catalog expected "$id")"
}

ensure_real_data() {
    local download_args=(--suite count)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (count/$COUNT_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    local time_ids check_ids id filename
    mapfile -t time_ids < <(bench_catalog ids --suite count --set "$COUNT_SET" --role time --no-expand)
    mapfile -t check_ids < <(bench_catalog ids --suite count --set "$COUNT_SET" --role check --no-expand)
    if [[ ${#time_ids[@]} -ne 3 ]]; then
        echo "error: count/$COUNT_SET/time must list exactly 3 ids" >&2
        exit 1
    fi
    DENSE_ID="${time_ids[0]}"
    VARIABLE_ID="${time_ids[1]}"
    LONG_ID="${time_ids[2]}"
    bind_dataset Dense "$DENSE_ID"
    bind_dataset Variable "$VARIABLE_ID"
    bind_dataset Long "$LONG_ID"

    local timed=""
    timed=" ${DENSE_ID} ${VARIABLE_ID} ${LONG_ID} "
    EXTRA_ORDER=()
    for id in "${check_ids[@]}"; do
        [[ "$timed" == *" $id "* ]] && continue
        filename="$(bench_catalog field "$id" filename)"
        EXTRA_ORDER+=("$id")
        EXTRA_GZ["$id"]="$DATA_DIR/$filename"
        EXTRA_PLAIN["$id"]="$PLAIN_DIR/${filename%.gz}"
        EXTRA_EXPECTED["$id"]="$(bench_catalog expected "$id")"
    done

    for slot in "${REAL_ORDER[@]}"; do
        [[ -f "${REAL_GZ[$slot]}" ]] || {
            echo "error: missing gzip dataset ${REAL_GZ[$slot]}" >&2
            exit 1
        }
        [[ -f "${REAL_PLAIN[$slot]}" ]] || {
            echo "error: missing plain cache ${REAL_PLAIN[$slot]}" >&2
            exit 1
        }
        REAL_DECODED["$slot"]="$(file_size_bytes "${REAL_PLAIN[$slot]}")"
    done
    for id in "${EXTRA_ORDER[@]}"; do
        [[ -f "${EXTRA_GZ[$id]}" && -f "${EXTRA_PLAIN[$id]}" ]] || {
            echo "error: missing extra check files for $id" >&2
            exit 1
        }
    done

    echo "  Dense gzip: ${REAL_GZ[Dense]}  (manifest id $DENSE_ID, expected ${REAL_EXPECTED[Dense]})"
    echo "  Variable gzip: ${REAL_GZ[Variable]}  (manifest id $VARIABLE_ID, expected ${REAL_EXPECTED[Variable]})"
    echo "  Long gzip:  ${REAL_GZ[Long]}  (manifest id $LONG_ID, expected ${REAL_EXPECTED[Long]})"
    for id in "${EXTRA_ORDER[@]}"; do
        echo "  extra check $id: ${EXTRA_GZ[$id]}  (expected ${EXTRA_EXPECTED[$id]})"
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

    local name id
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

    for id in "${EXTRA_ORDER[@]}"; do
        log_verify "--- extra check $id ---"
        check_same_count "${EXTRA_PLAIN[$id]}" "${EXTRA_EXPECTED[$id]}" false \
            "$(should_run_oracle "${EXTRA_EXPECTED[$id]}" && echo true || echo false)"
        check_same_count "${EXTRA_GZ[$id]}" "${EXTRA_EXPECTED[$id]}" true \
            "$(should_run_oracle "${EXTRA_EXPECTED[$id]}" && echo true || echo false)"
    done

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
    local name decoded
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
        printf '  "schema_version": "count-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "count",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$real_set")"
        printf '  "datasets": %s,\n' "$COUNT_DATASETS_JSON"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${COUNT_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${COUNT_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${COUNT_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${COUNT_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${COUNT_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$COUNT_TOOLS_JSON"
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
echo "z-fastq count bench  $TIMESTAMP"
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

COUNT_CHECK_EXTRA="$(IFS=,; printf '%s' "${EXTRA_ORDER[*]}")"
COUNT_DATASETS_JSON="$(count_datasets_json "$DENSE_ID" "$VARIABLE_ID" "$LONG_ID" "$COUNT_CHECK_EXTRA")"
COUNT_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
COUNT_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
COUNT_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
COUNT_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
COUNT_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
COUNT_SEQTK_VER="$(bench_tool_version seqtk || true)"
COUNT_FQTOOLS_VER="$(bench_tool_version fqtools || true)"
COUNT_NEEDLETAIL_VER="$(bench_tool_version needletail || true)"
COUNT_HELICASE_VER="$(bench_tool_version helicase || true)"
COUNT_SEQFU_VER="$(bench_tool_version seqfu || true)"
COUNT_TOOLS_JSON="$(
    printf '{'
    printf '"needletail":%s,' "$(zebrac_json_string "$COUNT_NEEDLETAIL_VER")"
    printf '"helicase":%s,' "$(zebrac_json_string "$COUNT_HELICASE_VER")"
    printf '"seqtk":%s,' "$(zebrac_json_string "$COUNT_SEQTK_VER")"
    printf '"seqfu":%s,' "$(zebrac_json_string "$COUNT_SEQFU_VER")"
    printf '"fqtools":%s' "$(zebrac_json_string "$COUNT_FQTOOLS_VER")"
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

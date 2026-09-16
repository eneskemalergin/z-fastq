#!/usr/bin/env bash
# Check benchmark runner: contract checks, capability probes, then zebrac.
#
# Workloads come from features.tsv:
#   check/<publication|small>/<se|paired|casava|slash|interleaved|gzip>
#
# Validator semantics are not identical on every edge case. z-fastq and both
# z-fastq binaries are hard-checked; external validators are timed only after
# they accept the positive workload and are recorded as descriptive peers.
#
# Usage:
#   bash bench/check/run.sh
#   bash bench/check/run.sh --small-real --runs 5 --warmup 3
#   bash bench/check/run.sh --skip-tests --skip-report
#   CHECK_RUN_TIMESTAMP=<ts> bash bench/check/run.sh --skip-tests

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

if [[ -n "${CHECK_RUN_TIMESTAMP:-}" ]]; then
    TIMESTAMP="$CHECK_RUN_TIMESTAMP"
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

CHECK_SET="publication"
if $SMALL_REAL; then
    CHECK_SET="small"
fi

declare -a SE_IDS=()
declare -a PAIRED_IDS=()
declare -a CASAVA_IDS=()
declare -a SLASH_IDS=()
declare -a INTERLEAVED_IDS=()
declare -a GZIP_IDS=()
declare -a WORKLOAD_KEYS=()
declare -a WORKLOAD_ROWS=()
declare -A DATA_GZ=()
declare -A DATA_PLAIN=()
declare -A DATA_EXPECTED=()
declare -A DATA_DECODED=()
declare -A LANE_ENABLED=()
declare -A LANE_REASON=()
declare -a PEER_FIXTURE_ROWS=()

CHECK_PEER_ORDER=(fq fastqvalidator seqfu fqtools)

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
    local key plain=true gzip=true id
    key="$(workload_key "$role" "${ids[@]}")"
    [[ "$role" == gzip ]] && plain=false
    WORKLOAD_KEYS+=("$key")
    local -a accessions=()
    local -a records=()
    for id in "${ids[@]}"; do
        accessions+=("${CATALOG_ACCESSION[$id]}")
        records+=("${CATALOG_EXPECTED[$id]}")
    done
    WORKLOAD_ROWS+=("$(printf \
        '{"key":%s,"role":%s,"ids":%s,"accessions":%s,"records":%s,"plain":%s,"gzip":%s}' \
        "$(zebrac_json_string "$key")" \
        "$(zebrac_json_string "$role")" \
        "$(json_string_array "${ids[@]}")" \
        "$(json_string_array "${accessions[@]}")" \
        "$(json_number_array "${records[@]}")" \
        "$plain" \
        "$gzip")")
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
        bench_catalog ids --suite check --set "$CHECK_SET" --role "$role" --no-expand
    )
    ((${#destination[@]} > 0)) || {
        echo "error: check/$CHECK_SET/$role has no datasets" >&2
        exit 1
    }
}

ensure_real_data() {
    local download_args=(--suite check)
    if $SMALL_REAL; then
        download_args+=(--small)
    fi
    echo "Ensuring REAL gzip FASTQ under $DATA_DIR (check/$CHECK_SET) ..."
    bash "$BENCH_ROOT/shared/download_data.sh" "${download_args[@]}"

    load_role_ids se SE_IDS
    load_role_ids paired PAIRED_IDS
    load_role_ids casava CASAVA_IDS
    load_role_ids slash SLASH_IDS
    load_role_ids interleaved INTERLEAVED_IDS
    load_role_ids gzip GZIP_IDS

    local -a all_ids=(
        "${SE_IDS[@]}"
        "${PAIRED_IDS[@]}"
        "${CASAVA_IDS[@]}"
        "${SLASH_IDS[@]}"
        "${INTERLEAVED_IDS[@]}"
        "${GZIP_IDS[@]}"
    )
    local id
    for id in "${all_ids[@]}"; do
        bind_dataset "$id"
    done

    local se_id
    for se_id in "${SE_IDS[@]}"; do
        register_workload se "$se_id"
    done
    register_workload paired "${PAIRED_IDS[@]}"
    register_workload casava "${CASAVA_IDS[@]}"
    register_workload slash "${SLASH_IDS[@]}"
    register_workload interleaved "${INTERLEAVED_IDS[@]}"
    register_workload gzip "${GZIP_IDS[@]}"

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

check_fail() {
    local label="$1" expected="$2" actual="$3"
    {
        echo "CHECK CONTRACT FAIL"
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

build_check_args() {
    local -n result="$1"
    local tool="$2"
    local role="$3"
    shift 3
    local -a paths=("$@")
    result=()

    case "$tool" in
        z-fastq|z-fastq-native)
            local binary="$ZFASTQ"
            [[ "$tool" == z-fastq-native ]] && binary="$ZFASTQ_NATIVE"
            result=("$binary" check)
            case "$role" in
                paired|casava|slash)
                    result+=(--paired --pair-names illumina)
                    ;;
                interleaved)
                    result+=(--interleaved --pair-names illumina)
                    ;;
            esac
            result+=("${paths[@]}")
            ;;
        fq)
            result=("$FQ" lint --lint-mode panic "${paths[@]}")
            ;;
        fastqvalidator)
            result=("$FASTQVALIDATOR" --file "${paths[0]}" --minReadLen 0)
            if [[ "$role" == interleaved ]]; then
                result+=(--interleaved)
            else
                result+=(--disableSeqIDCheck)
            fi
            result+=(--quiet)
            ;;
        seqfu)
            result=("$SEQFU" check --deep --quiet)
            if [[ "${#paths[@]}" == 1 ]]; then
                result+=(--no-paired)
            fi
            result+=("${paths[@]}")
            ;;
        fqtools)
            result=("$FQTOOLS" -d -a -u -l -q s)
            [[ "$role" == interleaved ]] && result+=(-i)
            result+=(validate "${paths[@]}")
            ;;
        *)
            echo "error: unknown check tool $tool" >&2
            return 1
            ;;
    esac
}

tool_supports_role() {
    local tool="$1" role="$2"
    case "$tool:$role" in
        fq:se|fq:paired|fq:casava|fq:slash|fq:gzip)
            return 0
            ;;
        fastqvalidator:se|fastqvalidator:interleaved|fastqvalidator:gzip)
            return 0
            ;;
        seqfu:se|seqfu:paired|seqfu:casava|seqfu:slash|seqfu:gzip)
            return 0
            ;;
        fqtools:se|fqtools:paired|fqtools:casava|fqtools:slash|fqtools:interleaved|fqtools:gzip)
            return 0
            ;;
        *)
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
    local fixture="$1" tool="$2" role="$3"
    shift 3
    local -a paths=("$@")
    if ! tool_supports_role "$tool" "$role"; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "not supported for role"
        log_verify "  fixture $fixture $tool unsupported"
        return 0
    fi
    if ! bench_has_tool "$tool"; then
        record_peer_fixture "$fixture" "$tool" "unsupported" "tool unavailable"
        log_verify "  fixture $fixture $tool unavailable"
        return 0
    fi
    local -a args=()
    build_check_args args "$tool" "$role" "${paths[@]}"
    local stem="fixture_${fixture}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local out="$CHECK_DIR/${stem}.out" err="$CHECK_DIR/${stem}.err"
    local status
    status="$(capture_command "$out" "$err" "${args[@]}")"
    if [[ "$status" == 0 ]]; then
        record_peer_fixture "$fixture" "$tool" "pass" "exit 0"
        log_verify "  fixture PASS $fixture $tool"
    else
        record_peer_fixture "$fixture" "$tool" "fail" "exit $status"
        log_verify "  fixture FAIL $fixture $tool (exit $status)"
    fi
}

run_peer_fixture_probes() {
    local tool gzip_id gzip_path
    gzip_id="${GZIP_IDS[0]}"
    gzip_path="${DATA_GZ[$gzip_id]}"

    log_verify "--- peer fixture probes (descriptive; do not fail the run) ---"
    for tool in "${CHECK_PEER_ORDER[@]}"; do
        probe_peer_fixture iupac "$tool" se "$FIXTURE_DIR/iupac_valid.fastq"
        probe_peer_fixture empty "$tool" se "$FIXTURE_DIR/empty_valid.fastq"
        probe_peer_fixture crlf "$tool" se "$FIXTURE_DIR/crlf.fastq"
        probe_peer_fixture missing_newline "$tool" se "$FIXTURE_DIR/missing_final_newline.fastq"
        probe_peer_fixture slash_pair "$tool" paired "$CHECK_DIR/pair_r1.fastq" "$CHECK_DIR/pair_r2.fastq"
        probe_peer_fixture casava_pair "$tool" casava "$CHECK_DIR/casava_r1.fastq" "$CHECK_DIR/casava_r2.fastq"
        probe_peer_fixture odd_interleaved "$tool" interleaved "$CHECK_DIR/odd.fastq"
        probe_peer_fixture concat_gzip "$tool" gzip "$gzip_path"
    done
}

expect_ok_binary() {
    local binary="$1" label="$2"
    shift 2
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    [[ "$status" == 0 && ! -s "$out" && ! -s "$err" ]] ||
        check_fail "$label" "exit 0 and empty output" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
}

expect_code_binary() {
    local binary="$1" label="$2" code="$3" message="$4"
    shift 4
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    if [[ "$status" != 1 ]] ||
        ! rg -Fq -- "$code:" "$err" ||
        ! rg -Fq -- "$message" "$err"; then
        check_fail "$label" "exit 1, $code, $message" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
    fi
}

expect_pair_code() {
    local binary="$1" label="$2" code="$3"
    shift 3
    local out="$CHECK_DIR/contract.out" err="$CHECK_DIR/contract.err"
    local status
    status="$(capture_command "$out" "$err" "$binary" "$@")"
    if [[ "$status" != 1 ]] || ! rg -Fq -- "$code:" "$err"; then
        check_fail "$label" "exit 1 and $code" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
    fi
}

run_contract_tests() {
    : >"$VERIFY_LOG"
    log_verify "=== check contract ${TIMESTAMP} ==="

    local fixture fixture_gz bad code message binary
    fixture_gz="$CHECK_DIR/basic_valid.fastq.gz"
    gzip -c -- "$FIXTURE_DIR/basic_valid.fastq" >"$fixture_gz"

    log_verify "--- valid fixtures ---"
    for fixture in \
        acgtn_valid.fastq \
        basic_valid.fastq \
        crlf.fastq \
        missing_final_newline.fastq \
        empty_valid.fastq \
        iupac_valid.fastq; do
        for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
            expect_ok_binary "$binary" "$(basename "$binary") $fixture" \
                check "$FIXTURE_DIR/$fixture"
        done
    done
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_ok_binary "$binary" "$(basename "$binary") basic_valid.fastq.gz" \
            check "$fixture_gz"
    done

    log_verify "--- S001-S006 rejection contract ---"
    while IFS=$'\t' read -r bad code message; do
        [[ -n "$bad" ]] || continue
        for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
            expect_code_binary "$binary" "$(basename "$binary") $bad" "$code" "$message" \
                check "$FIXTURE_DIR/$bad"
        done
    done <<'EOF'
bad_plus.fastq	S001	plus line must start with '+'
bad_alphabet.fastq	S002	sequence byte is outside the selected alphabet
bad_header.fastq	S003	header line must start with '@' and contain a nonempty identifier
truncated_record.fastq	S004	unexpected end of file in quality line
bad_qual_length.fastq	S005	sequence and quality lengths differ
bad_quality_range.fastq	S006	quality byte must be ASCII 33 through 126
EOF

    log_verify "--- alphabet and precedence contract ---"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_ok_binary "$binary" "$(basename "$binary") explicit iupac" \
            check --alphabet iupac "$FIXTURE_DIR/iupac_valid.fastq"
        expect_code_binary "$binary" "$(basename "$binary") narrow alphabet" S002 \
            "sequence byte is outside the selected alphabet" \
            check --alphabet acgtn "$FIXTURE_DIR/iupac_valid.fastq"
    done

    log_verify "--- pair-name and pair-count contract ---"
    local pair_r1="$CHECK_DIR/pair_r1.fastq"
    local pair_r2="$CHECK_DIR/pair_r2.fastq"
    local odd="$CHECK_DIR/odd.fastq"
    local casava_r1="$CHECK_DIR/casava_r1.fastq"
    local casava_r2="$CHECK_DIR/casava_r2.fastq"
    printf '@cluster/1\nA\n+\n!\n' >"$pair_r1"
    printf '@cluster/2\nA\n+\n!\n' >"$pair_r2"
    printf '@cluster/1\nA\n+\n!\n' >"$odd"
    printf '@cluster 1:N:0:ATCACG\nA\n+\n!\n' >"$casava_r1"
    printf '@cluster 2:N:0:ATCACG\nA\n+\n!\n' >"$casava_r2"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_ok_binary "$binary" "$(basename "$binary") slash pair" \
            check --paired "$pair_r1" "$pair_r2"
        expect_pair_code "$binary" "$(basename "$binary") exact slash pair" P001 \
            check --paired --pair-names exact "$pair_r1" "$pair_r2"
        expect_ok_binary "$binary" "$(basename "$binary") casava pair" \
            check --paired "$casava_r1" "$casava_r2"
        expect_ok_binary "$binary" "$(basename "$binary") exact casava pair" \
            check --paired --pair-names exact "$casava_r1" "$casava_r2"
        printf '@other/2\nA\n+\n!\n' >"$pair_r2"
        expect_pair_code "$binary" "$(basename "$binary") name mismatch" P001 \
            check --paired "$pair_r1" "$pair_r2"
        : >"$pair_r2"
        expect_pair_code "$binary" "$(basename "$binary") count mismatch" P002 \
            check --paired "$pair_r1" "$pair_r2"
        expect_pair_code "$binary" "$(basename "$binary") odd interleaved" P002 \
            check --interleaved "$odd"
        printf '@cluster/2\nA\n+\n!\n' >"$pair_r2"
    done

    log_verify "--- concat gzip member traversal ---"
    local gzip_id gzip_path expected_records binary out err status got
    gzip_id="${GZIP_IDS[0]}"
    gzip_path="${DATA_GZ[$gzip_id]}"
    expected_records="${DATA_EXPECTED[$gzip_id]}"
    for binary in "$ZFASTQ" "$ZFASTQ_NATIVE"; do
        expect_ok_binary "$binary" "$(basename "$binary") ConcatGzip check" \
            check "$gzip_path"
        out="$CHECK_DIR/concat.count.out"
        err="$CHECK_DIR/concat.count.err"
        status="$(capture_command "$out" "$err" "$binary" count "$gzip_path")"
        got="$(parse_single_int "$(cat "$out")" || true)"
        [[ "$status" == 0 && "$got" == "$expected_records" ]] ||
            check_fail "$(basename "$binary") ConcatGzip count" \
                "exit 0 and $expected_records records" \
                "exit=$status count=${got:-?} stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
        log_verify "  $(basename "$binary") ConcatGzip count $got"
    done

    run_peer_fixture_probes

    log_verify "--- real workloads and external validator capability ---"
    prepare_real_lanes
    log_verify "ALL PASSED"
}

preflight_peer() {
    local section="$1" compression="$2" role="$3" key="$4" tool="$5"
    shift 5
    if ! tool_supports_role "$tool" "$role"; then
        LANE_REASON["$section|$key|$tool"]="not supported for role"
        log_verify "  skip $section $key $tool (role not supported)"
        return 0
    fi
    if ! bench_has_tool "$tool"; then
        LANE_REASON["$section|$key|$tool"]="tool unavailable"
        log_verify "  skip $section $key $tool (tool unavailable)"
        return 0
    fi
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "$@"
    build_check_args args "$tool" "$role" "${paths[@]}"
    local stem="${section}_${key}_${tool}"
    stem="${stem//[^A-Za-z0-9_.+-]/_}"
    local out="$CHECK_DIR/${stem}.out" err="$CHECK_DIR/${stem}.err"
    local status
    status="$(capture_command "$out" "$err" "${args[@]}")"
    if [[ "$status" == 0 ]]; then
        if [[ "$role" == gzip ]]; then
            LANE_REASON["$section|$key|$tool"]="probed, not timed: member traversal unproven"
            log_verify "  peer PROBE $section $key $tool (not timed)"
        else
            LANE_ENABLED["$section|$key|$tool"]=1
            LANE_REASON["$section|$key|$tool"]="accepted positive workload"
            log_verify "  peer PASS $section $key $tool"
        fi
    else
        LANE_REASON["$section|$key|$tool"]="exit $status on positive workload"
        log_verify "  peer SKIP $section $key $tool (exit $status)"
    fi
}

prepare_workload_lanes() {
    local section="$1" compression="$2" role="$3"
    shift 3
    local -a ids=("$@")
    local key
    key="$(workload_key "$role" "${ids[@]}")"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"
    build_check_args args z-fastq "$role" "${paths[@]}"
    local out="$CHECK_DIR/zfastq_${section}_${key}.out"
    local err="$CHECK_DIR/zfastq_${section}_${key}.err"
    local status
    status="$(capture_command "$out" "$err" "${args[@]}")"
    [[ "$status" == 0 ]] ||
        check_fail "$section $key z-fastq" "exit 0" \
            "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"

    if [[ "$compression" == gzip ]]; then
        build_check_args args z-fastq-native "$role" "${paths[@]}"
        out="$CHECK_DIR/native_${section}_${key}.out"
        err="$CHECK_DIR/native_${section}_${key}.err"
        status="$(capture_command "$out" "$err" "${args[@]}")"
        [[ "$status" == 0 ]] ||
            check_fail "$section $key z-fastq-native" "exit 0" \
                "exit=$status stdout=$(summarize_output "$out") stderr=$(summarize_output "$err")"
    fi

    local tool
    for tool in "${CHECK_PEER_ORDER[@]}"; do
        preflight_peer "$section" "$compression" "$role" "$key" "$tool" "${ids[@]}"
    done
}

prepare_real_lanes() {
    local id
    for id in "${SE_IDS[@]}"; do
        prepare_workload_lanes perf_plain plain se "$id"
        prepare_workload_lanes perf_gzip gzip se "$id"
    done
    prepare_workload_lanes perf_plain plain paired "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip paired "${PAIRED_IDS[@]}"
    prepare_workload_lanes perf_plain plain casava "${CASAVA_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip casava "${CASAVA_IDS[@]}"
    prepare_workload_lanes perf_plain plain slash "${SLASH_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip slash "${SLASH_IDS[@]}"
    prepare_workload_lanes perf_plain plain interleaved "${INTERLEAVED_IDS[@]}"
    prepare_workload_lanes perf_gzip gzip interleaved "${INTERLEAVED_IDS[@]}"
    # ConcatGzip is intentionally a gzip-only robustness workload.
    prepare_workload_lanes perf_gzip gzip gzip "${GZIP_IDS[@]}"
}

check_add_command() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    zebrac_add_command "check" "$section" "$workload" "$tool" "$family" \
        "$input_bytes" "$decoded_bytes" "$json_out" "$command"
}

run_zebrac_tool() {
    local section="$1" workload="$2" tool="$3" family="$4"
    local json_out="$5" command="$6" input_bytes="$7" decoded_bytes="$8"
    local t0 t1 elapsed
    echo "  >> $section $workload $tool  (runs=$RUNS warmup=$WARMUP)"
    t0="$(date +%s)"
    zebrac_clear_commands
    check_add_command "$section" "$workload" "$tool" "$family" \
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
    local section="$1" compression="$2" role="$3"
    shift 3
    local -a ids=("$@")
    local key
    key="$(workload_key "$role" "${ids[@]}")"
    local input_bytes decoded_bytes json command tool
    input_bytes="$(workload_input_bytes "$compression" "${ids[@]}")"
    decoded_bytes="$(workload_decoded_bytes "${ids[@]}")"
    local -a args=() paths=()
    resolve_dataset_paths paths "$compression" "${ids[@]}"

    build_check_args args z-fastq "$role" "${paths[@]}"
    command="$(zebrac_command "${args[@]}")"
    json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq.json"
    run_zebrac_tool "$section" "$key" z-fastq z-fastq "$json" \
        "$command" "$input_bytes" "$decoded_bytes"

    if [[ "$compression" == gzip ]]; then
        build_check_args args z-fastq-native "$role" "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__z-fastq-native.json"
        run_zebrac_tool "$section" "$key" z-fastq-native z-fastq "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    fi

    for tool in "${CHECK_PEER_ORDER[@]}"; do
        [[ "${LANE_ENABLED[$section|$key|$tool]:-0}" == 1 ]] || continue
        build_check_args args "$tool" "$role" "${paths[@]}"
        command="$(zebrac_command "${args[@]}")"
        json="$RESULTS_DIR/${section}_${TIMESTAMP}/${key}__${tool}.json"
        run_zebrac_tool "$section" "$key" "$tool" "$tool" "$json" \
            "$command" "$input_bytes" "$decoded_bytes"
    done
}

run_perf() {
    mkdir -p "$RESULTS_DIR/perf_plain_${TIMESTAMP}" "$RESULTS_DIR/perf_gzip_${TIMESTAMP}"
    echo "=== perf_plain ==="
    local id
    for id in "${SE_IDS[@]}"; do
        run_timed_workload perf_plain plain se "$id"
    done
    run_timed_workload perf_plain plain paired "${PAIRED_IDS[@]}"
    run_timed_workload perf_plain plain casava "${CASAVA_IDS[@]}"
    run_timed_workload perf_plain plain slash "${SLASH_IDS[@]}"
    run_timed_workload perf_plain plain interleaved "${INTERLEAVED_IDS[@]}"

    echo "=== perf_gzip ==="
    for id in "${SE_IDS[@]}"; do
        run_timed_workload perf_gzip gzip se "$id"
    done
    run_timed_workload perf_gzip gzip paired "${PAIRED_IDS[@]}"
    run_timed_workload perf_gzip gzip casava "${CASAVA_IDS[@]}"
    run_timed_workload perf_gzip gzip slash "${SLASH_IDS[@]}"
    run_timed_workload perf_gzip gzip interleaved "${INTERLEAVED_IDS[@]}"
    run_timed_workload perf_gzip gzip gzip "${GZIP_IDS[@]}"
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
        printf '  "schema_version": "check-run.v1",\n'
        printf '  "timestamp": %s,\n' "$(zebrac_json_string "$TIMESTAMP")"
        printf '  "runner": "zebrac",\n'
        printf '  "mode": "warm",\n'
        printf '  "suite": "check",\n'
        printf '  "real_set": %s,\n' "$(zebrac_json_string "$CHECK_SET")"
        printf '  "workloads": %s,\n' "$(workloads_json)"
        printf '  "zebrac": %s,\n' "$(zebrac_json_string "${CHECK_ZEBRAC_VER}")"
        printf '  "z_fastq": %s,\n' "$(zebrac_json_string "${CHECK_ZFASTQ_VER}")"
        printf '  "z_fastq_native": %s,\n' "$(zebrac_json_string "${CHECK_ZFASTQ_NATIVE_VER}")"
        printf '  "z_fastq_bytes": %s,\n' "$(zebrac_json_number_or_null "${CHECK_ZFASTQ_BYTES}")"
        printf '  "z_fastq_native_bytes": %s,\n' "$(zebrac_json_number_or_null "${CHECK_ZFASTQ_NATIVE_BYTES}")"
        printf '  "runs": %s,\n' "$(zebrac_json_number_or_null "$RUNS")"
        printf '  "warmup": %s,\n' "$(zebrac_json_number_or_null "$WARMUP")"
        printf '  "duration_ms": %s,\n' "$(zebrac_json_number_or_null "$ZEBRAC_DURATION_MS")"
        printf '  "metadata": %s,\n' "$(zebrac_json_string "metadata_${TIMESTAMP}.jsonl")"
        printf '  "verify_log": %s,\n' "$(zebrac_json_string "verify_${TIMESTAMP}.log")"
        printf '  "verify_skipped": %s,\n' "$verify_skipped_json"
        printf '  "verify_pass": %s,\n' "$verify_pass_json"
        printf '  "tools": %s,\n' "$CHECK_TOOLS_JSON"
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

# --- run ---
echo "z-fastq check bench  $TIMESTAMP"
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
    bench_require_tool fq
    bench_require_tool fastqvalidator
    bench_require_tool seqfu
    bench_require_tool fqtools
    run_contract_tests
else
    : >"$VERIFY_LOG"
    log_verify "WARNING: skipping the check contract preflight (--skip-tests); not for a published report"
    if $DO_BENCHMARKS; then
        prepare_real_lanes
    fi
fi

if $DO_BENCHMARKS; then
    bench_require_tool zebrac
    bench_require_tool fq
    bench_require_tool fastqvalidator
    bench_require_tool seqfu
    bench_require_tool fqtools
    if ! $DO_TESTS; then
        # prepare_real_lanes already ran above so the skip path remains explicit.
        :
    fi
    if $DO_FULL; then
        run_perf
    fi
fi

CHECK_WORKLOADS_JSON="$(workloads_json)"
CHECK_ZEBRAC_VER="$(bench_tool_version zebrac || true)"
CHECK_ZFASTQ_VER="$(bench_tool_version z-fastq || true)"
CHECK_ZFASTQ_NATIVE_VER="$(bench_tool_version z-fastq-native || true)"
CHECK_ZFASTQ_BYTES="$(file_size_bytes "$ZFASTQ")"
CHECK_ZFASTQ_NATIVE_BYTES="$(file_size_bytes "$ZFASTQ_NATIVE")"
CHECK_FQ_VER="$(bench_tool_version fq || true)"
CHECK_FASTQVALIDATOR_VER="$(bench_tool_version fastqvalidator || true)"
CHECK_SEQFU_VER="$(bench_tool_version seqfu || true)"
CHECK_FQTOOLS_VER="$(bench_tool_version fqtools || true)"
CHECK_TOOLS_JSON="$(
    printf '{'
    printf '"fq":%s,' "$(zebrac_json_string "$CHECK_FQ_VER")"
    printf '"fastqvalidator":%s,' "$(zebrac_json_string "$CHECK_FASTQVALIDATOR_VER")"
    printf '"seqfu":%s,' "$(zebrac_json_string "$CHECK_SEQFU_VER")"
    printf '"fqtools":%s' "$(zebrac_json_string "$CHECK_FQTOOLS_VER")"
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

#!/usr/bin/env bash
# Path resolution and zebrac helpers for bench/*/run.sh.
#
# Source this file. It does not install tools: peers, adapters, zebrac, and
# report Python live under tools/.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "tools.sh is a library; source it from a benchmark script." >&2
    exit 1
fi

ulimit -S -c 0 2>/dev/null || true

BENCH_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$BENCH_SHARED_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_ROOT/.." && pwd)"
TOOLS_DIR="$PROJECT_ROOT/tools"
TOOLS_BIN_DIR="$TOOLS_DIR/bin"
TOOLS_VENV_DIR="$TOOLS_DIR/venv"

# shellcheck disable=SC1091
source "$BENCH_SHARED_DIR/catalog.sh"

# shellcheck disable=SC1091
source "$TOOLS_DIR/versions.sh"

ZFASTQ="${ZFASTQ:-$PROJECT_ROOT/zig-out/bin/z-fastq}"
ZFASTQ_NATIVE="${ZFASTQ_NATIVE:-$PROJECT_ROOT/zig-out/bin/z-fastq-native}"
ZEBRAC="${ZEBRAC:-$TOOLS_DIR/zebrac}"
SEQTK="${SEQTK:-$TOOLS_BIN_DIR/seqtk}"
FQTOOLS="${FQTOOLS:-$TOOLS_BIN_DIR/fqtools}"
NEEDLETAIL="${NEEDLETAIL:-$TOOLS_BIN_DIR/needletail-adapter}"
HELICASE="${HELICASE:-$TOOLS_BIN_DIR/helicase-adapter}"
SEQFU="${SEQFU:-$TOOLS_BIN_DIR/seqfu}"
SEQKIT="${SEQKIT:-$TOOLS_BIN_DIR/seqkit}"
FQ="${FQ:-$TOOLS_BIN_DIR/fq}"
FASTQVALIDATOR="${FASTQVALIDATOR:-$TOOLS_BIN_DIR/fastQValidator}"

ZEBRAC_DURATION_MS="${ZEBRAC_DURATION_MS:-5000}"
ZEBRAC_MIN_SAMPLES="${ZEBRAC_MIN_SAMPLES:-25}"
ZEBRAC_MAX_SAMPLES="${ZEBRAC_MAX_SAMPLES:-}"
ZEBRAC_WARMUP="${ZEBRAC_WARMUP:-5}"
ZEBRAC_ALLOW_FAILURES="${ZEBRAC_ALLOW_FAILURES:-false}"

declare -a ZEBRAC_BENCH_COMMANDS=()
declare -a ZEBRAC_BENCH_METADATA=()

bench_tool_path() {
    local name="$1"
    case "$name" in
        z-fastq) echo "$ZFASTQ" ;;
        z-fastq-native) echo "$ZFASTQ_NATIVE" ;;
        zebrac) echo "$ZEBRAC" ;;
        seqtk) echo "$SEQTK" ;;
        fqtools) echo "$FQTOOLS" ;;
        needletail) echo "$NEEDLETAIL" ;;
        helicase) echo "$HELICASE" ;;
        seqfu) echo "$SEQFU" ;;
        seqkit) echo "$SEQKIT" ;;
        fq) echo "$FQ" ;;
        fastqvalidator) echo "$FASTQVALIDATOR" ;;
        *) return 1 ;;
    esac
}

bench_has_tool() {
    local path
    path="$(bench_tool_path "$1")" || return 1
    [[ -n "$path" && -x "$path" ]]
}

bench_require_tool() {
    local name="$1"
    if ! bench_has_tool "$name"; then
        echo "error: required benchmark tool not found or not executable: $name" >&2
        echo "       resolved path: $(bench_tool_path "$name" 2>/dev/null || echo '<unknown>')" >&2
        return 1
    fi
}

bench_tool_version() {
    local name="$1"
    local path
    path="$(bench_tool_path "$name")" || return 1
    case "$name" in
        z-fastq|z-fastq-native|zebrac)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        seqtk)
            [[ -x "$path" ]] || return 1
            { "$path" 2>&1 || true; } | awk '/^Version:/{print "seqtk " $2; exit}'
            ;;
        fqtools)
            echo "fqtools $FQTOOLS_VERSION (HTSlib $HTSLIB_VERSION)"
            ;;
        needletail|helicase)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        seqfu)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        seqkit)
            [[ -x "$path" ]] && "$path" version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        fq)
            [[ -x "$path" ]] && "$path" --version 2>&1 | awk 'NR==1{print; exit}'
            ;;
        fastqvalidator)
            [[ -x "$path" ]] && printf 'fastQValidator %s\n' "$FASTQ_VALIDATOR_VERSION"
            ;;
        *)
            return 1
            ;;
    esac
}

bench_catalog() {
    catalog_main "$@"
}

report_python() {
    if [[ -x "$TOOLS_VENV_DIR/bin/python" ]]; then
        echo "$TOOLS_VENV_DIR/bin/python"
        return 0
    fi
    echo "error: tools/venv/bin/python is missing; run tools/install.sh venv" >&2
    return 1
}

file_size_bytes() {
    local path="$1"
    stat --printf='%s' "$path" 2>/dev/null || stat -f '%z' "$path" 2>/dev/null || echo 0
}

zebrac_json_escape() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    value="${value//$'\n'/\\n}"
    value="${value//$'\r'/\\r}"
    value="${value//$'\t'/\\t}"
    printf '%s' "$value"
}

zebrac_json_string() {
    printf '"%s"' "$(zebrac_json_escape "$1")"
}

zebrac_json_number_or_null() {
    local value="$1"
    if [[ -z "$value" ]]; then
        printf 'null'
    else
        printf '%s' "$value"
    fi
}

zebrac_clear_commands() {
    ZEBRAC_BENCH_COMMANDS=()
    ZEBRAC_BENCH_METADATA=()
}

zebrac_add_command() {
    if [[ $# -ne 9 ]]; then
        echo "usage: zebrac_add_command <suite> <section> <workload> <tool> <tool_family> <input_bytes> <decoded_bytes> <raw_json> <command>" >&2
        return 1
    fi

    local suite="$1"
    local section="$2"
    local workload="$3"
    local tool="$4"
    local tool_family="$5"
    local input_bytes="$6"
    local decoded_bytes="$7"
    local raw_json="$8"
    local command="$9"

    ZEBRAC_BENCH_COMMANDS+=("$command")
    ZEBRAC_BENCH_METADATA+=("$(
        printf '{'
        printf '"schema_version":"bench.meta.v1",'
        printf '"raw_json":%s,' "$(zebrac_json_string "$raw_json")"
        printf '"suite":%s,' "$(zebrac_json_string "$suite")"
        printf '"section":%s,' "$(zebrac_json_string "$section")"
        printf '"workload":%s,' "$(zebrac_json_string "$workload")"
        printf '"tool":%s,' "$(zebrac_json_string "$tool")"
        printf '"tool_family":%s,' "$(zebrac_json_string "$tool_family")"
        printf '"command":%s,' "$(zebrac_json_string "$command")"
        printf '"input_bytes":%s,' "$(zebrac_json_number_or_null "$input_bytes")"
        printf '"decoded_bytes":%s' "$(zebrac_json_number_or_null "$decoded_bytes")"
        printf '}'
    )")
}

zebrac_write_metadata() {
    local metadata_jsonl="$1"
    mkdir -p "$(dirname "$metadata_jsonl")"
    for row in "${ZEBRAC_BENCH_METADATA[@]}"; do
        printf '%s\n' "$row" >> "$metadata_jsonl"
    done
}

zebrac_run_current_group() {
    if [[ $# -ne 2 ]]; then
        echo "usage: zebrac_run_current_group <raw_json> <metadata_jsonl>" >&2
        return 1
    fi

    local raw_json="$1"
    local metadata_jsonl="$2"

    if [[ "${#ZEBRAC_BENCH_COMMANDS[@]}" -eq 0 ]]; then
        echo "error: no zebrac commands queued" >&2
        return 1
    fi

    bench_require_tool zebrac
    mkdir -p "$(dirname "$raw_json")"

    local args=(
        "$ZEBRAC"
        --quiet
        --duration "$ZEBRAC_DURATION_MS"
        --min-samples "$ZEBRAC_MIN_SAMPLES"
        --warmup "$ZEBRAC_WARMUP"
    )

    if [[ -n "$ZEBRAC_MAX_SAMPLES" ]]; then
        args+=(--max-samples "$ZEBRAC_MAX_SAMPLES")
    fi
    if [[ "$ZEBRAC_ALLOW_FAILURES" == "true" ]]; then
        args+=(--allow-failures)
    fi

    args+=(--json "$raw_json" --)
    args+=("${ZEBRAC_BENCH_COMMANDS[@]}")

    "${args[@]}"
    zebrac_write_metadata "$metadata_jsonl"
}

bench_group() {
    local json_out="$1"
    # shellcheck disable=SC2153
    zebrac_run_current_group "$json_out" "$METADATA_JSONL"
    zebrac_clear_commands
}

quote_arg() {
    printf '%q' "$1"
}

# Zebrac splits this string itself. It does not start a shell, and pipes or
# redirects stay literal. Wrapping the tool in `bash -c` would report bash RSS
# (~4 MB) for every lane and add bash spawn to wall time.
zebrac_quote() {
    local arg="$1"
    case "$arg" in
        *\'*)
            if [[ "$arg" == *'"'* ]]; then
                echo "error: cannot quote argument for zebrac: $arg" >&2
                return 1
            fi
            printf '"%s"' "$arg"
            ;;
        *[!A-Za-z0-9._+/=@:-]*)
            printf "'%s'" "$arg"
            ;;
        *)
            printf '%s' "$arg"
            ;;
    esac
}

zebrac_command() {
    local out="" sep="" arg
    for arg in "$@"; do
        out+="${sep}$(zebrac_quote "$arg")"
        sep=" "
    done
    printf '%s' "$out"
}

#!/usr/bin/env bash
# Prepare repository-local benchmark commands and parser adapters.

set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$TOOLS_DIR/.." && pwd)"
BIN_DIR="$TOOLS_DIR/bin"
VENV_DIR="$TOOLS_DIR/venv"
LOCAL_DIR="$TOOLS_DIR/.local"
INSTALLS_DIR="$LOCAL_DIR/installs"
WRAPPERS_DIR="$TOOLS_DIR/wrappers"
WRAPPER_INSTALLS_DIR="$INSTALLS_DIR/wrappers"
WRAPPER_CURRENT="$LOCAL_DIR/wrappers-current"
TOOL_JOBS="${TOOL_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}"
KEEP_TOOL_WORK="${KEEP_TOOL_WORK:-0}"
ACTIVE_WORK=""
ACTIVE_STAGE=""
BUILD_DESCRIPTION=""
PEER_NAMES=(
    seqtk
    fqtools
    fq
    fqkit
    seqkit
    rasusa
    fasten
    seqfu
    irma-core
    fastqvalidator
    fastp
    bbtools
)

# shellcheck source=tools/versions.sh
source "$TOOLS_DIR/versions.sh"

usage() {
    cat <<'EOF'
Usage:
  tools/install.sh                     create or check tools/venv
  tools/install.sh <name>              install or reuse one target
  tools/install.sh peers               install every enabled command-line peer
  tools/install.sh all                 install venv, wrappers, and peers
  tools/install.sh --check [name]      verify one target without changing it
  tools/install.sh --list              show targets and local state
  tools/install.sh --help              show this help

Names: venv, wrappers, seqtk, fqtools, fq, fqkit, seqkit, rasusa,
       fasten, seqfu, irma-core, fastqvalidator, fastp, bbtools, peers, all

The installer uses commands already on PATH. It does not use the repository
root .venv, install system packages, edit PATH, or activate an environment.
Downloads, package caches, source, and build files stay under /tmp. Only the
checked runtime is published below tools/.local and linked from tools/bin.
EOF
}

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "error: required command not found: $1" >&2
        return 1
    fi
}

validate_settings() {
    if [[ ! "$TOOL_JOBS" =~ ^[1-9][0-9]*$ ]]; then
        echo "error: TOOL_JOBS must be a positive integer" >&2
        return 1
    fi
    case "$KEEP_TOOL_WORK" in
        0|1) ;;
        *)
            echo "error: KEEP_TOOL_WORK must be 0 or 1" >&2
            return 1
            ;;
    esac
}

cleanup() {
    if [[ -n "$ACTIVE_STAGE" && -d "$ACTIVE_STAGE" ]]; then
        case "$ACTIVE_STAGE" in
            "$LOCAL_DIR"/stage/*)
                if [[ -f "$ACTIVE_STAGE/.z-fastq-tool" ]]; then
                    rm -rf -- "$ACTIVE_STAGE"
                fi
                ;;
        esac
    fi

    if [[ -n "$ACTIVE_WORK" && -d "$ACTIVE_WORK" ]]; then
        case "$ACTIVE_WORK" in
            /tmp/z-fastq-tools.*)
                if [[ "$KEEP_TOOL_WORK" == 1 ]]; then
                    echo "keep: $ACTIVE_WORK"
                else
                    rm -rf -- "$ACTIVE_WORK"
                fi
                ;;
        esac
    fi
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

finish_work() {
    local work="$ACTIVE_WORK"
    ACTIVE_WORK=""
    [[ -n "$work" && -d "$work" ]] || return

    case "$work" in
        /tmp/z-fastq-tools.*)
            if [[ "$KEEP_TOOL_WORK" == 1 ]]; then
                echo "keep: $work"
            else
                rm -rf -- "$work"
            fi
            ;;
        *)
            echo "warning: keep unexpected work path: $work" >&2
            ;;
    esac
}

check_platform() {
    if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
        echo "error: benchmark peers currently support Linux x86-64 only" >&2
        return 1
    fi
}

acquire_install_lock() {
    require_command flock
    mkdir -p "$LOCAL_DIR"
    exec 9>"$LOCAL_DIR/install.lock"
    flock 9
}

switch_relative_link() {
    local link="$1"
    local target="$2"
    local temporary="$link.tmp.$$"

    if [[ -e "$link" && ! -L "$link" ]]; then
        echo "error: installer will not replace non-link path: $link" >&2
        return 1
    fi
    if [[ -e "$temporary" || -L "$temporary" ]]; then
        echo "error: temporary link already exists: $temporary" >&2
        return 1
    fi
    ln -s "$target" "$temporary"
    mv -Tf -- "$temporary" "$link"
}

receipt_has() {
    local receipt="$1"
    local field="$2"
    local value="$3"
    grep -Fqx "$field"$'\t'"$value" "$receipt"
}

download() {
    local url="$1"
    local output="$2"
    curl --fail --location --retry 3 --show-error --silent "$url" --output "$output"
}

find_binary() {
    local root="$1"
    local name="$2"
    local path
    path="$(find "$root" -type f -name "$name" -perm -u+x -print -quit)"
    if [[ -z "$path" ]]; then
        echo "error: downloaded archive does not contain executable $name" >&2
        return 1
    fi
    printf '%s\n' "$path"
}

check_venv() {
    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        echo "error: tools/venv is missing; run tools/install.sh" >&2
        return 1
    fi

    if ! "$VENV_DIR/bin/python" -I -c '
import pathlib
import sys

expected = pathlib.Path(sys.argv[1]).resolve()
actual = pathlib.Path(sys.prefix).resolve()
raise SystemExit(0 if actual == expected and sys.prefix != sys.base_prefix else 1)
' "$VENV_DIR"; then
        echo "error: tools/venv does not contain its expected Python environment" >&2
        return 1
    fi

    echo "ok: tools/venv ($("$VENV_DIR/bin/python" --version 2>&1))"
}

install_venv() {
    if [[ -e "$VENV_DIR" ]]; then
        check_venv
        return
    fi

    require_command uv
    require_command python3

    local python_path
    python_path="$(command -v python3)"

    echo "create: tools/venv"
    if ! uv venv \
        --no-project \
        --no-config \
        --no-python-downloads \
        --no-cache \
        --python "$python_path" \
        "$VENV_DIR"; then
        rm -rf -- "$VENV_DIR"
        return 1
    fi
    if ! check_venv; then
        rm -rf -- "$VENV_DIR"
        return 1
    fi
}

rust_host() {
    rustc -vV | sed -n 's/^host: //p'
}

cpu_name() {
    local value=""
    if [[ -r /proc/cpuinfo ]]; then
        value="$(sed -n -E '/^(model name|Hardware|Processor)[[:space:]]*:/{s/^[^:]*:[[:space:]]*//;p;q;}' /proc/cpuinfo)"
    fi
    if [[ -z "$value" ]]; then
        value="$(uname -m)"
    fi
    printf '%s\n' "$value"
}

check_helicase_host() {
    local architecture
    local features=""
    architecture="$(uname -m)"
    if [[ -r /proc/cpuinfo ]]; then
        features="$(sed -n -E '/^(flags|Features)[[:space:]]*:/{s/^[^:]*:[[:space:]]*//;p;q;}' /proc/cpuinfo)"
    fi

    case "$architecture" in
        x86_64)
            if [[ " $features " != *' avx2 '* && " $features " != *' ssse3 '* ]]; then
                echo "error: Helicase needs AVX2 or SSSE3 for this engine comparison" >&2
                return 1
            fi
            ;;
        aarch64)
            if [[ " $features " != *' asimd '* && " $features " != *' neon '* ]]; then
                echo "error: Helicase needs NEON for this engine comparison" >&2
                return 1
            fi
            ;;
        *)
            echo "error: unsupported Helicase benchmark host: $architecture" >&2
            return 1
            ;;
    esac
}

wrapper_receipt() {
    if [[ ! -L "$WRAPPER_CURRENT" ]]; then
        return 1
    fi
    printf '%s/receipt.tsv\n' "$WRAPPER_CURRENT"
}

wrappers_current() {
    local receipt
    receipt="$(wrapper_receipt)" || return 1
    [[ -f "$receipt" ]] || return 1
    [[ -x "$BIN_DIR/needletail-adapter" ]] || return 1
    [[ -x "$BIN_DIR/helicase-adapter" ]] || return 1

    check_helicase_host >/dev/null 2>&1 || return 1
    receipt_has "$receipt" recipe "$WRAPPERS_RECIPE" || return 1
    receipt_has "$receipt" needletail "$NEEDLETAIL_VERSION" || return 1
    receipt_has "$receipt" helicase "$HELICASE_VERSION" || return 1
    receipt_has "$receipt" cpu "$(cpu_name)" || return 1
    receipt_has "$receipt" target_cpu native || return 1

    if find "$WRAPPERS_DIR" -type f \
        \( -name Cargo.toml -o -name Cargo.lock -o -name '*.rs' \) \
        -newer "$receipt" -print -quit | grep -q .; then
        return 1
    fi
}

check_version_line() {
    local command="$1"
    local engine="$2"
    local version="$3"
    local output
    output="$("$command" --version)"
    [[ "$output" == *"engine=$engine $version;"* ]] || {
        echo "error: unexpected $engine adapter version: $output" >&2
        return 1
    }
    [[ "$output" == *"domain=screened-four-line-lf-fastq;"* ]] || {
        echo "error: $engine adapter does not report its input domain" >&2
        return 1
    }
    [[ "$output" == *"target-cpu=native; threads=1" ]] || {
        echo "error: $engine adapter does not report the expected build mode" >&2
        return 1
    }
}

check_stripped_binary() {
    local binary="$1"
    local output
    require_command file
    output="$(file -L "$binary")"
    if [[ "$output" == *'not stripped'* || "$output" != *stripped* ]]; then
        echo "error: retained native executable is not stripped: $binary" >&2
        return 1
    fi
}

check_wrapper_smoke() {
    local needletail="$1"
    local helicase="$2"
    local fixture="$ROOT_DIR/tests/data/synthetic/basic_valid.fastq"
    local smoke_dir
    smoke_dir="$(mktemp -d /tmp/z-fastq-wrapper-smoke.XXXXXX)"

    (
        trap 'rm -rf -- "$smoke_dir"' EXIT
        local gzip_fixture="$smoke_dir/basic_valid.fastq.gz"
        local needletail_stats="$smoke_dir/needletail.stats"
        local helicase_stats="$smoke_dir/helicase.stats"
        local needletail_gzip_stats="$smoke_dir/needletail-gzip.stats"
        local helicase_gzip_stats="$smoke_dir/helicase-gzip.stats"

        gzip -c "$fixture" >"$gzip_fixture"

        [[ "$("$needletail" count "$fixture")" == 5 ]]
        [[ "$("$helicase" count "$fixture")" == 5 ]]
        [[ "$("$needletail" count "$gzip_fixture")" == 5 ]]
        [[ "$("$helicase" count "$gzip_fixture")" == 5 ]]

        "$needletail" stats "$fixture" >"$needletail_stats"
        "$helicase" stats "$fixture" >"$helicase_stats"
        "$needletail" stats "$gzip_fixture" >"$needletail_gzip_stats"
        "$helicase" stats "$gzip_fixture" >"$helicase_gzip_stats"

        sed "1s|.*|input: $fixture|" "$needletail_gzip_stats" | cmp -s "$needletail_stats" -
        sed "1s|.*|input: $fixture|" "$helicase_gzip_stats" | cmp -s "$helicase_stats" -
        cmp -s "$needletail_stats" "$helicase_stats"
        grep -Fqx 'reads: 5' "$needletail_stats"
        grep -Fqx 'bases: 21' "$needletail_stats"
        grep -Fqx 'quality_sum: 644' "$needletail_stats"
        grep -Fqx 'q20_bases: 16' "$needletail_stats"
        grep -Fqx 'q30_bases: 16' "$needletail_stats"
    ) || {
        echo "error: wrapper smoke check failed" >&2
        return 1
    }
}

check_wrappers() {
    require_command gzip
    if ! wrappers_current; then
        echo "error: parser adapters are missing or stale; run tools/install.sh wrappers" >&2
        return 1
    fi

    check_stripped_binary "$BIN_DIR/needletail-adapter"
    check_stripped_binary "$BIN_DIR/helicase-adapter"
    check_version_line "$BIN_DIR/needletail-adapter" needletail "$NEEDLETAIL_VERSION" || return
    check_version_line "$BIN_DIR/helicase-adapter" helicase "$HELICASE_VERSION" || return
    check_wrapper_smoke "$BIN_DIR/needletail-adapter" "$BIN_DIR/helicase-adapter" || return
    echo "ok: wrappers (Needletail $NEEDLETAIL_VERSION, Helicase $HELICASE_VERSION)"
}

remove_stale_stages() {
    local name="$1"
    local candidate
    while IFS= read -r -d '' candidate; do
        if [[ -f "$candidate/.z-fastq-tool" ]] &&
            [[ "$(<"$candidate/.z-fastq-tool")" == "$name" ]]; then
            rm -rf -- "$candidate"
        else
            echo "warning: keep unmarked stage path: $candidate" >&2
        fi
    done < <(find "$LOCAL_DIR/stage" -mindepth 1 -maxdepth 1 -type d -name "$name.*" -print0)
}

write_wrapper_receipt() {
    local receipt="$1"
    local git_commit="unknown"
    local git_dirty="unknown"

    if git -C "$ROOT_DIR" rev-parse --verify HEAD >/dev/null 2>&1; then
        git_commit="$(git -C "$ROOT_DIR" rev-parse HEAD)"
        if [[ -n "$(git -C "$ROOT_DIR" status --porcelain -- tools/wrappers tools/install.sh tools/versions.sh)" ]]; then
            git_dirty="yes"
        else
            git_dirty="no"
        fi
    fi

    {
        printf 'recipe\t%s\n' "$WRAPPERS_RECIPE"
        printf 'needletail\t%s\n' "$NEEDLETAIL_VERSION"
        printf 'helicase\t%s\n' "$HELICASE_VERSION"
        printf 'rustc\t%s\n' "$(rustc --version)"
        printf 'rust_host\t%s\n' "$(rust_host)"
        printf 'cpu\t%s\n' "$(cpu_name)"
        printf 'target_cpu\tnative\n'
        printf 'features\tnative-cpu\n'
        printf 'compression\tplain,gzip\n'
        printf 'runtime\tsystem C and compiler runtime\n'
        printf 'jobs\t%s\n' "$TOOL_JOBS"
        printf 'build\tcargo build --locked --release --bins --features native-cpu\n'
        printf 'git_commit\t%s\n' "$git_commit"
        printf 'git_dirty\t%s\n' "$git_dirty"
    } >"$receipt"
}

remove_old_installs() {
    local name="$1"
    local current="$2"
    local directory="$INSTALLS_DIR/$name"
    local candidate
    while IFS= read -r -d '' candidate; do
        if [[ "$candidate" == "$current" ]]; then
            continue
        fi
        if [[ -f "$candidate/.z-fastq-tool" ]] &&
            [[ "$(<"$candidate/.z-fastq-tool")" == "$name" ]]; then
            rm -rf -- "$candidate"
        else
            echo "warning: keep unmarked install path: $candidate" >&2
        fi
    done < <(find "$directory" -mindepth 1 -maxdepth 1 -type d -print0)
}

install_wrappers() {
    validate_settings
    require_command gzip
    acquire_install_lock
    check_helicase_host

    if wrappers_current && check_wrappers; then
        echo "reuse: wrappers"
        return
    fi

    require_command cargo
    require_command rustc
    require_command git
    require_command install

    mkdir -p "$LOCAL_DIR/stage" "$WRAPPER_INSTALLS_DIR" "$BIN_DIR"
    remove_stale_stages wrappers
    ACTIVE_WORK="$(mktemp -d /tmp/z-fastq-tools.XXXXXX)"
    ACTIVE_STAGE="$(mktemp -d "$LOCAL_DIR/stage/wrappers.XXXXXX")"
    printf 'wrappers\n' >"$ACTIVE_STAGE/.z-fastq-tool"
    mkdir -p "$ACTIVE_STAGE/bin"

    echo "build: wrappers (Needletail $NEEDLETAIL_VERSION, Helicase $HELICASE_VERSION)"
    env \
        CARGO_HOME="$ACTIVE_WORK/cargo-home" \
        CARGO_TARGET_DIR="$ACTIVE_WORK/cargo-target" \
        RUSTFLAGS='-C target-cpu=native' \
        cargo build \
            --locked \
            --release \
            --bins \
            --jobs "$TOOL_JOBS" \
            --manifest-path "$WRAPPERS_DIR/Cargo.toml" \
            --features native-cpu

    install -m 755 \
        "$ACTIVE_WORK/cargo-target/release/needletail-adapter" \
        "$ACTIVE_STAGE/bin/needletail-adapter"
    install -m 755 \
        "$ACTIVE_WORK/cargo-target/release/helicase-adapter" \
        "$ACTIVE_STAGE/bin/helicase-adapter"
    write_wrapper_receipt "$ACTIVE_STAGE/receipt.tsv"

    check_version_line "$ACTIVE_STAGE/bin/needletail-adapter" needletail "$NEEDLETAIL_VERSION"
    check_version_line "$ACTIVE_STAGE/bin/helicase-adapter" helicase "$HELICASE_VERSION"
    check_stripped_binary "$ACTIVE_STAGE/bin/needletail-adapter"
    check_stripped_binary "$ACTIVE_STAGE/bin/helicase-adapter"
    check_wrapper_smoke \
        "$ACTIVE_STAGE/bin/needletail-adapter" \
        "$ACTIVE_STAGE/bin/helicase-adapter"

    local install_id
    install_id="install-$(date -u +%Y%m%dT%H%M%SZ)-$$"
    local published="$WRAPPER_INSTALLS_DIR/$install_id"
    mv -- "$ACTIVE_STAGE" "$published"
    ACTIVE_STAGE=""

    switch_relative_link "$WRAPPER_CURRENT" "installs/wrappers/$install_id"
    switch_relative_link \
        "$BIN_DIR/needletail-adapter" \
        '../.local/wrappers-current/bin/needletail-adapter'
    switch_relative_link \
        "$BIN_DIR/helicase-adapter" \
        '../.local/wrappers-current/bin/helicase-adapter'
    remove_old_installs wrappers "$published"

    check_wrappers
    finish_work
}

peer_version() {
    case "$1" in
        seqtk) printf '%s\n' "$SEQTK_VERSION" ;;
        fqtools) printf '%s\n' "$FQTOOLS_VERSION" ;;
        fq) printf '%s\n' "$FQ_VERSION" ;;
        fqkit) printf '%s\n' "$FQKIT_VERSION" ;;
        seqkit) printf '%s\n' "$SEQKIT_VERSION" ;;
        rasusa) printf '%s\n' "$RASUSA_VERSION" ;;
        fasten) printf '%s\n' "$FASTEN_VERSION" ;;
        seqfu) printf '%s\n' "$SEQFU_VERSION" ;;
        irma-core) printf '%s\n' "$IRMA_CORE_VERSION" ;;
        fastqvalidator) printf '%s\n' "$FASTQ_VALIDATOR_VERSION" ;;
        fastp) printf '%s\n' "$FASTP_VERSION" ;;
        bbtools) printf '%s\n' "$BBTOOLS_VERSION" ;;
        *) return 1 ;;
    esac
}

peer_recipe() {
    case "$1" in
        seqtk) printf '%s\n' "$SEQTK_RECIPE" ;;
        fqtools) printf '%s\n' "$FQTOOLS_RECIPE" ;;
        fq) printf '%s\n' "$FQ_RECIPE" ;;
        fqkit) printf '%s\n' "$FQKIT_RECIPE" ;;
        seqkit) printf '%s\n' "$SEQKIT_RECIPE" ;;
        rasusa) printf '%s\n' "$RASUSA_RECIPE" ;;
        fasten) printf '%s\n' "$FASTEN_RECIPE" ;;
        seqfu) printf '%s\n' "$SEQFU_RECIPE" ;;
        irma-core) printf '%s\n' "$IRMA_CORE_RECIPE" ;;
        fastqvalidator) printf '%s\n' "$FASTQ_VALIDATOR_RECIPE" ;;
        fastp) printf '%s\n' "$FASTP_RECIPE" ;;
        bbtools) printf '%s\n' "$BBTOOLS_RECIPE" ;;
        *) return 1 ;;
    esac
}

peer_binaries() {
    case "$1" in
        seqtk) printf 'seqtk\n' ;;
        fqtools) printf 'fqtools\n' ;;
        fq) printf 'fq\n' ;;
        fqkit) printf 'fqkit\n' ;;
        seqkit) printf 'seqkit\n' ;;
        rasusa) printf 'rasusa\n' ;;
        fasten) printf 'fasten_sample\n' ;;
        seqfu) printf 'seqfu\n' ;;
        irma-core) printf 'irma-core\n' ;;
        fastqvalidator) printf 'fastQValidator\n' ;;
        fastp) printf 'fastp\n' ;;
        bbtools) printf 'bbversion.sh reformat.sh repair.sh\n' ;;
        *) return 1 ;;
    esac
}

peer_source() {
    case "$1" in
        seqtk) printf 'tag v%s\n' "$SEQTK_VERSION" ;;
        fqtools) printf 'commit %s; htslib %s\n' "$FQTOOLS_COMMIT" "$HTSLIB_VERSION" ;;
        fq) printf 'release v%s\n' "$FQ_VERSION" ;;
        fqkit) printf 'tag version%s\n' "$FQKIT_VERSION" ;;
        seqkit) printf 'release v%s\n' "$SEQKIT_VERSION" ;;
        rasusa) printf 'release %s\n' "$RASUSA_VERSION" ;;
        fasten) printf 'published crate %s with Cargo.lock\n' "$FASTEN_VERSION" ;;
        seqfu) printf 'tag v%s; Nim %s; locked Nimble dependencies\n' "$SEQFU_VERSION" "$NIM_VERSION" ;;
        irma-core) printf 'release v%s\n' "$IRMA_CORE_VERSION" ;;
        fastqvalidator)
            printf 'commit %s; libStatGen v%s\n' \
                "$FASTQ_VALIDATOR_COMMIT" \
                "$LIBSTATGEN_VERSION"
            ;;
        fastp) printf 'release v%s\n' "$FASTP_VERSION" ;;
        bbtools) printf 'tag v%s\n' "$BBTOOLS_VERSION" ;;
        *) return 1 ;;
    esac
}

peer_uses_build_jobs() {
    case "$1" in
        seqtk|fqtools|fqkit|fasten|seqfu|fastqvalidator) return 0 ;;
        *) return 1 ;;
    esac
}

peer_build_environment() {
    case "$1" in
        seqtk|fqtools)
            require_command cc >/dev/null
            cc --version | sed -n '1p'
            ;;
        fqkit|fasten)
            require_command rustc >/dev/null
            rustc --version
            ;;
        seqfu)
            require_command cc >/dev/null
            printf 'Nim %s; %s\n' "$NIM_VERSION" "$(cc --version | sed -n '1p')"
            ;;
        fastqvalidator)
            require_command c++ >/dev/null
            c++ --version | sed -n '1p'
            ;;
        bbtools)
            require_command java >/dev/null
            java -version 2>&1 | sed -n '1p'
            ;;
        *)
            printf 'upstream release binary\n'
            ;;
    esac
}

peer_current_link() {
    printf '%s/%s-current\n' "$LOCAL_DIR" "$1"
}

peer_current() {
    local name="$1"
    local current
    current="$(peer_current_link "$name")"
    [[ -L "$current" ]] || return 1

    local receipt="$current/receipt.tsv"
    [[ -f "$receipt" ]] || return 1
    [[ -f "$current/.z-fastq-tool" ]] || return 1
    [[ "$(<"$current/.z-fastq-tool")" == "$name" ]] || return 1
    receipt_has "$receipt" version "$(peer_version "$name")" || return 1
    receipt_has "$receipt" recipe "$(peer_recipe "$name")" || return 1
    receipt_has "$receipt" platform x86_64-linux || return 1
    receipt_has "$receipt" source "$(peer_source "$name")" || return 1
    if ! peer_uses_build_jobs "$name" && grep -q '^jobs'$'\t' "$receipt"; then
        return 1
    fi

    local binary
    for binary in $(peer_binaries "$name"); do
        [[ -x "$current/bin/$binary" ]] || return 1
        [[ -x "$BIN_DIR/$binary" ]] || return 1
    done

    if [[ "$name" == fqtools ]] &&
        [[ "$TOOLS_DIR/patches/fqtools-gzfile-casts.patch" -nt "$receipt" ]]; then
        return 1
    fi
    if [[ "$name" == seqfu ]] &&
        [[ "$TOOLS_DIR/seqfu.nimble.lock" -nt "$receipt" ]]; then
        return 1
    fi
}

check_peer_version() {
    local name="$1"
    local root="$2"
    local output
    case "$name" in
        seqtk)
            output="$("$root/bin/seqtk" 2>&1 || true)"
            [[ "$output" == *"Version: $SEQTK_VERSION"* ]]
            ;;
        fqtools)
            output="$("$root/bin/fqtools" -v 2>&1)"
            [[ "$output" == *"$FQTOOLS_VERSION"* ]]
            ;;
        fq)
            output="$("$root/bin/fq" --version 2>&1)"
            [[ "$output" == *"$FQ_VERSION"* ]]
            ;;
        fqkit)
            output="$("$root/bin/fqkit" --version 2>&1)"
            [[ "$output" == *"$FQKIT_VERSION"* ]]
            ;;
        seqkit)
            output="$("$root/bin/seqkit" version 2>&1)"
            [[ "$output" == *"$SEQKIT_VERSION"* ]]
            ;;
        rasusa)
            output="$("$root/bin/rasusa" --version 2>&1)"
            [[ "$output" == *"$RASUSA_VERSION"* ]]
            ;;
        fasten)
            output="$("$root/bin/fasten_sample" --version 2>&1)"
            [[ "$output" == *"$FASTEN_VERSION"* ]]
            ;;
        seqfu)
            output="$("$root/bin/seqfu" version 2>&1)"
            [[ "$output" == "$SEQFU_VERSION" ]]
            ;;
        irma-core)
            output="$("$root/bin/irma-core" --version 2>&1)"
            [[ "$output" == *"$IRMA_CORE_VERSION"* ]]
            ;;
        fastqvalidator)
            return 0
            ;;
        fastp)
            output="$("$root/bin/fastp" --version 2>&1)"
            [[ "$output" == *"$FASTP_VERSION"* ]]
            ;;
        bbtools)
            output="$("$root/bin/bbversion.sh" 2>&1)"
            [[ "$output" == *"$BBTOOLS_VERSION"* ]]
            ;;
    esac || {
        echo "error: $name did not report version $(peer_version "$name")" >&2
        return 1
    }
}

check_peer_smoke() {
    local name="$1"
    local root="$2"
    local fixture="$ROOT_DIR/tests/data/synthetic/basic_valid.fastq"
    local short_fixture="$ROOT_DIR/tests/data/synthetic/acgtn_valid.fastq"
    local smoke_dir
    smoke_dir="$(mktemp -d /tmp/z-fastq-peer-smoke.XXXXXX)"

    (
        trap 'rm -rf -- "$smoke_dir"' EXIT
        local output
        case "$name" in
            seqtk)
                output="$("$root/bin/seqtk" size "$fixture")"
                [[ "$output" == 5$'\t'* ]]
                ;;
            fqtools)
                output="$("$root/bin/fqtools" count "$fixture")"
                [[ "$output" == 5 ]]
                ;;
            fq)
                "$root/bin/fq" lint "$short_fixture" >/dev/null
                ;;
            fqkit)
                output="$("$root/bin/fqkit" size "$fixture" 2>/dev/null)"
                [[ "$output" == reads:5$'\t'bases:21$'\t'* ]]
                ;;
            seqkit)
                output="$("$root/bin/seqkit" stats -T -j 1 "$fixture")"
                [[ "$output" == *$'\t5\t21\t'* ]]
                ;;
            rasusa)
                "$root/bin/rasusa" reads \
                    --num 2 \
                    --seed 1 \
                    --output "$smoke_dir/rasusa.fastq" \
                    "$fixture" >/dev/null 2>&1
                [[ "$(wc -l <"$smoke_dir/rasusa.fastq")" == 8 ]]
                ;;
            fasten)
                printf '%s\n' \
                    '@pair/1' 'ACGTN' '+' '#####' \
                    '@pair/2' 'TGCAN' '+' '#####' \
                    >"$smoke_dir/interleaved.fastq"
                "$root/bin/fasten_sample" --paired-end --frequency 1 \
                    <"$smoke_dir/interleaved.fastq" \
                    >"$smoke_dir/fasten.fastq"
                cmp "$smoke_dir/interleaved.fastq" "$smoke_dir/fasten.fastq"
                ;;
            seqfu)
                output="$("$root/bin/seqfu" count "$fixture")"
                [[ "$output" == "$fixture"$'\t5\t'* ]]
                gzip -c "$fixture" >"$smoke_dir/basic.fastq.gz"
                output="$("$root/bin/seqfu" count "$smoke_dir/basic.fastq.gz")"
                [[ "$output" == "$smoke_dir/basic.fastq.gz"$'\t5\t'* ]]
                ;;
            irma-core)
                "$root/bin/irma-core" sampler \
                    --subsample-target 1 \
                    --rng-seed 1 \
                    --output "$smoke_dir/sample.fastq" \
                    "$short_fixture"
                [[ "$(wc -l <"$smoke_dir/sample.fastq")" == 4 ]]

                sed -n '1,4p' "$fixture" | sed '1s/$/\/1/' \
                    >"$smoke_dir/source-r1.fastq"
                sed -n '1,4p' "$fixture" | sed '1s/$/\/2/' \
                    >"$smoke_dir/source-r2.fastq"
                "$root/bin/irma-core" xleave \
                    --output "$smoke_dir/interleaved.fastq" \
                    "$smoke_dir/source-r1.fastq" \
                    "$smoke_dir/source-r2.fastq"
                "$root/bin/irma-core" xleave \
                    --output "$smoke_dir/r1.fastq" \
                    --output2 "$smoke_dir/r2.fastq" \
                    "$smoke_dir/interleaved.fastq"
                cmp "$smoke_dir/source-r1.fastq" "$smoke_dir/r1.fastq"
                cmp "$smoke_dir/source-r2.fastq" "$smoke_dir/r2.fastq"
                ;;
            fastqvalidator)
                "$root/bin/fastQValidator" \
                    --file "$short_fixture" \
                    --minReadLen 0 \
                    --disableSeqIDCheck \
                    --quiet
                ;;
            fastp)
                "$root/bin/fastp" \
                    --in1 "$short_fixture" \
                    --out1 "$smoke_dir/fastp.fastq" \
                    --json "$smoke_dir/fastp.json" \
                    --html "$smoke_dir/fastp.html" \
                    --thread 1 \
                    --disable_adapter_trimming >/dev/null 2>&1
                [[ -s "$smoke_dir/fastp.json" ]]
                ;;
            bbtools)
                printf '%s\n' \
                    '@pair/1' 'ACGTA' '+' '#####' \
                    '@pair/2' 'TGCAC' '+' '#####' \
                    >"$smoke_dir/interleaved.fastq"
                "$root/bin/reformat.sh" \
                    -Xmx200m \
                    in="$smoke_dir/interleaved.fastq" \
                    out="$smoke_dir/reformat.fastq" \
                    int=t \
                    verifyinterleaved=t \
                    samplerate=1 \
                    sampleseed=11 \
                    overwrite=t \
                    threads=1 \
                    qin=33 \
                    qout=33 >/dev/null 2>&1
                cmp "$smoke_dir/interleaved.fastq" "$smoke_dir/reformat.fastq"
                sed -n '1,4p' "$smoke_dir/interleaved.fastq" \
                    >"$smoke_dir/source-r1.fastq"
                sed -n '5,8p' "$smoke_dir/interleaved.fastq" \
                    >"$smoke_dir/source-r2.fastq"
                "$root/bin/repair.sh" \
                    in1="$smoke_dir/source-r1.fastq" \
                    in2="$smoke_dir/source-r2.fastq" \
                    out1="$smoke_dir/repaired-r1.fastq" \
                    out2="$smoke_dir/repaired-r2.fastq" \
                    overwrite=t \
                    threads=1 \
                    qin=33 >/dev/null 2>&1
                cmp "$smoke_dir/source-r1.fastq" "$smoke_dir/repaired-r1.fastq"
                cmp "$smoke_dir/source-r2.fastq" "$smoke_dir/repaired-r2.fastq"
                ;;
        esac
    ) || {
        echo "error: $name smoke check failed" >&2
        return 1
    }
}

check_peer_at() {
    local name="$1"
    local root="$2"
    if [[ "$name" == bbtools ]]; then
        require_command java
    else
        check_stripped_binary "$root/bin/$(peer_binaries "$name")"
    fi
    check_peer_version "$name" "$root" || return
    check_peer_smoke "$name" "$root" || return
}

check_peer() {
    local name="$1"
    if ! peer_current "$name"; then
        echo "error: $name is missing or stale; run tools/install.sh $name" >&2
        return 1
    fi
    local current
    current="$(peer_current_link "$name")"
    check_peer_at "$name" "$current"
    echo "ok: $name $(peer_version "$name")"
}

write_peer_receipt() {
    local name="$1"
    local receipt="$2"
    local build="$3"
    {
        printf 'tool\t%s\n' "$name"
        printf 'version\t%s\n' "$(peer_version "$name")"
        printf 'recipe\t%s\n' "$(peer_recipe "$name")"
        printf 'source\t%s\n' "$(peer_source "$name")"
        printf 'platform\tx86_64-linux\n'
        printf 'build_environment\t%s\n' "$(peer_build_environment "$name")"
        if peer_uses_build_jobs "$name"; then
            printf 'jobs\t%s\n' "$TOOL_JOBS"
        fi
        printf 'build\t%s\n' "$build"
        printf 'runtime\t%s\n' "$(peer_binaries "$name")"
    } >"$receipt"
}

build_seqtk() {
    local stage="$1"
    local archive="$ACTIVE_WORK/seqtk.tar.gz"
    download \
        "https://github.com/lh3/seqtk/archive/refs/tags/v$SEQTK_VERSION.tar.gz" \
        "$archive"
    tar -xzf "$archive" -C "$ACTIVE_WORK"
    make --no-print-directory -s -C "$ACTIVE_WORK/seqtk-$SEQTK_VERSION" \
        -j "$TOOL_JOBS" \
        CC=cc
    install -m 755 "$ACTIVE_WORK/seqtk-$SEQTK_VERSION/seqtk" "$stage/bin/seqtk"
    printf -v BUILD_DESCRIPTION 'make -j %s CC=cc' "$TOOL_JOBS"
}

build_fqtools() {
    local stage="$1"
    local source_archive="$ACTIVE_WORK/fqtools.tar.gz"
    local htslib_archive="$ACTIVE_WORK/htslib.tar.bz2"
    local source="$ACTIVE_WORK/fqtools-$FQTOOLS_COMMIT"
    local htslib="$ACTIVE_WORK/htslib-$HTSLIB_VERSION"
    local build_log="$ACTIVE_WORK/fqtools-build.log"

    download \
        "https://github.com/alastair-droop/fqtools/archive/$FQTOOLS_COMMIT.tar.gz" \
        "$source_archive"
    download \
        "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2" \
        "$htslib_archive"
    tar -xzf "$source_archive" -C "$ACTIVE_WORK"
    tar -xjf "$htslib_archive" -C "$ACTIVE_WORK"
    patch -d "$source" -p1 <"$TOOLS_DIR/patches/fqtools-gzfile-casts.patch"

    echo "build: fqtools HTSlib $HTSLIB_VERSION"
    if ! (
        cd "$htslib"
        ./configure \
            --quiet \
            --disable-bz2 \
            --disable-lzma \
            --disable-s3 \
            --disable-gcs \
            --disable-libcurl \
            --disable-plugins \
            --without-libdeflate
        make --no-print-directory -s -j "$TOOL_JOBS" lib-static
    ) >"$build_log" 2>&1; then
        echo "error: HTSlib build failed for fqtools" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi
    echo "build: fqtools $FQTOOLS_VERSION"
    if ! make --no-print-directory -s -C "$source" \
        -j "$TOOL_JOBS" \
        HTSDIR="$htslib" \
        CFLAGS="-O3 -fcommon" \
        CPPFLAGS="-Wall -Wextra -Wno-unused-parameter -I$htslib -I$htslib/htslib" \
        LIBS="$htslib/libhts.a -lz -lm -lpthread" \
        >>"$build_log" 2>&1; then
        echo "error: fqtools build failed" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi
    install -m 755 "$source/bin/fqtools" "$stage/bin/fqtools"
    printf -v BUILD_DESCRIPTION 'patched fqtools; HTSlib %s static; make -j %s' \
        "$HTSLIB_VERSION" \
        "$TOOL_JOBS"
}

install_release_binary() {
    local stage="$1"
    local url="$2"
    local archive_name="$3"
    local binary_name="$4"
    local archive="$ACTIVE_WORK/$archive_name"
    download "$url" "$archive"
    mkdir -p "$ACTIVE_WORK/extract"
    tar -xf "$archive" -C "$ACTIVE_WORK/extract"
    local source
    source="$(find_binary "$ACTIVE_WORK/extract" "$binary_name")"
    install -m 755 "$source" "$stage/bin/$binary_name"
}

build_fq() {
    local stage="$1"
    local asset="fq-$FQ_VERSION-x86_64-unknown-linux-gnu.tar.gz"
    install_release_binary \
        "$stage" \
        "https://github.com/stjude-rust-labs/fq/releases/download/v$FQ_VERSION/$asset" \
        "$asset" \
        fq
    BUILD_DESCRIPTION='upstream Linux x86-64 release binary'
}

build_fqkit() {
    local stage="$1"
    local archive="$ACTIVE_WORK/fqkit.tar.gz"
    local source="$ACTIVE_WORK/fqkit-version$FQKIT_VERSION"
    download \
        "https://github.com/BioinfoToolbox/fqkit/archive/refs/tags/version$FQKIT_VERSION.tar.gz" \
        "$archive"
    tar -xzf "$archive" -C "$ACTIVE_WORK"
    env \
        CARGO_HOME="$ACTIVE_WORK/cargo-home" \
        CARGO_TARGET_DIR="$ACTIVE_WORK/cargo-target" \
        cargo build \
            --quiet \
            --locked \
            --release \
            --jobs "$TOOL_JOBS" \
            --manifest-path "$source/Cargo.toml"
    install -m 755 "$ACTIVE_WORK/cargo-target/release/fqkit" "$stage/bin/fqkit"
    printf -v BUILD_DESCRIPTION 'cargo build --locked --release -j %s' "$TOOL_JOBS"
}

build_seqkit() {
    local stage="$1"
    local asset="seqkit_linux_amd64.tar.gz"
    install_release_binary \
        "$stage" \
        "https://github.com/shenwei356/seqkit/releases/download/v$SEQKIT_VERSION/$asset" \
        "$asset" \
        seqkit
    BUILD_DESCRIPTION='upstream Linux amd64 release binary'
}

build_rasusa() {
    local stage="$1"
    local asset="rasusa-$RASUSA_VERSION-x86_64-unknown-linux-gnu.tar.gz"
    install_release_binary \
        "$stage" \
        "https://github.com/mbhall88/rasusa/releases/download/$RASUSA_VERSION/$asset" \
        "$asset" \
        rasusa
    BUILD_DESCRIPTION='upstream Linux x86-64 GNU release binary'
}

build_fasten() {
    local stage="$1"
    local archive="$ACTIVE_WORK/fasten-$FASTEN_VERSION.crate"
    local source="$ACTIVE_WORK/fasten-$FASTEN_VERSION"
    local cargo_home="$ACTIVE_WORK/cargo-home"
    local target="$ACTIVE_WORK/target"
    local build_log="$ACTIVE_WORK/fasten-build.log"
    download \
        "https://static.crates.io/crates/fasten/fasten-$FASTEN_VERSION.crate" \
        "$archive"
    tar -xzf "$archive" -C "$ACTIVE_WORK"
    echo "build: fasten $FASTEN_VERSION"
    if ! CARGO_HOME="$cargo_home" CARGO_TARGET_DIR="$target" \
        cargo build --locked --release --bin fasten_sample \
        --manifest-path "$source/Cargo.toml" -j "$TOOL_JOBS" \
        >"$build_log" 2>&1; then
        echo "error: Fasten build failed" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi
    install -m 755 "$target/release/fasten_sample" "$stage/bin/fasten_sample"
    printf -v BUILD_DESCRIPTION \
        'cargo build --locked --release --bin fasten_sample -j %s' \
        "$TOOL_JOBS"
}

build_seqfu() {
    local stage="$1"
    local source_archive="$ACTIVE_WORK/seqfu.tar.gz"
    local nim_archive="$ACTIVE_WORK/nim.tar.xz"
    local source="$ACTIVE_WORK/seqfu2-$SEQFU_VERSION"
    local nim_root="$ACTIVE_WORK/nim-$NIM_VERSION"
    local nimble_dir="$ACTIVE_WORK/nimble"
    local nimcache="$ACTIVE_WORK/nimcache"
    local output="$ACTIVE_WORK/seqfu"
    local build_log="$ACTIVE_WORK/seqfu-build.log"

    download \
        "https://github.com/telatin/seqfu2/archive/refs/tags/v$SEQFU_VERSION.tar.gz" \
        "$source_archive"
    download \
        "https://nim-lang.org/download/nim-$NIM_VERSION-linux_x64.tar.xz" \
        "$nim_archive"
    printf '%s  %s\n' "$NIM_LINUX_X64_SHA256" "$nim_archive" | sha256sum -c - >/dev/null

    tar -xzf "$source_archive" -C "$ACTIVE_WORK"
    tar -xJf "$nim_archive" -C "$ACTIVE_WORK"
    cp "$TOOLS_DIR/seqfu.nimble.lock" "$source/nimble.lock"
    mkdir -p "$nimble_dir" "$nimcache"

    echo "build: seqfu $SEQFU_VERSION with temporary Nim $NIM_VERSION"
    if ! (
        cd "$source"
        env \
            GIT_CONFIG_GLOBAL=/dev/null \
            GIT_CONFIG_NOSYSTEM=1 \
            NIMBLE_DIR="$nimble_dir" \
            PATH="$nim_root/bin:$PATH" \
            "$nim_root/bin/nimble" \
                --nimbleDir:"$nimble_dir" \
                --nim:"$nim_root/bin/nim" \
                --useSystemNim \
                --disableNimBinaries \
                --accept \
                --noColor \
                install \
                --depsOnly
    ) >"$build_log" 2>&1; then
        echo "error: SeqFu dependency installation failed" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi
    cmp "$TOOLS_DIR/seqfu.nimble.lock" "$source/nimble.lock"

    local -a nim_paths=()
    local package_dir
    while IFS= read -r -d '' package_dir; do
        nim_paths+=("--path:$package_dir")
    done < <(find "$nimble_dir/pkgs2" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

    if ! (
        cd "$source"
        env \
            NIMBLE_DIR="$nimble_dir" \
            PATH="$nim_root/bin:$PATH" \
            "$nim_root/bin/nim" c \
                --noNimblePath \
                "${nim_paths[@]}" \
                --mm:orc \
                -d:NimblePkgVersion="$SEQFU_VERSION" \
                -d:release \
                --opt:speed \
                --passC:-Wno-error=incompatible-pointer-types \
                --parallelBuild:"$TOOL_JOBS" \
                --nimcache:"$nimcache" \
                --out:"$output" \
                src/sfu.nim
    ) >>"$build_log" 2>&1; then
        echo "error: SeqFu build failed" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi

    install -m 755 "$output" "$stage/bin/seqfu"
    printf -v BUILD_DESCRIPTION \
        'source build with temporary Nim %s, locked Nimble dependencies, parallel build %s' \
        "$NIM_VERSION" \
        "$TOOL_JOBS"
}

build_irma_core() {
    local stage="$1"
    local asset="irma-core-linux-x86_64-v$IRMA_CORE_VERSION.tar.gz"
    install_release_binary \
        "$stage" \
        "https://github.com/CDCgov/irma-core/releases/download/v$IRMA_CORE_VERSION/$asset" \
        "$asset" \
        irma-core
    BUILD_DESCRIPTION='upstream Linux x86-64 release binary'
}

build_fastqvalidator() {
    local stage="$1"
    local source_archive="$ACTIVE_WORK/fastQValidator.tar.gz"
    local library_archive="$ACTIVE_WORK/libStatGen.tar.gz"
    local source="$ACTIVE_WORK/fastQValidator-$FASTQ_VALIDATOR_COMMIT"
    local library="$ACTIVE_WORK/libStatGen-$LIBSTATGEN_VERSION"
    local build_log="$ACTIVE_WORK/fastqvalidator-build.log"

    download \
        "https://github.com/statgen/fastQValidator/archive/$FASTQ_VALIDATOR_COMMIT.tar.gz" \
        "$source_archive"
    download \
        "https://github.com/statgen/libStatGen/archive/refs/tags/v$LIBSTATGEN_VERSION.tar.gz" \
        "$library_archive"
    tar -xzf "$source_archive" -C "$ACTIVE_WORK"
    tar -xzf "$library_archive" -C "$ACTIVE_WORK"
    echo "build: fastqvalidator $FASTQ_VALIDATOR_VERSION"
    if ! make --no-print-directory -s -C "$source" \
        -j "$TOOL_JOBS" \
        LIB_PATH_FASTQ_VALIDATOR="$library" \
        >"$build_log" 2>&1; then
        echo "error: FastQValidator build failed" >&2
        tail -n 100 "$build_log" >&2
        return 1
    fi
    local binary
    binary="$(find_binary "$source" fastQValidator)"
    install -m 755 "$binary" "$stage/bin/fastQValidator"
    printf -v BUILD_DESCRIPTION \
        'make -j %s with libStatGen v%s' \
        "$TOOL_JOBS" \
        "$LIBSTATGEN_VERSION"
}

build_fastp() {
    local stage="$1"
    download \
        "https://opengene.org/fastp/fastp.$FASTP_VERSION" \
        "$ACTIVE_WORK/fastp"
    install -m 755 "$ACTIVE_WORK/fastp" "$stage/bin/fastp"
    BUILD_DESCRIPTION='upstream Linux release binary'
}

build_bbtools() {
    local stage="$1"
    local archive="$ACTIVE_WORK/bbtools.tar.gz"
    local source="$ACTIVE_WORK/BBTools-$BBTOOLS_VERSION"
    local script
    download \
        "https://github.com/bbushnell/BBTools/archive/refs/tags/v$BBTOOLS_VERSION.tar.gz" \
        "$archive"
    tar -xzf "$archive" -C "$ACTIVE_WORK"

    mkdir -p "$stage/runtime"
    cp -a "$source/current" "$stage/runtime/current"
    find "$stage/runtime/current" -type f ! -name '*.class' -delete
    for script in bbversion.sh reformat.sh repair.sh javasetup.sh memdetect.sh; do
        install -m 755 "$source/$script" "$stage/runtime/$script"
    done
    for script in bbversion.sh reformat.sh repair.sh; do
        ln -s "../runtime/$script" "$stage/bin/$script"
    done
    BUILD_DESCRIPTION='upstream compiled Java classes; source removed'
}

build_peer() {
    local name="$1"
    local stage="$2"
    case "$name" in
        seqtk) build_seqtk "$stage" ;;
        fqtools) build_fqtools "$stage" ;;
        fq) build_fq "$stage" ;;
        fqkit) build_fqkit "$stage" ;;
        seqkit) build_seqkit "$stage" ;;
        rasusa) build_rasusa "$stage" ;;
        fasten) build_fasten "$stage" ;;
        seqfu) build_seqfu "$stage" ;;
        irma-core) build_irma_core "$stage" ;;
        fastqvalidator) build_fastqvalidator "$stage" ;;
        fastp) build_fastp "$stage" ;;
        bbtools) build_bbtools "$stage" ;;
        *) return 1 ;;
    esac
}

strip_peer() {
    local name="$1"
    local stage="$2"
    [[ "$name" != bbtools ]] || return 0

    strip --strip-all "$stage/bin/$(peer_binaries "$name")"
    BUILD_DESCRIPTION="$BUILD_DESCRIPTION; stripped"
}

require_peer_commands() {
    local name="$1"
    require_command curl
    require_command install
    if [[ "$name" != bbtools ]]; then
        require_command strip
    fi
    if [[ "$name" != fastp ]]; then
        require_command tar
    fi
    case "$name" in
        seqtk)
            require_command make
            require_command cc
            ;;
        fqtools)
            require_command make
            require_command cc
            require_command patch
            ;;
        fqkit|fasten)
            require_command cargo
            require_command rustc
            ;;
        seqfu)
            require_command cc
            require_command git
            require_command gzip
            require_command sha256sum
            require_command xz
            ;;
        fastqvalidator)
            require_command make
            require_command c++
            ;;
        bbtools)
            require_command java
            ;;
    esac
}

install_peer() {
    local name="$1"
    validate_settings
    check_platform
    acquire_install_lock

    if peer_current "$name" && check_peer "$name"; then
        echo "reuse: $name"
        return
    fi

    require_peer_commands "$name"

    local installs="$INSTALLS_DIR/$name"
    mkdir -p "$LOCAL_DIR/stage" "$installs" "$BIN_DIR"
    remove_stale_stages "$name"
    ACTIVE_WORK="$(mktemp -d /tmp/z-fastq-tools.XXXXXX)"
    ACTIVE_STAGE="$(mktemp -d "$LOCAL_DIR/stage/$name.XXXXXX")"
    printf '%s\n' "$name" >"$ACTIVE_STAGE/.z-fastq-tool"
    mkdir -p "$ACTIVE_STAGE/bin"

    echo "install: $name $(peer_version "$name")"
    BUILD_DESCRIPTION=""
    build_peer "$name" "$ACTIVE_STAGE"
    strip_peer "$name" "$ACTIVE_STAGE"
    write_peer_receipt "$name" "$ACTIVE_STAGE/receipt.tsv" "$BUILD_DESCRIPTION"
    check_peer_at "$name" "$ACTIVE_STAGE"

    local install_id
    install_id="install-$(date -u +%Y%m%dT%H%M%SZ)-$$"
    local published="$installs/$install_id"
    mv -- "$ACTIVE_STAGE" "$published"
    ACTIVE_STAGE=""

    switch_relative_link "$(peer_current_link "$name")" "installs/$name/$install_id"
    local binary
    for binary in $(peer_binaries "$name"); do
        switch_relative_link \
            "$BIN_DIR/$binary" \
            "../.local/$name-current/bin/$binary"
    done
    remove_old_installs "$name" "$published"
    check_peer "$name"
    finish_work
}

install_peers() {
    local name
    for name in "${PEER_NAMES[@]}"; do
        install_peer "$name"
    done
}

check_peers() {
    local name
    for name in "${PEER_NAMES[@]}"; do
        check_peer "$name"
    done
}

is_peer_name() {
    local requested="$1"
    local name
    for name in "${PEER_NAMES[@]}"; do
        if [[ "$requested" == "$name" ]]; then
            return 0
        fi
    done
    return 1
}

list_targets() {
    local state="missing"
    local detail="run tools/install.sh"

    if [[ -x "$VENV_DIR/bin/python" ]]; then
        state="ready"
        detail="$("$VENV_DIR/bin/python" --version 2>&1)"
    fi
    printf '%-16s %-8s %s\n' venv "$state" "$detail"

    state="missing"
    detail="run tools/install.sh wrappers"
    if wrappers_current; then
        state="ready"
        detail="Needletail $NEEDLETAIL_VERSION, Helicase $HELICASE_VERSION"
    elif [[ -e "$WRAPPER_CURRENT" || -L "$WRAPPER_CURRENT" ]]; then
        state="stale"
    fi
    printf '%-16s %-8s %s\n' wrappers "$state" "$detail"

    local name
    for name in "${PEER_NAMES[@]}"; do
        state="missing"
        detail="run tools/install.sh $name"
        if peer_current "$name"; then
            state="ready"
            detail="$(peer_version "$name")"
        elif [[ -e "$(peer_current_link "$name")" || -L "$(peer_current_link "$name")" ]]; then
            state="stale"
        fi
        printf '%-16s %-8s %s\n' "$name" "$state" "$detail"
    done
}

run_install_target() {
    case "$1" in
        venv) install_venv ;;
        wrappers) install_wrappers ;;
        peers) install_peers ;;
        all)
            install_venv
            install_wrappers
            install_peers
            ;;
        seqtk|fqtools|fq|fqkit|seqkit|rasusa|fasten|seqfu|irma-core|fastqvalidator|fastp|bbtools)
            install_peer "$1"
            ;;
        *) return 1 ;;
    esac
}

run_check_target() {
    case "$1" in
        venv) check_venv ;;
        wrappers) check_wrappers ;;
        peers) check_peers ;;
        all)
            check_venv
            check_wrappers
            check_peers
            ;;
        seqtk|fqtools|fq|fqkit|seqkit|rasusa|fasten|seqfu|irma-core|fastqvalidator|fastp|bbtools)
            check_peer "$1"
            ;;
        *) return 1 ;;
    esac
}

is_named_target() {
    case "$1" in
        venv|wrappers|peers|all) return 0 ;;
        *) is_peer_name "$1" ;;
    esac
}

case "${1:-}" in
    "")
        if [[ $# -ne 0 ]]; then
            usage >&2
            exit 2
        fi
        install_venv
        ;;
    --check)
        if [[ $# -gt 2 ]]; then
            usage >&2
            exit 2
        fi
        if ! is_named_target "${2:-venv}"; then
            echo "error: unknown check target: ${2:-venv}" >&2
            usage >&2
            exit 2
        fi
        run_check_target "${2:-venv}"
        ;;
    --list)
        if [[ $# -ne 1 ]]; then
            usage >&2
            exit 2
        fi
        list_targets
        ;;
    --help|-h)
        usage
        ;;
    *)
        if [[ $# -ne 1 ]]; then
            usage >&2
            exit 2
        fi
        if ! is_named_target "$1"; then
            echo "error: unknown tool or option: $1" >&2
            usage >&2
            exit 2
        fi
        run_install_target "$1"
        ;;
esac

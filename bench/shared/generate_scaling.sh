#!/usr/bin/env bash

set -euo pipefail

SCALING_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCALING_DATA_DIR="${SCALING_SHARED_DIR}/cache/scaling"
SCALING_STAMP="${SCALING_DATA_DIR}/.stamp"
SCALING_SIZE_COUNT=100000
SCALING_READ_LENGTH=150
SCALING_HEADER_WIDTH=17
SCALING_RECORD_OVERHEAD=4
SCALING_DNA_SEQUENCE=""

for ((SCALING_DNA_INDEX = 0; SCALING_DNA_INDEX < 16; SCALING_DNA_INDEX++)); do
    SCALING_DNA_SEQUENCE+="ACGTACGTAC"
done
SCALING_DNA_HASH="$(printf '%s' "${SCALING_DNA_SEQUENCE}" | sha256sum | awk '{ print substr($1, 1, 16) }')"

scaling_csv() {
    local joined=""
    local value
    for value in "$@"; do
        joined+="${value},"
    done
    printf '%s' "${joined%,}"
}

SCALING_SIZE_MBS=(1 5 10 25 50 100 250 500)
SCALING_READ_COUNTS=(100000 250000 500000 1000000)
SCALING_GZIP_LEVEL=6

scaling_stamp_payload() {
    printf 'schema=fastq-scaling.v1\n'
    printf 'size_mbs=%s\n' "$(scaling_csv "${SCALING_SIZE_MBS[@]}")"
    printf 'reads_fixed_counts=%s\n' "$(scaling_csv "${SCALING_READ_COUNTS[@]}")"
    printf 'size_read_count=%s\n' "${SCALING_SIZE_COUNT}"
    printf 'fixed_read_len=%s\n' "${SCALING_READ_LENGTH}"
    printf 'gzip_level=%s\n' "${SCALING_GZIP_LEVEL}"
    printf 'header_width=%s\n' "${SCALING_HEADER_WIDTH}"
    printf 'dna_hash=%s\n' "${SCALING_DNA_HASH}"
    printf 'qual_byte=I\n'
}

scaling_output_paths() {
    local mode="$1"
    local size_mb
    if [[ "${mode}" == all || "${mode}" == size ]]; then
        for size_mb in "${SCALING_SIZE_MBS[@]}"; do
            printf '%s\n' "${SCALING_DATA_DIR}/size_${size_mb}mb.fastq"
        done
    fi
    local read_count
    if [[ "${mode}" == all || "${mode}" == reads || "${mode}" == seq ]]; then
        for read_count in "${SCALING_READ_COUNTS[@]}"; do
            printf '%s\n' "${SCALING_DATA_DIR}/reads_fixed_${read_count}.fastq"
        done
    fi
    if [[ "${mode}" == all || "${mode}" == gzip ]]; then
        for read_count in "${SCALING_READ_COUNTS[@]}"; do
            printf '%s\n' "${SCALING_DATA_DIR}/reads_fixed_${read_count}.fastq.gz"
        done
    fi
}

scaling_ready() {
    [[ -f "${SCALING_STAMP}" ]] || return 1
    local expected_stamp
    expected_stamp="$(scaling_stamp_payload)"
    [[ "$(<"${SCALING_STAMP}")" == "${expected_stamp}" ]] || return 1

    local output_path
    while IFS= read -r output_path; do
        [[ -s "${output_path}" ]] || return 1
    done < <(scaling_output_paths "${SCALING_READY_MODE}")
}

scaling_write_plain() {
    local output_path="$1"
    local record_count="$2"
    local read_length="$3"
    local temporary_path="$4"

    awk -v count="${record_count}" -v read_len="${read_length}" '
        function make_sequence(wanted_length,    base, result, i) {
            base = ""
            for (i = 0; i < 16; i++) {
                base = base "ACGTACGTAC"
            }
            result = ""
            while (length(result) < wanted_length) {
                result = result base
            }
            return substr(result, 1, wanted_length)
        }
        BEGIN {
            sequence = make_sequence(read_len)
            quality = ""
            for (i = 0; i < read_len; i++) {
                quality = quality "I"
            }
            for (i = 1; i <= count; i++) {
                printf "@r%015d\n%s\n+\n%s\n", i, sequence, quality
            }
        }
    ' >"${temporary_path}"
    mv -- "${temporary_path}" "${output_path}"
}

scaling_write_gzip() {
    local plain_path="$1"
    local gzip_path="$2"
    local temporary_path="$3"
    gzip -n -"${SCALING_GZIP_LEVEL}" -c "${plain_path}" >"${temporary_path}"
    mv -- "${temporary_path}" "${gzip_path}"
}

scaling_size_read_length() {
    local size_mb="$1"
    local target_bytes=$((size_mb * 1024 * 1024))
    local bytes_per_record=$((target_bytes / SCALING_SIZE_COUNT))
    local read_length=$(((bytes_per_record - SCALING_HEADER_WIDTH - SCALING_RECORD_OVERHEAD) / 2))
    ((read_length < 1)) && read_length=1
    printf '%s\n' "${read_length}"
}

scaling_write_stamp() {
    local temporary_path="$1"
    scaling_stamp_payload >"${temporary_path}"
    mv -- "${temporary_path}" "${SCALING_STAMP}"
}

scaling_main() {
    local mode="all"
    local force=0
    while (($# > 0)); do
        case "$1" in
            --mode)
                (($# >= 2)) || {
                    printf 'generate_scaling.sh: --mode needs a value\n' >&2
                    return 1
                }
                mode="$2"
                shift 2
                ;;
            --force)
                force=1
                shift
                ;;
            -h | --help)
                printf 'usage: %s [--mode all|size|reads|gzip] [--force]\n' "${BASH_SOURCE[0]}"
                return 0
                ;;
            *)
                printf 'generate_scaling.sh: unknown option: %s\n' "$1" >&2
                return 1
                ;;
        esac
    done

    case "${mode}" in
        all | size | reads | gzip) ;;
        seq) mode=reads ;;
        *)
            printf 'generate_scaling.sh: invalid mode: %s\n' "${mode}" >&2
            return 1
            ;;
    esac

    mkdir -p "${SCALING_DATA_DIR}"
    SCALING_TEMP_DIR="$(mktemp -d "${SCALING_DATA_DIR}/.scaling.XXXXXX")"
    trap 'rm -rf -- "${SCALING_TEMP_DIR}"' EXIT

    SCALING_READY_MODE="${mode}"
    if ((force == 0)) && scaling_ready; then
        printf 'scaling: ready\n'
        return 0
    fi

    local size_mb read_count read_length plain_path gzip_path
    if [[ "${mode}" == all || "${mode}" == size ]]; then
        for size_mb in "${SCALING_SIZE_MBS[@]}"; do
            read_length="$(scaling_size_read_length "${size_mb}")"
            plain_path="${SCALING_DATA_DIR}/size_${size_mb}mb.fastq"
            scaling_write_plain "${plain_path}" "${SCALING_SIZE_COUNT}" "${read_length}" "${SCALING_TEMP_DIR}/size_${size_mb}mb.fastq"
        done
    fi

    if [[ "${mode}" == all || "${mode}" == reads ]]; then
        for read_count in "${SCALING_READ_COUNTS[@]}"; do
            plain_path="${SCALING_DATA_DIR}/reads_fixed_${read_count}.fastq"
            scaling_write_plain "${plain_path}" "${read_count}" "${SCALING_READ_LENGTH}" "${SCALING_TEMP_DIR}/reads_fixed_${read_count}.fastq"
        done
    fi

    if [[ "${mode}" == all || "${mode}" == gzip ]]; then
        for read_count in "${SCALING_READ_COUNTS[@]}"; do
            plain_path="${SCALING_DATA_DIR}/reads_fixed_${read_count}.fastq"
            [[ -s "${plain_path}" ]] || {
                printf 'scaling: missing %s; run --mode reads first\n' "${plain_path}" >&2
                return 1
            }
            gzip_path="${SCALING_DATA_DIR}/reads_fixed_${read_count}.fastq.gz"
            scaling_write_gzip "${plain_path}" "${gzip_path}" "${SCALING_TEMP_DIR}/reads_fixed_${read_count}.fastq.gz"
        done
    fi

    if [[ "${mode}" == all || "${mode}" == size || "${mode}" == reads || "${mode}" == gzip ]]; then
        scaling_write_stamp "${SCALING_TEMP_DIR}/.scaling"
        printf 'scaling: generated\n'
    fi
}

scaling_main "$@"

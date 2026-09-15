#!/usr/bin/env bash

set -euo pipefail

DERIVED_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DERIVED_DATA_DIR="${DERIVED_SHARED_DIR}/data"
DERIVED_PLAIN_DIR="${DERIVED_SHARED_DIR}/cache/plain"
# shellcheck disable=SC1091
source "${DERIVED_SHARED_DIR}/catalog.sh"

derived_sha256() {
    sha256sum "$1" | awk '{ print $1 }'
}

derived_plain_path() {
    local id="$1"
    local filename="${CATALOG_FILENAME[$id]}"
    printf '%s/%s\n' "${DERIVED_PLAIN_DIR}" "${filename%.gz}"
}

derived_gzip_path() {
    local id="$1"
    printf '%s/%s\n' "${DERIVED_DATA_DIR}" "${CATALOG_FILENAME[$id]}"
}

derived_stamp_path() {
    local id="$1"
    printf '%s/%s.derived\n' "${DERIVED_DATA_DIR}" "${CATALOG_FILENAME[$id]}"
}

derived_kind() {
    local id="$1"
    local rest="${CATALOG_URL[$id]#derived:}"
    printf '%s\n' "${rest%%:*}"
}

derived_dependency_list() {
    local id="$1"
    catalog_derived_deps "${id}"
}

derived_stamp_payload() {
    local id="$1"
    local kind="$2"
    local expected="$3"
    shift 3
    printf 'schema=fastq-derived.v1\n'
    printf 'id=%s\n' "${id}"
    printf 'kind=%s\n' "${kind}"
    printf 'expected=%s\n' "${expected}"
    local IFS=,
    printf 'deps=%s\n' "$*"

    local dep dep_plain dep_hash
    for dep in "$@"; do
        dep_plain="$(derived_plain_path "${dep}")"
        [[ -f "${dep_plain}" ]] || return 1
        dep_hash="$(derived_sha256 "${dep_plain}")"
        printf 'dep_sha=%s:%s\n' "$(basename "${dep_plain}")" "${dep_hash}"
    done
}

derived_ready() {
    local id="$1"
    local kind="$2"
    local expected="$3"
    shift 3
    local gzip_path plain_path stamp_path expected_stamp
    gzip_path="$(derived_gzip_path "${id}")"
    plain_path="$(derived_plain_path "${id}")"
    stamp_path="$(derived_stamp_path "${id}")"
    [[ -s "${gzip_path}" && -s "${plain_path}" && -f "${stamp_path}" ]] || return 1
    expected_stamp="$(derived_stamp_payload "${id}" "${kind}" "${expected}" "$@")" || return 1
    [[ "$(<"${stamp_path}")" == "${expected_stamp}" ]]
}

derived_concat_gzip() {
    local input_path="$1"
    local split_records="$2"
    local first_path="$3"
    local second_path="$4"
    awk -v split_records="${split_records}" -v first_path="${first_path}" -v second_path="${second_path}" '
        {
            sub(/\r$/, "")
            if (NR <= split_records * 4) {
                print > first_path
            } else {
                print > second_path
            }
        }
        END {
            close(first_path)
            close(second_path)
        }
    ' "${input_path}"
    cat -- "${first_path}" "${second_path}"
}

derived_interleave() {
    local left_path="$1"
    local right_path="$2"
    awk -v left_path="${left_path}" -v right_path="${right_path}" '
        function clean(line) {
            sub(/\r$/, "", line)
            return line
        }
        function read_record(path,    i, status) {
            status = getline line < path
            if (status != 1) {
                return status
            }
            print clean(line)
            for (i = 1; i < 4; i++) {
                status = getline line < path
                if (status != 1) {
                    exit 2
                }
                print clean(line)
            }
            return 1
        }
        BEGIN {
            while (1) {
                left_status = read_record(left_path)
                if (left_status == 0) {
                    right_status = getline line < right_path
                    if (right_status != 0) {
                        exit 3
                    }
                    break
                }
                if (left_status < 0) {
                    exit 4
                }
                right_status = read_record(right_path)
                if (right_status != 1) {
                    exit 5
                }
            }
            close(left_path)
            close(right_path)
        }
    '
}

derived_rewrite_headers() {
    local input_path="$1"
    local output_path="$2"
    local kind="$3"
    local mate="$4"
    awk -v kind="${kind}" -v mate="${mate}" '
        NR % 4 == 1 {
            header = $0
            sub(/\r$/, "", header)
            sub(/^@/, "", header)
            sub(/[[:space:]].*$/, "", header)
            if (kind == "casava") {
                print "@" header " " mate ":N:0:ATCACG"
            } else {
                sub(/\/[12]$/, "", header)
                print "@" header "/" mate
            }
            next
        }
        NR % 4 == 3 {
            if (kind == "casava") {
                print "+"
            } else {
                line = $0
                sub(/\r$/, "", line)
                print line
            }
            next
        }
        {
            line = $0
            sub(/\r$/, "", line)
            print line
        }
    ' "${input_path}" >"${output_path}"
}

derived_record_count() {
    awk 'END { if (NR % 4 != 0) exit 1; print NR / 4 }' "$1"
}

derived_build_one() {
    local id="$1"
    catalog_has_id "${id}" || catalog_error "unknown dataset id ${id}"
    [[ "${CATALOG_URL[$id]}" == derived:* ]] || catalog_error "${id} is not derived"

    local kind expected
    kind="$(derived_kind "${id}")"
    expected="${CATALOG_EXPECTED[$id]}"

    local dependencies_text
    dependencies_text="$(derived_dependency_list "${id}")" || return 1
    [[ -n "${dependencies_text}" ]] || catalog_error "${id} has no dependencies"
    local -a dependencies=()
    mapfile -t dependencies <<<"${dependencies_text}"

    local dep dep_plain
    for dep in "${dependencies[@]}"; do
        [[ -n "${dep}" ]] || continue
        catalog_has_id "${dep}" || catalog_error "unknown dependency ${dep} for ${id}"
        dep_plain="$(derived_plain_path "${dep}")"
        [[ -s "${dep_plain}" ]] || catalog_error "missing plain dependency ${dep}: ${dep_plain}"
    done

    if ((DERIVED_FORCE == 0)) && derived_ready "${id}" "${kind}" "${expected}" "${dependencies[@]}"; then
        printf 'derived: %s ready\n' "${id}"
        return 0
    fi

    local output_gzip output_plain output_stamp
    output_gzip="$(derived_gzip_path "${id}")"
    output_plain="$(derived_plain_path "${id}")"
    output_stamp="$(derived_stamp_path "${id}")"
    local work_plain="${DERIVED_TEMP_DIR}/${id}.fastq"
    local concat_first_plain=""
    local concat_second_plain=""

    case "${kind}" in
        concat-gzip)
            ((${#dependencies[@]} == 1)) || catalog_error "${id}: concat-gzip needs one dependency"
            local source_id="${dependencies[0]}"
            local source_expected="${CATALOG_EXPECTED[$source_id]}"
            ((source_expected % 2 == 0)) || catalog_error "${id}: source record count is odd"
            concat_first_plain="${DERIVED_TEMP_DIR}/${id}.member1.fastq"
            concat_second_plain="${DERIVED_TEMP_DIR}/${id}.member2.fastq"
            derived_concat_gzip \
                "$(derived_plain_path "${dependencies[0]}")" \
                "$((source_expected / 2))" \
                "${concat_first_plain}" \
                "${concat_second_plain}" >"${work_plain}"
            ;;
        interleave)
            ((${#dependencies[@]} == 2)) || catalog_error "${id}: interleave needs two dependencies"
            [[ "${CATALOG_EXPECTED[${dependencies[0]}]}" == "${CATALOG_EXPECTED[${dependencies[1]}]}" ]] ||
                catalog_error "${id}: interleave dependencies have different record counts"
            derived_interleave \
                "$(derived_plain_path "${dependencies[0]}")" \
                "$(derived_plain_path "${dependencies[1]}")" >"${work_plain}"
            ;;
        casava1 | casava2 | slash1 | slash2)
            ((${#dependencies[@]} == 1)) || catalog_error "${id}: ${kind} needs one dependency"
            local mate rewrite_kind
            [[ "${kind}" == *2 ]] && mate=2 || mate=1
            [[ "${kind}" == casava* ]] && rewrite_kind=casava || rewrite_kind=slash
            derived_rewrite_headers \
                "$(derived_plain_path "${dependencies[0]}")" \
                "${work_plain}" \
                "${rewrite_kind}" \
                "${mate}"
            ;;
        *) catalog_error "${id}: unsupported derived kind ${kind}" ;;
    esac

    local actual_records
    actual_records="$(derived_record_count "${work_plain}")" ||
        catalog_error "${id}: generated FASTQ does not contain complete four-line records"
    [[ "${actual_records}" == "${expected}" ]] ||
        catalog_error "${id}: generated ${actual_records} records, expected ${expected}"

    local temporary_gzip="${DERIVED_TEMP_DIR}/${id}.fastq.gz"
    local temporary_plain="${DERIVED_TEMP_DIR}/${id}.plain.fastq"
    local temporary_stamp="${DERIVED_TEMP_DIR}/${id}.derived"
    if [[ "${kind}" == concat-gzip ]]; then
        gzip -n -6 -c "${concat_first_plain}" >"${DERIVED_TEMP_DIR}/${id}.member1.gz"
        gzip -n -6 -c "${concat_second_plain}" >"${DERIVED_TEMP_DIR}/${id}.member2.gz"
        cat -- "${DERIVED_TEMP_DIR}/${id}.member1.gz" "${DERIVED_TEMP_DIR}/${id}.member2.gz" >"${temporary_gzip}"
    else
        gzip -n -6 -c "${work_plain}" >"${temporary_gzip}"
    fi
    gzip -dc "${temporary_gzip}" >"${temporary_plain}"
    derived_stamp_payload "${id}" "${kind}" "${expected}" "${dependencies[@]}" >"${temporary_stamp}"

    mv -- "${temporary_gzip}" "${output_gzip}"
    mv -- "${temporary_plain}" "${output_plain}"
    mv -- "${temporary_stamp}" "${output_stamp}"
    printf 'derived: %s generated\n' "${id}"
}

derived_main() {
    DERIVED_FORCE=0
    local -a requested=()
    while (($# > 0)); do
        case "$1" in
            --force)
                DERIVED_FORCE=1
                shift
                ;;
            -h | --help)
                printf 'usage: %s [--force] DATASET_ID...\n' "${BASH_SOURCE[0]}"
                return 0
                ;;
            *)
                requested+=("$1")
                shift
                ;;
        esac
    done
    ((${#requested[@]} > 0)) || {
        printf 'generate_derived.sh: no dataset ids supplied\n' >&2
        return 1
    }

    catalog_load
    local expanded_text
    expanded_text="$(catalog_expand_ids "${requested[@]}")" || return 1
    local -a expanded=()
    mapfile -t expanded <<<"${expanded_text}"
    ((${#expanded[@]} > 0)) || {
        printf 'generate_derived.sh: selection is empty\n' >&2
        return 1
    }

    mkdir -p "${DERIVED_DATA_DIR}" "${DERIVED_PLAIN_DIR}"
    DERIVED_TEMP_DIR="$(mktemp -d "${DERIVED_DATA_DIR}/.derived.XXXXXX")"
    trap 'rm -rf -- "${DERIVED_TEMP_DIR}"' EXIT

    local id
    for id in "${expanded[@]}"; do
        [[ "${CATALOG_URL[$id]}" == derived:* ]] || continue
        derived_build_one "${id}"
    done
}

derived_main "$@"

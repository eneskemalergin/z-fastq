#!/usr/bin/env bash

# Shared benchmark registry. This file can be sourced by benchmark scripts or
# invoked directly for simple registry queries.

CATALOG_SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CATALOG_MANIFEST="${CATALOG_SHARED_DIR}/datasets.manifest"
CATALOG_FEATURES="${CATALOG_SHARED_DIR}/features.tsv"

declare -g CATALOG_LOADED="${CATALOG_LOADED:-0}"
declare -gA CATALOG_FILENAME
declare -gA CATALOG_EXPECTED
declare -gA CATALOG_ACCESSION
declare -gA CATALOG_URL
declare -gA CATALOG_LOCAL_PATH
declare -ga CATALOG_ID_ORDER
declare -ga CATALOG_FEATURE_SUITE
declare -ga CATALOG_FEATURE_SET
declare -ga CATALOG_FEATURE_ROLE
declare -ga CATALOG_FEATURE_IDS

catalog_error() {
    printf 'catalog: %s\n' "$*" >&2
    return 1
}

catalog_load() {
    [[ "${CATALOG_LOADED}" == 1 ]] && return 0

    [[ -f "${CATALOG_MANIFEST}" ]] || catalog_error "missing ${CATALOG_MANIFEST}"
    [[ -f "${CATALOG_FEATURES}" ]] || catalog_error "missing ${CATALOG_FEATURES}"

    local id filename expected accession url local_path extra
    while IFS=$'\t' read -r id filename expected accession url local_path extra || [[ -n "${id:-}" ]]; do
        [[ -z "${id:-}" || "${id}" == \#* ]] && continue
        [[ -n "${filename:-}" && -n "${expected:-}" && -n "${url:-}" ]] ||
            catalog_error "invalid manifest row for ${id}"
        [[ -z "${extra:-}" ]] ||
            catalog_error "too many fields in manifest row for ${id}"

        [[ -z "${CATALOG_FILENAME[$id]+present}" ]] ||
            catalog_error "duplicate dataset id ${id}"
        [[ "${expected}" =~ ^[0-9]+$ ]] ||
            catalog_error "invalid expected record count for ${id}: ${expected}"

        CATALOG_FILENAME["${id}"]="${filename}"
        CATALOG_EXPECTED["${id}"]="${expected}"
        CATALOG_ACCESSION["${id}"]="${accession:-}"
        CATALOG_URL["${id}"]="${url}"
        CATALOG_LOCAL_PATH["${id}"]="${local_path:-}"
        CATALOG_ID_ORDER+=("${id}")
    done <"${CATALOG_MANIFEST}"

    local suite set_name role ids feature_extra
    while IFS=$'\t' read -r suite set_name role ids feature_extra || [[ -n "${suite:-}" ]]; do
        [[ -z "${suite:-}" || "${suite}" == \#* ]] && continue
        [[ -n "${set_name:-}" && -n "${role:-}" && -n "${ids:-}" ]] ||
            catalog_error "invalid feature row for ${suite}"
        [[ -z "${feature_extra:-}" ]] ||
            catalog_error "too many fields in feature row for ${suite}/${set_name}/${role}"

        CATALOG_FEATURE_SUITE+=("${suite}")
        CATALOG_FEATURE_SET+=("${set_name}")
        CATALOG_FEATURE_ROLE+=("${role}")
        CATALOG_FEATURE_IDS+=("${ids}")
    done <"${CATALOG_FEATURES}"

    CATALOG_LOADED=1
}

catalog_has_id() {
    catalog_load
    [[ -n "${CATALOG_FILENAME[$1]+present}" ]]
}

catalog_derived_deps() {
    catalog_load
    local id="$1"
    local spec="${CATALOG_URL[$id]}"
    [[ "${spec}" == derived:* ]] || return 0

    local rest="${spec#derived:}"
    [[ "${rest}" == *:* ]] || catalog_error "invalid derived specification for ${id}: ${spec}"
    local deps="${rest#*:}"
    [[ -n "${deps}" ]] || return 0

    local -a dep_list=()
    local dep
    IFS=',' read -r -a dep_list <<<"${deps}"
    for dep in "${dep_list[@]}"; do
        [[ -n "${dep}" ]] && printf '%s\n' "${dep}"
    done
}

catalog_emit_with_deps() {
    local id="$1"
    catalog_has_id "${id}" || catalog_error "unknown dataset id ${id}"

    local state="${CATALOG_VISIT[${id}]:-0}"
    [[ "${state}" != 2 ]] || return 0
    [[ "${state}" != 1 ]] || catalog_error "cycle involving dataset ${id}"
    CATALOG_VISIT["${id}"]=1

    local dep_list dep
    dep_list="$(catalog_derived_deps "${id}")" || return 1
    while IFS= read -r dep; do
        [[ -n "${dep}" ]] || continue
        catalog_emit_with_deps "${dep}" || return 1
    done <<<"${dep_list}"

    CATALOG_VISIT["${id}"]=2
    printf '%s\n' "${id}"
}

catalog_expand_ids() {
    catalog_load
    (($# > 0)) || catalog_error "no dataset ids supplied"

    declare -gA CATALOG_VISIT=()
    local id
    for id in "$@"; do
        catalog_emit_with_deps "${id}" || return 1
    done
}

catalog_suite_ids() {
    catalog_load
    local wanted_suite="$1"
    local wanted_set="$2"
    local wanted_role="${3:-}"
    local -A seen=()
    local found=0
    local i id
    local -a feature_id_list=()

    for ((i = 0; i < ${#CATALOG_FEATURE_SUITE[@]}; i++)); do
        [[ "${CATALOG_FEATURE_SUITE[$i]}" == "${wanted_suite}" ]] || continue
        [[ "${CATALOG_FEATURE_SET[$i]}" == "${wanted_set}" ]] || continue
        [[ -z "${wanted_role}" || "${CATALOG_FEATURE_ROLE[$i]}" == "${wanted_role}" ]] || continue

        found=1
        IFS=',' read -r -a feature_id_list <<<"${CATALOG_FEATURE_IDS[$i]}"
        for id in "${feature_id_list[@]}"; do
            [[ -n "${id}" ]] || continue
            [[ -n "${seen[$id]+present}" ]] && continue
            seen["${id}"]=1
            printf '%s\n' "${id}"
        done
    done

    ((found == 1)) || catalog_error "no feature row for ${wanted_suite}/${wanted_set}/${wanted_role}"
}

catalog_ids_command() {
    catalog_load
    local suite=""
    local set_name=""
    local role=""
    local ids_arg=""
    local select_all=0
    local no_expand=0

    while (($# > 0)); do
        case "$1" in
            --suite)
                (($# >= 2)) || catalog_error "--suite needs a value"
                suite="$2"
                shift 2
                ;;
            --set)
                (($# >= 2)) || catalog_error "--set needs a value"
                set_name="$2"
                shift 2
                ;;
            --role)
                (($# >= 2)) || catalog_error "--role needs a value"
                role="$2"
                shift 2
                ;;
            --ids)
                (($# >= 2)) || catalog_error "--ids needs a value"
                ids_arg="$2"
                shift 2
                ;;
            --all)
                select_all=1
                shift
                ;;
            --no-expand)
                no_expand=1
                shift
                ;;
            *)
                catalog_error "unknown ids option: $1"
                ;;
        esac
    done

    local -a requested=()
    local id suite_ids_text
    if ((select_all == 1)); then
        requested=("${CATALOG_ID_ORDER[@]}")
    elif [[ -n "${ids_arg}" ]]; then
        local -a requested_ids=()
        IFS=',' read -r -a requested_ids <<<"${ids_arg}"
        requested=("${requested_ids[@]}")
    elif [[ -n "${suite}" && -n "${set_name}" ]]; then
        suite_ids_text="$(catalog_suite_ids "${suite}" "${set_name}" "${role}")" || return 1
        mapfile -t requested <<<"${suite_ids_text}"
    else
        catalog_error "choose --all, --ids, or --suite with --set"
    fi

    ((${#requested[@]} > 0)) || catalog_error "selection is empty"
    if ((no_expand == 1)); then
        printf '%s\n' "${requested[@]}"
    else
        catalog_expand_ids "${requested[@]}"
    fi
}

catalog_field_command() {
    catalog_load
    (($# == 2)) || catalog_error "field needs an id and a field name"
    local id="$1"
    local field="$2"
    catalog_has_id "${id}" || catalog_error "unknown dataset id ${id}"

    case "${field}" in
        id) printf '%s\n' "${id}" ;;
        filename) printf '%s\n' "${CATALOG_FILENAME[$id]}" ;;
        expected_records) printf '%s\n' "${CATALOG_EXPECTED[$id]}" ;;
        accession) printf '%s\n' "${CATALOG_ACCESSION[$id]}" ;;
        url) printf '%s\n' "${CATALOG_URL[$id]}" ;;
        local_path) printf '%s\n' "${CATALOG_LOCAL_PATH[$id]}" ;;
        derived)
            [[ "${CATALOG_URL[$id]}" == derived:* ]] && printf 'true\n' || printf 'false\n'
            ;;
        derived_kind)
            if [[ "${CATALOG_URL[$id]}" == derived:* ]]; then
                local rest="${CATALOG_URL[$id]#derived:}"
                printf '%s\n' "${rest%%:*}"
            else
                printf '\n'
            fi
            ;;
        derived_deps)
            local deps_text
            deps_text="$(catalog_derived_deps "${id}")" || return 1
            if [[ -n "${deps_text}" ]]; then
                printf '%s\n' "${deps_text}" | paste -sd, -
            else
                printf '\n'
            fi
            ;;
        *) catalog_error "unknown field ${field}" ;;
    esac
}

catalog_expected_command() {
    (($# == 1)) || catalog_error "expected needs a dataset id"
    catalog_field_command "$1" expected_records
}

catalog_check_command() {
    catalog_load
    ((${#CATALOG_ID_ORDER[@]} > 0)) || catalog_error "manifest has no datasets"
    ((${#CATALOG_FEATURE_SUITE[@]} > 0)) || catalog_error "features file has no rows"

    local i id dep feature_key
    local -A feature_rows=()
    for ((i = 0; i < ${#CATALOG_FEATURE_SUITE[@]}; i++)); do
        case "${CATALOG_FEATURE_SUITE[$i]}" in
            count | stats | check | sample | interleave | deinterleave) ;;
            *) catalog_error "unknown suite ${CATALOG_FEATURE_SUITE[$i]}" ;;
        esac
        feature_key="${CATALOG_FEATURE_SUITE[$i]}/${CATALOG_FEATURE_SET[$i]}/${CATALOG_FEATURE_ROLE[$i]}"
        [[ -z "${feature_rows[$feature_key]+present}" ]] ||
            catalog_error "duplicate feature row ${feature_key}"
        feature_rows["${feature_key}"]=1

        local -a feature_id_list=()
        IFS=',' read -r -a feature_id_list <<<"${CATALOG_FEATURE_IDS[$i]}"
        ((${#feature_id_list[@]} > 0)) || catalog_error "empty dataset list for ${feature_key}"
        for id in "${feature_id_list[@]}"; do
            [[ -n "${id}" ]] || catalog_error "empty dataset id in ${feature_key}"
            catalog_has_id "${id}" || catalog_error "unknown dataset ${id} in ${feature_key}"
        done
    done

    for id in "${CATALOG_ID_ORDER[@]}"; do
        local derived_spec="${CATALOG_URL[$id]}"
        [[ "${derived_spec}" == derived:* ]] || continue
        local deps_text
        deps_text="$(catalog_derived_deps "${id}")" || return 1
        [[ -n "${deps_text}" ]] || catalog_error "derived dataset ${id} has no dependencies"
        while IFS= read -r dep; do
            [[ -n "${dep}" ]] || continue
            catalog_has_id "${dep}" || catalog_error "unknown dependency ${dep} for ${id}"
        done <<<"${deps_text}"
    done

    local required
    for required in \
        count/publication/time \
        count/publication/check \
        count/small/time \
        count/small/check \
        stats/publication/time \
        stats/small/time; do
        [[ -n "${feature_rows[$required]+present}" ]] ||
            catalog_error "missing feature row ${required}"
    done

    local set_name time_text time_count
    for set_name in publication small; do
        time_text="$(catalog_suite_ids count "${set_name}" time)" || return 1
        if [[ -n "${time_text}" ]]; then
            time_count="$(printf '%s\n' "${time_text}" | awk 'NF { count += 1 } END { print count + 0 }')"
        else
            time_count=0
        fi
        [[ "${time_count}" == 3 ]] ||
            catalog_error "count/${set_name}/time must contain exactly 3 datasets"
    done

    for id in "${CATALOG_ID_ORDER[@]}"; do
        catalog_expand_ids "${id}" >/dev/null || return 1
    done
}

catalog_main() {
    (($# > 0)) || catalog_error "use ids, expected, field, or check"
    case "$1" in
        ids)
            shift
            catalog_ids_command "$@"
            ;;
        expected)
            shift
            catalog_expected_command "$@"
            ;;
        field)
            shift
            catalog_field_command "$@"
            ;;
        check)
            shift
            (($# == 0)) || catalog_error "check does not take options"
            catalog_check_command
            ;;
        *) catalog_error "unknown command $1" ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail
    catalog_main "$@"
fi

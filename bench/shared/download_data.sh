#!/usr/bin/env bash
# Copy or download REAL gzip FASTQ, then build derived fixtures.
#
# Default: feature set count/publication (Dense, Variable, Long, ConcatGzip).
#   --small / --suite NAME --small
#   --suite NAME [--small]
#   --pair-small   MiniSeq R1 and R2
#   --pair         trimmed MiSeq R1 and R2
#   --all          every manifest row
#   --ids X,Y      specific manifest ids (derived deps included)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/tools.sh"

DATA_DIR="$SCRIPT_DIR/data"
CACHE_PLAIN="$SCRIPT_DIR/cache/plain"
SELECT="suite:count:publication"

usage() {
    cat <<'EOF'
Usage: bash bench/shared/download_data.sh [options]

  --suite NAME     count|stats|check|sample|interleave|deinterleave
  --small          feature set "small" (with --suite, or count/small alone)
  --pair-small     MiniSeq R1 and R2
  --pair           trimmed MiSeq R1 and R2
  --all            every row in datasets.manifest
  --ids            comma-separated manifest ids
EOF
}

SUITE=""
SET_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --small)
            SET_NAME="small"
            if [[ "$SELECT" == suite:count:publication ]]; then
                SELECT="suite:count:small"
            fi
            shift
            ;;
        --suite)
            SUITE="$2"
            shift 2
            ;;
        --pair-small) SELECT="ids:DenseSmall,PairSmallR2"; shift ;;
        --pair) SELECT="ids:PairR1,PairR2"; shift ;;
        --all) SELECT="all"; shift ;;
        --ids)
            SELECT="ids:$2"
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "error: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -n "$SUITE" ]]; then
    [[ -n "$SET_NAME" ]] || SET_NAME="publication"
    SELECT="suite:${SUITE}:${SET_NAME}"
fi

catalog_check_command

selected_catalog_ids() {
    case "$SELECT" in
        all) catalog_ids_command --all ;;
        suite:*)
            local rest="${SELECT#suite:}"
            local suite="${rest%%:*}"
            local set_name="${rest#*:}"
            catalog_ids_command --suite "$suite" --set "$set_name"
            ;;
        ids:*)
            catalog_ids_command --ids "${SELECT#ids:}"
            ;;
        *)
            echo "error: unknown selection $SELECT" >&2
            exit 1
            ;;
    esac
}

file_sha256() {
    local path="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum -- "$path" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 -- "$path" | awk '{print $1}'
    else
        echo "error: need sha256sum or shasum" >&2
        exit 1
    fi
}

meta_field() {
    local meta="$1" key="$2"
    [[ -f "$meta" ]] || return 0
    awk -F'\t' -v k="$key" '$1==k{print $2; exit}' "$meta"
}

ensure_gzip() {
    local id="$1" filename="$2" expected_records="$3" accession="$4" url="$5" local_rel="$6"
    local dest="$DATA_DIR/$filename"
    local local_src="$PROJECT_ROOT/$local_rel"
    local meta="$DATA_DIR/${filename}.meta"
    local stem="${filename%.gz}"
    local plain="$CACHE_PLAIN/$stem"

    mkdir -p "$DATA_DIR"
    echo "[$id] $filename  $accession  (expected $expected_records records)"

    if [[ -f "$dest" ]]; then
        local old_source old_expected
        old_source="$(meta_field "$meta" source)"
        old_expected="$(meta_field "$meta" expected_records)"
        if [[ ! -f "$meta" || "$old_source" != "$url" || "$old_expected" != "$expected_records" ]]; then
            echo "  stale gzip (source or record count changed); replacing"
            rm -f -- "$dest" "$meta" "$plain"
        fi
    fi

    if [[ -f "$dest" ]]; then
        echo "  present: $dest ($(file_size_bytes "$dest") bytes)"
    elif [[ -n "$local_rel" && -f "$local_src" ]]; then
        echo "  copy from $local_rel"
        local incoming
        incoming="$(mktemp "$DATA_DIR/.incoming.XXXXXX")"
        if ! cp -a -- "$local_src" "$incoming"; then
            rm -f -- "$incoming"
            echo "error: failed to copy $local_src" >&2
            exit 1
        fi
        mv -- "$incoming" "$dest"
    else
        echo "  download $url"
        local incoming
        incoming="$(mktemp "$DATA_DIR/.incoming.XXXXXX")"
        if ! curl -fL --retry 3 --retry-delay 2 -o "$incoming" "$url"; then
            rm -f -- "$incoming"
            echo "error: download failed for $url" >&2
            exit 1
        fi
        mv -- "$incoming" "$dest"
    fi

    if [[ ! -s "$dest" ]]; then
        echo "error: missing or empty $dest" >&2
        exit 1
    fi

    {
        printf 'id\t%s\n' "$id"
        printf 'filename\t%s\n' "$filename"
        printf 'expected_records\t%s\n' "$expected_records"
        printf 'accession\t%s\n' "$accession"
        printf 'bytes\t%s\n' "$(file_size_bytes "$dest")"
        printf 'sha256\t%s\n' "$(file_sha256 "$dest")"
        printf 'source\t%s\n' "$url"
    } >"$meta"
    echo "  meta: $meta"
}

ensure_plain() {
    local filename="$1"
    local gz="$DATA_DIR/$filename"
    local stem="${filename%.gz}"
    local plain="$CACHE_PLAIN/$stem"

    mkdir -p "$CACHE_PLAIN"
    if [[ -f "$plain" ]]; then
        echo "  plain cache present: $plain"
        return 0
    fi
    echo "  gunzip -> $plain"
    local incoming
    incoming="$(mktemp "$CACHE_PLAIN/.plain.XXXXXX")"
    if ! gzip -dc -- "$gz" >"$incoming"; then
        rm -f -- "$incoming"
        echo "error: cannot decompress $gz" >&2
        exit 1
    fi
    mv -- "$incoming" "$plain"
    if [[ ! -s "$plain" ]]; then
        echo "error: gunzip produced an empty $plain" >&2
        exit 1
    fi
}

mapfile -t SELECTED_IDS < <(selected_catalog_ids)
if [[ "${#SELECTED_IDS[@]}" -eq 0 ]]; then
    echo "error: no matching rows for selection $SELECT" >&2
    exit 1
fi

echo "=== FASTQ datasets ==="
echo "Target: $DATA_DIR"
echo "Selection: $SELECT"
echo "Ids: ${SELECTED_IDS[*]}"
echo

entry_count=0
derived_ids=()
for id in "${SELECTED_IDS[@]}"; do
    filename="${CATALOG_FILENAME[$id]}"
    expected_records="${CATALOG_EXPECTED[$id]}"
    accession="${CATALOG_ACCESSION[$id]}"
    url="${CATALOG_URL[$id]}"
    local_path="${CATALOG_LOCAL_PATH[$id]}"
    if [[ "$url" == derived:* ]]; then
        derived_ids+=("$id")
        continue
    fi
    ensure_gzip "$id" "$filename" "$expected_records" "$accession" "$url" "${local_path:-}"
    ensure_plain "$filename"
    entry_count=$((entry_count + 1))
    echo
done

if [[ "${#derived_ids[@]}" -gt 0 ]]; then
    echo "=== derived FASTQ ==="
    bash "$SCRIPT_DIR/generate_derived.sh" "${derived_ids[@]}"
    for id in "${derived_ids[@]}"; do
        entry_count=$((entry_count + 1))
    done
    echo
fi

if [[ "$entry_count" -eq 0 ]]; then
    echo "error: no matching rows for selection $SELECT" >&2
    exit 1
fi

echo "=== done ($entry_count files) ==="

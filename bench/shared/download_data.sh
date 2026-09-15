#!/usr/bin/env bash
# Copy or download REAL gzip FASTQ into bench/shared/data/, then gunzip a plain cache.
#
# Default: Dense, Variable, Long.
#   --small   DenseSmall, Variable, Long (MiniSeq R1 instead of SRR1810900).
#   --all     every manifest row.
#   --ids X,Y specific manifest ids.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools.sh
source "$SCRIPT_DIR/tools.sh"

DATA_DIR="$SCRIPT_DIR/data"
CACHE_PLAIN="$SCRIPT_DIR/cache/plain"
MANIFEST="$SCRIPT_DIR/datasets.manifest"
SELECT="default"

usage() {
    cat <<'EOF'
Usage: bash bench/shared/download_data.sh [--small|--all] [--ids id,id]

  --small   MiniSeq Dense plus Variable and Long (bring-up).
  --all     every row in datasets.manifest.
  --ids     comma-separated manifest ids.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --small) SELECT="small"; shift ;;
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

if [[ ! -f "$MANIFEST" ]]; then
    echo "error: manifest missing: $MANIFEST" >&2
    exit 1
fi

want_id() {
    local id="$1"
    case "$SELECT" in
        all) return 0 ;;
        small)
            [[ "$id" == "DenseSmall" || "$id" == "Variable" || "$id" == "Long" ]]
            ;;
        default)
            [[ "$id" == "Dense" || "$id" == "Variable" || "$id" == "Long" ]]
            ;;
        ids:*)
            local list="${SELECT#ids:}"
            [[ ",$list," == *",$id,"* ]]
            ;;
        *) return 1 ;;
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

ensure_gzip() {
    local id="$1" filename="$2" expected_records="$3" url="$4" local_rel="$5"
    local dest="$DATA_DIR/$filename"
    local local_src="$PROJECT_ROOT/$local_rel"
    local meta="$DATA_DIR/${filename}.meta"

    mkdir -p "$DATA_DIR"
    echo "[$id] $filename  (expected $expected_records records)"

    if [[ -f "$dest" ]]; then
        echo "  present: $dest ($(file_size_bytes "$dest") bytes)"
    elif [[ -n "$local_rel" && -f "$local_src" ]]; then
        echo "  copy from $local_rel"
        cp -a -- "$local_src" "$dest"
    else
        echo "  download $url"
        curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"
    fi

    if [[ ! -s "$dest" ]]; then
        echo "error: missing or empty $dest" >&2
        exit 1
    fi

    {
        printf 'id\t%s\n' "$id"
        printf 'filename\t%s\n' "$filename"
        printf 'expected_records\t%s\n' "$expected_records"
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
    gzip -dc -- "$gz" >"$plain"
    if [[ ! -s "$plain" ]]; then
        echo "error: gunzip produced an empty $plain" >&2
        exit 1
    fi
}

echo "=== FASTQ datasets ==="
echo "Target: $DATA_DIR"
echo

entry_count=0
while IFS=$'\t' read -r id filename expected_records url local_path || [[ -n "${id:-}" ]]; do
    [[ -z "${id:-}" || "$id" == \#* ]] && continue
    if [[ -z "${filename:-}" || -z "${expected_records:-}" || -z "${url:-}" ]]; then
        echo "error: incomplete manifest row for id=${id:-?}" >&2
        exit 1
    fi
    want_id "$id" || continue
    ensure_gzip "$id" "$filename" "$expected_records" "$url" "${local_path:-}"
    ensure_plain "$filename"
    entry_count=$((entry_count + 1))
    echo
done <"$MANIFEST"

if [[ "$entry_count" -eq 0 ]]; then
    echo "error: no matching rows in $MANIFEST for selection $SELECT" >&2
    exit 1
fi

echo "=== done ($entry_count files) ==="

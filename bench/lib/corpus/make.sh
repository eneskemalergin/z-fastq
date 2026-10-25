#!/usr/bin/env bash
# Build the real-data corpus: fetch pinned originals, then cut nested prefix
# tiers as plain FASTQ and gzip, with one manifest.tsv per source.
set -euo pipefail
export LC_ALL=C
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
catalog="${CORPUS_CATALOG:-$here/real-data.json}"
root="${CORPUS_ROOT:-$here/../../data/real}"
threads="${CORPUS_THREADS:-4}"
# pigz 2.8 writes a different stream with one thread; two or more agree.
((threads >= 2)) || threads=2
usage() {
    cat <<'HELP'
Usage: bash bench/lib/corpus/make.sh [--force] SOURCE...
SOURCE: a catalogue ID, or "all"
  --force   Recut tiers even when the manifest is complete
  --list    Print catalogue IDs, layouts, tiers and original sizes

Output per source, under bench/data/real/ID/:
  single:  TIER.fastq  TIER.fastq.gz
  paired:  TIER_R1.fastq  TIER_R2.fastq  TIER_interleaved.fastq  (+ .gz each)
  manifest.tsv: source tier file records bytes min_len max_len
Missing originals (or pinned prefixes) are downloaded and resumed; verified\noriginals are never fetched again.
A source is fresh, and is not recut, when manifest.tsv lists every catalogue
tier and each listed file exists with its recorded size.
Tiers are record-aligned prefixes of the original, sized by plain bytes with
both mates counted; each smaller tier is a prefix of the larger ones.
Environment: CORPUS_ROOT, CORPUS_CATALOG, CORPUS_THREADS (pigz threads, default 4, minimum 2).
HELP
}
die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}
force=0 ids=()
while (($#)); do
    case "$1" in
        -h | --help) usage; exit 0 ;;
        --force) force=1 ;;
        --list)
            jq -r '.sources[] | [.id, .layout, (.tiers | join(",")), ([.files[].bytes] | add)] | @tsv' "$catalog"
            exit 0
            ;;
        -*) die "unknown option: $1" ;;
        all) mapfile -t -O "${#ids[@]}" ids < <(jq -r '.sources[].id' "$catalog") ;;
        *) ids+=("$1") ;;
    esac
    shift
done
((${#ids[@]})) || { usage >&2; exit 2; }
command -v pigz >/dev/null || die 'pigz is required'

# Fetch the first BYTES of URL into NAME.part, resuming from its current size,
# then accept it only on size and MD5. BYTES below the published size pins a
# prefix: tiers need only the leading records of a large run.
fetch_file() {
    local dest="$1" url="$2" bytes="$3" md5="$4" have last=-1 stalls=0 sum
    printf 'fetch: %s (%s bytes)\n' "${dest##*/}" "$bytes"
    # Servers drop long transfers; retry while each attempt makes progress.
    while ((stalls < 5)); do
        have="$(stat -c%s -- "$dest.part" 2>/dev/null || echo 0)"
        ((have < bytes)) || break
        if ((have > last)); then stalls=0; else stalls=$((stalls + 1)); fi
        last="$have"
        curl --fail --location --silent --show-error --proto '=https' --speed-limit 1024 --speed-time 60 \
            --range "$have-$((bytes - 1))" "$url" >>"$dest.part" || printf 'retry: %s at %s bytes\n' "${dest##*/}" "$have" >&2
    done
    [[ "$(stat -c%s -- "$dest.part" 2>/dev/null || echo 0)" == "$bytes" ]] || die "size mismatch: $dest.part"
    sum="$(md5sum <"$dest.part")"
    [[ "${sum%% *}" == "$md5" ]] || die "md5 mismatch: $dest.part has ${sum%% *}"
    mv -- "$dest.part" "$dest"
}

# The manifest covers exactly the catalogue tiers and every file it lists is intact.
fresh() {
    local dir="$1" tiers="$2" listed file bytes
    [[ -f "$dir/manifest.tsv" ]] || return 1
    listed="$(awk -F'\t' 'NR > 1 && !seen[$2]++ { printf "%s%s", sep, $2; sep = " " }' "$dir/manifest.tsv")"
    [[ "$listed" == "$tiers" ]] || return 1
    while IFS=$'\t' read -r file bytes; do
        [[ -f "$dir/$file" && "$(stat -c%s -- "$dir/$file")" == "$bytes" ]] || return 1
    done < <(awk -F'\t' 'NR > 1 { print $3 "\t" $5 }' "$dir/manifest.tsv")
}

build() {
    local id="$1" src layout dir name url bytes md5 tiers suffixes s tier count total low high recs file
    src="$(jq -ce --arg id "$id" '.tiers as $t | .sources[] | select(.id == $id) |
        .caps = [.tiers[] | "\(.)=\($t[.])"]' "$catalog")" || die "unknown source: $id"
    dir="$root/$id"
    mkdir -p -- "$dir/source"
    if ((force == 0)) && fresh "$dir" "$(jq -r '.tiers | join(" ")' <<<"$src")"; then
        printf 'fresh: %s\n' "$id"
        return
    fi
    local inputs=()
    while IFS=$'\t' read -r name url bytes md5; do
        [[ -f "$dir/source/$name" ]] || fetch_file "$dir/source/$name" "$url" "$bytes" "$md5"
        inputs+=("pigz -dc '$dir/source/$name'")
    done < <(jq -r '.files[] | [.name, .url, .bytes, .md5] | @tsv' <<<"$src")
    layout="$(jq -r .layout <<<"$src")"
    tiers="$(jq -r '.caps | join(" ")' <<<"$src")"
    suffixes=(.fastq)
    [[ "$layout" == single ]] || suffixes=(_R1.fastq _R2.fastq _interleaved.fastq)

    printf 'cut: %s [%s]\n' "$id" "$tiers"
    R1="${inputs[0]}" R2="${inputs[1]:-}" TIERS="$tiers" OUT="$dir/" \
        awk -f "$here/cut.awk" >"$dir/manifest.part"
    printf 'source\ttier\tfile\trecords\tbytes\tmin_len\tmax_len\n' >"$dir/manifest.tsv.part"
    while read -r tier count total low high; do
        for s in "${suffixes[@]}"; do
            file="$tier$s" recs="$count"
            [[ "$s" != _interleaved.fastq ]] || recs=$((count * 2))
            pigz -p "$threads" -6 -n -c -- "$dir/$file.part" >"$dir/$file.gz.part"
            mv -- "$dir/$file.part" "$dir/$file"
            mv -- "$dir/$file.gz.part" "$dir/$file.gz"
            for f in "$file" "$file.gz"; do
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$tier" "$f" "$recs" \
                    "$(stat -c%s -- "$dir/$f")" "$low" "$high" >>"$dir/manifest.tsv.part"
            done
        done
        printf '  %s: %s %s, %s plain bytes\n' "$tier" "$count" \
            "$([[ "$layout" == single ]] && echo records || echo pairs)" "$total"
    done <"$dir/manifest.part"
    rm -- "$dir/manifest.part"
    mv -- "$dir/manifest.tsv.part" "$dir/manifest.tsv"
}

for id in "${ids[@]}"; do build "$id"; done

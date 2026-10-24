#!/usr/bin/env bash
# Check tier cutting on tiny local sources; never contacts dataset servers.
set -euo pipefail
export LC_ALL=C
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
make=(bash "$here/../lib/corpus/make.sh")
work="$(mktemp -d /var/tmp/fastq-corpus-check.XXXXXX)"
export CORPUS_ROOT="$work/data" CORPUS_CATALOG="$work/catalog.json" CORPUS_THREADS=1
cleanup() {
    local d
    rm -f -- "$work"/catalog.json "$work"/log "$work"/data/*/source/* "$work"/data/*/*.*
    for d in "$work"/data/*/source "$work"/data/* "$work/data" "$work"; do [[ ! -d "$d" ]] || rmdir -- "$d"; done
}
trap cleanup EXIT
fail() {
    printf 'fail: %s\n' "$*" >&2
    exit 1
}
# Each mate record is 18 bytes, so a pair is 36 bytes.
records() { for i in $(seq "$2"); do printf '@r%s/%s\nACGT\n+\nIIII\n' "$i" "$1"; done; }
src() {
    mkdir -p -- "$work/data/$1/source"
    gzip -n <<<"$3" >"$work/data/$1/source/$2"
}
src pair R1.gz "$(records 1 3)"
src pair R2.gz "$(records 2 3)"
src crlf reads.gz "$(printf '@a\r\nAC\r\n+\r\nII\r\n@b\nAC\n+\nII\n@c\nAC\n+\nII\n@d\nAC\n+\nII')"
src short R1.gz "$(records 1 2)"
src short R2.gz "$(records 2 2)"
src swap R1.gz "$(records 1 2)"
src swap R2.gz "$(records 2 2 | sed 's/@r1/@rx/')"
src missing x.gz ''
src tech R1.gz "$(records 1 3 | sed 's#/1$#/3#')"
src tech R2.gz "$(records 2 3 | sed 's#/2$#/4#')"
rm -- "$work/data/missing/source/x.gz"
f='{"name":"N","url":"https://fixture.invalid/N","bytes":1,"md5":"0"}'
jq -n --argjson f "$f" '{tiers:{one:40,two:80},sources:[
    {id:"pair",layout:"paired",tiers:["one","two"],files:[($f|.name="R1.gz"),($f|.name="R2.gz")]},
    {id:"crlf",layout:"single",tiers:["one"],files:[($f|.name="reads.gz")]},
    {id:"short",layout:"paired",tiers:["two"],files:[($f|.name="R1.gz"),($f|.name="R2.gz")]},
    {id:"swap",layout:"paired",tiers:["one"],files:[($f|.name="R1.gz"),($f|.name="R2.gz")]},
    {id:"tech",layout:"paired",tiers:["one"],files:[($f|.name="R1.gz"),($f|.name="R2.gz")]},
    {id:"missing",layout:"single",tiers:["one"],files:[($f|.name="x.gz")]}]}' >"$work/catalog.json"

"${make[@]}" pair crlf tech >/dev/null
grep -q '^@r1/3$' "$work/data/tech/one_R1.fastq" || fail 'technical read numbers /3 and /4 were not accepted as mates'
d="$work/data/pair"
cmp -s "$d/one_R1.fastq" <(records 1 1) || fail 'one_R1 is not the first pair'
cmp -s "$d/two_R2.fastq" <(records 2 2) || fail 'two_R2 is not the first two pairs'
cmp -s "$d/two_interleaved.fastq" <(paste -d '\n' <(records 1 2 | paste - - - -) <(records 2 2 | paste - - - -) | tr '\t' '\n') ||
    fail 'interleaved order'
for f in "$d"/*.fastq; do gzip -dc "$f.gz" | cmp -s "$f" - || fail "gzip differs: $f"; done
[[ "$(awk -F'\t' '$3=="two_interleaved.fastq" {print $4, $5, $6, $7}' "$d/manifest.tsv")" == '4 72 4 4' ]] || fail 'manifest row'
[[ "$(wc -l <"$d/manifest.tsv")" == 13 ]] || fail 'manifest rows'
printf '@a\r\nAC\r\n+\r\nII\r\n@b\nAC\n+\nII\n@c\nAC\n+\nII\n' | cmp -s "$work/data/crlf/one.fastq" - || fail 'CRLF bytes changed'
grep -qx 'fresh: pair' <<<"$("${make[@]}" pair)" || fail 'complete manifest was recut'
rm -- "$d/one_R2.fastq.gz"
grep -qx 'cut: pair \[one=40 two=80\]' <<<"$("${make[@]}" pair)" || fail 'missing output was not recut'
[[ -f "$d/one_R2.fastq.gz" ]] || fail 'missing output was not restored'
jq '.tiers.three = 1000 | .sources[0].tiers += ["three"]' "$work/catalog.json" >"$work/log" && mv -- "$work/log" "$work/catalog.json"
grep -q '^cut: pair' <<<"$("${make[@]}" pair 2>/dev/null)" || fail 'added tier did not trigger a recut'
if "${make[@]}" short 2>"$work/log" >/dev/null; then fail 'accepted an underfilled tier'; fi
grep -q 'source ends before tier two' "$work/log" || fail 'underfilled message'
if "${make[@]}" swap 2>"$work/log" >/dev/null; then fail 'accepted mismatched mates'; fi
grep -q 'mate names differ' "$work/log" || fail 'mate message'
if "${make[@]}" missing >"$work/log" 2>&1; then fail 'built a source without its original'; fi
grep -q '^fetch: x.gz' "$work/log" || fail 'missing original was not fetched'
printf 'pass: tiers, pairs, interleave, gzip, CRLF, manifest, fresh, repair, new tier, underfill, mates, /3 /4 mates, fetch\n'

#!/usr/bin/env bash
# Catalogue-driven FASTQ generation. JSON is data; Bash manages files and workers.
set -euo pipefail
export LC_ALL=C
unset GZIP PIGZ
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$here/../../.." && pwd)"
die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}
usage() {
    cat <<'HELP'
Usage: make.sh [OPTIONS] SHAPE TIER [ROOT]
       make.sh --count N [--start N] [--length N] [OPTIONS] SHAPE LABEL [ROOT]
       make.sh --list [--catalog FILE]
Options:
  --plan                    Print the selected build as JSON without writes
  --catalog FILE            Use another validated recipe catalogue
  --root DIR                Output root (default: catalogue data_root)
  --formats LIST            plain,gzip,multi_gzip or all (default: catalogue)
  --paired-layout LAYOUT    mates, interleaved or both (default: mates)
  --jobs N                  Generation workers, 1..256 (default: 1)
  --compression-jobs N      pigz threads, 1..256 (default: generation workers)
  --gzip-level N            Compression level, 1..9 (default: catalogue value)
  --allow-large             Permit tiers larger than small
  --seed N                  Replicate content from a pool, 1..999999; files are named LABEL-sN
Count mode uses a custom LABEL instead of a tier and has no byte cap.
A recipe with a pool takes sequence and quality bytes from a real corpus tier
(bench/lib/corpus/make.sh builds it); file sizes match the same recipe without one.
HELP
}
integer() { [[ "$2" =~ ^[1-9][0-9]{0,8}$ ]] || die "$1 must be an integer in 1..999999999"; }
catalog="$here/synthetic.json" root='' mode=build formats="${MAKE_FORMAT:-}" layout=''
jobs="${MAKE_JOBS:-1}" compression_jobs='' level='' count='' start=1 length='' seed=0
allow_large="${ALLOW_BIG_TIERS:-0}" args=()
while (($#)); do
    case "$1" in
        --help | -h)
            usage
            exit 0
            ;;
        --plan | --list)
            [[ "$mode" == build ]] || die 'choose one of --plan and --list'
            mode="${1#--}"
            shift
            ;;
        --allow-large)
            allow_large=1
            shift
            ;;
        --catalog | --root | --formats | --paired-layout | --jobs | --compression-jobs | --gzip-level | --count | --start | --length | --seed)
            (($# >= 2)) && [[ -n "$2" ]] || die "missing value for $1"
            case "$1" in
                --catalog) catalog="$2" ;; --root) root="$2" ;; --formats) formats="$2" ;;
                --paired-layout) layout="$2" ;; --jobs) jobs="$2" ;; --compression-jobs) compression_jobs="$2" ;;
                --gzip-level) level="$2" ;; --count) count="$2" ;; --start) start="$2" ;; --length) length="$2" ;; --seed) seed="$2" ;;
            esac
            shift 2
            ;;
        --)
            shift
            args+=("$@")
            break
            ;;
        -*) die "unknown option: $1" ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done
# Reject unknown fields rather than silently accepting a misspelled setting.
config="$(jq -ces '
    def keys_are($allowed): type == "object" and ((keys - $allowed) | length == 0);
    def uint($max): type == "number" and . == floor and . >= 1 and . <= $max;
    def id: type == "string" and length <= 64 and test("^[a-z][a-z0-9_-]*\\z");
    def recipe:
        keys_are(["id","engine","layout","description","lengths","pool","sequence_unit","quality_unit","header","line_endings","nominal_record_bytes","multi_gzip_members","outliers"])
        and (.id | id) and (.description | type == "string")
        and (.engine | IN("fixed","variable","paired","longread"))
        and (.layout == (if .engine == "paired" then "paired" else "single" end))
        and (.lengths | type == "array" and length >= 1 and length <= 1024 and all(.[]; uint(999999999)))
        and (if .engine == "fixed" or .engine == "paired" then (.lengths | length == 1) else true end)
        and (if has("pool") then
            (has("sequence_unit") or has("quality_unit") | not)
            and (.pool | type == "string" and test("^[a-z0-9]+(-[a-z0-9]+)*/[A-Za-z0-9_-]+\\.fastq\\z"))
            and (if .engine == "paired" then (.pool | endswith("_R1.fastq")) else true end)
        else
            (.sequence_unit | type == "string" and length <= 65536 and test("^[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+\\z"))
            and (.quality_unit | type == "string" and length <= 65536 and test("^[!-~]+\\z"))
        end)
        and (.line_endings | IN("lf","alternating_crlf_lf"))
        and (.nominal_record_bytes | uint(9007199254740991))
        and (if has("multi_gzip_members") then .layout == "single" and (.multi_gzip_members | uint(256)) else true end)
        and (if .engine == "longread" then
            (.header | keys_are(["runid_offset","channels","barcodes","unclassified_every","start_time"])
                and (.runid_offset | uint(999999999)) and (.channels | uint(65536))
                and (.barcodes | uint(99)) and (.unclassified_every | uint(999999999))
                and (.start_time | type == "string" and length <= 64 and test("^[0-9T:+.-]+\\z")))
            and (.outliers | keys_are(["every","base_length","index_multiplier","modulus"])
                and (.every | uint(999999999)) and (.base_length | uint(999999999))
                and (.index_multiplier | uint(9007199)) and (.modulus | uint(999999999))
                and .base_length + .modulus - 1 <= 999999999)
        else
            (has("outliers") | not)
            and (.header | keys_are(["prefix","index_width","barcode"])
                and (.prefix | type == "string" and length <= 64 and test("^[A-Za-z][A-Za-z0-9_-]*\\z"))
                and (.index_width | uint(9)))
            and (if .engine == "paired" then (.header.barcode | type == "string" and length <= 64 and test("^[ACGTN]+\\z"))
                else (.header | has("barcode") | not) end)
        end);
    select(length == 1) | .[0] |
    select(keys_are(["schema_version","revision","generator_revision","data_root","pool_root","tiers","defaults","sizing","shapes"])
        and .schema_version == 1 and (.revision | uint(999999999)) and (.generator_revision | uint(999999999))
        and (.data_root | type == "string" and test("^[A-Za-z0-9_-]+(/[A-Za-z0-9_-]+)*\\z"))
        and (.pool_root | type == "string" and test("^[A-Za-z0-9_-]+(/[A-Za-z0-9_-]+)*\\z"))
        and (.tiers | keys == ["large","medium","small","tiny","xsmall"] and all(.[]; uint(1099511627776)))
        and .tiers.tiny < .tiers.xsmall and .tiers.xsmall < .tiers.small
        and .tiers.small < .tiers.medium and .tiers.medium < .tiers.large
        and (.defaults | keys_are(["formats","paired_layout","gzip_level"])
            and (.formats | type == "array" and length >= 1 and length == (unique | length) and all(.[]; IN("plain","gzip")))
            and (.paired_layout | IN("mates","interleaved","both")) and (.gzip_level | uint(9)))
        and .sizing == {budget:"decompressed_bytes",paired_budget:"sum_of_mates",count:"floor_cap_over_nominal_then_exact_prefix_fit",record_index_max:999999999}
        and (.shapes | type == "array" and length > 0 and all(.[]; recipe) and (map(.id) | length == (unique | length))))
' "$catalog")" || die "invalid synthetic catalogue: $catalog"
if [[ "$mode" == list ]]; then
    ((${#args[@]} == 0)) || die '--list takes no shape or tier'
    jq . <<<"$config"
    exit 0
fi
((${#args[@]} >= 2 && ${#args[@]} <= 3)) || {
    usage >&2
    exit 1
}
shape="${args[0]}" label="${args[1]}"
[[ "$label" =~ ^[a-z][a-z0-9_-]*$ && ${#label} -le 64 ]] || die "invalid tier or custom label: $label"
if ((${#args[@]} == 3)); then
    [[ -z "$root" ]] || die 'specify the root only once'
    root="${args[2]}"
fi
integer jobs "$jobs"
((jobs <= 256)) || die 'jobs must not exceed 256'
compression_jobs="${compression_jobs:-$jobs}"
integer compression-jobs "$compression_jobs"
((compression_jobs <= 256)) || die 'compression-jobs must not exceed 256'
integer start "$start"
[[ "$seed" =~ ^(0|[1-9][0-9]{0,5})$ ]] || die 'seed must be an integer in 0..999999'
[[ -z "$length" ]] || integer length "$length"
settings="$(jq -re --arg id "$shape" --arg tier "$label" --argjson len "${length:-0}" '
    (.shapes[] | select(.id == $id) | if $len > 0 then .lengths=[$len] else . end) as $r |
    [.data_root, (.defaults.formats | join(",")), .defaults.paired_layout, .defaults.gzip_level,
    $r.engine, ($r.lengths | map(tostring) | join(" ")), ($r.sequence_unit // ""), ($r.quality_unit // ""), $r.line_endings,
    ($r.header.prefix // ""), ($r.header.index_width // 1), ($r.header.barcode // ""),
    ($r.header.runid_offset // 0), ($r.header.channels // 0), ($r.header.barcodes // 0),
    ($r.header.unclassified_every // 0), ($r.header.start_time // ""), ($r.outliers.every // 0),
    ($r.outliers.base_length // 0), ($r.outliers.index_multiplier // 0), ($r.outliers.modulus // 0),
    $r.nominal_record_bytes, ($r.multi_gzip_members // 0), (.tiers[$tier] // 0), .tiers.small, ($r | tojson), ($r.pool // ""), .pool_root] | .[]' <<<"$config")" || die "unknown shape: $shape"
mapfile -t values <<<"$settings"
cap="${values[23]}" recipe="${values[25]}" pool="${values[26]}"
((seed == 0)) || [[ -n "$pool" ]] || die '--seed requires a recipe with a pool'
pools=() pool_bytes=null
if [[ -n "$pool" ]]; then
    pools=("$repo/${values[27]}/$pool")
    [[ "${values[4]}" != paired ]] || pools+=("${pools[0]%_R1.fastq}_R2.fastq")
    for file in "${pools[@]}"; do
        [[ -f "$file" ]] || die "missing pool $file; build it with: bash bench/lib/corpus/make.sh ${pool%%/*}"
    done
    pool_bytes="$(stat -c%s -- "${pools[@]}" | jq -cs .)"
fi
export S_POOL="${pools[0]:-}" S_POOL2="${pools[1]:-}" S_SEED="$seed"
[[ -z "$length" || "${values[4]}" == fixed || "${values[4]}" == paired ]] || die '--length requires a fixed or paired engine'
if [[ -n "$count" ]]; then
    integer count "$count"
    ((cap == 0)) || die 'use a custom label with --count, not a tier name'
    ((start + count - 1 <= 999999999)) || die 'final record index exceeds 999999999'
else
    ((cap > 0)) || die "unknown tier: $label"
    [[ "$start" == 1 && -z "$length" ]] || die '--start and --length require --count'
    if [[ "$mode" == build && "$allow_large" != 1 ]]; then
        ((cap <= values[24])) || die 'use --allow-large for tiers larger than small'
    fi
fi
((seed == 0)) || label+="-s$seed"
root="${root:-$repo/${values[0]}}" formats="${formats:-${values[1]}}" layout="${layout:-${values[2]}}" level="${level:-${values[3]}}"
[[ "$layout" == mates || "$layout" == interleaved || "$layout" == both ]] || die 'invalid paired layout'
[[ "$level" =~ ^[1-9]$ ]] || die 'gzip-level must be in 1..9'
export S_ENGINE="${values[4]}" S_LENGTHS="${values[5]}" S_SEQUENCE="${values[6]}" S_QUALITY="${values[7]}" S_ENDINGS="${values[8]}"
export S_PREFIX="${values[9]}" S_WIDTH="${values[10]}" S_BARCODE="${values[11]}" S_OFFSET="${values[12]}" S_CHANNELS="${values[13]}"
export S_BARCODES="${values[14]}" S_UNCLASSIFIED="${values[15]}" S_TIME="${values[16]}" S_EVERY="${values[17]}" S_BASE="${values[18]}"
export S_MULTIPLIER="${values[19]}" S_MODULUS="${values[20]}" S_NOMINAL="${values[21]}" S_START="$start" S_COUNT="${count:-0}" S_CAP="$cap" S_LAYOUT="$layout"
members="${values[22]}"
if [[ "$formats" == all ]]; then
    formats=plain,gzip
    ((members == 0)) || formats+=,multi_gzip
fi
[[ "$formats" =~ ^(plain|gzip|multi_gzip)(,(plain|gzip|multi_gzip))*$ ]] || die 'invalid formats'
IFS=, read -r -a formats <<<"$formats"
declare -A seen=()
for format in "${formats[@]}"; do
    [[ ! -v "seen[$format]" ]] || die "duplicate format: $format"
    seen[$format]=1
    [[ "$format" != multi_gzip ]] || ((members > 0)) || die 'this recipe has no multi_gzip_members'
done
planned="$(S_MODE=plan awk -f "$here/records.awk")"
read -r count expected_bytes <<<"$planned"
((count > 0 && expected_bytes > 0 && expected_bytes <= 9007199254740991)) || die 'no complete record fits, or size exceeds exact integer arithmetic'
((jobs <= count)) || jobs="$count"
((members <= count)) || members="$count"
export S_COUNT="$count"
roles=(single) suffixes=(.fastq) role_bytes=("$expected_bytes")
if [[ "$S_ENGINE" == paired ]]; then
    roles=() suffixes=() role_bytes=()
    if [[ "$layout" != interleaved ]]; then
        roles+=(R1 R2)
        suffixes+=(_R1.fastq _R2.fastq)
        role_bytes+=("$((expected_bytes / 2))" "$((expected_bytes / 2))")
    fi
    if [[ "$layout" != mates ]]; then
        roles+=(interleaved)
        suffixes+=(_interleaved.fastq)
        role_bytes+=("$expected_bytes")
    fi
fi
products=() product_formats=() product_roles=()
for ((r = 0; r < ${#roles[@]}; r++)); do
    for format in "${formats[@]}"; do
        name="$label${suffixes[r]}"
        case "$format" in gzip) name+=.gz ;; multi_gzip) name="$label.multi.fastq.gz" ;; esac
        products+=("$name")
        product_formats+=("$format")
        product_roles+=("$r")
    done
done
plan="$(jq -cn --arg shape "$shape" --arg label "$label" --arg root "$root" --argjson count "$count" --argjson start "$start" \
    --argjson cap "$cap" --argjson bytes "$expected_bytes" --argjson jobs "$jobs" --arg layout "$layout" \
    --arg unit "$([[ "$S_ENGINE" == paired ]] && printf pairs || printf records)" --args \
    '{shape:$shape,label:$label,root:$root,count:$count,unit:$unit,start:$start,byte_cap:$cap,decompressed_bytes:$bytes,jobs:$jobs,paired_layout:$layout,files:$ARGS.positional}' "${products[@]}")"
if [[ "$mode" == plan ]]; then
    jq . <<<"$plan"
    exit 0
fi
# Check every destination before mkdir, locking or generation; reject symlinks.
root="$(realpath -m -- "$root")" dest="$root/$shape"
[[ ! -L "$dest" && (! -e "$dest" || -d "$dest") ]] || die "output must be a directory: $dest"
for name in state.json "${products[@]}"; do
    [[ ! -L "$dest/$name" && (! -e "$dest/$name" || -f "$dest/$name") ]] || die "output must be a regular file: $dest/$name"
done
mkdir -p -- "$dest"
# Lock the directory inode so replacing state.json never changes the lock.
exec {lock_fd}<"$dest"
flock -x "$lock_fd"
state='{"schema_version":1,"entries":{}}'
if [[ -s "$dest/state.json" ]]; then
    state="$(jq -ce 'select(.schema_version == 1 and (.entries | type == "object"))' "$dest/state.json")" || state='{"schema_version":1,"entries":{}}'
fi
identity="$(jq -cS --argjson r "$recipe" --argjson cap "$cap" --argjson count "$count" --argjson start "$start" \
    --argjson pool_bytes "$pool_bytes" --argjson seed "$seed" \
    '{generator_revision,recipe:($r | del(.description,.nominal_record_bytes,.multi_gzip_members)),cap:$cap,count:$count,start:$start,pool_bytes:$pool_bytes,seed:$seed}' <<<"$config")"
compressor=() compression_identity=null
if [[ -v seen[gzip] || -v seen[multi_gzip] ]]; then
    # pigz 2.8 writes a different stream with one thread; two or more agree,
    # so cached gzip bytes do not depend on the worker count.
    if command -v pigz >/dev/null; then compressor=(pigz -p "$((compression_jobs < 2 ? 2 : compression_jobs))"); else compressor=(gzip); fi
    version="$("${compressor[0]}" --version 2>&1)"
    version="${version%%$'\n'*}"
    compression_identity="$(jq -cn --arg version "$version" --argjson level "$level" '{version:$version,level:$level}')"
    compressor+=(-n "-$level" -c)
fi
declare -A saved_ids=() saved_stats=()
saved="$(jq -r --arg label "$label" '.entries[$label].products // {} | to_entries[] |
    [.key,(.value.identity | if has("compression") then .compression |= del(.threads) else . end | tojson),.value.stat] | join("\t")' <<<"$state")"
while IFS=$'\t' read -r name id stamp; do
    [[ -n "$name" && "$name" =~ ^[a-z][A-Za-z0-9_.-]*$ ]] || continue
    saved_ids[$name]="$id"
    saved_stats[$name]="$stamp"
done <<<"$saved"
stamp() { [[ -f "$1" && ! -L "$1" && -s "$1" ]] && stat -c '%d:%i:%s:%y:%z' -- "$1"; }
fresh() {
    local current
    [[ "${saved_ids[$1]:-}" == "$2" ]] || return 1
    current="$(stamp "$dest/$1")" || return 1
    [[ "$current" == "${saved_stats[$1]:-}" ]]
}
plain_id="{\"data\":$identity,\"format\":\"plain\"}"
ids=() needed=() sources=() staged=() all_fresh=true
for ((i = 0; i < ${#products[@]}; i++)); do
    format="${product_formats[i]}"
    id="$plain_id"
    if [[ "$format" != plain ]]; then
        id="{\"data\":$identity,\"format\":\"$format\",\"compression\":$compression_identity"
        [[ "$format" != multi_gzip ]] || id+=",\"members\":$members"
        id+='}'
    fi
    ids+=("$id")
    staged+=("")
    if fresh "${products[i]}" "$id"; then needed+=(false); else
        needed+=(true)
        all_fresh=false
    fi
done
if $all_fresh; then
    printf 'fresh: %s/%s (%s decompressed bytes)\n' "$shape" "$label" "$expected_bytes"
    exit 0
fi
work="$(mktemp -d /var/tmp/fastq-shapes.XXXXXX)"
files=() pids=()
cleanup() {
    local rc=$? pid file
    trap - EXIT HUP INT TERM
    for pid in "${pids[@]}"; do kill -- "-$pid" 2>/dev/null || true; done
    for pid in "${pids[@]}"; do wait "$pid" 2>/dev/null || true; done
    for file in "${files[@]}"; do [[ ! -e "$file" && ! -L "$file" ]] || rm -- "$file" || rc=1; done
    rmdir -- "$work" || rc=1
    exit "$rc"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
run() (
    exec {lock_fd}>&-
    exec setsid "$@"
)
wait_workers() {
    local pid failed=0
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    ((failed == 0)) || die 'worker failed; completed outputs were preserved'
    pids=()
}
# Reuse verified plain inputs for added or damaged compressed formats.
generate=false
for ((r = 0; r < ${#roles[@]}; r++)); do
    sources+=("$dest/$label${suffixes[r]}")
    for ((i = 0; i < ${#products[@]}; i++)); do
        if [[ "${product_roles[i]}" == "$r" ]] && ${needed[i]} && ! fresh "$label${suffixes[r]}" "$plain_id"; then generate=true; fi
    done
done
if $generate; then
    base=$((count / jobs)) rem=$((count % jobs)) next="$start"
    for ((k = 0; k < jobs; k++)); do
        n=$base
        ((k >= rem)) || n=$((n + 1))
        for suffix in "${suffixes[@]}"; do files+=("$work/part-$k$suffix"); done
        output="$work/part-$k"
        [[ "$S_ENGINE" == paired ]] || output+=.fastq
        stdout="$output"
        [[ "$S_ENGINE" != paired ]] || stdout=/dev/null
        S_MODE=generate S_COUNT="$n" S_START="$next" S_OUT="$output" run awk -f "$here/records.awk" >"$stdout" &
        pids+=($!)
        next=$((next + n))
    done
    wait_workers
    for ((r = 0; r < ${#roles[@]}; r++)); do
        source="$work/part-0${suffixes[r]}"
        for ((k = 1; k < jobs; k++)); do
            part="$work/part-$k${suffixes[r]}"
            cat -- "$part" >>"$source"
            rm -- "$part"
        done
        [[ "$(stat -c%s -- "$source")" == "${role_bytes[r]}" ]] || die "unexpected byte count for ${roles[r]}"
        sources[r]="$source"
    done
fi
for ((i = 0; i < ${#products[@]}; i++)); do
    ${needed[i]} || continue
    r="${product_roles[i]}" source="${sources[${product_roles[i]}]}" format="${product_formats[i]}"
    if [[ "$format" == plain ]]; then
        staged[i]="$source"
        continue
    fi
    target="$work/product-$i"
    files+=("$target")
    staged[i]="$target"
    if [[ "$format" == gzip ]]; then
        run "${compressor[@]}" -- "$source" >"$target" &
    else
        # GNU head restores the offset after read-ahead on regular files.
        # Each member contains whole records; no second plain copy is made.
        run bash -o pipefail -c '
            source=$1 count=$2 members=$3; shift 3
            exec 3<"$source"
            for ((m=0; m<members; m++)); do
                n=$((count/members + (m < count%members)))
                head -n "$((n*4))" <&3 | "$@" || exit
            done
        ' bash "$source" "$count" "$members" "${compressor[@]}" >"$target" &
    fi
    pids=($!)
    wait_workers
done
# All workers succeeded. File stats make interrupted publication a cache miss.
# Moving across filesystems may copy; consumers must wait for successful exit.
for ((i = 0; i < ${#products[@]}; i++)); do
    [[ -z "${staged[i]}" ]] || mv -f -- "${staged[i]}" "$dest/${products[i]}"
done
files+=("$work/products.jsonl" "$work/state.json")
for ((i = 0; i < ${#products[@]}; i++)); do
    name="${products[i]}"
    jq -cn --arg name "$name" --argjson identity "${ids[i]}" --arg stat "$(stamp "$dest/$name")" \
        --argjson bytes "$(stat -c%s -- "$dest/$name")" '{key:$name,value:{identity:$identity,stat:$stat,bytes:$bytes}}'
done >"$work/products.jsonl"
jq --arg shape "$shape" --arg label "$label" --argjson plan "$plan" --slurpfile products "$work/products.jsonl" \
    '.shape=$shape | .entries[$label] = (($plan | del(.root,.files,.jobs)) + {products:((.entries[$label].products // {}) + ($products | from_entries))})' \
    <<<"$state" >"$work/state.json"
mv -f -- "$work/state.json" "$dest/state.json"
printf 'made: %s/%s (%s decompressed bytes)\n' "$shape" "$label" "$expected_bytes"

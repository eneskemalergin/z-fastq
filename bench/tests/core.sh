#!/usr/bin/env bash
# Check the benchmark core without timing anything: syntax, help, declaration
# validation, the batch-reuse rule, the lane identity, case parameters, cores
# per lane, output directories, and the per-CPU timing lock.
# Runs in seconds on the corpus already cut under bench/data.
set -euo pipefail
export LC_ALL=C
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bench="$(dirname -- "$here")"
work="$(mktemp -d /var/tmp/bench-core.XXXXXX)"
cleanup() {
    rm -f -- "$work"/feat/results/batches/*/* "$work"/feat/*.sh "$work"/*.log "$work"/lib "$work"/data "$work"/cache
    rmdir -- "$work"/feat/results/batches/* "$work"/feat/results/batches "$work"/feat/results "$work"/feat "$work" 2>/dev/null || true
}
trap cleanup EXIT
fail() {
    printf 'fail: %s\n' "$*" >&2
    exit 1
}

# Syntax of every script and of plot.py.
for f in "$bench"/lib/bench.sh "$bench"/*/run.sh; do bash -n "$f" || fail "syntax: $f"; done
python3 -I -c 'import ast, sys; ast.parse(open(sys.argv[1]).read())' "$bench/lib/plot.py" || fail 'syntax: plot.py'

# A scratch feature over the real core, data, and cache.
mkdir -p -- "$work/feat"
ln -s -- "$bench/lib" "$work/lib"
ln -s -- "$bench/data" "$work/data"
ln -s -- "$bench/cache" "$work/cache"
feature() { # extra declaration lines -> feat/run.sh
    cat >"$work/feat/run.sh" <<EOF
source "\$(dirname -- "\$0")/../lib/bench.sh"
job "records"
truth "manifest counts"
headline small
scaling ""
lane z-fastq subject "{zfastq} count {in}"
lane seqtk partial "{bin} size {in}"
family file "one file"
lanes all
input real illumina-drr1001190 small interleaved
$1
fact() { cat; }
[[ "\${1:-}" == __probe ]] && return 0
bench_main "\$@"
EOF
}

# Help and an unknown mode.
feature ''
bash "$work/feat/run.sh" help | grep -q '^  peer ' || fail 'help lacks peer'
bash "$work/feat/run.sh" nonsense >/dev/null 2>&1 && fail 'unknown mode accepted'
bash "$work/feat/run.sh" list >"$work/list.log" || fail 'list failed on a valid feature'
[[ "$(wc -l <"$work/list.log")" == 2 ]] || fail 'list: expected plain and gzip cases'

# Validation names each mistake and stops before work.
expect_invalid() { # declarations, message part
    feature "$1"
    if bash "$work/feat/run.sh" list >/dev/null 2>"$work/err.log"; then fail "accepted: $1"; fi
    grep -qF -- "$2" "$work/err.log" || fail "message lacks '$2': $(cat "$work/err.log")"
}
expect_invalid "unset 'DESCRIBE[illumina-drr147610]'
input real illumina-drr147610 small" 'source illumina-drr147610 has no describe line'
expect_invalid "unset 'GROUP_OF[illumina-drr147610]'
input real illumina-drr147610 small" 'source illumina-drr147610 is in no group'
expect_invalid 'pick only nosuch:small:plain' 'pick: no case nosuch:small:plain'
expect_invalid 'scaling illumina-drr1001190' 'scaling source illumina-drr1001190 is declared at 1 size'
expect_invalid 'remark summary "x"' 'remark: unknown section summary'

# The batch-reuse rule: a cached batch stands for a planned one when form,
# inputs, runs, CPU, and zexec match, it has at least as many warmups, and
# every planned lane is in it with the same identity.
feature ''
# shellcheck disable=SC2034,SC2154  # read by the sourced core's functions
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    ZEXEC_ID=Z
    LANE_IDS=(["z-fastq|"]=A ["seqtk|"]=B)
    REG_PATH[seqtk]=/usr/bin/true
    pre="$work/feat/results/batches/c-"
    stored() { # name warm lanes... (identity lines as batch_meta writes them)
        batch_meta file IN 11 "$2" "${@:3}"
        mkdir -p -- "$pre$1"
        printf '%s' "$META" >"$pre$1/meta.tsv"
        : >"$pre$1/summary.tsv"
    }
    stored both 3 z-fastq seqtk
    batch_meta file IN 11 1 z-fastq
    find_batch "$pre" && [[ "$BATCH" == "${pre}both" ]] || fail 'a superset batch with more warmups was not reused'
    batch_meta file IN 11 5 z-fastq
    find_batch "$pre" && fail 'a batch with fewer warmups was reused'
    batch_meta file IN 21 1 z-fastq
    find_batch "$pre" && fail 'a batch with other runs was reused'
    LANE_IDS["z-fastq|"]=A2
    batch_meta file IN 11 1 z-fastq
    find_batch "$pre" && fail 'a batch with another z-fastq identity was reused'
    LANE_IDS["z-fastq|"]=A
    batch_meta stdin IN 11 1 z-fastq
    find_batch "$pre" && fail 'a file batch was reused for stdin'
    # A lane's identity covers its command as well as its binary.
    a="$(lane_id_now seqtk)"
    LANE_TMPL[seqtk]="{bin} size -x {in}"
    [[ "$(lane_id_now seqtk)" != "$a" ]] || fail 'changing a command kept the lane identity'
    exit 0
) || fail 'batch-reuse checks'

# Case parameters: a short:A,long:B value follows the source's read class, a
# varied parameter names its value in the case id, every placeholder expands,
# and two parameter sets never share a check key.
feature 'family k "by count"
lanes all
param k short:100,long:7
param seed 11
input real illumina-drr1001190,ont-err14758576 small interleaved
family kv "by K"
lanes all
vary k 5,50
input real illumina-drr1001190 small interleaved'
# shellcheck disable=SC2034,SC2154  # read by the sourced core's functions
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    list="$(cases)"
    [[ "$(awk -F'\t' '$1 == "k" && $7 == "illumina-drr1001190" { print $10; exit }' <<<"$list")" == "k=100;seed=11" ]] || fail 'a short-read source did not get the short value'
    [[ "$(awk -F'\t' '$1 == "k" && $7 == "ont-err14758576" { print $10; exit }' <<<"$list")" == "k=7;seed=11" ]] || fail 'a long-read source did not get the long value'
    [[ "$(awk -F'\t' '$1 == "kv" { print $2 }' <<<"$list" | head -2 | tr '\n' ' ')" == "illumina-drr1001190:small:k=5:plain illumina-drr1001190:small:k=50:plain " ]] ||
        fail "a varied parameter is not in the case id: $(awk -F'\t' '$1 == "kv" { print $2 }' <<<"$list" | head -2)"
    [[ "$(expand '{bin} -n {k} -s {seed} {in} {out1}' seqtk IN 'k=3;seed=9')" == "${REG_PATH[seqtk]} -n 3 -s 9 IN /dev/null" ]] || fail 'expand left a placeholder or used the wrong value'
    [[ "$(OUT1=/x expand '{out1}' seqtk IN -)" == /x ]] || fail 'expand ignored OUT1'
    LANE_IDS=(["seqtk|k"]=B)
    [[ "$(status_key seqtk k file F -)" == "seqtk|B|file|F" ]] || fail 'a case without parameters changed its check key'
    [[ "$(status_key seqtk k file F "k=5")" != "$(status_key seqtk k file F "k=50")" ]] || fail 'two parameter sets share a check key'
    exit 0
) || fail 'case-parameter checks'
expect_invalid 'command seqtk "{bin} size -n {nope} {in}"' 'uses {nope}, which no param line declares'

# Cores per lane: CPUs on distinct physical cores, the CPU set in the batch
# identity, and the declaration checks.
feature ''
# shellcheck disable=SC2034,SC2154  # read by the sourced core's functions
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    set2="$(cpu_set 2)"
    [[ "$(tr ',' '\n' <<<"$set2" | wc -l)" == 2 ]] || fail "cpu_set 2 gave $set2"
    for c in ${set2//,/ }; do
        first="$(<"/sys/devices/system/cpu/cpu$c/topology/thread_siblings_list")"
        [[ "${first%%[,-]*}" == "$c" ]] || fail "cpu_set picked CPU $c, a second thread of its core"
    done
    [[ -z "$(cpu_set 100000)" ]] || fail 'cpu_set gave CPUs for more cores than the machine has'
    [[ "$(cpu_set 1 count)" -ge 2 ]] || fail 'cpu_set count is below 2'
    ZEXEC_ID=Z BATCH_CPUS="$set2"
    batch_meta file IN 3 1
    [[ "$META" == *$'\n'"cpu $set2"$'\n'* ]] || fail 'the CPU set is not in the batch identity'
    exit 0
) || fail 'cores checks'
expect_invalid 'lane fastp partial "{bin} -w 2 {in}" cores=two' 'cores=two is not a positive count'
expect_invalid 'lane fastp partial "{bin} -w 2 {in}" of=nosuch' 'of=nosuch names no declared lane'
expect_invalid 'lane fastp partial "{bin} -o {stem}.d {in}" outdir=.d' 'defines no outdir_entries'

# An output directory goes file by file: its named entries, then the
# directory; a file nobody named stops the removal and stays.
feature 'lane fastp partial "{bin} -o {stem}.d {in}" outdir=.d
facts piped
outdir_entries() { printf "a\\nsub/b\\nsub/\\n"; }'
# shellcheck disable=SC2034,SC2154  # read by the sourced core's functions
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    [[ "$(outdir_of fastp file)" == .d && -z "$(outdir_of z-fastq file)" ]] || fail 'outdir_of named the wrong lanes'
    d="$work/od"
    mkdir -p -- "$d/sub"
    : >"$d/a"; : >"$d/sub/b"
    clear_outdir "$d" fastp file IN || fail 'clear_outdir failed on named entries'
    [[ ! -e "$d" ]] || fail 'clear_outdir left the directory'
    mkdir -- "$d"
    : >"$d/a"; : >"$d/extra"
    clear_outdir "$d" fastp file IN && fail 'clear_outdir removed a directory holding an unnamed file'
    [[ -f "$d/extra" ]] || fail 'clear_outdir removed an unnamed file'
    rm -f -- "$d/extra"
    rmdir -- "$d"
    exit 0
) || fail 'output directory checks'

# Only one timing run at a time on a CPU.
cpu="${BENCH_CPU:-3}"
exec 8>"$bench/cache/timing-cpu$cpu.lock"
flock -n 8 || fail "the timing lock of CPU $cpu is already held; is a benchmark running?"
feature ''
BENCH_CPU="$cpu" bash "$work/feat/run.sh" peer >/dev/null 2>"$work/lock.log" && fail 'a timing mode ran while the lock was held'
grep -q "timing-cpu$cpu.lock" "$work/lock.log" || fail "lock message missing: $(cat "$work/lock.log")"
exec 8>&-

printf 'bench core: ok\n'

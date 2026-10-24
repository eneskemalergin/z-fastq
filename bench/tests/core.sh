#!/usr/bin/env bash
# Check the benchmark core without timing anything: syntax, help, declaration
# validation, the batch-reuse rule, the lane identity, and the per-CPU timing lock.
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

# Only one timing run at a time on a CPU.
cpu="${BENCH_CPU:-3}"
exec 8>"$bench/cache/timing-cpu$cpu.lock"
flock -n 8 || fail "the timing lock of CPU $cpu is already held; is a benchmark running?"
feature ''
BENCH_CPU="$cpu" bash "$work/feat/run.sh" peer >/dev/null 2>"$work/lock.log" && fail 'a timing mode ran while the lock was held'
grep -q "timing-cpu$cpu.lock" "$work/lock.log" || fail "lock message missing: $(cat "$work/lock.log")"
exec 8>&-

printf 'bench core: ok\n'

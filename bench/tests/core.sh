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
    local name
    if [[ -f "$work/restarts.log" ]]; then
        while IFS= read -r name; do
            [[ "$name" == "$work/feat/results/batches/"* ]] || continue
            rm -f -- "$name/meta.tsv" "$name/summary.tsv" "$name/partial.tsv"
            rmdir -- "$name" 2>/dev/null || true
        done <"$work/restarts.log"
    fi
    for name in c-both file-test-current; do
        rm -f -- "$work/feat/results/batches/$name/meta.tsv" "$work/feat/results/batches/$name/summary.tsv"
        rmdir -- "$work/feat/results/batches/$name" 2>/dev/null || true
    done
    rm -f -- "$work/feat/run.sh" "$work/list.log" "$work/err.log" "$work/lock.log" \
        "$work/checks.tsv" "$work/collect.log" "$work/build.log" "$work/producer.sh" "$work/input" "$work/subject" "$work/external" \
        "$work/restarts.log" "$work/worker.log" "$work/reference.tsv" "$work/reference.tsv.part" \
        "$work/strict.log" "$work/reference-error.log" \
        "$work/feat/results/reference.key" "$work/feat/results/current.tsv.part" \
        "$work/od/a" "$work/od/sub/b" "$work/od/extra" "$work/lib" "$work/data" "$work/cache"
    rmdir -- "$work/od/sub" "$work/od" "$work/feat/results/batches" "$work/feat/results" "$work/feat" "$work" 2>/dev/null || true
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
    CHECK_ID=C
    [[ "$(status_key seqtk k file F -)" == "seqtk|B|file|F|check=C" ]] || fail 'the verification identity is missing from the check key'
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

# Producer and checker exits are independent, for captured stdout, a pipe,
# a FIFO, and named output files.
feature ''
printf 'input\n' >"$work/input"
cat >"$work/producer.sh" <<'EOF'
#!/usr/bin/env bash
if (($#)); then printf 'read\n' >"$1"; else printf 'read\n'; fi
EOF
# shellcheck disable=SC2034,SC2154  # Variables consumed by the sourced core.
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    REG_PATH[z-fastq]=/usr/bin/true
    LANE_IDS=(["z-fastq|file"]=A)
    CHECK_ID=C EXPECT=accept AGREE_LANE=""
    fact() {
        if [[ -n "${OUT1:-}" ]]; then cat -- "$OUT1" >/dev/null; else cat >/dev/null; fi
        printf 'accept\n'
        return "$checker_exit"
    }
    for form in captured pipe fifo named; do
        FACT_PIPE=1
        LANE_OUTS[z-fastq]=""
        LANE_TMPL[z-fastq]="/usr/bin/bash $work/producer.sh"
        case "$form" in
            captured) FACT_PIPE=0 ;;
            fifo) LANE_TMPL[z-fastq]+=" {out1}" ;;
            named) LANE_TMPL[z-fastq]+=" {stem}.fq"; LANE_OUTS[z-fastq]=.fq ;;
        esac
        for checker_exit in 0 7; do
            check_exit=0
            row="$(check_one z-fastq file file test plain "$work/input" accept 6 - accept)" || check_exit=$?
            fields="$(cut -f2-4 <<<"$row")"
            if ((checker_exit == 0)); then
                [[ "$fields" == $'pass\t0\taccept' ]] || fail "$form checker success: $fields"
                ((check_exit == 0)) || fail "$form successful checker failed its worker"
            else
                [[ "$fields" == $'fail\t0\tchecker failed (exit 7)' ]] || fail "$form checker failure: $fields"
                ((check_exit != 0)) || fail "$form failed checker left a successful worker"
            fi
        done
    done
) || fail 'checker exit checks'

# A checker failure must remain a failed cache row and stop dispatch.
if bash -s -- "$work/feat/run.sh" "$work/input" "$work/checks.tsv" > /dev/null 2>"$work/strict.log" <<'EOF'
set -euo pipefail
BASH_ARGV0="$1"
source "$1" __probe
LANE_TMPL[z-fastq]='/usr/bin/printf ignored'
REG_PATH[z-fastq]=/usr/bin/printf
LANE_IDS=(["z-fastq|file"]=A)
CHECK_ID=C EXPECT=accept FACT_PIPE=1 CHECKS="$3"
fact() { cat >/dev/null; printf 'accept\n'; return 7; }
printf 'z-fastq\tfile\tfile\ttest\tplain\t%s\taccept\t6\t-\taccept\n' "$2" | dispatch_checks
EOF
then fail 'dispatch accepted a failed checker'; fi
[[ "$(cut -f2-4 "$work/checks.tsv")" == $'fail\t0\tchecker failed (exit 7)' ]] || fail 'a failed checker produced a passing or missing cache row'

# A command failure inside fact must stop it before a later success print.
if bash -s -- "$work/feat/run.sh" "$work/input" >"$work/strict.log" <<'EOF'
set -euo pipefail
BASH_ARGV0="$1"
source "$1" __probe
LANE_TMPL[z-fastq]='/usr/bin/printf ignored'
REG_PATH[z-fastq]=/usr/bin/printf
LANE_IDS=(["z-fastq|file"]=A)
CHECK_ID=C EXPECT=accept FACT_PIPE=1
fact() { cat >/dev/null; false; printf 'accept\n'; }
check_one z-fastq file file test plain "$2" accept 6 - accept
EOF
then fail 'checker continued after a failed command'; fi
[[ "$(cut -f2-4 "$work/strict.log")" == $'fail\t0\tchecker failed (exit 1)' ]] || fail 'checker failure did not reach the result row'

# Cached verdicts follow the checker and expected results. Reports keep only
# matching timings after a successful check under the current rules.
# shellcheck disable=SC2034,SC2154
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    LANES=(z-fastq) FAMILIES=(file) SCALING="" VERIFY=0
    FAM_LANES[file]=z-fastq
    REG_PATH[z-fastq]=/usr/bin/true
    ZEXEC=/usr/bin/true CHECKS="$work/checks.tsv"
    expected=accept case_parameters=-
    all_cases() { printf 'file\ttest\tplain\t%s\t%s\t6\tsource\tsmall\treal\t%s\n' "$work/input" "$expected" "$case_parameters"; }
    cases() { all_cases; }
    verify_cases() { :; }
    current_key() { status_key z-fastq file file "$(fileid "$work/input")" "$case_parameters"; }
    save_check() {
        CHECK_STATUS=() CHECK_VALUE=()
        cache_lane_ids
        printf '%s\tpass\t0\taccept\t%s\n' "$(current_key)" "$work/input" >"$CHECKS"
        load_checks
    }
    save_check
    [[ -z "$(check_jobs lanes)" ]] || fail 'an unchanged passing check was queued again'
    expected=reject
    cache_lane_ids
    [[ -n "$(check_jobs lanes)" ]] || fail 'changed expected output reused an old verdict'
    expected=accept
    save_check
    fact() { printf 'a changed checker\n'; }
    cache_lane_ids
    [[ -n "$(check_jobs lanes)" ]] || fail 'changed checker reused an old verdict'
    REFERENCE="$work/reference.tsv"
    printf 'accept\n' >"$REFERENCE"
    touch -d @1700000000 "$REFERENCE"
    save_check
    reference_identity="$(fileid "$REFERENCE")"
    reference_lane="$(lane_id z-fastq file)"
    [[ -z "$(check_jobs lanes)" ]] || fail 'unchanged reference results were not reused'
    printf 'reject\n' >"$REFERENCE"
    touch -d @1700000000 "$REFERENCE"
    [[ "$(fileid "$REFERENCE")" == "$reference_identity" ]] || fail 'reference fixture metadata changed'
    cache_lane_ids
    [[ -n "$(check_jobs lanes)" ]] || fail 'changed reference results reused an old verdict'
    [[ "$(lane_id z-fastq file)" == "$reference_lane" ]] || fail 'reference results changed the timing identity'
    REFERENCE=""
    save_check
    edge_expectation=accept
    edge_list() { printf 'probe\t%s\ttest\n' "$edge_expectation"; }
    edge_path() { printf '%s\n' "$work/input"; }
    edge_before="$(edge_key)"
    edge_expectation=reject
    [[ "$(edge_key)" != "$edge_before" ]] || fail 'changed edge expectations kept the old key'
    ZEXEC_ID="$(fileid "$ZEXEC")"
    plan_batch file "$work/input" small real "$case_parameters"
    dir="$RESULTS/batches/file-test-current"
    mkdir -p -- "$dir"
    printf '%s' "$META" >"$dir/meta.tsv"
    printf 'file\ttest\tz-fastq\tsubject\t12345\n' >"$dir/summary.tsv"
    collect >"$work/collect.log" || fail 'report rejected a verified matching batch'
    [[ "$(sed -n '2p' "$work/collect.log")" == $'file\ttest\tz-fastq\tsubject\t12345' ]] || fail 'report lost its matching timing row'
    for change in input binary command params checker; do
        saved_command="${LANE_TMPL[z-fastq]}"
        saved_binary="${REG_PATH[z-fastq]}"
        case "$change" in
            input) printf 'changed input\n' >"$work/input" ;;
            binary) REG_PATH[z-fastq]=/usr/bin/false ;;
            command) LANE_TMPL[z-fastq]+=" --changed" ;;
            params) case_parameters=k=2 ;;
            checker) fact() { printf 'another checker\n'; } ;;
        esac
        if [[ "$change" != checker ]]; then save_check; fi
        collect >"$work/collect.log" 2>/dev/null && fail "report accepted stale $change data"
        [[ "$(wc -l <"$work/collect.log")" == 1 ]] || fail "report emitted a stale $change row"
        printf 'input\n' >"$work/input"
        LANE_TMPL[z-fastq]="$saved_command" case_parameters=-
        REG_PATH[z-fastq]="$saved_binary"
        save_check
        plan_batch file "$work/input" small real "$case_parameters"
        printf '%s' "$META" >"$dir/meta.tsv"
        collect >"$work/collect.log" || fail "report did not recover after $change was checked"
    done
    touch -d @1700000000.100000000 "$work/input"
    save_check
    plan_batch file "$work/input" small real -
    printf '%s' "$META" >"$dir/meta.tsv"
    printf 'other\n' >"$work/input"
    touch -d @1700000000.200000000 "$work/input"
    save_check
    collect >"$work/collect.log" 2>/dev/null && fail 'report reused timing after a same-size edit within one second'
    [[ "$(cases | with_input_ids | cut -f11)" == "$(fileid "$work/input")" ]] || fail 'batched and individual input identities differ'
    # Forced repeats get new directories, preserving every earlier sample.
    save_check
    fail_timing=0
    time_lanes() {
        mkdir -p -- "$1"
        printf '%s\n' "$1" >>"$work/restarts.log"
        printf '%s' "$META" >"$1/meta.tsv"
        if ((fail_timing)); then
            printf 'incomplete sample\n' >"$1/partial.tsv"
            exit 7
        fi
        printf 'fresh\n' >"$1/summary.tsv"
    }
    RETIME=1
    for repeat in 1 2; do
        run_batch file test plain "$work/input" accept 6 source small real -
    done
    [[ "$(sort -u "$work/restarts.log" | wc -l)" == 2 ]] || fail 'forced repeats share one directory'
    [[ "$(cat "$dir/summary.tsv")" == $'file\ttest\tz-fastq\tsubject\t12345' ]] || fail 'forced repeat overwrote the earlier measurements'
    mapfile -t completed <"$work/restarts.log"
    fail_timing=1
    rerun_status=0
    (run_batch file test plain "$work/input" accept 6 source small real -) || rerun_status=$?
    [[ "$rerun_status" == 7 ]] || fail 'failed timing did not stop its run'
    incomplete="$(tail -n 1 "$work/restarts.log")"
    [[ -f "$incomplete/partial.tsv" && ! -e "$incomplete/summary.tsv" ]] || fail 'failed rerun published a completed summary'
    for finished in "${completed[@]}"; do
        [[ "$finished" != "$incomplete" && "$(cat "$finished/summary.tsv")" == fresh ]] || fail 'failed rerun damaged a completed batch'
    done
    unset RETIME
    run_batch file test plain "$work/input" accept 6 source small real -
    [[ "$BATCH" == "${completed[0]}" || "$BATCH" == "${completed[1]}" ]] || fail 'cache selected an incomplete rerun'
    [[ "$(wc -l <"$work/restarts.log")" == 3 ]] || fail 'cache reuse unexpectedly measured again'
    : >"$CHECKS"
    CHECK_STATUS=() CHECK_VALUE=()
    collect >"$work/collect.log" 2>/dev/null && fail 'report accepted a batch without verification'
    exit 0
) || fail 'verdict and report cache checks'

# Source timestamps do not identify the installed optimization mode.
# shellcheck disable=SC2034
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    REPO="$work"
    ZFASTQ_DEFAULT="$work/subject" ZFASTQ="$work/subject"
    compiler_exit=0
    zig() {
        [[ "$PWD" == "$REPO" ]] || return 9
        printf '%s\n' "$*" >>"$work/build.log"
        [[ "$*" == 'build -Dstatic=true -Doptimize=ReleaseFast' ]] || return 8
        ((compiler_exit == 0)) || return "$compiler_exit"
        printf 'ReleaseFast\n' >"$ZFASTQ"
    }
    for mode in Debug ReleaseSafe ReleaseFast absent; do
        if [[ "$mode" == absent ]]; then
            rm -- "$ZFASTQ"
        else
            printf '%s\n' "$mode" >"$ZFASTQ"
            touch -d @2000000000 "$ZFASTQ"
        fi
        build_subject || fail "build failed for $mode fixture"
        [[ "$(cat "$ZFASTQ")" == ReleaseFast ]] || fail "$mode subject was reused"
    done
    [[ "$(wc -l <"$work/build.log")" == 4 ]] || fail 'a default subject bypassed the build check'
    compiler_exit=7
    printf 'Debug\n' >"$ZFASTQ"
    status=0
    build_subject || status=$?
    [[ "$status" == 7 && "$(cat "$ZFASTQ")" == Debug ]] || fail 'compiler failure was hidden'
    ZFASTQ="$work/external"
    printf 'external\n' >"$ZFASTQ"
    build_subject
    [[ "$(wc -l <"$work/build.log")" == 5 && "$(cat "$ZFASTQ")" == external ]] || fail 'an external subject was rebuilt'
) || fail 'subject build checks'

# Failed reference jobs cannot replace or certify the previously saved data.
# shellcheck disable=SC2034
(
    BASH_ARGV0="$work/feat/run.sh"
    source "$work/feat/run.sh" __probe
    BENCH_JOBS=2
    jobs() { printf 'first\nmiddle\nlater\nlast\n'; }
    failed_job=middle
    worker() {
        if [[ "$1" == "$failed_job" ]]; then printf 'partial %s\n' "$1"; return 7; fi
        printf '%s\n' "$1"
    }
    jobs | pool worker >"$work/worker.log" && fail 'pool discarded a worker failure'
    printf 'previous\n' >"$work/reference.tsv"
    printf 'old-key\n' >"$RESULTS/reference.key"
    for failed_job in first middle last; do
        cached_step reference new-key jobs worker "$work/reference.tsv" && fail "failed $failed_job reference job was cached"
        [[ "$(cat "$work/reference.tsv")" == previous ]] || fail 'failed reference jobs replaced the earlier output'
        [[ "$(cat "$RESULTS/reference.key")" == old-key ]] || fail 'failed reference jobs replaced the earlier key'
        grep -q '^partial ' "$work/reference.tsv.part" || fail 'the failure fixture did not produce partial output'
    done
    failed_job=""
    cached_step reference new-key jobs worker "$work/reference.tsv" || fail 'successful reference jobs were rejected'
    [[ "$(wc -l <"$work/reference.tsv")" == 4 ]] || fail 'successful reference jobs lost output'
    [[ "$(cat "$RESULTS/reference.key")" == new-key ]] || fail 'successful reference jobs did not update the key'
    jobs() { return 9; }
    worker() { return 9; }
    cached_step reference new-key jobs worker "$work/reference.tsv" || fail 'unchanged successful reference results were not reused'
    [[ "$(wc -l <"$work/reference.tsv")" == 4 ]] || fail 'cache reuse changed the successful results'
) || fail 'reference worker checks'

# Process reference producers must report failure to the shared worker pool.
(
    BASH_ARGV0="$bench/process/run.sh"
    # shellcheck disable=SC1090  # Load declarations without invoking bench_main.
    source <(sed '$d' "$BASH_ARGV0")
    lane_cmd() { printf '/usr/bin/false\n'; }
    for reference_family in trim-truseq qc-tables; do
        reference_one "$reference_family"$'\t'"$work/input $work/input" 2>"$work/reference-error.log" && fail 'failed reference producer returned success'
    done
    exit 0
) || fail 'process reference producer checks'

# Only one timing run at a time on a CPU.
cpu="${BENCH_CPU:-3}"
exec 8>"$bench/cache/timing-cpu$cpu.lock"
flock -n 8 || fail "the timing lock of CPU $cpu is already held; is a benchmark running?"
feature ''
BENCH_CPU="$cpu" bash "$work/feat/run.sh" peer >/dev/null 2>"$work/lock.log" && fail 'a timing mode ran while the lock was held'
grep -q "timing-cpu$cpu.lock" "$work/lock.log" || fail "lock message missing: $(cat "$work/lock.log")"
exec 8>&-

printf 'bench core: ok\n'

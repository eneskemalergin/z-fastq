#!/usr/bin/env bash
# Check zexec's redirects, per-run token, argument and environment passing,
# exit status, and failure messages. Builds bench/cache/zexec when stale.
set -euo pipefail
export LC_ALL=C
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bench="$(dirname -- "$here")"
src="$bench/lib/zexec.zig" x="$bench/cache/zexec"
if [[ ! -x "$x" || "$src" -nt "$x" ]]; then
    mkdir -p -- "$bench/cache"
    zig build-exe "$src" -OReleaseSmall -fstrip -fsingle-threaded -target x86_64-linux-none -femit-bin="$x"
fi
work="$(mktemp -d /var/tmp/zexec-check.XXXXXX)"
cleanup() {
    rm -f -- "$work"/*
    rmdir -- "$work"
}
trap cleanup EXIT
fail() {
    printf 'fail: %s\n' "$*" >&2
    exit 1
}
# expect STATUS MESSAGE-PART -- zexec args: exit status and the stderr line.
expect() {
    local status="$1" part="$2" got=0
    shift 3
    "$x" "$@" 2>"$work/err" >/dev/null || got=$?
    [[ "$got" == "$status" ]] || fail "status $got, wanted $status: zexec $*"
    [[ -z "$part" ]] || grep -qF -- "$part" "$work/err" || fail "stderr lacks '$part': $(cat "$work/err")"
}

# Arguments reach the program unchanged, including spaces and empty strings.
"$x" -- /usr/bin/printf '[%s]\n' 'a b' '' '-o' '--' >"$work/args"
printf '[a b]\n[]\n[-o]\n[--]\n' | cmp -s - "$work/args" || fail 'arguments changed'

# stdin and stdout are the named files, freshly opened for each run.
printf '@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n' >"$work/in.fq"
for run in 1 2; do "$x" -i "$work/in.fq" -o "$work/out$run" -- /usr/bin/wc -c; done
size="$(wc -c <"$work/in.fq")"
[[ "$(cat "$work/out1")" == "$size" && "$(cat "$work/out2")" == "$size" ]] || fail 'stdin was not read in full on each run'
printf 'stale and longer than the new output\n' >"$work/trunc"
"$x" -o "$work/trunc" -- /usr/bin/printf x
[[ "$(cat "$work/trunc")" == x ]] || fail '-o did not truncate'

# @RUN@ becomes the same pid in -o paths and arguments, and differs per run.
"$x" -o "$work/run.@RUN@" -- /usr/bin/printf '%s' 'pid=@RUN@ again=@RUN@'
file="$(ls "$work" | grep '^run\.')"
[[ "$(cat "$work/$file")" == "pid=${file#run.} again=${file#run.}" ]] || fail 'token differs between path and argument'
"$x" -o "$work/run.@RUN@" -- /usr/bin/printf x
[[ "$(ls "$work" | grep -c '^run\.')" == 2 ]] || fail 'two runs shared one @RUN@ value'

# Environment passes through; the exit status is the program's own.
ZEXEC_PROBE=yes "$x" -- /usr/bin/env | grep -qx 'ZEXEC_PROBE=yes' || fail 'environment not passed'
expect 1 '' -- -- /usr/bin/false
expect 143 '' -- -- /usr/bin/bash -c 'kill -TERM $$'

# A closed stdin still ends up as the named file, not closed by close-on-exec.
"$x" -i "$work/in.fq" -o "$work/closed" -- /usr/bin/wc -c 0<&-
[[ "$(cat "$work/closed")" == "$size" ]] || fail 'closed stdin was not replaced'

# Long arguments and many arguments are passed whole.
long="$(head -c 100000 /dev/zero | tr '\0' a)@RUN@"
"$x" -o "$work/long" -- /usr/bin/printf '%s' "$long"
[[ "$(wc -c <"$work/long")" -gt 100000 ]] || fail 'long argument was cut'
mapfile -t many < <(seq 3000)
"$x" -o "$work/many" -- /usr/bin/printf '%s\n' "${many[@]}"
[[ "$(wc -l <"$work/many")" == 3000 ]] || fail 'many arguments were not all passed'

# -r runs the program as a child and appends one usage line per run; zexec ends
# with the program's exit status or signal, and the line records it.
"$x" -r "$work/usage" -i "$work/in.fq" -o "$work/rout" -- /usr/bin/wc -c
[[ "$(cat "$work/rout")" == "$size" ]] || fail '-r lost the redirects'
"$x" -r "$work/usage" -- /usr/bin/head -c 4000000 /dev/zero >/dev/null
expect 1 '' -- -r "$work/usage" -- /usr/bin/false
expect 143 '' -- -r "$work/usage" -- /usr/bin/bash -c 'kill -TERM $$'
[[ "$(wc -l <"$work/usage")" == 4 ]] || fail '-r did not write one line per run'
awk -F'\t' 'NF != 6 || $0 ~ /[^0-9\t]/ { exit 1 }' "$work/usage" || fail '-r line is not six integers'
[[ "$(awk -F'\t' 'NR == 3 { print $5 }' "$work/usage")" == 256 ]] || fail '-r did not record exit status 1'
[[ "$(awk -F'\t' 'NR == 4 { print $5 }' "$work/usage")" == 15 ]] || fail '-r did not record SIGTERM'
expect 127 'open /nonexistent/usage: NOENT' -- -r /nonexistent/usage -- /usr/bin/true

# Every failure exits 127 with one line naming the failed step.
expect 127 'usage:' --
expect 127 'usage:' -- -i "$work/in.fq"
expect 127 'unknown option -x' -- -x y -- /usr/bin/true
expect 127 'missing value for -i' -- -i
expect 127 'missing PROGRAM' -- --
expect 127 "needs a path with '/'" -- -- true
expect 127 'exec /nonexistent/prog: NOENT' -- -- /nonexistent/prog
expect 127 'open /nonexistent/in: NOENT' -- -i /nonexistent/in -- /usr/bin/true
expect 127 'open /nonexistent/dir/out: NOENT' -- -o /nonexistent/dir/out -- /usr/bin/true
expect 127 'exec' -- -- "$work/in.fq"
printf 'pass: arguments, stdin, stdout, truncate, @RUN@, environment, exit status, closed stdin, long and many arguments, -r usage lines and status, failures\n'

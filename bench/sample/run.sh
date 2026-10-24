#!/usr/bin/env bash
# Benchmark: sampling records from FASTQ files. Modes: init, check, release, report, peer, ab REV.
source "$(dirname -- "$0")/../lib/bench.sh"

job "sampling records from FASTQ files"
truth "every output must be an ordered subset of the input's records, unchanged, mates together, with the count asked for (exact) or a binomial count (fraction); the subset checker (bench/lib/subset.zig) decides"
expect accept
facts piped
# count and check already run z-fastq on every corpus file; here every source is in a family.
verify off
SUBSET_CHECKER=1
# Speed tables name what each tool lets through, since the work is not equal.
THROUGH=1
order z-fastq,rasusa,seqtk,fq,fqkit,seqkit,irma-core,fasten

word title "Sampling records"
word subject "z-fastq sample"
word others "samplers"
word timed "its sample was right"
word right "the right sample"
word wrong "a wrong sample"
word missed "gave a wrong sample"
word unseen "the sample is wrong"
word legend "wrong sample, exit 0"
word verb "samples"
word past "sampled"
word payer "a sampler"
word same "sampling"
word many "Sampling"
word truth "A sample is right when every kept record is an input record, unchanged in header, sequence, and quality, in input order, with mates kept together, and the count is the one asked for (exact count) or passes a two-sided binomial test at 1e-6 (fraction); the subset checker, which shares no code with z-fastq, decides."
word cells "A cell is the tool's verdict when asked for every record (K = 1,000, above every edge file's count, or p = 1): \`accept\` (exit 0, every record back unchanged and in order), \`reject S003\` (with the rule it named), \`reject\`, \`wrong: accept\` (exit 0 on a malformed input), \`wrong: N kept\` or \`wrong: not a subset\` (exit 0 on a valid input with another output: records lost, changed, or reordered), \`refused\` (rejected a valid input), or \`crash\`."
word probes "The probe count is how many of four malformed inputs the tool rejects: a sequence and quality length mismatch, a truncated last record, \`-\` instead of \`+\`, and a quality byte below 33."

lane z-fastq subject "{zfastq} sample --fraction {p} --seed {seed} {in}"
lane fq partial "{bin} subsample -p {p} -s {seed} --r1-dst {out1} {in}" note="writes only to named files"
lane fqkit partial "{bin} sample -n {k} -s {seed} -2 -@ 1 {in}" note="exact count only; 4 threads unless given -@ 1"
lane rasusa partial "{bin} reads -f {p} -1 -s {seed} {in}" note="-1 is its one-pass fraction; without it, two passes and exactly round(pN)"
lane irma-core complete "{bin} sampler -p {pct} -s {seed} {in}" note="whole percents only; an exact count (-t) on gzip input comes out of input order, so exact count is timed on plain input"
lane fasten partial "{bin} -f {p}" note="stdin only, plain only, no seed"
lane seqkit partial "{bin} sample -p {p} -s {seed} -j 1 {in}" note="exact counts need sample2 -2 -w 0; 4 threads unless given -j 1"
lane seqtk partial "{bin} sample -s {seed} {in} {p}" note="exact counts need -2 for input order"

family fraction "a share of the records, p = 0.1"
lanes z-fastq,seqtk,seqkit,rasusa,fq,irma-core
look full
param p 0.1
param pct 10
param seed 11
corpus_inputs

# The exact-count commands, shared by the two exact families.
exact_commands() {
    command z-fastq "{zfastq} sample --count {k} --seed {seed} {in}"
    command seqtk "{bin} sample -2 -s {seed} {in} {k}"
    command seqkit "{bin} sample2 -n {k} -2 -w 0 -s {seed} -j 1 {in}"
    command rasusa "{bin} reads -n {k} -s {seed} {in}"
    command fq "{bin} subsample -n {k} -s {seed} --r1-dst {out1} {in}"
    command irma-core "{bin} sampler -t {k} -s {seed} {in}"
}

family exact "an exact number of records: K = 100,000, or 1,000 for long reads"
lanes z-fastq,seqtk,seqkit,rasusa,fq,fqkit,irma-core
look full
formats plain
param k short:100000,long:1000
param seed 11
exact_commands
input real "$MEDIUM_SOURCES" medium interleaved
input shape "$SHAPE_SOURCES" medium

family exact-k "an exact number of records against K, on one 1 GiB short-read file"
lanes z-fastq,seqtk,seqkit,rasusa,fq,fqkit,irma-core
look lines
formats plain
param seed 11
vary k 1000,10000,100000,1000000
exact_commands
input real illumina-drr1001190 medium interleaved

family paired-fraction "mates in two files, p = 0.1"
lanes z-fastq,rasusa,fq,irma-core
look table
formats plain
param p 0.1
param pct 10
param seed 11
command z-fastq "{zfastq} sample --paired --fraction {p} --seed {seed} {in}"
command rasusa "{bin} reads -f {p} -1 -s {seed} -o {out1} -o {out2} {in}"
command fq "{bin} subsample -p {p} -s {seed} --r1-dst {out1} --r2-dst {out2} {in}"
input pair "$PAIRED_SOURCES,paired" medium

family paired-exact "mates in two files, K = 100,000 pairs"
lanes z-fastq,rasusa,fq,irma-core
look table
formats plain
param k 100000
param seed 11
command z-fastq "{zfastq} sample --paired --count {k} --seed {seed} {in}"
command rasusa "{bin} reads -n {k} -s {seed} -o {out1} -o {out2} {in}"
command fq "{bin} subsample -n {k} -s {seed} --r1-dst {out1} --r2-dst {out2} {in}"
command irma-core "{bin} sampler -t {k} -s {seed} {in}"
input pair "$PAIRED_SOURCES,paired" medium

family interleaved-fraction "mates interleaved in one file, on stdin, p = 0.1"
lanes z-fastq,fasten stdin
look table
formats plain
param p 0.1
param seed 11
command z-fastq "{zfastq} sample --interleaved --fraction {p} --seed {seed} {in}"
command fasten "{bin} -p -f {p}"
input real "$PAIRED_SOURCES" medium interleaved
input shape paired medium interleaved

# Edge inputs ask for every record back: an exact count above any edge file's
# record count, or p = 1 for a tool with no exact count (Fasten, single-end).
edges sample exact fraction "k=1000;p=1;pct=100;seed=11"

# The rule a failing tool named, as in check.
reject_code() {
    local code
    code="$(sed 's/\x1b\[[0-9;]*m//g' "$ERR" | grep -o -m1 -E ': \[?[SP]00[1-6]\]?:?' | grep -o -E '[SP]00[1-6]' | head -1)" || code=""
    printf 'reject%s\n' "${code:+ $code}"
}
# A plain path for an input: the plain copy beside a gzip file, or the
# decompressed stream when there is none.
plain_of() { if [[ "$1" != *.gz ]]; then printf '%s\n' "$1"; elif [[ -f "${1%.gz}" ]]; then printf '%s\n' "${1%.gz}"; else printf '%s\n' "<(gzip -dc '$1')"; fi; }
# Facts of the subset checker on one output: the stream on stdin, or the
# named outputs OUT1 and OUT2, each against its own mate file.
subset_facts() { # paths
    local -a in
    read -r -a in <<<"$1"
    local a b
    a="$(plain_of "${in[0]}")"
    if ((${#in[@]} == 2)); then b="$(plain_of "${in[1]}")"; fi
    if [[ -n "${OUT2:-}" ]]; then
        # Two named outputs: one checker per mate file, side by side, then
        # mates are in step when both kept the same index sequence.
        local f1 f2
        f1="$(mktemp /var/tmp/bench-sample-facts.XXXXXX)" f2="$f1.2"
        eval "\"\$SUBSET\" $a" <"$OUT1" >"$f1" &
        eval "\"\$SUBSET\" $b" <"$OUT2" >"$f2"
        wait
        awk 'FNR == NR { a[$1] = $2; next } { b[$1] = $2 }
            END { print "kept", a["kept"]
                  print "order", (a["order"] == "yes" && b["order"] == "yes") ? "yes" : "no"
                  print "pairs", (a["kept"] == b["kept"] && a["fingerprint"] == b["fingerprint"]) ? "yes" : "no"
                  print "plus", a["plus"] + b["plus"]; print "bare", a["bare"] + b["bare"]
                  print "cr", a["cr"] + b["cr"]; print "crlines", a["crlines"] + b["crlines"]; print "fingerprint", a["fingerprint"] }' "$f1" "$f2"
        rm -f -- "$f1" "$f2"
    elif [[ -n "${OUT1:-}" ]]; then
        eval "\"\$SUBSET\" $a" <"$OUT1"
    elif [[ -n "${b:-}" ]]; then
        eval "\"\$SUBSET\" $a $b"
    elif [[ "${FAMILY:-}" == interleaved-* && -z "${EDGE:-}" ]]; then
        eval "\"\$SUBSET\" --pairs $a"
    else
        eval "\"\$SUBSET\" $a"
    fi
}
# The verdict, then the checker's facts: accept, or wrong with what was wrong.
fact() {
    if [[ "${RC:-0}" != 0 ]]; then
        reject_code
        return
    fi
    local facts n k="" p="" kv
    facts="$(subset_facts "$IN")"
    for kv in ${PARAMS//;/ }; do
        case "$kv" in k=*) k="${kv#k=}" ;; p=*) p="${kv#p=}" ;; esac
    done
    # The record (or pair) count of the input: the case truth, or counted for an edge file.
    n="${TRUTH:-}"
    if [[ -n "${EDGE:-}" || -z "$n" || "$n" == - ]]; then n="$(eval "cat $(plain_of "$IN")" | awk 'END { print int((NR + 3) / 4) }')"; fi
    [[ "${FAMILY:-}" != interleaved-* || -n "${EDGE:-}" ]] || n=$((n / 2))
    awk -v n="$n" -v k="$k" -v p="$p" '{ f[$1] = $2 } END {
            kept = f["kept"] + 0
            if (f["order"] != "yes") print "wrong: not a subset"
            else if (("pairs" in f) && f["pairs"] != "yes") print "wrong: mates out of step"
            else if (k != "") print (kept == (k + 0 < n + 0 ? k + 0 : n + 0)) ? "accept" : "wrong: " kept " kept"
            else if (p + 0 >= 1) print (kept == n + 0) ? "accept" : "wrong: " kept " kept"
            else { mu = n * p; sd = sqrt(n * p * (1 - p)); z = sd > 0 ? (kept - mu) / sd : 0
                   print (z <= 4.89 && z >= -4.89) ? "accept" : "wrong: " kept " kept" }
        }' <<<"$facts"
    printf '%s\nn %s\n' "$facts" "$n"
}

# --- feature steps: fairness and abilities ------------------------------------
FAIR="$RESULTS/fair"
FAIR_SEEDS=200
# Fixtures whose headers are the record index (@u0000001; mates share their
# pair's index), so kept indexes are read straight from a tool's output.
make_fixtures() {
    mkdir -p -- "$FAIR"
    [[ -f "$FAIR/ends.fq" && -f "$FAIR/pairs_interleaved.fq" ]] && return 0
    awk 'NR % 4 == 1 { printf "@u%07d\n", (NR + 3) / 4; next } { print }' "$REAL/illumina-drr147610/xsmall.fastq" >"$FAIR/body.fq"
    head -n 800 "$FAIR/body.fq" >"$FAIR/ends.fq"
    awk -v m=1 'NR % 4 == 1 { printf "@u%07d/%d\n", (NR + 3) / 4, m; next } { print }' "$REAL/illumina-drr1001190/xsmall_R1.fastq" >"$FAIR/pairs_R1.fq"
    awk -v m=2 'NR % 4 == 1 { printf "@u%07d/%d\n", (NR + 3) / 4, m; next } { print }' "$REAL/illumina-drr1001190/xsmall_R2.fastq" >"$FAIR/pairs_R2.fq"
    paste -d '\n' <(paste - - - - <"$FAIR/pairs_R1.fq") <(paste - - - - <"$FAIR/pairs_R2.fq") | tr '\t' '\n' >"$FAIR/pairs_interleaved.fq"
}
# The fixture and parameters of each family's fairness runs.
fair_setup() { # family fixture -> paths, records, params (no seed)
    case "$1:$2" in
        fraction:body) printf '%s\t%s\t%s\n' "$FAIR/body.fq" 79100 "p=0.1;pct=10" ;;
        fraction:ends) printf '%s\t%s\t%s\n' "$FAIR/ends.fq" 200 "p=0.1;pct=10" ;;
        exact:body) printf '%s\t%s\t%s\n' "$FAIR/body.fq" 79100 "k=1000" ;;
        exact:ends) printf '%s\t%s\t%s\n' "$FAIR/ends.fq" 200 "k=20" ;;
        paired-fraction:pairs) printf '%s\t%s\t%s\n' "$FAIR/pairs_R1.fq $FAIR/pairs_R2.fq" 14413 "p=0.1;pct=10" ;;
        paired-exact:pairs) printf '%s\t%s\t%s\n' "$FAIR/pairs_R1.fq $FAIR/pairs_R2.fq" 14413 "k=1000" ;;
        interleaved-fraction:pairs) printf '%s\t%s\t%s\n' "$FAIR/pairs_interleaved.fq" 14413 "p=0.1" ;;
    esac
}
fair_jobs() {
    local fam fixture name seed
    for fam in fraction exact paired-fraction paired-exact interleaved-fraction; do
        for fixture in body ends pairs; do
            [[ -n "$(fair_setup "$fam" "$fixture")" ]] || continue
            for name in ${FAM_LANES[$fam]}; do
                for ((seed = 1; seed <= FAIR_SEEDS; seed++)); do printf '%s %s %s %s 1\n' "$name" "$fam" "$fixture" "$seed"; done
                printf '%s %s %s 11 2\n' "$name" "$fam" "$fixture"
            done
        done
    done
}
# One run: lane family fixture seed rep -> the kept count, a hash of the kept
# indexes, the kept count per tenth of the fixture, whether the first and the
# last record were kept, and whether mates came together.
fair_one() {
    local name fam fixture seed rep paths n params cmd out o1="" o2="" tmpl
    read -r name fam fixture seed rep <<<"$1"
    IFS=$'\t' read -r paths n params < <(fair_setup "$fam" "$fixture")
    out="$(mktemp /var/tmp/bench-sample-fair.XXXXXX)"
    tmpl="$(lane_tmpl "$name" "$fam")"
    [[ "$tmpl" != *"{out1}"* ]] || o1="$out.1"
    [[ "$tmpl" != *"{out2}"* ]] || o2="$out.2"
    if [[ "${FAM_STDIN[$fam]}" == 1 ]]; then
        cmd="$(OUT1="$o1" OUT2="$o2" lane_cmd "$name" - "$fam" "$params;seed=$seed")"
        $cmd <"$paths" >"$out" 2>/dev/null || true
    else
        cmd="$(OUT1="$o1" OUT2="$o2" lane_cmd "$name" "$paths" "$fam" "$params;seed=$seed")"
        $cmd >"$out" 2>/dev/null || true
    fi
    # Indexes in output order: the named outputs, else stdout.
    awk -v n="$n" -v lane="$name" -v fam="$fam" -v fx="$fixture" -v seed="$seed" -v rep="$rep" -v two="$([[ -n "$o2" ]] && echo 1 || echo 0)" '
        FNR % 4 == 1 { i = substr($1, 3) + 0; if (two) { if (FILENAME == ARGV[1]) a[++na] = i; else b[++nb] = i } else s[++ns] = i }
        END {
            mates = 1
            if (two) { if (na != nb) mates = 0; for (j = 1; j <= na; j++) { if (a[j] != b[j]) mates = 0; s[j] = a[j] }; ns = na }
            else if (fam ~ /^paired|^interleaved/) {
                # One stream of pairs: mates follow each other with one index.
                if (ns % 2) mates = 0
                for (j = 1; j < ns; j += 2) if (s[j] != s[j + 1]) mates = 0
                m = 0; for (j = 1; j <= ns; j += 2) t[++m] = s[j]; ns = m; for (j = 1; j <= m; j++) s[j] = t[j]
            }
            h = 0; for (j = 1; j <= ns; j++) h = (h * 31 + s[j]) % 2147483647
            for (d = 0; d < 10; d++) tenth[d] = 0
            first = last = 0
            for (j = 1; j <= ns; j++) { d = int((s[j] - 1) * 10 / n); if (d > 9) d = 9; tenth[d]++; if (s[j] == 1) first = 1; if (s[j] == n) last = 1 }
            line = lane "\t" fam "\t" fx "\t" seed "\t" rep "\t" ns "\t" h
            for (d = 0; d < 10; d++) line = line "\t" tenth[d]
            print line "\t" first "\t" last "\t" mates
        }' ${o1:+"$o1"} ${o2:+"$o2"} $([[ -z "$o1" ]] && echo "$out")
    rm -f -- "$out" "$o1" "$o2"
}
# Abilities measured on tiny files: reading stdin and gzip, each with the
# fraction command (exact for a lane with none) and the subset checker.
ABILITY="$RESULTS/abilities.tsv"
ability_jobs() {
    local name
    for name in "${LANES[@]}"; do printf '%s stdin\n%s gzip\n' "$name" "$name"; done
}
ability_one() {
    local name what fam params file cmd out rc=0 o1="" verdict facts
    read -r name what <<<"$1"
    fam=fraction params="p=0.1;pct=10;seed=11"
    [[ " ${FAM_LANES[fraction]} " == *" $name "* ]] || { fam=exact params="k=1000;seed=11"; }
    [[ " ${FAM_LANES[$fam]} " == *" $name "* ]] || { fam=interleaved-fraction params="p=0.1;seed=11"; }
    file="$REAL/illumina-drr1001190/tiny_interleaved.fastq"
    out="$(mktemp /var/tmp/bench-sample-ability.XXXXXX)"
    [[ "$(lane_tmpl "$name" "$fam")" != *"{out1}"* ]] || o1="$out.1"
    if [[ "$what" == stdin ]]; then
        cmd="$(OUT1="$o1" lane_cmd "$name" - "$fam" "$params")"
        $cmd <"$file" >"$out" 2>/dev/null || rc=$?
    else
        cmd="$(OUT1="$o1" lane_cmd "$name" "$file.gz" "$fam" "$params")"
        if [[ "${FAM_STDIN[$fam]}" == 1 ]]; then $cmd <"$file.gz" >"$out" 2>/dev/null || rc=$?; else $cmd >"$out" 2>/dev/null || rc=$?; fi
    fi
    verdict=no
    if ((rc == 0)); then
        facts="$(OUT1="$o1" FAMILY="$fam" subset_facts "$file" <"$out")"
        local kept
        kept="$(awk '$1 == "kept" { print $2 }' <<<"$facts")"
        [[ "$(awk '$1 == "order" { print $2 }' <<<"$facts")" != yes || "${kept:-0}" -le 0 ]] || verdict=yes
    fi
    printf '%s\t%s\t%s\n' "$name" "$what" "$verdict"
    rm -f -- "$out" "$o1"
}
feature_check() {
    make_fixtures
    cached_step fairness "$(lane_keys; printf '%s\n' "$FAIR_SEEDS"; fileid "$FAIR"/*.fq; declare -f fair_setup fair_jobs fair_one | cksum)" \
        fair_jobs fair_one "$FAIR/runs.tsv"
    cached_step abilities "$(lane_keys; fileid "$REAL/illumina-drr1001190/tiny_interleaved.fastq"{,.gz}; declare -f ability_one subset_facts plain_of | cksum)" \
        ability_jobs ability_one "$ABILITY"
}

# --- report sections ----------------------------------------------------------
# Per lane and family: the fairness tests, one row each. Columns: family, lane,
# count bias z (or "exact"), spread over the binomial's (or "fixed"), worst
# tenth z, first kept, last kept (of 200, "-" for pairs), same twice, seeds
# differ, mates, fair.
fair_table_rows() {
    awk -F'\t' -v S="$FAIR_SEEDS" '
        BEGIN {
            # The two-sided 1e-6 acceptance range of Binomial(200, 0.1).
            q = 0.1; for (x = 0; x <= S; x++) pmf[x] = exp(lgam(S + 1) - lgam(x + 1) - lgam(S - x + 1) + x * log(q) + (S - x) * log(1 - q))
            c = 0; for (x = 0; x <= S; x++) { c += pmf[x]; if (c > 0.5e-6) { lo = x; break } }
            c = 0; for (x = S; x >= 0; x--) { c += pmf[x]; if (c > 0.5e-6) { hi = x; break } }
        }
        function lgam(x,   s) { if (x < 2) return 0; s = 0; for (i = 2; i < x; i++) s += log(i); return s }
        {
            key = $1 "\t" $2
            if ($3 == "ends") { if ($5 == 1) { ends[key] = 1; first[key] += $18; last[key] += $19 } next }
            if ($5 == 2) { again[key] = $7; next }
            keys[key] = 1; fx[key] = $3; runs[key]++; sum[key] += $6; sq[key] += $6 * $6
            if ($4 == 11) h11[key] = $7; if ($4 == 12) h12[key] = $7
            for (d = 0; d < 10; d++) tenth[key, d] += $(8 + d)
            if (!$20) badmates[key] = 1
        }
        END {
            for (key in keys) {
                split(key, kf, "\t"); fam = kf[2]; R = runs[key]
                N = (fx[key] == "body") ? 79100 : 14413
                exact = fam ~ /exact/
                mean = sum[key] / R; var = sq[key] / R - mean * mean; if (var < 0) var = 0
                if (exact) { K = 1000; qq = K / N; bias = (mean == K && var == 0) ? "exact" : sprintf("%.1f", mean - K); biasok = (mean == K && var == 0) }
                else { qq = 0.1; sdb = sqrt(N * qq * (1 - qq)); z = (mean - N * qq) / (sdb / sqrt(R)); bias = sprintf("%+.1f", z); biasok = (z <= 4.89 && z >= -4.89) }
                if (exact || var == 0) { spread = exact ? "-" : "fixed"; spreadok = 1 }
                else { r = sqrt(var) / sdb; spread = sprintf("%.2f", r); spreadok = (r >= 0.5 && r <= 1.5) }
                worst = 0
                for (d = 0; d < 10; d++) {
                    nt = int(N * (d + 1) / 10) - int(N * d / 10)
                    e = R * qq * nt
                    v = exact ? R * K * (nt / N) * (1 - nt / N) * (N - K) / (N - 1) : R * nt * qq * (1 - qq)
                    zt = (tenth[key, d] - e) / sqrt(v); if (zt < 0) zt = -zt; if (zt > worst) worst = zt
                }
                tenthok = worst <= 5.3
                if (key in ends) { fk = first[key]; lk = last[key]; endsok = (fk >= lo && fk <= hi && lk >= lo && lk <= hi) } else { fk = lk = "-"; endsok = 1 }
                twice = (key in again && again[key] == h11[key]) ? "yes" : "no"
                differ = (h11[key] != h12[key]) ? "yes" : "no"
                mates = (fam ~ /^paired|^interleaved/) ? ((key in badmates) ? "no" : "yes") : "-"
                fair = (biasok && spreadok && tenthok && endsok && differ == "yes" && mates != "no") ? "yes" : "no"
                printf "%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\n", fam, kf[1], bias, spread, worst, fk, lk, twice, differ, mates, fair
            }
            printf "#range\t%d\t%d\n", lo, hi
        }' "$FAIR/runs.tsv"
}
# The fairness table, families in declaration order, lanes in registry order.
fair_table() {
    local rows fam name
    rows="$(fair_table_rows)"
    {
        printf 'Mode\tTool\tCount bias (z)\tSpread\tWorst tenth (z)\tFirst kept\tLast kept\tSame twice\tSeeds differ\tMates\tFair\n'
        printf 'l\tl\tr\tr\tr\tr\tr\tl\tl\tl\tl\n'
        for fam in fraction exact paired-fraction paired-exact interleaved-fraction; do
            for name in "${LANES[@]}"; do
                awk -F'\t' -v OFS='\t' -v f="$fam" -v l="$name" -v label="${REG_LABEL[$name]}" -v mode="$(mode_name "$fam")" \
                    '$1 == f && $2 == l { $1 = mode; $2 = label; print }' <<<"$rows"
            done
        done
    } | mdtable
}
# A family's name for text.
mode_name() { case "$1" in fraction) echo Fraction ;; exact) echo "Exact count" ;; exact-k) echo "Exact count against K" ;;
    paired-fraction) echo "Paired fraction" ;; paired-exact) echo "Paired exact count" ;; interleaved-fraction) echo "Interleaved fraction" ;; *) echo "$1" ;; esac; }
# fair_verdicts: lane family fair, for the findings.
fair_verdicts() { fair_table_rows | awk -F'\t' '$1 !~ /^#/ { print $2, $1, $11 }'; }
# Per lane: one row of what it can do and what it changes in kept records.
# Modes come from the checks of the families it is declared in; stdin and
# gzip from the ability runs; seed from the fairness runs; plus line and line
# endings from its checks on the variable (annotated plus) and endings (CRLF)
# shapes at the headline tier.
tools_table() {
    local name fam fams cell stdin gz seed line
    {
        printf 'Tool\tClass\tFraction\tExact count\tPairs, two files\tPairs, interleaved\tStdin\tGzip\tSeed\n'
        printf 'l\tl\tl\tl\tl\tl\tl\tl\tl\n'
        for name in "${LANES[@]}"; do
            line="${REG_LABEL[$name]}"$'\t'"${LANE_CLASS[$name]}"
            for fams in fraction exact "paired-fraction paired-exact" interleaved-fraction; do
                cell="-"
                for fam in $fams; do [[ " ${FAM_LANES[$fam]} " != *" $name "* ]] || cell="$(mode_cell "$name" "$fam")"; done
                line+=$'\t'"$cell"
            done
            stdin="$(awk -F'\t' -v l="$name" '$1 == l && $2 == "stdin" { print $3 }' "$ABILITY")"
            gz="$(awk -F'\t' -v l="$name" '$1 == l && $2 == "gzip" { print $3 }' "$ABILITY")"
            seed="$(fair_table_rows | awk -F'\t' -v l="$name" '$2 == l { t = t $8 $9 } END { print (t ~ /no/) ? "no" : "yes" }')"
            printf '%s\t%s\t%s\t%s\n' "$line" "${stdin:--}" "${gz:--}" "$seed"
        done
    } | mdtable
}
# Per lane: what it rejects, by kind of malformed input, and what it changes
# in kept records. A kind is passed when the tool exits 0 on any of its edge
# inputs.
checks_table() {
    local name plus ends unfair
    {
        printf 'Tool\tLength mismatch\tTruncated record\tBad plus line\tBad header\tQuality range\tBlank line\tBroken gzip\tPlus line written\tLine endings written\tUnfair in\n'
        printf 'l\tl\tl\tl\tl\tl\tl\tl\tl\tl\tl\n'
        for name in "${LANES[@]}"; do
            plus="$(record_cell "$name" variable plus)" ends="$(record_cell "$name" endings cr)"
            unfair="$(fair_verdicts | awk -v l="$name" '$1 == l && $3 == "no" { print $2 }' | while read -r f; do mode_name "$f"; done | paste -sd, - | sed 's/,/, /g')"
            awk -F'\t' -v l="$name" -v label="${REG_LABEL[$name]}" -v tail="$plus"$'\t'"$ends"$'\t'"${unfair:--}" '
                BEGIN { split("lenmis trunc dash noat qlow,qhigh blank gztrunc,gzgarbage", kinds, " ") }
                FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next }
                $1 == l && ($2 in want) { cell[$2] = $3 }
                END {
                    line = label
                    for (i = 1; i <= 7; i++) { n = split(kinds[i], e, ","); passed = 0
                        for (j = 1; j <= n; j++) if (cell[e[j]] ~ /^wrong/) passed = 1
                        line = line "\t" (passed ? "**passes**" : "rejects") }
                    print line "\t" tail
                }' <(edge_list) "$RESULTS/edges.tsv"
        done
    } | mdtable
}
# yes when every check of LANE in FAMILY passed, partly when some did, no when none.
mode_cell() { # lane family
    awk -F'\t' -v l="$1" -v f="$2" '$1 == f && $2 == l { n++; if ($3 == "pass") p++ } END { print (n == 0) ? "-" : (p == n) ? "yes" : (p > 0) ? "partly" : "no" }' "$FAMILY_CHECKS"
}
# Every current check of every family into FAMILY_CHECKS: family, lane,
# status, value, case id. Written once per report.
FAMILY_CHECKS=""
family_checks() {
    local fam id fmt paths truth bytes src tier kind params inid name form key
    FAMILY_CHECKS="$REPORT_DIR/family-checks.tsv"
    load_checks
    cache_lane_ids
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params inid; do
        form=file
        [[ "${FAM_STDIN[$fam]}" == 0 ]] || form=stdin
        for name in ${FAM_LANES[$fam]}; do
            key="$(status_key "$name" "$fam" "$form" "$inid" "$params")"
            [[ ! -v "CHECK_STATUS[$key]" ]] || printf '%s\t%s\t%s\t%s\t%s\n' "$fam" "$name" "${CHECK_STATUS[$key]}" "${CHECK_VALUE[$key]}" "$id"
        done
    done < <(cases | with_input_ids) >"$FAMILY_CHECKS"
}
# What LANE's kept records on SHAPE's headline-tier plain case carry: the
# plus line (kept, bare +, changed) or the line endings (as input, LF, CRLF,
# mixed).
record_cell() { # lane shape fact
    local fam
    fam=fraction
    [[ " ${FAM_LANES[fraction]} " == *" $1 "* ]] || fam=exact
    [[ " ${FAM_LANES[$fam]} " == *" $1 "* ]] || { printf -- '-\n'; return; }
    awk -F'\t' -v f="$fam" -v l="$1" -v id="$2:$HEADLINE:plain" -v what="$3" '
        $1 == f && $2 == l && $5 == id { n = split($4, x, ";"); for (i = 1; i <= n; i++) { split(x[i], kv, " "); v[kv[1]] = kv[2] } found = 1 }
        END {
            if (!found || !("kept" in v)) { print "-"; exit }
            k = v["kept"] + 0
            if (what == "plus") print (v["plus"] == 0) ? "kept" : (v["bare"] + 0 >= k) ? "bare +" : "changed"
            else print (v["cr"] == 0) ? "as input" : (v["crlines"] == 0) ? "LF" : (v["crlines"] + 0 >= 4 * k) ? "CRLF" : "mixed"
        }' "$FAMILY_CHECKS"
}
# seqtk's selection against z-fastq's on every screened input both read
# (LF endings, bare plus lines): FAMILY matched total.
seqtk_match() {
    local fam
    for fam in fraction exact; do
        awk -F'\t' -v f="$fam" '$1 != f || $5 ~ /^(endings|variable):/ { next }
            $3 == "pass" && ($2 == "z-fastq" || $2 == "seqtk") { n = split($4, x, ";"); for (i = 1; i <= n; i++) if (x[i] ~ /^fingerprint /) fp[$2, $5] = substr(x[i], 13); cs[$5] = 1 }
            END { for (c in cs) if (("z-fastq", c) in fp && ("seqtk", c) in fp) { t++; if (fp["z-fastq", c] == fp["seqtk", c]) m++; else miss = miss " " c }
                  printf "%s\t%d\t%d\t%s\n", f, m, t, miss }' "$FAMILY_CHECKS"
    done
}
# Peak memory against K, from the exact-k family: per K, z-fastq's peak and
# the other tools' as multiples of it.
k_memory() {
    awk -F'\t' 'FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == "exact-k") { split($12, p, ";"); for (i in p) if (p[i] ~ /^k=/) k[$2] = substr(p[i], 3) } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
        FNR > 1 && ($2 in k) { m[k[$2], $3] = $8; ks[k[$2]] = 1; lanes[$3] = 1 }
        END { for (K in ks) if ((K, "z-fastq") in m) { lo = hi = 0
                for (l in lanes) if (l != "z-fastq" && (K, l) in m) { r = m[K, l] / m[K, "z-fastq"]; if (!lo || r < lo) { lo = r; ll = l } if (r > hi) { hi = r; hl = l } }
                printf "%d\t%.2f\t%.2f\t%s\t%.1f\t%s\n", K, m[K, "z-fastq"] / 1024, lo, label[ll], hi, label[hl] } }' \
        "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | sort -n
}
sample_findings() {
    local bad="" name fam fair
    while read -r name fam fair; do
        [[ "$fair" == no ]] && bad+="${bad:+, }${REG_LABEL[$name]} ($(mode_name "$fam" | tr 'A-Z' 'a-z'))"
    done < <(fair_verdicts)
    # Right answers over every tool and mode: the checks and the edge matrix.
    awk -F'\t' -v others="$((${#LANES[@]} - 1))" 'FILENAME ~ /checks.tsv$/ { if ($1 == "z-fastq") { if ($3 == "pass") zp++; else zf++ } next }
        FNR > 1 && $6 != "policy" { if ($4 == "z-fastq") { if ($6 != "right") zbad++ } else if ($6 == "wrong" || $6 == "crash") bad[$4] = 1 }
        END { line = zf ? sprintf("- **Right answers:** z-fastq gave a wrong sample on %d of %d checked inputs", zf, zp + zf) : sprintf("- **Right answers:** z-fastq returned the right sample on all %d checked inputs", zp)
              line = line (zbad ? "." : " and rejected every malformed edge case.")
              k = 0; for (b in bad) k++
              if (k) line = line sprintf(" %d of the %d other tools print a wrong sample or crash on at least one edge case.", k, others)
              print line }' "$REPORT_DIR/checks.tsv" "$REPORT_DIR/edges.tsv"
    if [[ -n "$bad" ]]; then
        printf -- '- **Fair samples:** every tool passes the fairness tests except %s.\n' "$bad"
    else
        printf -- '- **Fair samples:** every tool passes the fairness tests in every mode.\n'
    fi
    seqtk_match | awk -F'\t' '{ t += $3; m += $2; d = d sprintf("%s%s %d of %d", (NR > 1 ? ", " : ""), ($1 == "exact" ? "exact count" : $1), $2, $3) } END {
        printf "- **Matches seqtk:** z-fastq kept the same records as seqtk on %d of %d screened inputs (%s).\n", m, t, d }'
    ANY=1 findings fraction speed Fraction
    ANY=1 findings exact speed Exact-count
    k_memory | awk -F'\t' '$1 == 100000 || $1 == 1000000 { line = line sprintf("%sat K = %s, z-fastq peaks at %.2f MiB and the other tools at %.1fx to %.1fx that", (n++ ? "; " : ""), ($1 == 100000 ? "100,000" : "1,000,000"), $2, $3, $5) }
        END { if (line != "") print "- **Memory against K:** " line "." }'
    findings fraction where Fraction
    findings exact where "Exact count"
}
feature_readme() {
    family_checks
    readme_head
    printf '## Findings\n\n'
    sample_findings
    printf '\n## Compared tools\n> <sub>*complete* tools check records as much as z-fastq; *partial* tools check less, so they do less work per record. Modes: `yes` when the tool was right on every input of the mode, `partly` on some, `-` when the tool has no such mode. Stdin and gzip: measured on a 1 MiB file with the timed command (a two-pass exact count cannot read stdin). Seed: the same output twice with one seed, and different output with another.</sub>\n\n'
    tools_table
    printf '\n### What each tool checks and changes\n> <sub>Each kind of malformed input from the edge cases: `rejects`, or **passes** when the tool writes a sample from it with exit 0, so it does less work than z-fastq. Plus line and line endings written: on the `variable` (annotated plus line) and `endings` (CRLF and LF) shapes; `as input` and `kept` copy the input. Unfair in: the modes where its sample fails the fairness tests below.</sub>\n\n'
    checks_table
    printf '\n'
    local name
    for name in "${LANES[@]}"; do [[ -z "${LANE_NOTE[$name]}" ]] || printf -- '- `%s %s`: *%s*\n' "${REG_LABEL[$name]}" "$(lane_version "$name")" "${LANE_NOTE[$name]}"; done
    printf '\n> [!IMPORTANT]\n> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.\n\n'
    printf '## Matches seqtk\n> <sub>z-fastq promises seqtk'"'"'s selection on inputs with LF endings, nonempty reads, bare plus lines, printable headers, a seed that fits signed 64-bit, and a value seqtk reads as a fraction. Equivalent commands: `z-fastq sample --fraction P --seed S IN` and `seqtk sample -s S IN P`; `z-fastq sample --count K --seed S IN` and `seqtk sample -2 -s S IN K`.</sub>\n\n'
    seqtk_match | awk -F'\t' '{ printf "- `%s`: *the same records on %d of %d inputs%s*\n", ($1 == "fraction" ? "fraction" : "exact count"), $2, $3, ($4 == "" ? "" : "; different on" $4) }'
    printf '\n> [!NOTE]\n> seqtk without `-2` keeps the same records but writes them in reservoir order, reads `1.0` as a count of 1, and reads `0` as a fraction.\n'
    readme_right
    printf '## Is the sample fair?\n> <sub>%d seeds per tool and mode on files whose headers are the record index: 79,100 short reads (single-end modes) and 14,413 pairs (pair modes), p = 0.1 and K = 1,000; first and last record on a 200-record file with K = 20 or p = 0.1. Count bias: the mean kept count'"'"'s distance from N p in standard errors (exact: every run kept K). Spread: the standard deviation over the binomial'"'"'s (`fixed`: the same count every run). Worst tenth: the largest distance of a tenth of the file'"'"'s kept share from its expected share, in standard errors.</sub>\n\n' "$FAIR_SEEDS"
    fair_table
    printf '\n'
    fair_table_rows | awk -F'\t' '$1 == "#range" { printf "A tool fails at a count bias beyond 4.9, a spread outside 0.5 to 1.5, a tenth beyond 5.3 (each two-sided 1e-6), the first or last record kept outside %d to %d times of %d, or the same kept set from seeds 11 and 12. At these sizes a tenth'"'"'s share is tested to about 0.7%% (fraction) and 2%% (exact) of its value. Same twice is reproducibility, not fairness.\n\n", $2, $3, '"$FAIR_SEEDS"' }'
    readme_speed fraction "Fraction sampling, p = 0.1"
    readme_scaling fraction "###"
    readme_speed exact "Exact count, K = 100,000 (1,000 for long reads)"
    printf 'Plain input only: with gzip, most exact-count tools decompress the file twice, so the time is mostly zlib; fraction sampling shows the gzip cost. K is 100,000 reads for short-read inputs and 1,000 for long-read inputs (the Long reads group, `longread`, and `long`); the `long` shape holds 1,789 records, so its 1,000 is close to a copy.\n\n'
    printf '### Time and memory against K\n\n'
    picture exact-k-lines "Time and peak memory against K" && printf '\n'
    k_memory | awk -F'\t' 'BEGIN { print "K\tz-fastq peak\tLowest other\tHighest other"; print "r\tr\tl\tl" }
        function commas(n,   s) { s = n ""; while (s ~ /[0-9]{4}/) sub(/[0-9]{3}($|,)/, ",&", s); sub(/,,/, ",", s); return s }
        { printf "%s\t%.2f MiB\t%.1fx (%s)\t%.1fx (%s)\n", commas($1), $2, $3, $4, $5, $6 }' | mdtable
    printf '\n## Pairs\n> <sub>Each tool'"'"'s time and peak memory as a multiple of z-fastq'"'"'s, geometric mean over the 4 paired runs and the `paired` shape at 1 GiB. Plain input: paired fraction and paired exact count read two mate files, interleaved fraction one interleaved file on stdin.</sub>\n\n'
    summary_table "paired-fraction:Paired fraction" "paired-exact:Paired exact count" "interleaved-fraction:Interleaved fraction"
    printf '\n'
    readme_repro
}

bench_main "$@"

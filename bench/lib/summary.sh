# Benchmark reports, sourced by bench.sh: each feature's README.md, its report
# tables and figures, and the summary across features. Nothing here decides a
# verdict or a timing, so editing it keeps every cached check and batch.
# Aligns a TSV (header, alignment row of l/r, rows) into a Markdown table.
mdtable() {
    awk -F'\t' '
        NR == 1 { n = NF; for (i = 1; i <= NF; i++) { h[i] = $i; w[i] = length($i) } next }
        NR == 2 { for (i = 1; i <= NF; i++) a[i] = $i; next }
        { r++; for (i = 1; i <= n; i++) { c[r, i] = $i; if (length($i) > w[i]) w[i] = length($i) } }
        END {
            for (i = 1; i <= n; i++) if (w[i] < 3) w[i] = 3
            line = "|"; for (i = 1; i <= n; i++) line = line " " pad(h[i], w[i], a[i]) " |"; print line
            line = "|"; for (i = 1; i <= n; i++) line = line " " (a[i] == "r" ? dashes(w[i] - 1) ":" : dashes(w[i])) " |"; print line
            for (j = 1; j <= r; j++) { line = "|"; for (i = 1; i <= n; i++) line = line " " pad(c[j, i], w[i], a[i]) " |"; print line }
        }
        function dashes(k,   s) { s = ""; while (k-- > 0) s = s "-"; return s }
        function pad(s, k, al,   g) { g = k - length(s); while (g-- > 0) s = (al == "r") ? " " s : s " "; return s }'
}
picture() { # stem alt; nothing when the figure was not drawn
    [[ -e "$RESULTS/figures/$1-light.png" ]] || return 0
    cat <<EOF
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/figures/$1-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="results/figures/$1-light.png">
  <img src="results/figures/$1-light.png" alt="$2" width="100%">
</picture>
EOF
}
style_tsv() {
    printf 'lane\tlabel\tclass\tlight\tdark\tshape\torder\tlanguage\tcores\tof\n'
    local name
    for name in "${LANES[@]}"; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "${REG_LABEL[$name]}" "${LANE_CLASS[$name]}" "${REG_LIGHT[$name]}" "${REG_DARK[$name]}" "${REG_SHAPE[$name]}" "${REG_ORDER[$name]}" "${REG_LANG[$name]}" "${LANE_CORES[$name]:-1}" "${LANE_OF[$name]:-$name}"
    done
}
cases_tsv() {
    printf 'family\tcase\tformat\tsource\ttier\tkind\tbytes\tlabel\torder\tgroup\tgroup_order\tparams\n'
    local fam id fmt paths truth bytes src tier kind params n=0
    while IFS=$'\t' read -r fam id fmt paths truth bytes src tier kind params; do
        # A pair cut from a synthetic shape is a shape in the figures.
        [[ "$kind" != pair || ! -d "$SHAPES/$src" ]] || kind=shape
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$fam" "$id" "$fmt" "$src" "$tier" "$kind" "$bytes" \
            "$(label_of "$src")" $((n++)) "${GROUP_OF[$src]:-Other}" "${GROUP_RANK[$src]:-99}" "$params"
    done < <(cases)
}

lane_version() {
    local name="$1" bin receipt pin v
    bin="$(lane_bin "$name")" receipt="${REG_RECEIPT[$name]}" pin="${REG_PIN[$name]}"
    case "$name" in
        z-fastq) "$bin" --version | awk '{ print $2 }' ;;
        *)
            if [[ "$receipt" == adapter ]]; then
                "$bin" --version | sed -n 's/.*engine=[a-z]* \([^;]*\);.*/\1/p'
            else
                awk -F'\t' '$1 == "version" { print $2 }' "$REPO/tools/.local/$receipt-current/receipt.tsv"
            fi
            ;;
    esac
}
# Per lane, the geometric mean over the headline shapes of FAMILY of its time
# and peak RSS as multiples of z-fastq's, on the shapes both ran. With several
# FAMILY:MODE arguments, one table with a Mode column (MODE_HEADER names it),
# a row per mode and lane.
# With THROUGH=1, a last column counts the malformed inputs each tool lets through (lets_through).
# Thread variants are left out: a feature with them shows them in its own
# threads part. z-fastq's row gives its own time and memory; every other row
# gives the tool's as a multiple of z-fastq's, within each batch.
summary_table() { # family | family:mode...
    local fams
    fams="$(printf '%s|' "$@")"
    awk -F'\t' -v fams="${fams%|}" -v tier="$HEADLINE" -v through="${THROUGH:-0}" -v modehead="${MODE_HEADER:-Mode}" '
        function secs(ns) { return ns >= 1e9 ? sprintf("%.2f s", ns / 1e9) : ns >= 1e7 ? sprintf("%.0f ms", ns / 1e6) : sprintf("%.1f ms", ns / 1e6) }
        function mib(kib) { return kib >= 10240 ? sprintf("%.1f MiB", kib / 1024) : sprintf("%.2f MiB", kib / 1024) }
        BEGIN { nf = split(fams, fl, "|"); for (i = 1; i <= nf; i++) { split(fl[i], q, ":"); fam[i] = q[1]; mode[i] = q[2]; at[q[1]] = i } }
        FILENAME ~ /through.tsv$/ { lets[$2, $1] = $4; next }
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && ($1 in at) && $5 == tier && $6 != "many") { fmt[$1, $2] = $3; if ($3 == "gzip") hasgz = 1 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7; variant[$1] = ($10 != "" && $10 != $1) } next }
        FNR == 1 || !(($1, $2) in fmt) || variant[$3] { next }
        { w[$1, $2, $3] = $5; m[$1, $2, $3] = $8; lanes[$1, $3] = 1; cs[$1, $2] = 1 }
        END {
            moded = mode[1] != ""; nfmt = hasgz ? 2 : 1
            print (moded ? modehead "\t" : "") "Tool\tPlain time\tPlain memory" (hasgz ? "\tGzip time\tGzip memory" : "") (through ? "\tPasses malformed" : "")
            print (moded ? "l\t" : "") "l\tr\tr" (hasgz ? "\tr\tr" : "") (through ? "\tl" : "")
            for (x = 1; x <= nf; x++) {
                F = fam[x]; n = 0; delete list
                for (k in lanes) { split(k, q, SUBSEP); if (q[1] == F) list[++n] = q[2] }
                for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
                for (i = 1; i <= n; i++) {
                    l = list[i]; line = (moded ? mode[x] "\t" : "") label[l]; gsub("_", " ", line)
                    split("plain gzip", fs, " ")
                    for (f = 1; f <= nfmt; f++) {
                        k = sw = sm = 0
                        for (cc in cs) { split(cc, q, SUBSEP); if (q[1] != F) continue; c2 = q[2]
                            if (fmt[F, c2] == fs[f] && (F, c2, l) in w && (F, c2, "z-fastq") in w) {
                                k++
                                if (l == "z-fastq") { sw += log(w[F, c2, l]); sm += log(m[F, c2, l]) }
                                else { sw += log(w[F, c2, l] / w[F, c2, "z-fastq"]); sm += log(m[F, c2, l] / m[F, c2, "z-fastq"]) } } }
                        if (!k) line = line "\t-\t-"
                        else if (l == "z-fastq") line = line "\t" secs(exp(sw / k)) "\t" mib(exp(sm / k))
                        else line = line "\t" sprintf("%.2fx", exp(sw / k)) "\t" sprintf("%.2fx", exp(sm / k))
                    }
                    print line (through ? "\t" (((F, l) in lets) ? lets[F, l] : (("", l) in lets) ? lets["", l] : "-") : "")
                }
            }
        }' "$REPORT_DIR/through.tsv" "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# z-fastq's own geometric means over the headline inputs, as one sentence.
subject_sentence() { # family
    awk -F'\t' -v fam="$1" -v tier="$HEADLINE" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") fmt[$2] = $3; next }
        FNR == 1 || $1 != fam || !($2 in fmt) || $3 != "z-fastq" { next }
        { k[fmt[$2]]++; w[fmt[$2]] += log($5); m[fmt[$2]] += log($8)
          if ($11 > 0) { nc[fmt[$2]]++; c[fmt[$2]] += log($11) } }
        END {
            if (!k["gzip"]) {
                printf "z-fastq itself takes %.1f ms and %.2f MiB per plain input (geometric mean over %d inputs); the table gives every other tool as a multiple of that.\n",
                    exp(w["plain"] / k["plain"]) / 1e6, exp(m["plain"] / k["plain"]) / 1024, k["plain"]
                if (nc["plain"] == k["plain"]) printf "Geometric-mean process CPU time is %.1f ms. Elapsed time also includes launcher overhead and time waiting to run or for I/O.\n", exp(c["plain"] / nc["plain"]) / 1000
                exit
            }
            printf "z-fastq itself takes %.1f ms and %.2f MiB per plain input, %.1f ms and %.2f MiB per gzip input (geometric means over %s); the table gives every other tool as a multiple of that.\n",
                exp(w["plain"] / k["plain"]) / 1e6, exp(m["plain"] / k["plain"]) / 1024, exp(w["gzip"] / k["gzip"]) / 1e6, exp(m["gzip"] / k["gzip"]) / 1024,
                (k["plain"] == k["gzip"] ? k["plain"] " inputs each" : k["plain"] " and " k["gzip"] " inputs")
            if (nc["plain"] == k["plain"] && nc["gzip"] == k["gzip"]) printf "Geometric-mean process CPU time is %.1f ms for plain input and %.1f ms for gzip. Elapsed time also includes launcher overhead and time waiting to run or for I/O.\n", exp(c["plain"] / nc["plain"]) / 1000, exp(c["gzip"] / nc["gzip"]) / 1000
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/summary.tsv"
}
# Per lane, user-space counters per input byte, geometric mean over the plain headline inputs.
# The family's name in report text: its title when the feature gives one.
family_name() { printf '%s\n' "${FAMILY_TITLE[$1]:-$1}"; } # family
# A tool's thread variants in FAMILY (in any family when none is given), as
# "; also with -w 2, -w 4" after its command.
thread_options() { # lane [family]
    local name out="" pool
    if [[ -n "${2:-}" ]]; then pool="${FAM_LANES[$2]}"; else pool="${LANES[*]}"; fi
    for name in $pool; do
        [[ "${LANE_OF[$name]:-}" == "$1" ]] || continue
        out+="${out:+, }${REG_LABEL[$name]#"${REG_LABEL[$1]} "}"
    done
    [[ -z "$out" ]] || printf '; also with %s' "$out"
}
counters_table() { # family...
    local fams fam
    fams="$(for fam; do printf '%s\t%s\n' "$fam" "$(family_name "$fam")"; done)"
    awk -F'\t' -v tier="$HEADLINE" -v jobs="$#" '
        NR == FNR { want[$1] = 1; title[$1] = $2; fl[++nfam] = $1; next }
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && ($1 in want) && $5 == tier && $6 != "many" && $3 == "plain") bytes[$1, $2] = $7; next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7; variant[$1] = ($10 != "" && $10 != $1) } next }
        FNR == 1 || !(($1, $2) in bytes) || variant[$3] || $15 <= 0 || $16 <= 0 { next }
        { b = bytes[$1, $2]; mib = b / 1048576; key = $1 SUBSEP $3; k[key]++; seen[$1, $3] = 1
          ins[key] += log($15 / b); ipc[key] += log($15 / $16); cm[key] += log(($18 + 1) / mib); bm[key] += log(($19 + 1) / mib) }
        END {
            print (jobs > 1 ? "Job\t" : "") "Tool\tInstructions per byte\tInstructions per cycle\tCache misses per MiB\tBranch misses per MiB"
            print (jobs > 1 ? "l\t" : "") "l\tr\tr\tr\tr"
            for (f = 1; f <= nfam; f++) { F = fl[f]; n = 0; delete list
                for (s in seen) { split(s, q, SUBSEP); if (q[1] == F) list[++n] = q[2] }
                for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
                for (i = 1; i <= n; i++) { l = list[i]; key = F SUBSEP l; name = label[l]; gsub("_", " ", name)
                    printf "%s%s\t%.2f\t%.2f\t%.0f\t%.0f\n", (jobs > 1 ? title[F] "\t" : ""), name, exp(ins[key] / k[key]), exp(ipc[key] / k[key]), exp(cm[key] / k[key]), exp(bm[key] / k[key]) } }
        }' <(printf '%s\n' "$fams") "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# The headline findings, computed from the report tables so they stay true
# after every run: correctness, speed against the fastest full-check tool,
# memory against every tool, and the kinds of input where z-fastq is not
# fastest. PARTS picks the bullets (right speed memory where); LABEL prefixes
# each bullet's name; ANY=1 names the fastest tool of any class before the
# fastest full-check one.
findings() { # family [parts] [label]
    awk -F'\t' -v fam="$1" -v tier="$HEADLINE" -v right="${WORD[right]}" -v missed="${WORD[missed]}" -v wrong="${WORD[wrong]}" \
        -v parts=" ${2:-right speed memory where} " -v lbl="${3:-}" -v any="${ANY:-0}" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && $5 == tier && $6 != "many") { fmt[$2] = $3; grp[$2] = $10 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; cls[$1] = $3; ord[$1] = $7 } next }
        FILENAME ~ /checks.tsv$/ { if ($1 == "z-fastq") { if ($3 == "pass") zp++; else zf++ } next }
        FILENAME ~ /edges.tsv$/ {
            if (FNR == 1 || $6 == "policy") next
            if ($4 == "z-fastq") { if ($6 != "right") zbad++ } else if ($6 == "wrong" || $6 == "crash") bad[$4] = 1
            next
        }
        FNR == 1 || $1 != fam || !($2 in fmt) { next }
        { w[$2, $3] = $5; m[$2, $3] = $8; lanes[$3] = 1 }
        END {
            n = 0; for (l in lanes) if (l != "z-fastq") list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            # Correctness.
            line = sprintf("- **Right answers:** z-fastq returned %s on all %d checked inputs", right, zp)
            if (zf) line = sprintf("- **Right answers:** z-fastq %s on %d of %d checked inputs", missed, zf, zp + zf)
            line = line (zbad ? "." : " and rejected every malformed edge case.")
            k = 0; for (i = 1; i <= n; i++) if (label[list[i]] in bad) k++
            if (k) line = line sprintf(" %d of the %d other tools print %s or crash on at least one.", k, n, wrong)
            if (index(parts, " right ")) print line
            # Speed against the fastest complete lane, per format.
            split("plain gzip", fs, " ")
            for (f = 1; f <= 2; f++) {
                best[f] = ""; br[f] = 0
                for (i = 1; i <= n; i++) if (cls[list[i]] == "complete") { r = ratio(list[i], fs[f], w); if (r && (best[f] == "" || r < br[f])) { best[f] = list[i]; br[f] = r } }
            }
            if (index(parts, " speed ")) {
                head = "- **" (lbl == "" ? "Speed" : lbl " speed") ":** "
                if (any) {
                    # The fastest tool of any class first, per format.
                    for (f = 1; f <= 2; f++) { fa[f] = ""; fr[f] = 0
                        for (i = 1; i <= n; i++) { r = ratio(list[i], fs[f], w); if (r && (fa[f] == "" || r < fr[f])) { fa[f] = list[i]; fr[f] = r } } }
                    head = head sprintf("the fastest other tool takes %.2fx z-fastq%ss time on plain input (%s)", fr[1], "\047", label[fa[1]]) (fr[2] ? sprintf(" and %.2fx on gzip (%s)", fr[2], label[fa[2]]) : "") "; "
                }
                if (best[1] == "") print head "no other full-check tool runs it."
                else if (br[2] == 0)
                    printf "%s%s, the fastest other full-check tool, takes %.2fx z-fastq%ss time on plain input.\n", head, label[best[1]], br[1], "\047"
                else if (best[1] == best[2])
                    printf "%s%s, the fastest other full-check tool, takes %.2fx z-fastq%ss time on plain input and %.2fx on gzip.\n", head, label[best[1]], br[1], "\047", br[2]
                else
                    printf "%sthe fastest other full-check tool takes %.2fx z-fastq%ss time on plain input (%s) and %.2fx on gzip (%s).\n", head, br[1], "\047", label[best[1]], br[2], label[best[2]]
            }
            # Memory against every other lane, both formats.
            lo = hi = 0
            for (i = 1; i <= n; i++) for (f = 1; f <= 2; f++) { r = ratio(list[i], fs[f], m); if (r) { if (!lo || r < lo) lo = r; if (r > hi) hi = r } }
            for (c in fmt) if ((c, "z-fastq") in m) { zm += log(m[c, "z-fastq"]); zk++; if (fmt[c] == "gzip") zg = 1 }
            if (index(parts, " memory ")) {
                span = (hi < lo * 1.05) ? sprintf("%.1fx", lo) : sprintf("%.1fx to %.0fx", lo, hi)
                printf "- **Memory%s:** %s peaks at %s z-fastq%ss %.2f MiB.\n", (lbl == "" ? "" : ", " lc1(lbl)), (n == 1 ? label[list[1]] : "every other tool"), span, "\047", exp(zm / zk) / 1024
            }
            # Kinds of input where some tool is faster than z-fastq.
            out = ""
            for (c in fmt) gs[grp[c]] = 1
            for (g in gs) for (f = 1; f <= 2; f++) for (i = 1; i <= n; i++) {
                k = s = 0
                for (c in fmt) if (grp[c] == g && fmt[c] == fs[f] && (c, list[i]) in w && (c, "z-fastq") in w) { k++; s += log(w[c, list[i]] / w[c, "z-fastq"]) }
                if (k && exp(s / k) < 0.99) out = out sprintf("%s%s on %s %s (%.2fx)", (out == "" ? "" : "; "), label[list[i]], fs[f], tolower(g), exp(s / k))
            }
            where = "- **Where z-fastq is not fastest" (lbl == "" ? "" : ", " lc1(lbl)) ":** "
            if (!index(parts, " where ")) { }
            else if (out == "") print where "on no kind of input; it is the fastest tool on every kind of input (geometric mean over the group), " (zg ? "plain and gzip." : "plain.")
            else print where out "."
        }
        # A label inside a sentence: its first letter lowercase, unless it begins an acronym (QC).
        function lc1(s) { return substr(s, 2, 1) ~ /[A-Z]/ ? s : tolower(substr(s, 1, 1)) substr(s, 2) }
        function ratio(l, f, a,   c, k, s) { for (c in fmt) if (fmt[c] == f && (c, l) in a && (c, "z-fastq") in a) { k++; s += log(a[c, l] / a[c, "z-fastq"]) } return k ? exp(s / k) : 0 }
    ' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/checks.tsv" "$REPORT_DIR/edges.tsv" "$REPORT_DIR/summary.tsv"
}
# Two sentences under the scaling figures, computed from the scaling cases:
# z-fastq's throughput at the largest size, who passes it there, and whose
# peak memory grows with the input.
scaling_sentences() { # family
    local range
    range="$(scaling_range)"
    awk -F'\t' -v fam="$1" -v srcs=" $SCALING " -v toplabel="${range##* to }" -v verb="${WORD[verb]}" '
        FILENAME ~ /cases.tsv$/ { if (FNR > 1 && $1 == fam && index(srcs, " " $4 " ")) { fmt[$2] = $3; src[$2] = $4; size[$2] = $7 } next }
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7; variant[$1] = ($10 != "" && $10 != $1) } next }
        FNR == 1 || $1 != fam || !($2 in fmt) || variant[$3] { next }
        { tp = size[$2] / ($5 / 1e9); v[$2, $3] = tp; r[$2, $3] = $8; lanes[$3] = 1
          if (size[$2] > top[src[$2], fmt[$2]]) { top[src[$2], fmt[$2]] = size[$2]; big[src[$2], fmt[$2]] = $2 } }
        END {
            n = 0; for (l in lanes) if (l != "z-fastq") list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            split("plain gzip", fs, " ")
            line = "At " toplabel ", z-fastq " verb
            for (f = 1; f <= 2; f++) {
                lo = hi = 0
                for (k in big) { split(k, p, SUBSEP); if (p[2] != fs[f]) continue; c = big[k]; t = v[c, "z-fastq"] / 2^30
                    if (!lo || t < lo) lo = t; if (t > hi) hi = t }
                rate = sprintf("%.1f", lo) == sprintf("%.1f", hi) ? sprintf("%.1f", hi) : sprintf("%.1f to %.1f", lo, hi)
                line = line sprintf("%s %s GiB/s on %s input", (f > 1 ? " and" : ""), rate, fs[f])
            }
            out = ""
            for (i = 1; i <= n; i++) for (k in big) { c = big[k]; split(k, p, SUBSEP); l = list[i]
                if (v[c, l] > v[c, "z-fastq"] * 1.01)
                    out = out (out == "" ? "" : "; ") sprintf("%s on %s %s (%.2fx)", label[l], p[2], p[1], v[c, l] / v[c, "z-fastq"]) }
            print line (out == "" ? "; no other tool is faster there." : ". Faster on a single shape there, as throughput over z-fastq'"'"'s: " out ".")
            zlo = zhi = 0; grow = ""
            for (c in fmt) if ((c, "z-fastq") in r) { x = r[c, "z-fastq"]; if (!zlo || x < zlo) zlo = x; if (x > zhi) zhi = x }
            for (i = 1; i <= n; i++) { l = list[i]; g = 0
                for (k in big) { split(k, p, SUBSEP); small = 0
                    for (c in fmt) if (src[c] == p[1] && fmt[c] == p[2] && (c, l) in r && (!small || size[c] < size[sm])) { small = 1; sm = c }
                    if (small && (big[k], l) in r) { q = r[big[k], l] / r[sm, l]; if (q > g) { g = q; gk = p[2] " " p[1] } } }
                if (g >= 1.25) grow = grow (grow == "" ? "" : "; ") sprintf("%s by %.1fx on %s", label[l], g, gk) }
            # The shapes these lines cover; peaks closer than 512 KiB are one value, the precision of the RSS counter.
            ns = 0; for (c in fmt) if (!(src[c] in seen)) { seen[src[c]] = 1; shp[++ns] = src[c] }
            for (i = 1; i <= ns; i++) for (j = i + 1; j <= ns; j++) if (shp[j] < shp[i]) { t = shp[i]; shp[i] = shp[j]; shp[j] = t }
            names = ""; for (i = 1; i <= ns; i++) names = names (i == 1 ? "" : i == ns ? " and " : ", ") shp[i]
            of = ns == 1 ? "the " names " shape" : "the " names " shapes"
            if (zhi - zlo < 512) printf "z-fastq peaks at %.2f MiB at every size of %s.", zhi / 1024, of
            else printf "z-fastq peaks at %.2f to %.2f MiB across the sizes of %s.", zlo / 1024, zhi / 1024, of
            print (grow == "" ? " No other tool grows by more than a quarter either." : " Memory that grows with the input: " grow ".")
        }' "$REPORT_DIR/cases.tsv" "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
}
# A CAUTION callout naming the tools that print a wrong answer with exit 0 or
# crash on malformed edge cases; nothing when none does.
caution() {
    awk -F'\t' -v wrong="${WORD[wrong]}" -v unseen="${WORD[unseen]}" 'FILENAME ~ /style.tsv$/ { if (FNR > 1) ord[$2] = $7 + 0; next }
        FNR > 1 && $6 == "wrong" { w[$4] = 1 } FNR > 1 && $6 == "crash" { c[$4] = 1 }
        END {
            for (t in w) k1++
            for (t in c) k2++
            wl = names(w); cl = names(c)
            if (wl == "" && cl == "") exit
            printf "> [!CAUTION]\n> On the edge cases below,"
            if (wl != "") printf " %s %s %s and %s 0 on some malformed inputs, so a script cannot tell %s.", wl, (k1 > 1 ? "print" : "prints"), wrong, (k1 > 1 ? "exit" : "exits"), unseen
            if (cl != "") printf " %s %s on some malformed gzip input.", cl, (k2 > 1 ? "crash" : "crashes")
            printf "\n\n"
        }
        # The tools in a set, in registry order, joined with commas and "and".
        function names(set,   t, list, i, j, n, out) {
            n = 0; for (t in set) list[++n] = t
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= n; i++) out = out (i == 1 ? "" : i == n ? " and " : ", ") list[i]
            return out
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/edges.tsv"
}
# Per lane, its time and memory reading stdin as a multiple of its own on the
# same file by path, geometric mean over the shared shapes.
stdin_table() { # family
    awk -F'\t' -v fam="$1" -v base="${FAMILIES[0]}" '
        FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 } next }
        FNR == 1 { next }
        $1 == fam { sw[$2, $3] = $5; sm[$2, $3] = $8; lanes[$3] = 1; cs[$2] = 1; next }
        $1 == base { fw[$2, $3] = $5; fm[$2, $3] = $8 }
        END {
            print "Tool\tPlain time\tPlain memory\tGzip time\tGzip memory"
            print "l\tr\tr\tr\tr"
            n = 0; for (l in lanes) list[++n] = l
            for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (ord[list[j]] < ord[list[i]]) { t = list[i]; list[i] = list[j]; list[j] = t }
            for (i = 1; i <= n; i++) {
                l = list[i]; line = label[l]
                split("plain gzip", fs, " ")
                for (f = 1; f <= 2; f++) {
                    k = a = b = 0
                    for (c in cs) { p = split(c, q, ":"); if (q[p] == fs[f] && (c, l) in sw && (c, l) in fw) { k++; a += log(sw[c, l] / fw[c, l]); b += log(sm[c, l] / fm[c, l]) } }
                    line = line "\t" (k ? sprintf("%.2fx", exp(a / k)) : "-") "\t" (k ? sprintf("%.2fx", exp(b / k)) : "-")
                }
                print line
            }
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv" | mdtable
}
# One sentence for a family timed as one call over many files.
many_sentence() { # family
    local files bytes
    read -r files bytes < <(cases | awk -F'\t' -v f="$1" '$1 == f { print split($4, p, " "), $6; exit }')
    awk -F'\t' -v fam="$1" -v files="$files" -v bytes="$bytes" -v gerund="${WORD[many]}" -v verb="${WORD[verb]}" '
        FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
        FNR == 1 || $1 != fam { next }
        { n = split($2, p, ":"); f = p[n]; w[$3, f] = $5; m[$3, f] = $8; if ($3 != "z-fastq") other[$3] = 1 }
        END {
            if (!(("z-fastq", "plain") in w) || !(("z-fastq", "gzip") in w)) { print "Not timed yet."; exit }
            printf "%s all %d xsmall corpus files (%.0f MiB) in one call, z-fastq takes %.1f ms on plain input and %.1f ms on gzip, at %.2f MiB peak.", gerund, files, bytes / 1048576, w["z-fastq", "plain"] / 1e6, w["z-fastq", "gzip"] / 1e6, m["z-fastq", "plain"] / 1024
            k = 0; for (o in other) k++
            for (o in other) printf " %s%s takes %.1fx and %.1fx as long with %.0fx and %.0fx the memory.", label[o], (k == 1 ? ", the only other tool that " verb " many files in one call," : ""), w[o, "plain"] / w["z-fastq", "plain"], w[o, "gzip"] / w["z-fastq", "gzip"], m[o, "plain"] / m["z-fastq", "plain"], m[o, "gzip"] / m["z-fastq", "gzip"]
            printf "\n"
        }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/summary.tsv"
}
# Per lane (and per family, for a feature whose edge rows each run in one
# family), the malformed edge inputs it lets through with exit 0, as short
# names ("nothing" when it rejects them all), and their count: "N of M" or
# "none". One line per lane: lane, family (empty for all), names, count.
lets_through() {
    local name fam fams=""
    declare -F edge_family >/dev/null && fams="${FAMILIES[*]}"
    for name in "${LANES[@]}"; do
        for fam in "" $fams; do
            awk -F'\t' -v l="$name" -v f="$fam" -v rowfam="$(edge_list | while IFS=$'\t' read -r e _; do printf '%s=%s ' "$e" "$(declare -F edge_family >/dev/null && edge_family "$e")"; done)" 'BEGIN {
                    short["blank"] = "blank line"; short["lenmis"] = "length mismatch"; short["trunc"] = "truncated record"
                    short["dash"] = "bad plus line"; short["qlow"] = "quality range"; short["qhigh"] = "quality range"
                    short["noat"] = "bad header"; short["gztrunc"] = "broken gzip"; short["gzgarbage"] = "broken gzip"
                    n = split(rowfam, rf, " "); for (i = 1; i <= n; i++) { split(rf[i], kv, "="); of[kv[1]] = kv[2] } }
                FNR == NR { if ($2 ~ /^reject/) { want[$1] = 1; label[$1] = $3 } next }
                $1 == l && ($2 in want) && $3 != "-" && (f == "" || of[$2] == f) { total++
                    if ($3 ~ /^wrong/) { passed++; k = ($2 in short) ? short[$2] : label[$2]; if (!(k in seen)) { seen[k] = 1; out = out (out == "" ? "" : ", ") k } } }
                END { if (f == "" || total) printf "%s\t%s\t%s\t%s\n", l, f, (out == "" ? "nothing" : out), (passed ? passed " of " total : "none") }' <(edge_list) "$RESULTS/edges.tsv"
        done
    done
}
# One findings bullet: each valid edge input on which some tool exits 0 but
# writes the wrong output, with those tools (in table order) and how; nothing
# when none does.
valid_wrong_bullet() {
    awk -F'\t' 'FILENAME ~ /style.tsv$/ { if (FNR > 1) { label[$1] = $2; ord[$1] = $7 + 0 } next }
        FNR == NR { if ($2 !~ /^reject/ && $2 != "policy") { ok[$1] = 1; what[$1] = $3; row[$1] = ++rows; name[rows] = $1 } next }
        ($2 in ok) && $3 ~ /^wrong/ { s = $3; sub(/^wrong: /, "", s); if (s ~ /^[0-9]+$/) s = "printed " s; else if (s ~ /^[a-z_]+$/) s = s " wrong"
            n = ++cnt[$2]; tool[$2, n] = $1; how[$2, n] = s }
        END {
            for (r = 1; r <= rows; r++) { e = name[r]; if (!(e in cnt)) continue
                m = cnt[e]
                for (a = 1; a <= m; a++) for (b = a + 1; b <= m; b++) if (ord[tool[e, b]] < ord[tool[e, a]]) {
                    t = tool[e, a]; tool[e, a] = tool[e, b]; tool[e, b] = t; t = how[e, a]; how[e, a] = how[e, b]; how[e, b] = t }
                part = what[e] ": "
                for (a = 1; a <= m; a++) part = part (a > 1 ? ", " : "") label[tool[e, a]] ((a == m || how[e, a] != how[e, a + 1]) ? " (" how[e, a] ")" : "")
                out = out (out == "" ? "" : "; ") part }
            if (out != "") { gsub(/`/, "", out); print "- **Wrong on valid input, exit 0:** " out "." } }' \
        <(edge_list) "$REPORT_DIR/style.tsv" "$RESULTS/edges.tsv"
}
# Per lane, "N of M": the malformed edge inputs it rejects, over every one it ran.
rejects_all() { # lane
    awk -F'\t' -v l="$1" 'FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next } $1 == l && ($2 in want) && $3 != "-" { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' <(edge_list) "$RESULTS/edges.tsv"
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
# One row per lane: each kind of malformed edge input (SPEC is
# Label=edge,edge|...) as rejects, or **passes** when the lane exits 0 on any
# of its edge inputs; then the columns EXTRA (tab-separated headers) that
# FUNCTION prints for the lane, or none.
kinds_table() { # spec [extra-headers function]
    local name spec="$1" extra="${2:-}" fn="${3:-}" head="" align="" tail
    head="$(tr '|' '\n' <<<"$spec" | cut -d= -f1 | paste -sd$'\t' -)"
    align="$(tr '|' '\n' <<<"$spec" | sed 's/.*/l/' | paste -sd$'\t' -)"
    {
        printf 'Tool\t%s%s\n' "$head" "${extra:+$'\t'$extra}"
        printf 'l\t%s%s\n' "$align" "${extra:+$'\t'$(tr '\t' '\n' <<<"$extra" | sed 's/.*/l/' | paste -sd$'\t' -)}"
        for name in "${LANES[@]}"; do
            [[ ! -v "LANE_OF[$name]" ]] || continue
            tail=""
            [[ -z "$fn" ]] || tail="$("$fn" "$name")"
            awk -F'\t' -v l="$name" -v label="${REG_LABEL[$name]}" -v spec="$spec" -v tail="$tail" '
                BEGIN { nk = split(spec, kinds, "|"); for (i = 1; i <= nk; i++) sub(/^[^=]*=/, "", kinds[i]) }
                FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next }
                $1 == l && ($2 in want) { cell[$2] = $3 }
                END {
                    line = label
                    for (i = 1; i <= nk; i++) { n = split(kinds[i], e, ","); passed = 0; seen = 0
                        for (j = 1; j <= n; j++) { if (e[j] in cell) seen = 1; if (cell[e[j]] ~ /^wrong/) passed = 1 }
                        line = line "\t" (!seen || cell[e[1]] == "-" ? "-" : passed ? "**passes**" : "rejects") }
                    print line (tail == "" ? "" : "\t" tail)
                }' <(edge_list) "$RESULTS/edges.tsv"
        done
    } | mdtable
}
# Per lane, "N of M": the malformed edge inputs it rejects, over every one in the set.
rejects_all() { # lane
    awk -F'\t' -v l="$1" 'FNR == NR { if ($2 ~ /^reject/) want[$1] = 1; next } $1 == l && ($2 in want) { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' <(edge_list) "$RESULTS/edges.tsv"
}
probes() { # lane -> "N of M", over the edge set's probe rows
    awk -F'\t' -v l="$1" -v rows=" $(probe_rows | tr '\n' ' ')" '$1 == l && index(rows, " " $2 " ") { n++; if ($3 ~ /^reject/) r++ } END { printf "%d of %d", r, n }' "$RESULTS/edges.tsv"
}
# One sentence on z-fastq, then the other tools wrong on inputs they were
# given, as a table (tool, wrong, of, inputs), or a sentence when none was.
correctness_sentence() {
    local name pass fail rows="" zf="" where
    for name in "${LANES[@]}"; do
        read -r pass fail < <(awk -F'\t' -v l="$name" '$1 == l { if ($3 == "pass") a++; else b++ } END { print a + 0, b + 0 }' "$REPORT_DIR/checks.tsv")
        if [[ "$name" == z-fastq ]]; then
            local t top=""
            for t in tiny xsmall small medium large; do verify_cases | cut -f8 | grep -qx "$t" && top="${TIER_SIZE[$t]}"; done
            if ((pass + fail == 0)); then zf="z-fastq has no current validation results."
            elif ((VERIFY)); then zf="z-fastq returned ${WORD[right]} on all $pass inputs: every corpus file and synthetic shape on disk, ${TIER_SIZE[tiny]} to $top, plain and gzip."
            else zf="z-fastq returned ${WORD[right]} on all $pass inputs it was timed on."; fi
            ((fail == 0)) || zf="z-fastq was wrong on $fail of $((pass + fail)) inputs."
            continue
        fi
        ((fail > 0)) || continue
        where="$(awk -F'\t' -v l="$name" '$1 == l && $3 == "fail" { print $2 }' "$REPORT_DIR/checks.tsv" | sort -u |
            while read -r src; do if [[ -d "$SHAPES/$src" ]]; then printf '%s, ' "$src"; else printf '%s, ' "${DESCRIBE[$src]:-$src}"; fi; done)"
        rows+="${REG_LABEL[$name]}"$'\t'"$fail"$'\t'"$((pass + fail))"$'\t'"${where%, }"$'\n'
    done
    printf '%s ' "$zf"
    if [[ -z "$rows" ]]; then
        printf 'Every other tool was right on every input it was timed on.\n'
    else
        printf 'Every other tool was right on every input it was timed on except these, which are not timed where they %s:\n\n' "${WORD[missed]}"
        { printf 'Tool\tWrong\tOf\tInputs\n'; printf 'l\tr\tr\tl\n'; printf '%s' "$rows"; } | mdtable
    fi
}
# The edge-case matrix for plot.py: row, label, expected, tool, cell, verdict.
edges_tsv() {
    printf 'row\tlabel\texpected\ttool\tcell\tverdict\n'
    local row=0 edge expected label name cell verdict
    while IFS=$'\t' read -r edge expected label; do
        for name in "${LANES[@]}"; do
            [[ ! -v "LANE_OF[$name]" ]] || continue
            cell="$(awk -F'\t' -v l="$name" -v e="$edge" '$1 == l && $2 == e { print $3 }' "$RESULTS/edges.tsv")"
            case "$cell" in
                wrong*) verdict=wrong ;;
                crash*) verdict=crash ;;
                refused) verdict=refused ;;
                -) verdict=policy ;;
                *) verdict=right ;;
            esac
            [[ "$expected" != policy ]] || verdict=policy
            # Expected facts show as ok, as their right cells do.
            [[ -z "$AGREE_LANE" || "$expected" != *" "* ]] || expected=ok
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$row" "${label//\`/}" "${expected/policy/tool policy}" "${REG_LABEL[$name]}" "$cell" "$verdict"
        done
        row=$((row + 1))
    done < <(edge_list)
}
details() { printf '<details><summary>%s</summary>\n\n' "$1"; }
# The size span of the scaling inputs, as "1 MiB to 10 GiB".
scaling_range() {
    local tiers
    tiers="$(cases | awk -F'\t' -v s=" $SCALING " 'index(s, " " $7 " ") { print $8 }' | sort -u)"
    local t first="" last="" name
    for t in tiny xsmall small medium large; do
        grep -qx "$t" <<<"$tiers" || continue
        name="${TIER_SIZE[$t]}"
        [[ -n "$first" ]] || first="$name"
        last="$name"
    done
    printf '%s to %s\n' "$first" "$last"
}
# Lanes timed at the headline tier only, as "fqtools".
headline_only() { local name out=""; for name in "${LANES[@]}"; do [[ "${LANE_HEADLINE[$name]}" == 0 || -v "LANE_OF[$name]" ]] || out+="${REG_LABEL[$name]}, "; done; printf '%s\n' "${out%, }"; }
# The inputs behind every figure: one line per kind, then the scaling sources.
data_list() {
    local kind src reals=0 shapes=0 lines="" synth="" used
    used=" $(cases | awk -F'\t' -v f="${FAMILIES[0]}" -v t="$HEADLINE" '$1 == f && $8 == t { print $7 }' | sort -u | tr '\n' ' ')"
    for kind in "${GROUP_LIST[@]}"; do
        local real=""
        for src in ${GROUP_MEMBERS[$kind]}; do
            [[ "$used" == *" $src "* ]] || continue
            if [[ -d "$SHAPES/$src" ]]; then
                synth+="    - \`$src\`: *${DESCRIBE[$src]}*"$'\n' shapes=$((shapes + 1))
            else
                real+="${DESCRIBE[$src]}; " reals=$((reals + 1))
            fi
        done
        [[ -z "$real" ]] || lines+="    - \`$kind\`: *${real%; }*"$'\n'
    done
    printf -- '- `data`:\n'
    printf -- '  - **Real**: *%d public sequencing runs, %s each, pooled in the bars by the kind of reads they hold*\n%s' "$reals" "${TIER_SIZE[$HEADLINE]}" "$lines"
    printf -- '  - **Synthetic**: *%d %s generated from real reads, each changing one thing; one group in the bars at %s, and the scaling lines from %s*\n%s' "$shapes" "$( ((shapes == 1)) && echo shape || echo shapes)" "${TIER_SIZE[$HEADLINE]}" "$(scaling_range)" "$synth"
}
# The source revisions of the z-fastq binaries the current batches timed, as
# recorded when each was built, with batch counts when there are several.
zfastq_revision() {
    local dir
    awk -F'\t' '{ print $3 }' "$RESULTS/current.tsv" | sort -u | while IFS= read -r dir; do
        sed -n 's/^lane z-fastq \([^#]*\)#.*/\1/p' "$dir/meta.tsv"
    done | awk -F'\t' -v revs="$CACHE/z-fastq.revisions" 'BEGIN { while ((getline line < revs) > 0) { split(line, f, "\t"); rev[f[1]] = f[2] } }
        { n[($1 in rev) ? rev[$1] : "revision not recorded"]++ }
        END { for (r in n) { k++; list[k] = r }
            if (k == 1) { print (list[1] == "revision not recorded" ? "" : "at ") list[1]; exit }
            for (i = 1; i <= k; i++) printf "%s%s (%d batches)", (i > 1 ? "; " : "at "), list[i], n[list[i]]; print "" }'
}
# z-fastq as the current batches timed it: the version in build.zig.zon when
# every batch timed the current source, then the source revision. The binary
# on disk is not asked, since it may be a later build.
zfastq_label() {
    local rev ids dir version=""
    rev="$(zfastq_revision)"
    ids="$(awk -F'\t' '{ print $3 }' "$RESULTS/current.tsv" | sort -u | while IFS= read -r dir; do
        sed -n 's/^lane z-fastq \([^#]*\)#.*/\1/p' "$dir/meta.tsv"
    done | sort -u)"
    [[ "$ids" != "$(zfastq_id)" ]] || version="$(sed -n 's/^ *\.version = "\([^"]*\)",*$/\1/p' "$REPO/build.zig.zon")"
    printf '%s%s\n' "${version:+$version, }" "${rev:+source $rev}"
}
# The README's sections, each printing its Markdown; readme() puts them in
# the default order, and a feature with its own feature_readme orders them.
readme_head() {
    local first="${FAMILIES[0]}" name
    local date rev
    date="$(date -u +%Y-%m-%d)"
    rev="$(zfastq_label)"
    local counters reals shapes
    counters=0
    for name in "${LANES[@]}"; do
        # A thread variant is its tool again, not another tool.
        [[ -v "LANE_OF[$name]" ]] || case "${LANE_CLASS[$name]}" in complete | partial) counters=$((counters + 1)) ;; esac
    done
    # Sources of the first family at the headline tier: a shape has a synthetic directory.
    local used src fmts
    used="$(cases | awk -F'\t' -v f="$first" -v t="$HEADLINE" '$1 == f && $8 == t && $9 != "many" { print $7 }' | sort -u)"
    reals=0 shapes=0
    for src in $used; do if [[ -d "$SHAPES/$src" ]]; then shapes=$((shapes + 1)); else reals=$((reals + 1)); fi; done
    fmts="$(cases | awk -F'\t' '$3 == "gzip" { g = 1 } END { print g ? "plain and gzip" : "plain" }')"
    local cpu
    cpu="$(awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo | sed 's/^AMD //; s/^Intel(R) Core(TM) //; s/ [0-9]*-Core Processor$//')"
    local per="one CPU per run" on="on one CPU"
    if ((${#LANE_CORES[@]})); then per="one CPU per run, threaded variants on as many cores as threads" on="on one CPU, a threaded variant on as many cores as its threads"; fi
    cat <<EOF
<!-- markdownlint-disable MD013 MD033 MD041 -->

# ${WORD[title]}
> <sub>\`${WORD[subject]}\` against $counters other ${WORD[others]}, on $reals real sequencing run$( ((reals == 1)) || echo s) and $shapes synthetic shape$( ((shapes == 1)) || echo s), $fmts. A tool is timed only where ${WORD[timed]}; time and memory come from the same run, $on.</sub>

- \`z-fastq\`: *$rev*
- \`machine\`: *$cpu, Linux $(uname -r | cut -d- -f1); $per*
- \`runner\`: *$(zebrac --version) with zexec*
$(data_list)
- \`date\`: *$date*

EOF
}
readme_tools() {
    local name
    cat <<EOF
## Compared tools
> <sub>*complete* tools check as much as z-fastq; *partial* tools check less, so they do less work. ${WORD[probes]}</sub>

EOF
    for name in "${LANES[@]}"; do
        local what="${LANE_NOTE[$name]}" probe=""
        case "${LANE_CLASS[$name]}" in
            subject | complete) what="checks every record's structure${what:+; $what}" ;;
        esac
        probe="; **$(probes "$name") rejected**"
        printf -- '- `%s %s` **(%s)**; *%s*%s\n' "${REG_LABEL[$name]}" "$(lane_version "$name")" "${LANE_CLASS[$name]}" "$what" "$probe"
    done
    if [[ " ${LANE_CLASS[*]} " == *" partial "* ]]; then
        printf '\n> [!IMPORTANT]\n> Partial tools skip checks z-fastq makes, so they do less work per record; their times are not like-for-like with z-fastq or the complete tools.\n'
    fi
}
readme_right() {
    printf '\n## Right answers\n> <sub>%s The edge cases are small inputs where tools differ: each cell shows what the tool printed, and its color whether that was right.</sub>\n\n' "${WORD[truth]}"
    correctness_sentence
    printf '\n'
    caution
    picture "edges" "What each tool does on small edge-case inputs"
    printf '\n'
    printf '> [!NOTE]\n> %s "tool policy" rows have no single right answer.\n\n' "${WORD[cells]}"
    if awk -F'\t' '$3 == "fail" { found = 1 } END { exit !found }' "$REPORT_DIR/checks.tsv"; then
        details "Inputs where a tool ${WORD[missed]}"
        {
            printf 'Tool\tInput\tExit\tPrinted\n'
            printf 'l\tl\tr\tl\n'
            # A verdict feature prints its verdict, the first fact; a row seen in
            # two families shows once.
            awk -F'\t' -v repo="$REPO/" -v verdict="${EXPECT:+1}" 'FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
                $3 == "fail" { sub(repo, "", $6); v = $5; if (verdict) sub(/;.*/, "", v); printf "%s\t%s\t%s\t%s\n", label[$1], $6, $4, (v == "" ? "nothing" : v) }' "$REPORT_DIR/style.tsv" "$REPORT_DIR/checks.tsv" | sort -u
        } | mdtable
        printf '\n</details>\n\n'
    fi
}
readme_speed() { # family [title]
    local fam="$1"
    printf '## %s\n' "${2:-Speed and memory}"
    local cpu="Every tool runs on one CPU"
    ((${#LANE_CORES[@]} == 0)) || cpu="Every tool runs on one CPU, a threaded variant on as many physical cores as its threads, in the same batch"
    printf '> <sub>The real runs and synthetic shapes at the %s tier, about %s of FASTQ each, pooled into kinds of input. %s, with the file in the page cache.</sub>\n\n' "$HEADLINE" "${TIER_SIZE[$HEADLINE]}" "$cpu"
    picture "$fam-overview" "Time and peak RSS relative to z-fastq"
    printf '\n'
    subject_sentence "$fam"
    printf '\n'
    summary_table "$fam"
    printf '\n'
    picture "$fam-speed" "Wall time by kind of input" && printf '\n'
    picture "$fam-memory" "Peak memory by kind of input" && printf '\n'
    if [[ -e "$RESULTS/figures/$fam-faults-light.png" ]]; then
        details "Minor page faults by kind of input"
        picture "$fam-faults" "Minor page faults by kind of input"
        printf '\n</details>\n\n'
    fi
}
readme_scaling() { # family [heading level]
    local fam="$1" slow
    slow="$(headline_only)"
    if [[ "$SCALING" == *" "* ]]; then
        printf '%s Scaling\n> <sub>Each synthetic shape from %s of FASTQ: throughput and peak memory against input size. A shape changes one thing %s pays for, so its line shows what that one thing costs each tool.</sub>\n\n' "${2:-##}" "$(scaling_range)" "${WORD[payer]}"
    else
        printf '%s Scaling\n> <sub>The `%s` shape (%s) from %s of FASTQ: throughput and peak memory against input size.</sub>\n\n' "${2:-##}" "$SCALING" "${DESCRIBE[$SCALING]}" "$(scaling_range)"
    fi
    picture "$fam-scaling" "Throughput and peak memory by input size" && printf '\n'
    picture "$fam-scaling-speed" "Throughput by input size" && printf '\n'
    picture "$fam-scaling-memory" "Peak memory by input size" && printf '\n'
    scaling_sentences "$fam"
    printf '\n'
    if [[ -n "$slow" || -n "${REMARK[scaling]:-}" ]]; then
        printf '> [!NOTE]\n> %s\n\n' "${slow:+$slow is timed at the $HEADLINE tier only: it is slow enough that sweeping sizes with it would take most of the run. }${REMARK[scaling]:-}"
    fi
}
readme_forms() {
    local fam
    printf '## %s\n> <sub>%s</sub>\n\n' "${WORD[forms]}" "${WORD[forms_sub]:-The same ${WORD[same]}, with the input on stdin or as many files in one call.}"
    for fam in "${FAMILIES[@]:1}"; do
        printf '### %s\n\n%s.\n\n' "$(tr '-' ' ' <<<"${fam^}")" "${FAM_DESC[$fam]^}"
        if [[ "${FAM_KIND[$fam]}" == many ]]; then
            many_sentence "$fam"
            printf '\n'
        elif [[ "${FAM_STDIN[$fam]}" == 1 ]]; then
            printf "Each tool's time and peak memory reading the file from a pipe, as a multiple of its own with the path (geometric mean over the shared shapes):\n\n"
            stdin_table "$fam"
            printf '\n'
        else
            subject_sentence "$fam"
            printf '\n'
            summary_table "$fam"
            printf '\n'
        fi
    done
}
# The families with their own memory cap, as ", 8G for deinterleave".
mem_caps() { local fam out=""; for fam in "${FAMILIES[@]}"; do [[ -z "${FAM_MEM[$fam]:-}" ]] || out+=", ${FAM_MEM[$fam]} for $fam"; done; printf '%s' "$out"; }
# The memory cap the current batches ran under: the one bounded execution
# recorded, else the cap this run would use.
batch_memory() {
    local m
    m="$(awk -F'\t' '{ print $3 }' "$RESULTS/current.tsv" | sort -u | while IFS= read -r d; do sed -n 's/^execution .* memory=\([^ ]*\).*/\1/p' "$d/meta.tsv"; done | sort -u | paste -sd/ -)"
    printf '%s\n' "${m:-$BENCH_MEM}"
}
# The day the last release run that ended well finished.
last_release() { awk -F'\t' '$3 == "release" && $9 == "ok" { d = substr($2, 1, 10) } END { print (d == "" ? "an unrecorded day" : d) }' "$RESULTS/runs.tsv"; }
readme_repro() {
    local name fam date rev
    date="$(date -u +%Y-%m-%d)"
    rev="$(zfastq_label)"
    printf '## Reproducibility\n\n'
    details "How it was measured"
    cat <<EOF
Every tool runs through \`zexec -r\` inside one zebrac batch per case, held to one CPU core by a systemd scope with a CPU quota and a $(batch_memory) memory cap$(mem_caps)$( ((${#LANE_CORES[@]} == 0)) || printf '; a batch with threaded variants holds as many physical cores as its largest lane, one logical CPU each, with each lane pinned to its own cores and z-fastq to one'), after the input was read once so every run reads from the page cache. zebrac gives wall time; zexec records each run's peak memory, CPU time, and minor faults for the tool alone. Peak RSS is the kernel's peak counter that \`wait4\` reports, which can lag the true peak by up to about 0.5 MiB, so memory differences smaller than that are inconclusive. Runs per tier, each after one warmup: tiny ${RUNS[tiny]}, xsmall ${RUNS[xsmall]}, small ${RUNS[small]}, medium ${RUNS[medium]}; sizes 1 MiB to 1 GiB. Timed by \`bench/$FEATURE/run.sh release\`, last on $(last_release); this page written on $date.

</details>

EOF
    details "Identity"
    printf -- '- `z-fastq`: *%s*\n' "$rev"
    printf -- '- `zebrac`: *%s*\n' "$(zebrac --version | awk '{ print $2 }')"
    printf -- '- `zexec`: *%s bytes, static, no libc*\n' "$(stat -c %s "$ZEXEC")"
    printf -- '- `kernel`: *%s*\n' "$(uname -r)"
    printf -- '- `cpu`: *%s*\n' "$(awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo)"
    printf -- '- `load at batch start`: *%s to %s*\n' "$(awk -F'\t' 'NR > 1 { print $21 }' "$REPORT_DIR/summary.tsv" | sort -n | head -1)" "$(awk -F'\t' 'NR > 1 { print $21 }' "$REPORT_DIR/summary.tsv" | sort -n | tail -1)"
    printf '\n</details>\n\n'
    details "Commands"
    local any=""
    for fam in "${FAMILIES[@]}"; do any+="${FAM_PARAMS[$fam]}${FAM_VARY[$fam]}"; done
    if [[ -z "$any" ]]; then
        for name in "${LANES[@]}"; do
            [[ ! -v "LANE_OF[$name]" ]] || continue
            printf -- '- `%s`: *%s*%s\n' "${REG_LABEL[$name]}" "$(lane_cmd "$name" IN | sed "s#$REPO/##g")" "$(thread_options "$name")"
        done
        for fam in "${FAMILIES[@]}"; do
            for name in ${FAM_LANES[$fam]}; do
                [[ -v "FAM_TMPL[$fam|$name]" && ! -v "LANE_OF[$name]" ]] || continue
                printf -- '- `%s`, %s: *%s*%s\n' "${REG_LABEL[$name]}" "$(family_name "$fam")" "$(lane_cmd "$name" IN "$fam" | sed "s#$REPO/##g")" "$(thread_options "$name" "$fam")"
            done
        done
    else
        # Every family's command, with the parameters of its first case; a
        # named output is the null device when timed.
        local params
        for fam in "${FAMILIES[@]}"; do
            params="$(cases | awk -F'\t' -v f="$fam" '$1 == f && !p { p = $10 } END { print p }')"
            for name in "${LANES[@]}"; do
                [[ " ${FAM_LANES[$fam]} " == *" $name "* && ! -v "LANE_OF[$name]" ]] || continue
                printf -- '- `%s`, %s: *%s*%s\n' "${REG_LABEL[$name]}" "$(family_name "$fam")" "$(lane_cmd "$name" IN "$fam" "$params" | sed "s#$REPO/##g")" "$(thread_options "$name" "$fam")"
            done
        done
    fi
    printf '\n</details>\n\n'
    if [[ -f "$RESULTS/picks.tsv" ]]; then
        details "Quick-check cases for peer and ab"
        if ((RIVAL_ANY)); then
            printf 'Picked from these results: per family, group, and format, the case where z-fastq is closest to, or furthest behind, its fastest rival of any class. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.\n\n'
        else
            printf 'Picked from these results: per group and format, the case where z-fastq leads its fastest full-check rival by the least. `run.sh peer` times z-fastq against a cached reference of the rival on each; `run.sh ab REV` times two z-fastq builds on each.\n\n'
        fi
        awk -F'\t' -v many="$(($(wc -w <<<"$(full_families)") > 1))" 'NR == FNR { title[$1] = $2; next }
            FILENAME ~ /style.tsv$/ { if (FNR > 1) label[$1] = $2; next }
            FILENAME ~ /cases.tsv$/ { if (FNR > 1) name[$1 "\t" $2] = $8; next }
            FNR > 1 { k = $1 "\t" $2
                printf "- `%s%s, %s`: *%s*", (many ? title[$1] ": " : ""), $3, $4, (k in name ? name[k] : $2)
                if ($5 != "-") printf "; rival %s at %.2fx z-fastq%ss time", label[$5], $7, "\047"
                if ($6 != "-") printf "; memory rival %s", label[$6]
                if ($8 == "declared") printf "; declared"
                printf "\n" }' <(for fam in "${FAMILIES[@]}"; do printf '%s\t%s\n' "$fam" "$(family_name "$fam")"; done) "$REPORT_DIR/style.tsv" "$REPORT_DIR/cases.tsv" "$RESULTS/picks.tsv"
        printf '\n</details>\n\n'
    fi
    local mains
    read -r -a mains <<<"$(full_families)"
    details "Hardware counters (user space, plain $HEADLINE-tier inputs, geometric means, tools at one thread)"
    counters_table "${mains[@]}"
    printf '\n</details>\n'
}
readme() {
    if declare -F feature_readme >/dev/null; then feature_readme; return; fi
    readme_head
    printf '## Findings\n\n%s\n' "$(findings "${FAMILIES[0]}")"
    valid_wrong_bullet
    printf '\n'
    readme_tools
    readme_right
    if declare -F feature_report >/dev/null; then feature_report; fi
    readme_speed "${FAMILIES[0]}"
    readme_scaling "${FAMILIES[0]}"
    readme_forms
    readme_repro
}
report() {
    STAGE=report
    REPORT_DIR="$RESULTS/report"
    local figures="$RESULTS/figures" stamp key f
    mkdir -p -- "$REPORT_DIR" "$figures"
    [[ -d "$RESULTS/batches" ]] || die "no timing run yet: run.sh release"
    collect >"$REPORT_DIR/summary.tsv.part" || die "no verified current timing batches; run.sh release is required before reporting"
    if [[ -n "$EDGE_SET" ]]; then
        [[ -f "$RESULTS/edges.tsv" && -f "$RESULTS/edges.key" && "$(<"$RESULTS/edges.key")" == "$(edge_key)" ]] ||
            die "edge checks are missing or stale; run.sh check is required before reporting"
    fi
    mv -- "$REPORT_DIR/summary.tsv.part" "$REPORT_DIR/summary.tsv"
    style_tsv >"$REPORT_DIR/style.tsv"
    cases_tsv >"$REPORT_DIR/cases.tsv"
    edges_tsv >"$REPORT_DIR/edges.tsv"
    current_checks >"$REPORT_DIR/checks.tsv"
    if [[ -f "$RESULTS/edges.tsv" ]]; then lets_through >"$REPORT_DIR/through.tsv"; else : >"$REPORT_DIR/through.tsv"; fi
    printf 'wrong\t%s\ndone\t%s\nsize\t%s\n' "${WORD[legend]}" "${WORD[past]}" "${TIER_SIZE[$HEADLINE]}" >"$REPORT_DIR/words.tsv"
    # Figures depend only on these files, plot.py, and the arguments; redraw
    # only when one of them changed or a figure is missing.
    local mains lines="" fam
    mains="$(full_families)" mains="${mains% }" mains="${mains// /,}"
    for fam in "${FAMILIES[@]}"; do [[ "${FAM_LOOK[$fam]}" != lines ]] || lines+="$fam:${FAM_VARY[$fam]%% *},"; done
    key="$(cat "$REPORT_DIR"/{summary,style,cases,edges,words}.tsv "$BENCH/lib/plot.py" | cksum) $HEADLINE $SCALING $mains${lines:+ ${lines%,}}${SECTIONS:+ sections $SECTIONS}"
    if [[ "$(cat "$REPORT_DIR/figures.key" 2>/dev/null)" != "$key" ]] || ! ls "$figures"/*.png >/dev/null 2>&1; then
        stamp="$REPORT_DIR/stamp"
        : >"$stamp"
        "$REPO/tools/venv/bin/python" -I "$BENCH/lib/plot.py" "$REPORT_DIR" "$figures" "$HEADLINE" "$SCALING" "$mains" "${lines%,}" "${SECTIONS:-}"
        # A figure the plot step did not just write belongs to nothing declared now.
        for f in "$figures"/*.svg "$figures"/*.png; do [[ "$f" -nt "$stamp" || ! -e "$f" ]] || rm -- "$f"; done
        printf '%s\n' "$key" >"$REPORT_DIR/figures.key"
    fi
    write_picks
    readme >"$FEATURE_DIR/README.md.part"
    mv -- "$FEATURE_DIR/README.md.part" "$FEATURE_DIR/README.md"
    log report DONE "${FEATURE_DIR#"$REPO"/}/README.md, $(ls "$figures"/*.png | wc -l) figures"
}


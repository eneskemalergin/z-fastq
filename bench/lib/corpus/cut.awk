# Cut nested prefix tiers from one FASTQ source in a single pass.
# Environment: R1 and optional R2 are decompression commands; TIERS is
# "name=cap ..." in ascending cap order; OUT is the output path prefix.
# A tier holds whole records (pairs count both mates) and stops before the
# first record that would exceed its plain-byte cap. Every output file is
# written as NAME.part; one line per tier is printed:
#   tier records bytes min_len max_len
function fail(msg) {
    printf "error: %s at record %d\n", msg, count + 1 > "/dev/stderr"
    failed = 1
    exit 1
}
function next_record(cmd, r,    i) {
    for (i = 1; i <= 4; i++)
        if ((cmd | getline r[i]) <= 0) {
            if (i == 1) return 0
            fail("truncated record")
        }
    if (substr(r[1], 1, 1) != "@" || substr(r[3], 1, 1) != "+" || length(r[2]) != length(r[4]))
        fail("malformed record")
    return 1
}
function size(r) { return length(r[1]) + length(r[2]) + length(r[3]) + length(r[4]) + 4 }
function name(h) { sub(/[ \t].*/, "", h); sub(/\/[0-9]+$/, "", h); return h }
function put(r, file) { printf "%s\n%s\n%s\n%s\n", r[1], r[2], r[3], r[4] > file }
function seen(len) {
    if (low == "" || len < low) low = len
    if (len > high) high = len
}
function finish(t,    s) {
    if (count == 0) fail("first record exceeds tier " tier[t])
    printf "%s %d %d %d %d\n", tier[t], count, total, low, high
    for (s = 1; s <= nsuffix; s++) close(OUT tier[t] suffix[s] ".part")
}
BEGIN {
    r1 = ENVIRON["R1"]; r2 = ENVIRON["R2"]; OUT = ENVIRON["OUT"]
    ntier = split(ENVIRON["TIERS"], spec, " ")
    for (t = 1; t <= ntier; t++) {
        split(spec[t], kv, "=")
        tier[t] = kv[1]; cap[t] = kv[2] + 0
    }
    if (r2 == "") { nsuffix = 1; suffix[1] = ".fastq" }
    else { nsuffix = 3; suffix[1] = "_R1.fastq"; suffix[2] = "_R2.fastq"; suffix[3] = "_interleaved.fastq" }
    open = 1
    while (open <= ntier && next_record(r1, a)) {
        bytes = size(a)
        if (r2 != "") {
            if (!next_record(r2, b)) fail("R2 has fewer records than R1")
            if (name(a[1]) != name(b[1])) fail("mate names differ")
            bytes += size(b)
        }
        while (open <= ntier && total + bytes > cap[open]) finish(open++)
        for (t = open; t <= ntier; t++) {
            if (r2 == "") put(a, OUT tier[t] ".fastq.part")
            else {
                put(a, OUT tier[t] "_R1.fastq.part"); put(b, OUT tier[t] "_R2.fastq.part")
                put(a, OUT tier[t] "_interleaved.fastq.part"); put(b, OUT tier[t] "_interleaved.fastq.part")
            }
        }
        seen(length(a[2]))
        if (r2 != "") seen(length(b[2]))
        total += bytes
        count++
    }
    if (open <= ntier) {
        if (r2 != "" && next_record(r2, b)) fail("R2 has more records than R1")
        fail("source ends before tier " tier[open] " is full")
    }
    exit 0
}
END { if (failed) exit 1 }

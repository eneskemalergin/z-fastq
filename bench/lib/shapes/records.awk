# Recipe values arrive through ENVIRON, without awk -v escape processing.
# Prefix arithmetic uses the same header and plus formats as emission.
# With a pool, sequence and quality bytes are slices of real reads; line
# lengths, and therefore all byte arithmetic, are the same as without one.
function digits(n, width,    total, p) {
    total = width * n
    for (p = 10 ^ width; p <= n; p *= 10) total += n - p + 1
    return total
}
function cycle(n, values, period,    j, total) {
    for (j = 1; j <= period; j++)
        total += (int(n / period) + ((j - 1 - (start - 1) % period + period) % period < n % period)) * values[j]
    return total
}
function bytes(n,    total, i, c, last) {
    if (n <= 0) return 0
    last = start + n - 1
    total = constant * n + cycle(n, costs, period)
    if (engine == "longread") {
        total += cycle(n, channel_widths, channels) + 3 * (int(last / unclassified) - int((start - 1) / unclassified))
        for (i = (int((start - 1) / every) + 1) * every; i <= last; i += every) {
            c = (i - 1) % period + 1
            total += 2 * (outlier_base + (i * multiplier) % modulus - lens[c])
        }
    } else {
        total += copies * (digits(last, width) - digits(start - 1, width))
        if (engine == "variable") total += digits(last, 1) - digits(start - 1, 1)
    }
    if (endings != "lf") total += 4 * mates * (int((last + 1) / 2) - int(start / 2))
    return total
}
function mulmod(seed, factor, modulus) {
    return (((int(seed / 65536) * factor) % modulus) * 65536 + (seed % 65536) * factor) % modulus
}
function uuid(seed,    a, b, c, low, high, product) {
    a = sprintf("%08x", mulmod(seed, 2654435761, 4294967295))
    b = sprintf("%04x", mulmod(seed, 40503, 65535))
    c = sprintf("%04x", mulmod(seed, 65599, 65535))
    product = (seed % 65536) * 1597334677
    low = product % 65536
    high = (int(product / 65536) + (int(seed / 65536) * 1597334677) % 4294967296) % 4294967296
    return a "-" b "-" c "-9f3a-" sprintf("%08x%04x", high, low)
}
# Concatenate the sequence and quality lines of a FASTQ pool into ps and pq.
function load(path,    line, n, len) {
    ps = pq = ""
    while ((getline line < path) > 0) {
        sub(/\r$/, "", line)
        if (++n % 4 == 2) {
            ps = ps line
            len = length(line)
            if (n == 2) read_length = len
            else if (len != read_length) read_length = 0
        } else if (n % 4 == 0) pq = pq line
    }
    if (n == 0 || n % 4 != 0 || length(ps) != length(pq)) {
        printf "error: unreadable pool %s\n", path > "/dev/stderr"
        exit 2
    }
    close(path)
}
# Choose a slice of length len for record i; mates share the offset so a
# pair of pools keeps real mates together. Fixed-length pools slice at read
# starts, which keeps the per-cycle quality profile.
function slice(i, len,    stride, slots, off) {
    stride = (read_length >= len) ? read_length : 1
    slots = int((pool_bytes - len) / stride) + 1
    off = mulmod(i, 2654435761, slots)
    if (seed) off = (off + mulmod(seed, 40503, slots)) % slots
    off = off * stride + 1
    s = substr(seq, off, len); q = substr(qual, off, len)
    if (pool2 != "") { s2 = substr(seq2, off, len); q2 = substr(qual2, off, len) }
}
BEGIN {
    engine = ENVIRON["S_ENGINE"]; endings = ENVIRON["S_ENDINGS"]
    period = split(ENVIRON["S_LENGTHS"], lens, " ")
    width = ENVIRON["S_WIDTH"] + 0
    start = ENVIRON["S_START"] + 0; n = ENVIRON["S_COUNT"] + 0
    cap = ENVIRON["S_CAP"] + 0
    mates = (engine == "paired") ? 2 : 1; copies = mates
    hfmt = "@" ENVIRON["S_PREFIX"] "%0" width "d"; pfmt = "+"
    if (engine == "variable") {
        hfmt = hfmt " sample=%d"
        pfmt = "+" ENVIRON["S_PREFIX"] "%0" width "d len=%d note"
        constant = length(sprintf(hfmt, 1, 1)) + length(sprintf(pfmt, 1, 1)) + 4 - 2 * width - 2
        copies = 2
    } else if (engine == "longread") {
        offset = ENVIRON["S_OFFSET"] + 0; channels = ENVIRON["S_CHANNELS"] + 0
        barcodes = ENVIRON["S_BARCODES"] + 0; unclassified = ENVIRON["S_UNCLASSIFIED"] + 0
        every = ENVIRON["S_EVERY"] + 0; outlier_base = ENVIRON["S_BASE"] + 0
        multiplier = ENVIRON["S_MULTIPLIER"] + 0; modulus = ENVIRON["S_MODULUS"] + 0
        hfmt = "@%s runid=%s ch=%d start_time=" ENVIRON["S_TIME"] " barcode=%s"
        constant = length(sprintf(hfmt, uuid(1), uuid(1 + offset), 1, "barcode01")) + 5 - 1
        for (i = 1; i <= channels; i++) channel_widths[i] = length(sprintf("%d", i % channels + 1))
    } else {
        if (engine == "paired") hfmt = hfmt "/%d %d:N:0:" ENVIRON["S_BARCODE"]
        constant = mates * (length(sprintf(hfmt, 1, 1, 1)) + 5 - width)
    }
    for (c = 1; c <= period; c++) {
        costs[c] = 2 * mates * lens[c]
        if (engine == "variable") costs[c] += length(sprintf("%d", lens[c]))
        if (lens[c] > longest) longest = lens[c]
    }
    if (ENVIRON["S_MODE"] == "plan") {
        if (cap > 0) {
            n = int(cap / ENVIRON["S_NOMINAL"])
            if (n < 1) n = 1
            if (n > 999999999) n = 999999999
            if (bytes(n) > cap) {
                lo = 0; hi = n
                while (lo < hi) {
                    mid = int((lo + hi + 1) / 2)
                    if (bytes(mid) <= cap) lo = mid; else hi = mid - 1
                }
                n = lo
            }
        }
        printf "%.0f %.0f\n", n, bytes(n)
        exit
    }
    if (engine == "longread" && outlier_base + modulus - 1 > longest) longest = outlier_base + modulus - 1
    pool = ENVIRON["S_POOL"]; pool2 = ENVIRON["S_POOL2"]; seed = ENVIRON["S_SEED"] + 0
    if (pool != "") {
        if (pool2 != "") {
            load(pool2); seq2 = ps; qual2 = pq; fixed2 = read_length
        }
        load(pool); seq = ps; qual = pq; ps = pq = ""
        pool_bytes = length(seq)
        if (pool2 != "" && (length(seq2) != pool_bytes || fixed2 != read_length)) {
            printf "error: mate pools differ in geometry\n" > "/dev/stderr"
            exit 2
        }
        if (pool_bytes < longest) {
            printf "error: pool holds %d bases; the recipe needs %d\n", pool_bytes, longest > "/dev/stderr"
            exit 2
        }
    } else {
        seq = ENVIRON["S_SEQUENCE"]; qual = ENVIRON["S_QUALITY"]
        while (length(seq) < longest) seq = seq seq
        while (length(qual) < longest) qual = qual qual
    }
    for (c = 1; c <= period; c++) {
        S[c] = substr(seq, 1, lens[c]); Q[c] = substr(qual, 1, lens[c])
    }
    fmt = hfmt "\n%s\n" pfmt "\n%s\n"
    crfmt = fmt; gsub(/\n/, "\r\n", crfmt)
    out = ENVIRON["S_OUT"]; layout = ENVIRON["S_LAYOUT"]
    s = s2 = S[1]; q = q2 = Q[1]
    if (engine == "paired") {
        r1 = out "_R1.fastq"; r2 = out "_R2.fastq"; inter = out "_interleaved.fastq"
        for (i = start; i < start + n; i++) {
            f = (endings != "lf" && i % 2) ? crfmt : fmt
            if (pool != "") slice(i, lens[1])
            if (layout != "interleaved") {
                printf f, i, 1, 1, s, q > r1
                printf f, i, 2, 2, s2, q2 > r2
            }
            if (layout != "mates") {
                printf f, i, 1, 1, s, q > inter
                printf f, i, 2, 2, s2, q2 > inter
            }
        }
        if (layout != "interleaved" && (close(r1) || close(r2))) exit 1
        if (layout != "mates" && close(inter)) exit 1
    } else if (engine == "fixed") {
        if (pool != "") {
            for (i = start; i < start + n; i++) {
                slice(i, lens[1])
                printf (endings != "lf" && i % 2 ? crfmt : fmt), i, s, q
            }
        } else if (endings == "lf") {
            for (i = start; i < start + n; i++) printf fmt, i, s, q
        } else {
            for (i = start; i < start + n; i++) printf (i % 2 ? crfmt : fmt), i, s, q
        }
    } else if (engine == "variable") {
        for (i = start; i < start + n; i++) {
            f = (endings != "lf" && i % 2) ? crfmt : fmt
            c = (i - 1) % period + 1
            if (pool != "") {
                slice(i, lens[c])
                printf f, i, i, s, i, lens[c], q
            } else printf f, i, i, S[c], i, lens[c], Q[c]
        }
    } else {
        for (i = start; i < start + n; i++) {
            f = (endings != "lf" && i % 2) ? crfmt : fmt
            c = (i - 1) % period + 1
            len = (i % every == 0) ? outlier_base + (i * multiplier) % modulus : lens[c]
            if (pool != "") slice(i, len)
            else if (i % every == 0) { s = substr(seq, 1, len); q = substr(qual, 1, len) }
            else { s = S[c]; q = Q[c] }
            barcode = (i % unclassified == 0) ? "unclassified" : sprintf("barcode%02d", i % barcodes + 1)
            printf f, uuid(i), uuid(i + offset), i % channels + 1, barcode, s, q
        }
    }
}

<!-- markdownlint-disable MD033 MD041 -->

# Benchmarks

Publication benches for z-fastq. Suites live in `bench/<name>/`. Each suite checks that timed commands agree on the compared fact, then [zebrac](https://github.com/eneskemalergin/zebrac) records wall, peak RSS, and minor page faults.

Zebrac is Linux-only. Peers, adapters, zebrac, and report Python are under [`tools/`](../tools/README.md). **This repo tests Linux x86-64 only.** Other OS cells below are from each tool's own docs or release assets, not from us.

## Suites

| Suite | Report | Compared fact |
| ----- | ------ | ------------- |
| [`count`](count/) | [REPORT.md](count/REPORT.md) | record count integer |

`count` builds both z-fastq binaries: native inflate (`-Disa-l=false`) as `zig-out/bin/z-fastq-native`, then the ISA-L product at `zig-out/bin/z-fastq`.

## Run

```bash
tools/install.sh all
tools/install.sh --check all
bash bench/shared/download_data.sh
bash bench/count/run.sh
```

Publication sampling (also the `run.sh` default): **5000 ms, 25 runs, 5 warmups**. Zebrac stops when both duration and min-samples are met.

Bring-up: `bash bench/count/run.sh --small-real --skip-scale --runs 2 --warmup 1 --duration 800`. MiniSeq is not the headline Dense file (SRR1810900).

Timed argv is the tool. Do not wrap in `bash -c` or a pipeline; that measures the shell. Plain vs gzip are separate sections. Gzip throughput is decoded FASTQ MiB/s.

## Binaries

Stripped Linux x86-64 sizes on 2026-09-15. Adapters are our wrappers, not upstream CLIs. BBTools is the retained class tree, not one ELF.

| Binary | Version | Bytes | This host |
| ------ | ------- | ----: | --------- |
| z-fastq (native) | 0.0.18 | 535,088 | static |
| z-fastq (ISA-L) | 0.0.18 | 754,560 | static |
| seqtk | 1.5-r133 | 77,600 | `libz`, `libm`, `libc` |
| FastQValidator | 0.1.1a | 169,120 | `libz`, `libstdc++`, `libgcc_s`, `libm`, `libc` |
| Needletail adapter | 0.7.3 | 425,896 | `libgcc_s`, `libc` |
| Helicase adapter | 0.2.0 | 428,944 | `libgcc_s`, `libc` |
| Fasten `fasten_sample` | 0.9.0 | 470,336 | `libgcc_s`, `libc` |
| fqtools | 2.3 + HTSlib 1.24 | 905,088 | `libz`, `libm`, `libc` (HTSlib in the binary) |
| SeqFu | 1.27.1 | 1,673,280 | `libz`, `libm`, `libc` |
| IRMA Core | 0.10.1 | 2,485,768 | `libgcc_s`, `libpthread`, `libm`, `libdl`, `libc` |
| fq | 0.12.0 | 2,740,832 | `libgcc_s`, `libpthread`, `libm`, `libc` |
| fastp | 1.3.7 | 3,007,424 | static |
| Rasusa | 5.1.0 | 5,640,352 | `libgcc_s`, `librt`, `libpthread`, `libdl`, `libc` |
| fqkit | 0.4.14 | 6,007,128 | `liblzma`, `libz`, `libbz2`, `libgcc_s`, plus fontconfig/freetype/harfbuzz/png/xml/glib |
| SeqKit | 2.13.0 | 20,076,696 | static |
| BBTools | 40.02 | 19,359,140 | Java class tree; needs a JVM |

seqtk is the smallest file (`libz` at runtime). z-fastq is larger because it is static. Native is 219,472 bytes smaller than ISA-L.

## Targets

`yes` = author ships or documents that OS/arch. `pkg` = Bioconda/conda only. `cargo` = `cargo install` / compile. `—` = not stated. `no` = author or this project says no. We have only run the Linux x86-64 column.

| Binary | linux amd64 | linux arm64 | macOS amd64 | macOS arm64 | Windows |
| ------ | ----------- | ----------- | ----------- | ----------- | ------- |
| z-fastq (ISA-L) | yes | no | no | no | no |
| z-fastq (native) | yes | no | no | no | no |
| seqtk | yes | pkg | — | pkg | — |
| FastQValidator | yes | — | — | — | — |
| Needletail | yes | — | CI wheels | CI wheels | — |
| Helicase | yes (AVX2/SSE3) | NEON | — | NEON | — |
| Fasten | yes | cargo | cargo | cargo | cargo |
| fqtools | yes | pkg | pkg | pkg | — |
| SeqFu | yes | — | yes | Darwin zip | no |
| IRMA Core | yes | — | — | — | — |
| fq | yes | — | yes | — | yes |
| fastp | yes | — | compile | compile | — |
| Rasusa | yes | yes | yes | yes | — |
| fqkit | cargo | cargo | cargo | cargo | cargo |
| SeqKit | yes | yes | yes | yes | yes |
| BBTools | JVM | JVM | JVM | JVM | JVM |

Needletail and Helicase are libraries; we time Linux x86-64 adapters. Helicase wants AVX2, SSE3, or NEON. fastp's Linux binary is the documented prebuilt; macOS is compile/conda. Intel ISA-L also documents Windows; that is ISA-L, not z-fastq.

## Dependencies

Runtime shared libraries on this install, then what the build actually pulls in.

| Binary | Runtime libs | Stack |
| ------ | ------------ | ----- |
| z-fastq (native) | none | Zig std inflate |
| z-fastq (ISA-L) | none | Zig std + vendored ISA-L inflate/CRC (NASM at build) |
| seqtk | `libz` | one C file |
| FastQValidator | `libz`, libstdc++ | libStatGen |
| Needletail adapter | libc, libgcc | 11 crates (`flate2`) |
| Helicase adapter | libc, libgcc | 12 crates (`deko`, `memmap2`, `flate2`) |
| Fasten | libc, libgcc | Rust crate graph |
| fqtools | `libz` | HTSlib 1.24 (bz2/lzma/curl off in our recipe) |
| SeqFu | `libz` | Nim 2.2.12 + 16 Nimble packages |
| IRMA Core | libc, libgcc, libpthread, libdl | upstream Rust release |
| fq | libc, libgcc, libpthread | upstream Rust release |
| fastp | none (static) | ISA-L + libdeflate + Highway when built from source |
| Rasusa | libc, libgcc, libpthread, libdl | upstream Rust release |
| fqkit | zlib, lzma, bz2, font stack | Rust + plotters/GUI libs |
| SeqKit | none (static) | Go, many codecs in-binary |
| BBTools | JVM | Java runtime + class tree |

Native z-fastq: no loader, no third-party lock. ISA-L stays static. seqtk is the shallowest dynamic C peer. fqtools is HTSlib. fqkit and SeqFu are the heavy package graphs. BBTools is a JVM.

## Files

`bench/shared/` has `tools.sh`, `datasets.manifest`, `download_data.sh`, `generate_scaling.py`. Cache and downloads are gitignored. Commit each suite's `REPORT.md` plus `results/figures/*.png`.

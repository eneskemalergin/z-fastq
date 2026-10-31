# For contributors

This page maps the code, the tests, and the wiki onto each other, so a change to one can find the others.

## Source map

| File                 | Holds                                                        |
| -------------------- | ------------------------------------------------------------ |
| `src/main.zig`       | the command line: options, commands, output, and messages    |
| `src/root.zig`       | the library's exports and version                            |
| `src/fastq.zig`      | records, the reader and writer, and validation               |
| `src/count_scan.zig` | the record counter                                           |
| `src/stats.zig`      | statistics                                                   |
| `src/sample.zig`     | sampling                                                     |
| `src/pair.zig`       | pair-name matching                                           |
| `src/preprocess.zig` | trimming, filtering, inserts, and pair duplication           |
| `src/qc.zig`         | the QC modules and the QC table writer                       |
| `src/simd.zig`       | vector kernels and their run-time selection                  |
| `src/io.zig`         | byte sources and sinks, limits, and gzip                     |
| `tests/`             | command-line and library tests, one file per command or part |

## Before changing documented behavior

1. Find the code and the tests behind the behavior.
2. Change the code and the tests together.
3. Update the one wiki page that owns the behavior, and link to it rather than repeat it elsewhere. The [Home](Home) page groups the pages by what they own.

## Tests

Run these from the repository root, one after another, with Zig 0.16.0 on Linux x86-64; on aarch64, the portable run repeats the first, since there is nothing to dispatch:

```bash
zig fmt --check src tests build.zig
zig build test --summary all
zig build test -Dkernel-backend=portable --summary all
zig build test -Dstatic=true -Doptimize=ReleaseSafe --summary all
zig build test -Dstatic=true -Doptimize=ReleaseFast --summary all
```

`zig build test-<suite>` runs one test binary while you work: `test-src` for the tests inside `src/`, or a file in `tests/` without its `test_` prefix, such as `test-reader`. `zig build -l` lists them.

- The runs share `zig-out/bin/z-fastq`, so run them in sequence; parallel runs need separate checkouts.
- The tests use small tracked files only, no downloads or other tools. They need `pidfd_open`, `/proc`, FIFOs, and hard and symbolic links.
- No test is expected to skip; investigate any skip.
- Debug builds target the host CPU; the release runs target generic x86-64 and are the ones that test the run-time AVX2 kernels against their scalar versions. A timing is never proof of correct output.

## Continuous integration

GitHub Actions checks every pull request and every push to `main`, and any branch when started by hand, on x86-64 and on an ARM runner. A change that touches only documentation, benchmark reports, or tool recipes skips the builds. The first failure stops the run. The runs above stay the local check before a commit; CI repeats only what a separate machine adds.

- **Source checks:** formatting, the syntax of every shell script, shellcheck on the CI scripts, the benchmark Python, and `.github/scripts/docs_check.sh`, which fails on a link, image, or heading anchor that does not exist in the README, the reports, or the wiki.
- **Tests:** one Debug run per architecture; on x86-64 it targets generic x86-64, so the AVX2 kernels run through run-time dispatch as in the release binary.
- **Release binaries:** each architecture's static ReleaseFast binary, packaged and smoke-tested on every command from the unpacked archive, with Zig off `PATH`.
- **Source package:** the files under `.paths` in `build.zig.zon`, fetched into a fresh Zig package and tested through the public module.
- **Workflows:** actionlint and zizmor, when `.github/` changes.

A `vX.Y.Z` tag runs the release. It goes ahead only for a tag on `main`'s head naming the version in `build.zig.zon` and its changelog entry, and only once `main`'s CI has passed on that commit; it then publishes the archives that run tested, with checksums, the changelog entry, and build attestations, without building again. A tag pushed before that CI run finishes fails, and the release is re-run after it passes. Started by hand on `main`, the release workflow rehearses the checks for `main`'s version without publishing. Pushes to `main` that change `wiki/` publish it to the GitHub Wiki.

## Benchmarks

`bench/` compares z-fastq with other tools on the data listed in [References](References). Each suite has a README with its method, findings, and figures; `tools/` installs the pinned tools they need.

## Wiki figures

The schematics in `wiki/assets` come from one script:

```bash
tools/venv/bin/python tools/wiki_figures.py
```

It writes a light and a dark PNG of each figure, plus an SVG for editing; the PNGs are tracked and the SVGs are ignored by git. Benchmark figures are not copied: pages link to them on `main`, so a committed bench rerun updates the wiki by itself.

Every figure on a page follows the same pattern: a centred `<picture>` with a dark source, a light source, and a light fallback image with alt text, then a small numbered legend:

```html
<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="assets/NAME-dark.png"
    />
    <source
      media="(prefers-color-scheme: light)"
      srcset="assets/NAME-light.png"
    />
    <img alt="What the figure shows." src="assets/NAME-light.png" width="680" />
  </picture>
</p>
<p align="center">
  <sub><b>Figure 1.</b> One or two sentences on what to see.</sub>
</p>
```

## Related pages

- [References](References): the tools and data the benchmarks use.
- [Zig library](Zig-Library): the public module.

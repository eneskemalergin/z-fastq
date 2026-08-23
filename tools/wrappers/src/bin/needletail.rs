//! Needletail adapter for valid-input engine measurements.

#![forbid(unsafe_code)]

use std::env;
use std::ffi::OsString;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;

use needletail::parse_fastx_file;

#[path = "../stats.rs"]
mod stats;

const ENGINE_VERSION: &str = "0.7.3";

#[derive(Clone, Copy)]
enum Job {
    Count,
    Stats,
}

fn main() {
    let (job, path) = arguments();
    match job {
        Job::Count => run_count(&path),
        Job::Stats => run_stats(&path),
    }
}

fn arguments() -> (Job, PathBuf) {
    let mut args = env::args_os().skip(1);
    let first = args.next().unwrap_or_else(|| usage_error(None));

    if first == "--help" || first == "-h" {
        usage(&mut io::stdout().lock()).unwrap_or_else(|error| fail("stdout", error));
        process::exit(0);
    }
    if first == "--version" || first == "-V" {
        version(&mut io::stdout().lock()).unwrap_or_else(|error| fail("stdout", error));
        process::exit(0);
    }

    let job = match first.to_str() {
        Some("count") => Job::Count,
        Some("stats") => Job::Stats,
        _ => usage_error(Some(first)),
    };
    let path = args.next().unwrap_or_else(|| usage_error(None));
    if let Some(extra) = args.next() {
        usage_error(Some(extra));
    }
    (job, PathBuf::from(path))
}

fn usage_error(argument: Option<OsString>) -> ! {
    if let Some(argument) = argument {
        eprintln!("error: unexpected argument: {}", argument.to_string_lossy());
    }
    usage(&mut io::stderr().lock()).unwrap_or_else(|error| fail("stderr", error));
    process::exit(2);
}

fn usage(output: &mut impl Write) -> io::Result<()> {
    writeln!(output, "Usage: needletail-adapter <count|stats> <in.fastq>")?;
    writeln!(output, "       needletail-adapter --version")
}

fn version(output: &mut impl Write) -> io::Result<()> {
    let target_cpu = if cfg!(feature = "native-cpu") {
        "native"
    } else {
        "portable"
    };
    writeln!(
        output,
        "needletail-adapter {}; engine=needletail {}; jobs=count,stats; domain=screened-four-line-lf-fastq; compression=plain,gzip; target={}-{}; target-cpu={}; threads=1",
        env!("CARGO_PKG_VERSION"),
        ENGINE_VERSION,
        env::consts::ARCH,
        env::consts::OS,
        target_cpu,
    )
}

fn run_count(path: &Path) {
    let mut reader = parse_fastx_file(path).unwrap_or_else(|error| fail_path(path, error));
    let mut reads = 0_u64;
    while let Some(record) = reader.next() {
        record.unwrap_or_else(|error| fail_path(path, error));
        reads = reads
            .checked_add(1)
            .unwrap_or_else(|| fail_path(path, "record count overflow"));
    }

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    writeln!(output, "{reads}").unwrap_or_else(|error| fail("stdout", error));
    output.flush().unwrap_or_else(|error| fail("stdout", error));
}

fn run_stats(path: &Path) {
    let mut reader = parse_fastx_file(path).unwrap_or_else(|error| fail_path(path, error));
    let mut aggregate = stats::Stats::default();
    while let Some(record) = reader.next() {
        let record = record.unwrap_or_else(|error| fail_path(path, error));
        let quality = record
            .qual()
            .unwrap_or_else(|| fail_path(path, "record has no quality"));
        aggregate
            .add_record(&record.raw_seq(), quality)
            .unwrap_or_else(|error| fail_path(path, error));
    }

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());
    aggregate
        .write_human(&mut output, path.display())
        .unwrap_or_else(|error| fail("stdout", error));
    output.flush().unwrap_or_else(|error| fail("stdout", error));
}

fn fail_path(path: &Path, error: impl std::fmt::Display) -> ! {
    eprintln!("error: {}: {error}", path.display());
    process::exit(1);
}

fn fail(name: &str, error: impl std::fmt::Display) -> ! {
    eprintln!("error: {name}: {error}");
    process::exit(1);
}

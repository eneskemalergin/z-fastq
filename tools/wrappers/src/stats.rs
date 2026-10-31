//! Shared aggregate statistics for the parser adapters.

use std::fmt;
use std::io::{self, Write};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StatsError {
    LengthMismatch { sequence: usize, quality: usize },
    InvalidQuality(u8),
    Overflow,
}

impl fmt::Display for StatsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch { sequence, quality } => {
                write!(
                    formatter,
                    "sequence length {sequence} differs from quality length {quality}"
                )
            }
            Self::InvalidQuality(byte) => write!(formatter, "invalid Phred+33 byte {byte}"),
            Self::Overflow => formatter.write_str("statistics overflow"),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Stats {
    reads: u64,
    bases: u64,
    min_length: u64,
    max_length: u64,
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other_bases: u64,
    quality_sum: u64,
    q20_bases: u64,
    q30_bases: u64,
}

impl Stats {
    #[inline]
    pub fn add_record(&mut self, sequence: &[u8], quality: &[u8]) -> Result<(), StatsError> {
        if sequence.len() != quality.len() {
            return Err(StatsError::LengthMismatch {
                sequence: sequence.len(),
                quality: quality.len(),
            });
        }

        if let Some(&encoded_quality) = quality.iter().find(|&&q| !(33..=126).contains(&q)) {
            return Err(StatsError::InvalidQuality(encoded_quality));
        }

        // One branch-free pass per counter, so each loop vectorizes; folding
        // bit 5 makes each base test case-insensitive.
        let length = u64::try_from(sequence.len()).map_err(|_| StatsError::Overflow)?;
        let count = |base: u8| sequence.iter().filter(|&&b| (b | 0x20) == base).count() as u64;
        let a = count(b'a');
        let c = count(b'c');
        let g = count(b'g');
        let t = count(b't');
        let n = count(b'n');
        let other_bases = length - a - c - g - t - n;
        // Scores are at most 93, so this bound keeps the plain sum exact.
        if length > u64::MAX / 93 {
            return Err(StatsError::Overflow);
        }
        let quality_sum = quality.iter().map(|&q| u64::from(q - 33)).sum::<u64>();
        let q20_bases = quality.iter().filter(|&&q| q >= 33 + 20).count() as u64;
        let q30_bases = quality.iter().filter(|&&q| q >= 33 + 30).count() as u64;
        let next = Self {
            reads: checked_add(self.reads, 1)?,
            bases: checked_add(self.bases, length)?,
            min_length: if self.reads == 0 {
                length
            } else {
                self.min_length.min(length)
            },
            max_length: if self.reads == 0 {
                length
            } else {
                self.max_length.max(length)
            },
            a: checked_add(self.a, a)?,
            c: checked_add(self.c, c)?,
            g: checked_add(self.g, g)?,
            t: checked_add(self.t, t)?,
            n: checked_add(self.n, n)?,
            other_bases: checked_add(self.other_bases, other_bases)?,
            quality_sum: checked_add(self.quality_sum, quality_sum)?,
            q20_bases: checked_add(self.q20_bases, q20_bases)?,
            q30_bases: checked_add(self.q30_bases, q30_bases)?,
        };
        *self = next;
        Ok(())
    }

    pub fn write_human(&self, writer: &mut impl Write, input: impl fmt::Display) -> io::Result<()> {
        writeln!(writer, "input: {input}")?;
        writeln!(writer, "reads: {}", self.reads)?;
        writeln!(writer, "bases: {}", self.bases)?;
        write_optional_integer(writer, "min_length", self.min_length, self.reads)?;
        write_optional_integer(writer, "max_length", self.max_length, self.reads)?;
        write_ratio(writer, "mean_length", self.bases, self.reads)?;
        writeln!(writer, "a: {}", self.a)?;
        writeln!(writer, "c: {}", self.c)?;
        writeln!(writer, "g: {}", self.g)?;
        writeln!(writer, "t: {}", self.t)?;
        writeln!(writer, "n: {}", self.n)?;
        writeln!(writer, "other_bases: {}", self.other_bases)?;
        let gc_population = self.a + self.c + self.g + self.t;
        write_ratio(writer, "gc_fraction", self.g + self.c, gc_population)?;
        writeln!(writer, "quality_sum: {}", self.quality_sum)?;
        write_ratio(writer, "mean_quality", self.quality_sum, self.bases)?;
        writeln!(writer, "q20_bases: {}", self.q20_bases)?;
        write_ratio(writer, "q20_fraction", self.q20_bases, self.bases)?;
        writeln!(writer, "q30_bases: {}", self.q30_bases)?;
        write_ratio(writer, "q30_fraction", self.q30_bases, self.bases)
    }
}

#[inline]
fn checked_add(left: u64, right: u64) -> Result<u64, StatsError> {
    left.checked_add(right).ok_or(StatsError::Overflow)
}

fn write_optional_integer(
    writer: &mut impl Write,
    name: &str,
    value: u64,
    population: u64,
) -> io::Result<()> {
    if population == 0 {
        writeln!(writer, "{name}: -")
    } else {
        writeln!(writer, "{name}: {value}")
    }
}

fn write_ratio(
    writer: &mut impl Write,
    name: &str,
    numerator: u64,
    denominator: u64,
) -> io::Result<()> {
    if denominator == 0 {
        writeln!(writer, "{name}: -")
    } else {
        writeln!(
            writer,
            "{name}: {:.6}",
            numerator as f64 / denominator as f64
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_exact_human_contract() {
        let mut stats = Stats::default();
        stats.add_record(b"AaCcGgTtNnR", b"!5?~!!!!!!!").unwrap();

        let mut output = Vec::new();
        stats.write_human(&mut output, "reads.fastq").unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            concat!(
                "input: reads.fastq\n",
                "reads: 1\n",
                "bases: 11\n",
                "min_length: 11\n",
                "max_length: 11\n",
                "mean_length: 11.000000\n",
                "a: 2\n",
                "c: 2\n",
                "g: 2\n",
                "t: 2\n",
                "n: 2\n",
                "other_bases: 1\n",
                "gc_fraction: 0.500000\n",
                "quality_sum: 143\n",
                "mean_quality: 13.000000\n",
                "q20_bases: 3\n",
                "q20_fraction: 0.272727\n",
                "q30_bases: 2\n",
                "q30_fraction: 0.181818\n",
            )
        );
    }

    #[test]
    fn failed_record_does_not_change_totals() {
        let mut stats = Stats::default();
        stats.add_record(b"AC", b"!!").unwrap();
        let expected = stats;

        assert_eq!(
            stats.add_record(b"A", b""),
            Err(StatsError::LengthMismatch {
                sequence: 1,
                quality: 0,
            })
        );
        assert_eq!(stats, expected);
        assert_eq!(
            stats.add_record(b"A", b"\x7f"),
            Err(StatsError::InvalidQuality(127))
        );
        assert_eq!(stats, expected);
    }
}

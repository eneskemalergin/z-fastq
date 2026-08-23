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

        let length = u64::try_from(sequence.len()).map_err(|_| StatsError::Overflow)?;
        let mut a = 0_u64;
        let mut c = 0_u64;
        let mut g = 0_u64;
        let mut t = 0_u64;
        let mut n = 0_u64;
        let mut other_bases = 0_u64;
        let mut quality_sum = 0_u128;
        let mut q20_bases = 0_u64;
        let mut q30_bases = 0_u64;

        for (&base, &encoded_quality) in sequence.iter().zip(quality) {
            match base {
                b'A' | b'a' => a += 1,
                b'C' | b'c' => c += 1,
                b'G' | b'g' => g += 1,
                b'T' | b't' => t += 1,
                b'N' | b'n' => n += 1,
                _ => other_bases += 1,
            }

            if !(33..=126).contains(&encoded_quality) {
                return Err(StatsError::InvalidQuality(encoded_quality));
            }
            let score = encoded_quality - 33;
            quality_sum += u128::from(score);
            q20_bases += (score >= 20) as u64;
            q30_bases += (score >= 30) as u64;
        }

        let quality_sum = u64::try_from(quality_sum).map_err(|_| StatsError::Overflow)?;
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

# Paired reads

Paired-end sequencing reads each DNA fragment from both ends, giving two reads per fragment: mate 1 and mate 2, usually stored in two files, R1 and R2. Everything downstream assumes the files stay in step, record for record. z-fastq treats a pair as one unit everywhere: it never keeps, drops, or moves one mate without the other.

A pair is one record from each of two files, or two neighbouring records in one interleaved file ([Interleaving](Interleaving)).

## Check pairs

```bash
z-fastq check --paired R1.fastq.gz R2.fastq.gz
z-fastq check --interleaved interleaved.fastq.gz
```

Each mate is first validated as a record ([Validation](Validation)), then the two names are compared. In an interleaved file the first two records form the first pair, the next two the second, and so on; an odd last record has no mate.

- `P001`: the mates' names do not match.
- `P002`: one input has more records than the other.

A record error comes before a pair error at the same place. Messages give the pair number, counted from 0, and shorten very long names.

## How names are matched

Sequencers mark the mate in the read name in one of two ways: a `/1` or `/2` at the end of the name, or a second word such as `1:N:0:1` (Illumina's Casava format). The default policy, `--pair-names illumina`, understands both:

| R1 header | R2 header | Pair? |
| --- | --- | --- |
| `@read42/1` | `@read42/2` | yes |
| `@read42 1:N:0:1` | `@read42 2:N:0:1` | yes |
| `@read42` | `@read42` | yes |
| `@read42/2` | `@read42/1` | no, mates reversed |
| `@read42/1` | `@read42` | no, marker on one side only |
| `@read42/3` | `@read42/4` | no, `/3` is part of the name |

Step by step:

1. Take the first word of each header, and drop a final `/1` or `/2` from it.
2. Read a mate marker from the second word, if there is one: a leading `1:` or `2:`, or a final `/1` or `/2`.
3. The identifiers left must be equal and not empty.
4. If both headers carry markers, mate 1 must say 1 and mate 2 must say 2, and markers within one header must agree. A marker on one side only is a mismatch.

Anything after the marker word is ignored. So `@read42 1:N:0/3` and `@read42 2:N:0/4` do pair: their first words agree, and their second words say 1 and 2.

**`--pair-names exact`** is the strict alternative: the whole first word must match byte for byte, and nothing is removed.

## Related pages

- [Sampling](Sampling#pairs): subsets that keep both mates or neither.
- [Interleaving](Interleaving): between two files and one interleaved file.
- [Processing paired reads](Processing): every pair is validated before it is trimmed.
- [Troubleshooting](Troubleshooting#pair-names-do-not-match): when names do not match.

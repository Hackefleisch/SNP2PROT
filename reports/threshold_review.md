# Are the E-score cutoffs right?

Written 2026-08-18, prompted by the Kock screening (`reports/kock2024_screening.md`): if two
labs measuring the same protein agree on only 46% of their positive calls, the cutoff that
produced those calls deserves a look.

**Short answer: `0.45` is the best absolute cutoff we could pick, and the absolute-cutoff
*rule* is what costs the agreement.** A rank-matched rule at the same stringency agrees
15-20 points better. Two sources are effectively label-dead at `0.45`, and that is a separate
defect the sweep turned up.

Nothing here has been changed. `configs/thresholds.yaml` still reads `positive: 0.45`,
`negative: 0.35`, `per_experiment: true`.

## Method

The 47 domains that two sources both measured — 49 source pairs, the same set
`reports/overlap.md` uses. For each pair, labels were recomputed from the stored `raw_score`
under each candidate rule and the positive sets compared by Jaccard.

One caveat on comparability: this sweep takes the median `raw_score` per
(domain, source, 8-mer), where `overlap.md` uses the stored labels, which are binarized per
experiment and then reconciled. So the absolute numbers run a little higher here — 0.520 at
`E >= 0.45` against the report's 0.458 — and only the *differences between rules* are the
result.

## Absolute cutoffs

| rule | median Jaccard | mean | median positives | median size ratio |
|---|---:|---:|---:|---:|
| `E >= 0.30` | 0.479 | 0.477 | 882 | 1.4 |
| `E >= 0.35` | 0.495 | 0.490 | 491 | 1.5 |
| `E >= 0.40` | 0.509 | 0.488 | 255 | 1.6 |
| `E >= 0.42` | 0.506 | 0.488 | 182 | 1.6 |
| **`E >= 0.45`** (live) | **0.520** | 0.487 | 88 | 1.6 |
| `E >= 0.47` | 0.500 | 0.476 | 42 | 1.6 |

The live cutoff is the maximum. Loosening it does not recover agreement — it adds 8-mers the
two labs disagree about at least as often. **The cutoff is not the problem.**

"Size ratio" is the larger positive count over the smaller, within a pair. At a fixed cutoff
the median pair has one lab calling 1.6x as many positives as the other, and the tail is far
worse: `Cell08`/`weirauch2014` on `Hoxa2` is 165 against 17, `SCI09`/`weirauch2014` on `Atf1`
is 13 against 96. Same statistic, same 8-mer set, same protein.

## Rank-matched rules

| rule | median Jaccard | mean | median positives | median size ratio |
|---|---:|---:|---:|---:|
| top 50 per experiment | **0.639** | 0.594 | 50 | 1.0 |
| top 100 per experiment | 0.600 | 0.583 | 100 | 1.0 |
| top 200 per experiment | 0.581 | 0.576 | 200 | 1.0 |
| top 500 per experiment | 0.555 | 0.531 | 500 | 1.0 |
| top 0.1% (33 8-mers) | 0.561 | 0.566 | 32 | 1.0 |
| top 0.3% (99 8-mers) | 0.581 | 0.580 | 98 | 1.0 |
| top 1% (329 8-mers) | 0.585 | 0.556 | 328 | 1.0 |

**Top-100 beats `E >= 0.45` by 8 points at the same stringency** — the live cutoff already
yields a median of 88 positives, so this is not a size artifact. Top-50 beats it by 12.

The mechanism is visible in the size-ratio column. A fixed cutoff lets each experiment's own
scale decide how many positives it contributes; a rank rule takes the same depth from every
experiment and compares like with like.

This is the same conclusion the Kock screening reached from the other direction: its
`affinityEstimate` has no transferable scale at all, and only a rate-calibrated rule lands
in the right place. Our E-scores are far better behaved — but the difference is of degree.

## The 0.45 cutoff is silently switching two sources off

Positive rate at the live cutoff, per source:

| source | positive rate | | source | positive rate |
|---|---:|---|---|---:|
| `EMBO10` | 0.4889% | | `NAR11` | 0.1702% |
| `Cell08` | 0.4439% | | `SCI09` | 0.1804% |
| `RAD13A` | 0.3496% | | `weirauch2014` | 0.1559% |
| `MBE14` | 0.3359% | | `PNAS08` | 0.1398% |
| `PNAS13` | 0.3147% | | `LIN14B` | 0.0859% |
| `GR09` | 0.3000% | | `CB11` | 0.0365% |
| `ROG18A` | 0.2807% | | **`Cell09`** | **0.0048%** |
| `BAR15A` | 0.2742% | | **`LIU18B`** | **0.0015%** |
| `DEV12`, `MAR17A`, `SHO18A` | 0.23-0.27% | | | |

**54 of 1,383 domain records (3.9%) have no positive 8-mer at all.** They are not spread
evenly:

| source | domains | zero-positive | per-domain max E (median) |
|---|---:|---:|---:|
| `Cell09` | 17 | **12 (71%)** | **0.428** |
| `BAR15A` | 90 | 18 (20%) | 0.495 |
| `LIU18B` | 2 | 1 | 0.468 |
| `LIN14B` | 12 | 1 | 0.484 |
| `weirauch2014` | 868 | 22 (2.5%) | 0.495 |
| `Cell08` | 147 | 0 | 0.498 |

The two rows to read are `Cell09` and `BAR15A`, and they fail differently.

- **`Cell09` is a scale problem.** Its per-domain maximum E-score has a median of **0.428** —
  the whole source sits below the cutoff. 12 of its 17 domains never reach `0.45` at their
  single best 8-mer, so `Cell09` contributes 559,232 rows that say only "does not bind". That
  is not a measurement of non-binding; it is a source whose E-score distribution is
  compressed relative to everyone else's, silenced by a cutoff calibrated on other arrays.
  `LIU18B` is the same story with 2 domains and 1 positive in the entire source.
- **`BAR15A`'s 18 are probably real.** It is the variant-series source, and a variant that
  abolishes binding *should* have no positive 8-mer. Its per-domain max E has a median of
  0.495 with a minimum of 0.342 — a healthy distribution with a dead tail, not a compressed
  one. This is exactly why a naive rank rule would be wrong: top-100 would manufacture 100
  positives for each dead variant, destroying the signal the corpus exists to carry.

## What follows

The two failure modes point in opposite directions, so no single rule fixes both:

| | fixed cutoff | rank rule |
|---|---|---|
| compressed source (`Cell09`) | silences it | rescues it |
| dead variant (`BAR15A`) | correctly gives no positives | invents ~100 |

A rule that survives both has to be a **rank cap under an absolute floor** — positives are the
experiment's top *N* by E-score *and* above some floor, so a dead variant stays dead and a
compressed source is still read at its own scale. Both parameters would need setting, and by
this project's convention that invalidates the dataset, every report and every provenance row.

Tracked as `TODO.md` `T20`; `Cell09` and `LIU18B` are `T21`. This report does not decide
either, and it changes nothing until one is decided.

**One thing is settled by these numbers, though:** if `T19` succeeds and Kock's arrays are
rescored to real E-scores, they can enter on the existing scale with the existing cutoff, and
none of this has to be reopened first.

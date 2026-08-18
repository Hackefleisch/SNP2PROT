# E-scores recomputed from raw scans — `T19`, validated against `BAR15A`

Generated 2026-08-18. **The pipeline is built, tested and measured. Nothing has been admitted
to the corpus and no threshold has been changed.**

`T19` asked whether a source that publishes no E-scores can still enter on the project's own
scale, by recomputing the statistic from the raw scanner output. The answer is a qualified
yes, and the qualification is the whole content of this report: **the ranking reproduces, the
values do not.**

Code: [`src/snp2prot/rawpbm.py`](../src/snp2prot/rawpbm.py),
[`scripts/validate_escores.py`](../scripts/validate_escores.py),
[`tests/test_rawpbm.py`](../tests/test_rawpbm.py) (9 tests, whole suite 171 green).

```bash
.venv/bin/python scripts/validate_escores.py --manifest <arrays.tsv> --source BAR15A
```

One array takes ~2 s end to end, so scoring all 122 Kock samples is a few minutes, not a job.

---

## 1. What was validated, and against what

GEO `GSE233827` deposits the raw scans behind **`BAR15A`** as well as Kock's own arrays — the
Barrera-era `.gpr` files, 2012 to 2015. That makes a closed loop possible: recompute E-scores
from Barrera's scans and compare them with the UniPROBE values already in
`data/interim/BAR15A/`.

**32 arrays, 14 alleles, 2 genes** — the whole ARX series (6 alleles) and the whole HOXD13
series (8) — each allele scored from all of its replicate arrays.

| | live alleles (11) | dead variants (3) |
|---|---|---|
| Spearman, all 32,896 8-mers | 0.505 – 0.635 (median **0.614**) | 0.59 – 0.61 |
| Spearman, 8-mers the stored source puts above 0.3 | 0.778 – 0.905 (median **0.851**) | n/a |
| rank-matched positive set, Jaccard vs stored | 0.565 – 0.766 (median **0.639**) | n/a |
| highest recomputed E-score | 0.448 – 0.493 | **0.325 – 0.427** |
| 8-mers at `E >= 0.45` | 0 – 32 (median 7) | 0 |
| stored positives at `E >= 0.45` | 129 – 204 (median 147) | 0, 0, 2 |

Read the third row against `reports/overlap.md`: two replicates **within** one source agree at
0.70–0.72, and two different sources measuring one protein agree at 0.458 median. A
recomputation that lands at 0.64 median, 0.77 at best, is inside the band the corpus already
lives in — it is as close to UniPROBE as UniPROBE's own sources are to each other.

Read the last two rows and the premise of `T19` breaks.

## 2. The finding: the ranking transfers, the cutoff does not

Our E-scores are the same statistic on the same [-0.5, +0.5] scale, and they reach the same
ceiling (0.493 against UniPROBE's 0.495). What differs is the **density just below the
ceiling**: where UniPROBE puts 129–204 8-mers at or above 0.45, we put 0–32.

Ask instead what cutoff on the recomputed scale selects as many 8-mers as the stored source
called positive, and the answer is unusually stable:

| allele | count-matching cutoff | | allele | count-matching cutoff |
|---|---|---|---|---|
| `ARX_REF` | 0.343 | | `HOXD13_REF` | 0.340 |
| `ARX_P353L` | 0.324 | | `HOXD13_S316C` | 0.324 |
| `ARX_T333N` | 0.327 | | `HOXD13_I322L` | 0.321 |
| `HOXD13_R306W` | 0.319 | | `HOXD13_Q325R` | 0.318 |
| `HOXD13_N298S` | 0.318 | | `HOXD13_Q325K` | 0.312 |
| `HOXD13_I297V` | 0.313 | | **median** | **0.321** |

**0.312 to 0.343 over 11 alleles and two genes.** A flat cutoff at 0.35 on the recomputed
scale gives 82–192 positives per allele and Jaccard 0.47–0.76 (median 0.565) against the
stored positive sets.

So admitting a source this way costs **one recalibrated number**, not a second score type, not
a second scale, and not a second `neg_provenance` vocabulary. That is a materially better
position than the `affinityQ` route screened in
[`kock2024_screening.md`](kock2024_screening.md) §6 — but it is **not** the "nothing
invalidated" that `T19` was written on. Carrying `pbm.positive = 0.45` over unchanged would
give this source a positive rate roughly a tenth of every other source's, which is the same
label-prior leak in the other direction.

**The dead variants are the good news.** `ARX_L343Q`, `ARX_P353R` and `ARX_R332H` are stored
with 0, 0 and 2 positives; recomputed they top out at 0.325, 0.427 and 0.394 and yield 0, 1
and 2 positives at a 0.35 cutoff. The absolute floor still does the work a rank rule cannot —
which is precisely the objection `T20` raises against rank-matched labelling. **A recomputed
E-score under a recalibrated absolute floor satisfies both halves of `T20` at once.**

## 3. Three implementation decisions, each found by measurement

**Join by `(Column, Row)`, never by probe ID.** Inherited from
`data/external/pbm_design/README.md` and confirmed here. The ID sets match exactly and name
different probes; joining by ID gives rho = -0.003.

**Keep saturated spots.** This was the largest single correction. Masking them looks obviously
right and is obviously wrong: on one ARX array 274 spots are saturated and **260 of them are
in the top 1,000 by signal**. Masking cost the top 8-mer `CTAATTAG` 12 of its 18 probes and
`GCTAATTA` 12 of its 37 — the E-score is a rank statistic, so a ceiling-clipped spot still
carries the only fact that matters. Keeping them lifted `ARX_REF` from 0.700 to 0.766 ranked
Jaccard and restored the top of the scale.

**Normalise Cy3 as a residual, not a ratio.** Dividing the protein channel by the raw Cy3
signal is what "Cy3 normalization" sounds like, and it destroys the result: 0.489 → 0.123 on
one array. The double-stranding signal is 53–74% explained by probe sequence composition, and
that part is not an artifact — A/T-rich probes carry both more dsDNA and more homeodomain
binding. Regressing Cy3 on mono- and dinucleotide composition and dividing only by the
unexplained residual removes the spot defect and leaves the biology: 0.595 → 0.621.

Spatial detrending over tiles was also tried and is **not** in the pipeline: it moved
agreement by at most ±0.01 once the other three were right, and it is one more knob.

## 4. What did not work, and what is still open

**More replicates do not fix the compression.** `HOXD13_REF` has five arrays; every subset was
scored, and the best 2-, 3-, 4- and 5-array combinations give 5, 8, 6 and 11 8-mers above
0.45 against the stored 144. The gap is not sampling noise.

The likeliest cause is **masliner**: Barrera's pipeline scans each array at several laser
powers and stitches them into one extended-dynamic-range image before scoring. GEO deposits
one scan per array, so that step cannot be reproduced from this deposit at all. This is worth
stating plainly because it bounds the method — no amount of engineering on our side recovers
it.

**Array quality control is unsolved.** Kock's own processing filters replicates on "upper-tail
probe intensity width", and it does separate good arrays from bad: on the HOXD13 series one
slide gave 1.1–1.7 for all eight chambers and reproduced badly (ranked Jaccard 0.10–0.39)
while the other gave 3.5–4.1 and reproduced well (0.61–0.70). **But the metric cannot be used
as a gate**, because a variant that has lost DNA binding has a narrow tail for the honest
reason: the three dead ARX variants measure 1.4–2.1 on arrays that are fine. Filtering on it
would delete exactly the negative evidence the corpus most needs. `rawpbm.tail_width` reports
it and nothing filters on it; a gate that distinguishes a bad slide from a dead protein — a
per-slide comparison across chambers, most likely — is still to be designed.

**Coverage of the validation is narrow.** 14 alleles, 2 genes, one lab, one array design. The
count-matching cutoff is tight across those 14, but it has not been tested on a second family
or a second era of arrays.

## 5. Where this leaves `T11`

Kock can enter as `pbm_escore`, on the existing scale, with the existing `neg_provenance`
vocabulary and the existing replicate machinery — at the cost of one source-specific cutoff
near 0.32–0.35 whose value is *derived*, not chosen, from the `BAR15A` overlap this report
measures. Against the alternative (a second score type, two new cutoffs, and a 20–40× positive
rate mismatch), that is the better deal, and it is now demonstrated rather than assumed.

What it needs before it happens:

1. a decision on `T20`, since the cutoff question is now the same question in both places;
2. the array-QC gate of §4, or an explicit decision to score every array and let the
   reconciliation step absorb the bad ones;
3. the downloads into `data/raw/` with `PROVENANCE.md` rows — 122 samples, roughly 660 files
   and about 2 GB;
4. a parser that maps GEO samples to alleles, pairs each protein scan with the Cy3 scan of
   its own slide and chamber, and joins the result to the 93 clone sequences of Supplementary
   Data 3.

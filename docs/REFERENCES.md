# References

The paper behind each dataset in [TFDNA_MERGE_BRIEF.md](TFDNA_MERGE_BRIEF.md). Every DOI
below was **resolved against the Crossref API**, not recalled — title, journal, volume and
page were checked to match the brief's citation. Anything unverified is marked as such.

PDFs live in [papers/](papers/) (git-ignored) under the naming scheme documented there;
run `scripts/check_papers.py` to see which are still missing.

Read the paper before writing its parser: the binarization rule, the negative's provenance,
and the DBD boundary convention all come from the methods section, not from the file format.

## Tier 1 — PBM 8-mer E-scores

| dataset | paper | DOI | notes |
|---|---|---|---|
| `BAR15A` | Barrera et al. 2016, *Science* 351:1450 — "Survey of variation in human transcription factors reveals prevalent DNA binding changes" | [`10.1126/science.aad2257`](https://doi.org/10.1126/science.aad2257) | **The priority read.** 41 reference + 117 variant TF alleles; the matched WT/mutant pairs are the backbone of the protein axis. |
| mouse homeodomains | Berger et al. 2008, *Cell* 133:1266 — "Variation in Homeodomain DNA Binding Revealed by High-Resolution Analysis of Sequence Preferences" | [`10.1016/j.cell.2008.05.024`](https://doi.org/10.1016/j.cell.2008.05.024) | ~168 mouse HDs. UniPROBE accession to be looked up on the site, not guessed. |
| *Drosophila* homeodomains | Noyes et al. 2008, *Cell* 133:1277 — "Analysis of Homeodomain Specificities Allows the Family-wide Prediction of Preferred Recognition Sites" | [`10.1016/j.cell.2008.05.023`](https://doi.org/10.1016/j.cell.2008.05.023) | Companion paper to Berger, same issue. |
| UniPROBE itself | Hume et al. 2015, *NAR* 43:D117 — "UniPROBE, update 2015" | [`10.1093/nar/gku1045`](https://doi.org/10.1093/nar/gku1045) | Current database paper — file formats, E-score definition, accession scheme. |
| UniPROBE, original | Newburger & Bulyk 2009, *NAR* 37:D77 — "UniPROBE: an online database of protein binding microarray data on protein-DNA interactions" | [`10.1093/nar/gkn660`](https://doi.org/10.1093/nar/gkn660) | Only if the 2015 update is unclear on something. |

## Tier 2 — C2H2 zinc-finger B1H

| dataset | paper | DOI | notes |
|---|---|---|---|
| `persikov2015` | Persikov et al. 2015, *NAR* 43(3):1965–1984 — "A systematic survey of the Cys2His2 zinc finger DNA-binding landscape" | [`10.1093/nar/gku1395`](https://doi.org/10.1093/nar/gku1395) | ⚠️ The brief cites this correctly by volume/page, but the DOI is **`gku1395`** — a plausible-looking `gkv030` is a different article. Needed to settle the provisional 3-AT binarization rule (open item #4). |
| `najafabadi2015` | Najafabadi et al. 2015, *Nat Biotechnol* 33:555 — "C2H2 zinc finger proteins greatly expand the human regulatory lexicon" | [`10.1038/nbt.3128`](https://doi.org/10.1038/nbt.3128) | Defines the B1H s-score; needed before `sscore_positive` / `sscore_negative` can be set. |

## Tier 3 — DNA-side depth

| dataset | paper | DOI | notes |
|---|---|---|---|
| SNP-SELEX / GVAT | Yan et al. 2021, *Nature* 591:147 — "Systematic analysis of binding of transcription factors to noncoding variants" | [`10.1038/s41586-021-03211-0`](https://doi.org/10.1038/s41586-021-03211-0) | Defines OBS and PBS. Needed before the Phase 4 OBS cutoff can be chosen (open item #5). |
| Codebook (optional) | Vorontsov et al. 2025, *Commun Biol* — "Cross-platform motif discovery and benchmarking to explore binding specificities of poorly studied human transcription factors" | [`10.1038/s42003-025-08909-9`](https://doi.org/10.1038/s42003-025-08909-9) | Peer-reviewed Codebook/GRECO-BIT consortium paper. Citation taken from the PDF's own metadata. |
| Codebook (optional) | Jolma et al. 2024, bioRxiv **preprint** (v2 posted Oct 2025) — "GHT-SELEX demonstrates unexpectedly high intrinsic sequence specificity and complex DNA binding of many human transcription factors" | [`10.1101/2024.11.11.618478`](https://doi.org/10.1101/2024.11.11.618478) | ⚠️ **Not peer reviewed.** Codebook is a consortium with several companion papers; the brief marks it lowest priority. Do not cite this where the Vorontsov paper covers the same ground. |

## Tier 4 — held-out quantitative test sets (never trained on)

| dataset | paper | DOI | notes |
|---|---|---|---|
| BET-seq | Le et al. 2018, *PNAS* 115:E3702 — "Comprehensive, high-resolution binding energy landscapes reveal context dependencies of transcription factor binding" | [`10.1073/pnas.1715888115`](https://doi.org/10.1073/pnas.1715888115) | Pho4 / Cbf1, ~10⁶ × 16 bp with real ΔG. (Preprint version is `10.1101/193904` — cite the PNAS one.) |
| STAMMP | Aditham et al. 2021, *Cell Systems* 12(2):112–127 — "High-Throughput Affinity Measurements of Transcription Factor and DNA Mutations Reveal Affinity and Specificity Determinants" | [`10.1016/j.cels.2020.11.012`](https://doi.org/10.1016/j.cels.2020.11.012) | ~210 Pho4 variants × 9 oligos. |
| MAX / k-STAMMP | Hastings et al. 2025, *Nat Commun* 16:636 — "Mutations to transcription factor MAX allosterically increase DNA selectivity by altering folding and binding pathways" | [`10.1038/s41467-024-55672-2`](https://doi.org/10.1038/s41467-024-55672-2) | Volume 16, article 636 — matches the brief. Note the `2024` in the DOI despite the 2025 issue; that is correct. |

## Supporting metadata

| resource | paper | DOI |
|---|---|---|
| CIS-BP | Weirauch et al. 2014, *Cell* 158:1431 — "Determination and Inference of Eukaryotic Transcription Factor Sequence Specificity" | [`10.1016/j.cell.2014.08.009`](https://doi.org/10.1016/j.cell.2014.08.009) |
| Pfam / InterPro / HMMER, UniProt | current database release papers | not resolved — these are tooling citations needed at write-up time, not for parsing decisions |

## Reading order for the parsing work

1. **Barrera 2016** — sets up the whole protein axis; read before Phase 1.
2. **Hume 2015** (UniPROBE) — file formats and the E-score scale, alongside Barrera.
3. **Persikov 2015** — settles open item #4 before Phase 3.
4. **Yan 2021** — settles open item #5 before Phase 4.
5. Berger / Noyes 2008 — before Phase 2.
6. Tier 4 papers — before Phase 6 only.

Codebook (Vorontsov, Jolma) is read only if Tier 3 turns out to need extra motif-centered
windows; the brief marks it lowest priority.

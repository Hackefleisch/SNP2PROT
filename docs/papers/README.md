# Papers

PDFs of the source papers. **Git-ignored** — publisher PDFs are copyrighted and would bloat
the repo. This README and [../REFERENCES.md](../REFERENCES.md) are tracked, so the manifest
travels with the project even though the files do not.

## Naming

```
<firstauthor><year>_<slug>.pdf              main article
<firstauthor><year>_<slug>_supp.pdf         its supplement, when there is one
<firstauthor><year>_<slug>_supp-<what>.pdf  when a paper has several supplement documents
```

All lowercase, hyphens inside the slug, underscore only before `supp`. The extension follows
the file — one Aditham supplement is a `.docx` and stays one. A supplement always shares the
main article's stem, so `ls` groups a paper with its own material.

⚠️ **Do not trust publisher filenames.** Elsevier's `mmc` numbering is not "supplement N":
`…S0092867414010368-mmc7.pdf` was the Weirauch **main article** and `…S2405471220304634-mmc4.pdf`
was the Aditham **main article**. Both were identified from the PDF's own title metadata.
Check contents before renaming, never the download name.

Check what is present, missing, or unlisted:

```bash
.venv/bin/python scripts/check_papers.py
.venv/bin/python scripts/check_papers.py --phase 1     # only what Phase 1 needs
```

## Manifest

Phase = when it is needed. `!` = blocking that phase. Supplement rows carry `↳`.

| filename | | paper | DOI |
|---|---|---|---|
| `barrera2016_human-tf-variation.pdf` | ! P1 | Barrera et al. 2016, *Science* 351:1450 | [10.1126/science.aad2257](https://doi.org/10.1126/science.aad2257) |
| `barrera2016_human-tf-variation_supp.pdf` | ! P1 | ↳ supplementary materials (35 pp) | |
| `hume2015_uniprobe-update.pdf` | ! P1 | Hume et al. 2015, *NAR* 43:D117 | [10.1093/nar/gku1045](https://doi.org/10.1093/nar/gku1045) |
| `newburger2009_uniprobe.pdf` | P1 | Newburger & Bulyk 2009, *NAR* 37:D77 | [10.1093/nar/gkn660](https://doi.org/10.1093/nar/gkn660) |
| `weirauch2014_cis-bp.pdf` | P1 | Weirauch et al. 2014, *Cell* 158:1431 | [10.1016/j.cell.2014.08.009](https://doi.org/10.1016/j.cell.2014.08.009) |
| `weirauch2014_cis-bp_supp-table-s1.pdf` | P1 | ↳ Table S1 legend, sources of other TF motifs | |
| `berger2008_homeodomain-pbm.pdf` | P2 | Berger et al. 2008, *Cell* 133:1266 | [10.1016/j.cell.2008.05.024](https://doi.org/10.1016/j.cell.2008.05.024) |
| `berger2008_homeodomain-pbm_supp.pdf` | P2 | ↳ supplemental data (59 pp) | |
| `noyes2008_homeodomain-specificities.pdf` | P2 | Noyes et al. 2008, *Cell* 133:1277 | [10.1016/j.cell.2008.05.023](https://doi.org/10.1016/j.cell.2008.05.023) |
| `noyes2008_homeodomain-specificities_supp.pdf` | P2 | ↳ supplemental data (64 pp) | |
| `persikov2015_c2h2-zf-landscape.pdf` | ! P3 | Persikov et al. 2015, *NAR* 43:1965 | [10.1093/nar/gku1395](https://doi.org/10.1093/nar/gku1395) |
| `najafabadi2015_c2h2-regulatory-lexicon.pdf` | P3 | Najafabadi et al. 2015, *Nat Biotechnol* 33:555 | [10.1038/nbt.3128](https://doi.org/10.1038/nbt.3128) |
| `najafabadi2015_c2h2-regulatory-lexicon_supp.pdf` | P3 | ↳ supplementary information (14 pp) | |
| `yan2021_snp-selex.pdf` | ! P4 | Yan et al. 2021, *Nature* 591:147 | [10.1038/s41586-021-03211-0](https://doi.org/10.1038/s41586-021-03211-0) |
| `yan2021_snp-selex_supp-motif-enrichment.pdf` | P4 | ↳ HOMER motif-enrichment output (679 pp) | |
| `vorontsov2025_codebook-benchmarking.pdf` | P4 | Vorontsov et al. 2025, *Commun Biol* | [10.1038/s42003-025-08909-9](https://doi.org/10.1038/s42003-025-08909-9) |
| `vorontsov2025_codebook-benchmarking_supp.pdf` | P4 | ↳ supplementary information (16 pp) | |
| `vorontsov2025_codebook-benchmarking_supp-data-index.pdf` | P4 | ↳ index of supplementary data files | |
| `jolma2024_ght-selex-preprint.pdf` | P4 | Jolma et al. 2024, bioRxiv **preprint** (v2, Oct 2025) | [10.1101/2024.11.11.618478](https://doi.org/10.1101/2024.11.11.618478) |
| `le2018_bet-seq.pdf` | P6 | Le et al. 2018, *PNAS* 115:E3702 | [10.1073/pnas.1715888115](https://doi.org/10.1073/pnas.1715888115) |
| `le2018_bet-seq_supp.pdf` | P6 | ↳ supporting information (32 pp) | |
| `aditham2021_stammp.pdf` | P6 | Aditham et al. 2021, *Cell Syst* 12:112 | [10.1016/j.cels.2020.11.012](https://doi.org/10.1016/j.cels.2020.11.012) |
| `aditham2021_stammp_supp-figures.pdf` | P6 | ↳ supplemental figures S1–S29 | |
| `aditham2021_stammp_supp-legends.docx` | P6 | ↳ figure legends / contents for the above | |
| `hastings2025_max-selectivity.pdf` | P6 | Hastings et al. 2025, *Nat Commun* 16:636 | [10.1038/s41467-024-55672-2](https://doi.org/10.1038/s41467-024-55672-2) |
| `hastings2025_max-selectivity_supp.pdf` | P6 | ↳ supplementary information (73 pp) | |
| `hastings2025_max-selectivity_supp-data-index.pdf` | P6 | ↳ index of supplementary data files | |

## Known gaps

Not asserted to exist — verify on the publisher's page before hunting.

- **Persikov 2015 supplement.** No supplementary document present, and this is the one that
  matters most: open item #4 (the provisional 2 mM / 10 mM 3-AT binarization rule) is likely
  settled in the supplementary methods rather than the 20-page main text. Worth checking the
  NAR page before Phase 3.
- **Yan 2021 supplementary information.** Only the 679-page HOMER motif-enrichment dump
  (`MOESM12`) is here. The narrative Supplementary Information — where the OBS/PBS definitions
  and cutoffs live — is a different MOESM and would settle open item #5.
- **UniPROBE database papers** (Hume, Newburger) have no supplement; nothing missing.
- **Jolma preprint** supplementary files are hosted separately on bioRxiv; low priority, since
  Codebook is the brief's lowest-priority source.

Tables and pure figure sets were deliberately skipped — only documents were collected.

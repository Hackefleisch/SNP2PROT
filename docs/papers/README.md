# Papers

PDFs of the source papers. **Git-ignored** — publisher PDFs are copyrighted and would bloat
the repo. This README and [../REFERENCES.md](../REFERENCES.md) are tracked, so the manifest
travels with the project even though the files do not.

## Naming

```
<firstauthor><year>_<slug>.pdf      all lowercase, hyphens in the slug
```

First author surname, publication year, then a short topic slug. No spaces, no journal names,
no "(1)" suffixes from the browser. Consistent names mean a parser's docstring can point at
its paper by filename and the reference stays findable.

Check what is present and what is missing:

```bash
.venv/bin/python scripts/check_papers.py
```

## Manifest

Phase column = when it is needed. Priority `!` = blocking the next phase.

| filename | | paper | DOI |
|---|---|---|---|
| `barrera2016_human-tf-variation.pdf` | ! P1 | Barrera et al. 2016, *Science* 351:1450 | [10.1126/science.aad2257](https://doi.org/10.1126/science.aad2257) |
| `hume2015_uniprobe-update.pdf` | ! P1 | Hume et al. 2015, *NAR* 43:D117 | [10.1093/nar/gku1045](https://doi.org/10.1093/nar/gku1045) |
| `newburger2009_uniprobe.pdf` | P1 | Newburger & Bulyk 2009, *NAR* 37:D77 | [10.1093/nar/gkn660](https://doi.org/10.1093/nar/gkn660) |
| `berger2008_homeodomain-pbm.pdf` | P2 | Berger et al. 2008, *Cell* 133:1266 | [10.1016/j.cell.2008.05.024](https://doi.org/10.1016/j.cell.2008.05.024) |
| `noyes2008_homeodomain-specificities.pdf` | P2 | Noyes et al. 2008, *Cell* 133:1277 | [10.1016/j.cell.2008.05.023](https://doi.org/10.1016/j.cell.2008.05.023) |
| `persikov2015_c2h2-zf-landscape.pdf` | ! P3 | Persikov et al. 2015, *NAR* 43:1965 | [10.1093/nar/gku1395](https://doi.org/10.1093/nar/gku1395) |
| `najafabadi2015_c2h2-regulatory-lexicon.pdf` | P3 | Najafabadi et al. 2015, *Nat Biotechnol* 33:555 | [10.1038/nbt.3128](https://doi.org/10.1038/nbt.3128) |
| `yan2021_snp-selex.pdf` | ! P4 | Yan et al. 2021, *Nature* 591:147 | [10.1038/s41586-021-03211-0](https://doi.org/10.1038/s41586-021-03211-0) |
| `weirauch2014_cis-bp.pdf` | P1+ | Weirauch et al. 2014, *Cell* 158:1431 | [10.1016/j.cell.2014.08.009](https://doi.org/10.1016/j.cell.2014.08.009) |
| `le2018_bet-seq.pdf` | P6 | Le et al. 2018, *PNAS* 115:E3702 | [10.1073/pnas.1715888115](https://doi.org/10.1073/pnas.1715888115) |
| `aditham2021_stammp.pdf` | P6 | Aditham et al. 2021, *Cell Syst* 12:112 | [10.1016/j.cels.2020.11.012](https://doi.org/10.1016/j.cels.2020.11.012) |
| `hastings2025_max-selectivity.pdf` | P6 | Hastings et al. 2025, *Nat Commun* 16:636 | [10.1038/s41467-024-55672-2](https://doi.org/10.1038/s41467-024-55672-2) |

Supplementary files, where a paper's supplement carries the actual data description, go
alongside as `<same-stem>_supp.pdf` (or `_supp01.xlsx` etc.).

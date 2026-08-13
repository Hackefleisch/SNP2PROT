# Provenance

One row per file under `data/raw/`. Append-only: never edit or delete a row, and never edit
the file it describes. This table is publication evidence, not a convenience log.

Regenerate the hash/size columns for a file with:

```bash
python scripts/record_provenance.py data/raw/<source>/<file>
```

## Downloaded files

| source_dataset | file (rel. to `data/raw/`) | source URL | accession | downloaded (UTC) | bytes | sha256 | contents |
|---|---|---|---|---|---|---|---|
| BAR15A | `BAR15A/BAR15A_contig8mers.zip` | https://thebrain.bwh.harvard.edu/uniprobe/downloads/BAR15A/BAR15A_contig8mers.zip | BAR15A | 2026-08-13T11:29:06Z | 104,707,619 | `916f2f266ba154b5c52a03dd414f203c63f1a6984fb22bfec274cf9707389cee` | UniPROBE normalized 8-mer data, Barrera 2016. 208 files, `GENE/GENE_ALLELE/REP/*_8mers.txt`; 161 alleles x 32,896 8-mers with E-score and Z-score. |
| BAR15A | `BAR15A/details/*.html` (41 files) | https://thebrain.bwh.harvard.edu/uniprobe/details34.php?id=584..624 | BAR15A | 2026-08-13T11:33Z | see manifest | see manifest | UniPROBE detail page per gene, carrying the clone insert sequence for every allele plus Pfam domain, Swiss-Prot accession and species. Per-file URL, size and sha256 are in `BAR15A/details/_manifest.tsv` (sha256 `5d4b101cfdde12d34edf2ea7527714cd89950c794aef630dfada8c878ee22c34`), one row per page, rather than 41 rows here. |

## Sources attempted and rejected

Record anything that turned out to be unusable under the project's constraints, with the
reason. A short honest catalogue beats a padded one.

| source | date | outcome | reason |
|---|---|---|---|
| Barrera 2016 Table S4 (Science supplement) | 2026-08-13 | not needed | science.org returns HTTP 403 to non-browser clients, so the Excel supplements are a manual download. Turned out to be unnecessary: UniPROBE's own detail pages carry the clone insert sequence for every allele, which is the authoritative record of what was on the array. |

## Manual-download queue

Sources behind a licence click-through, MTA, or an interactive portal, which the agent
cannot fetch unattended. Each has a `data/raw/<source>/HOWTO.md` with exact steps.

| source | blocker | HOWTO | status |
|---|---|---|---|
| Barrera 2016 Excel supplements (Tables S1-S4, S6, S7) | science.org blocks non-browser clients (HTTP 403) | not written | **not currently needed** — only fetch if Phase 5 needs the population-genetics annotations |

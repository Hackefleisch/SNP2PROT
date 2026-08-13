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
| Cell08 | `Cell08/Cell08_contig8mers.zip` | https://thebrain.bwh.harvard.edu/uniprobe/downloads/Cell08/Cell08_contig8mers.zip | Cell08 | 2026-08-13T11:52Z | 356,413,560 | `58d70b6f0962dbf418fee2d499bdf6ff56429882636f7c5026410bef21ec0876` | UniPROBE contiguous 8-mer E-scores, Berger et al., Cell 2008. ~168 mouse homeodomains; 178 experiments. |
| Cell08 | `Cell08/details/*.html` (168 files) | https://thebrain.bwh.harvard.edu/uniprobe/ (per-gene detail pages) | Cell08 | 2026-08-13T11:56Z | see manifest | see manifest | Detail page per gene: Pfam domain, UniProt accession, species, and the DBD and/or clone insert sequence. Per-file URL/size/sha256 in `Cell08/details/_manifest.tsv` (sha256 `2084d7b0c8546f7d798a14c8eb64fde30604cdb37ecaf9ec376b225a1b708942`). |
| EMBO10 | `EMBO10/EMBO10_contig8mers.zip` | https://thebrain.bwh.harvard.edu/uniprobe/downloads/EMBO10/EMBO10_contig8mers.zip | EMBO10 | 2026-08-13T11:52Z | 10,884,062 | `802714cc1c6004a14afd3a060f6b836052008e77d3e6e81d38b2de903bad0091` | UniPROBE contiguous 8-mer E-scores, Wei et al., EMBO J 2010. 22 mouse ETS-family TFs. |
| EMBO10 | `EMBO10/details/*.html` (22 files) | https://thebrain.bwh.harvard.edu/uniprobe/ (per-gene detail pages) | EMBO10 | 2026-08-13T11:56Z | see manifest | see manifest | Detail page per gene: Pfam domain, UniProt accession, species, and the DBD and/or clone insert sequence. Per-file URL/size/sha256 in `EMBO10/details/_manifest.tsv` (sha256 `546a86d335300cae6bc903f6e12e1ead268539bdf06f5122cd29238f16df4766`). |
| PNAS13 | `PNAS13/PNAS13_contig8mers.zip` | https://thebrain.bwh.harvard.edu/uniprobe/downloads/PNAS13/PNAS13_contig8mers.zip | PNAS13 | 2026-08-13T11:52Z | 10,599,624 | `bb2fd410feeb95260bf024d11de99709746d9630b4ef2827682b8c4249f3bacf` | UniPROBE contiguous 8-mer E-scores, Nakagawa et al., PNAS 2013. 20 forkhead TFs spanning human, fly and fungi. |
| PNAS13 | `PNAS13/details/*.html` (20 files) | https://thebrain.bwh.harvard.edu/uniprobe/ (per-gene detail pages) | PNAS13 | 2026-08-13T11:56Z | see manifest | see manifest | Detail page per gene: Pfam domain, UniProt accession, species, and the DBD and/or clone insert sequence. Per-file URL/size/sha256 in `PNAS13/details/_manifest.tsv` (sha256 `be00e413cd33e0b4171e7cfb1a9daf4b70d52e1f1fceee08e43cf427eb83213e`). |
| pfam | `../external/pfam/PF00010.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00010/?annotation=hmm | PF00010 | 2026-08-13T12:10Z | 25,851 | `8de4f7baac46ff96fab3591914cb525d...` | Pfam HMM for HLH, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00046.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00046/?annotation=hmm | PF00046 | 2026-08-13T12:10Z | 27,694 | `a366390981bc1cdcb14f8c6f72d65b2c...` | Pfam HMM for Homeodomain, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00096.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00096/?annotation=hmm | PF00096 | 2026-08-13T12:10Z | 11,893 | `7fe8d4f0a17470766133eaee65921401...` | Pfam HMM for zf-C2H2, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00105.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00105/?annotation=hmm | PF00105 | 2026-08-13T12:10Z | 33,753 | `d1a18903a2248909b825411371658612...` | Pfam HMM for zf-C4, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00157.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00157/?annotation=hmm | PF00157 | 2026-08-13T12:10Z | 34,683 | `baa8779de2025f7b777dc1ad4020e148...` | Pfam HMM for Pou, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00178.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00178/?annotation=hmm | PF00178 | 2026-08-13T12:10Z | 39,770 | `a0e1768113fb0b8ef35310340b9f7e33...` | Pfam HMM for Ets, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00250.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00250/?annotation=hmm | PF00250 | 2026-08-13T12:10Z | 42,102 | `eacba028c1423b307b4c34daae29c0b4...` | Pfam HMM for Forkhead, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF00292.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00292/?annotation=hmm | PF00292 | 2026-08-13T12:10Z | 59,305 | `99d81e8c6bb9116cae7512f63dd2d1be...` | Pfam HMM for PAX, used for DBD boundary annotation (HMMER3/f 3.3). |
| pfam | `../external/pfam/PF02178.hmm` | https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF02178/?annotation=hmm | PF02178 | 2026-08-13T12:10Z | 7,234 | `b584aa1b1ae57052ae5701c28fd6ac4f...` | Pfam HMM for AT_hook, used for DBD boundary annotation (HMMER3/f 3.3). |
| uniprot | `../external/uniprot/*.fasta` (204 files) | https://rest.uniprot.org/uniprotkb/<accession>.fasta | per-protein | 2026-08-13T12:15Z | varies | per-file | Canonical UniProt sequences for the protein table, one FASTA per accession, fetched on demand and cached. |

## Sources attempted and rejected

Record anything that turned out to be unusable under the project's constraints, with the
reason. A short honest catalogue beats a padded one.

| source | date | outcome | reason |
|---|---|---|---|
| Persikov 2015 B1H + Najafabadi 2015 C2H2 (Phase 3) | 2026-08-13 | **dropped** | Every entry is a C2H2 zinc-finger array: 2-6 separate ~23-residue folds on flexible linkers, each needing its own Zn(2+). Fails condition 2 of the domain policy (one continuous region). Persikov additionally varies one finger inside a fixed three-finger context, so the subunit that varies and the subunit that binds differ, failing conditions 1 and 3. ~8,000 domains excluded. See `docs/DOMAIN_POLICY.md`. |
| Noyes et al. 2008 (Drosophila homeodomains) | 2026-08-13 | not applicable to Phase 2 | Not in UniPROBE under any accession. Noyes used bacterial one-hybrid, not PBM, so it belongs to the Phase 3 B1H work rather than the UniPROBE panels. |
| Barrera 2016 Table S4 (Science supplement) | 2026-08-13 | not needed | science.org returns HTTP 403 to non-browser clients, so the Excel supplements are a manual download. Turned out to be unnecessary: UniPROBE's own detail pages carry the clone insert sequence for every allele, which is the authoritative record of what was on the array. |

## Manual-download queue

Sources behind a licence click-through, MTA, or an interactive portal, which the agent
cannot fetch unattended. Each has a `data/raw/<source>/HOWTO.md` with exact steps.

| source | blocker | HOWTO | status |
|---|---|---|---|
| Barrera 2016 Excel supplements (Tables S1-S4, S6, S7) | science.org blocks non-browser clients (HTTP 403) | not written | **not currently needed** — only fetch if Phase 5 needs the population-genetics annotations |

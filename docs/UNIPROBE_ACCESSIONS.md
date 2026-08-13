# UniPROBE accession survey

Every accession on <https://thebrain.bwh.harvard.edu/uniprobe/downloads.php>, with the
citation UniPROBE itself gives and a **heuristic** family count. Built in Phase 2 so that
the choice of panels is a decision on record rather than a guess.

The family columns count gene-name prefixes (`Fox*`, `Ets*/Elk*/Elf*`, `Atf*/Cebp*/Jun*`,
homeodomain stems). They are a **survey aid only** — the family actually written into
`dbd_family` comes from each detail page's Pfam `Domain` field at parse time.

Every accession exposes `<ACC>_contig8mers.zip`, the same contiguous-8-mer format used in
Phases 1-2. Downloads are not licence-gated (open item #6).

| accession | citation | proteins | forkhead | ETS | bZIP | homeodomain | parsed |
|---|---|---:|---:|---:|---:|---:|---|
| `Cell08` | Berger et al., Cell 2008 | 168 | 0 | 0 | 0 | 121 | Phase 2 |
| `SCI09` | Badis et al., Science 2009 | 104 | 5 | 5 | 4 | 3 | — |
| `GR09` | Zhu et al., Genome Res 2009 | 89 | 0 | 0 | 1 | 0 | — |
| `MAR17A` | Mariani et al., Cell Systems 2017 | 61 | 8 | 1 | 8 | 3 | — |
| `SHO18A` | Shokri et al., Cell Reports 2019 | 45 | 0 | 0 | 2 | 0 | — |
| `BAR15A` | Barrera et al., Science 2016 | 41 | 1 | 0 | 0 | 20 | Phase 1 |
| `GB11` | Gordan et al., Gen. Bio. 2011 | 27 | 0 | 0 | 0 | 0 | — |
| `EMBO10` | Wei et al., EMBO J 2010 | 22 | 0 | 20 | 0 | 0 | Phase 2 |
| `Cell09` | Grove et al., Cell 2009 | 21 | 0 | 0 | 0 | 0 | — |
| `PNAS13` | Nakagawa et al., PNAS 2013 | 20 | 11 | 0 | 0 | 0 | Phase 2 |
| `LAI20A` | Lai et al. 2020 | 20 | 3 | 0 | 1 | 3 | — |
| `Path10` | Campbell et al., PLoS Pathog 2010 | 19 | 0 | 0 | 0 | 0 | — |
| `LIN14B` | Lindemose et al, Nucleic Acids Res. 2014 | 13 | 0 | 0 | 0 | 0 | — |
| `DEV12` | Busser et al., Development 2012 | 10 | 0 | 0 | 0 | 1 | — |
| `ROG18A` | Rogers et al., Mol Cell 2019 | 7 | 7 | 0 | 0 | 0 | — |
| `PO10` | Del Bianco et al., PLoS ONE 2010 | 6 | 0 | 0 | 0 | 0 | — |
| `NBT06` | Berger et al., Nat Biotech 2006 | 5 | 0 | 0 | 0 | 0 | — |
| `PNAS08` | De Silva et al., PNAS 2008 | 3 | 0 | 0 | 0 | 0 | — |
| `NAR10` | Alibes et al., NAR 2010 | 3 | 0 | 0 | 2 | 0 | — |
| `GD12` | Peterson et al., Genes Dev 2012 | 3 | 0 | 0 | 0 | 0 | — |
| `LIU18B` | Liu et al., eLife 2018 | 2 | 0 | 0 | 0 | 0 | — |
| `LIU18A` | Liu et al., Cell 2018 | 2 | 0 | 0 | 0 | 0 | — |
| `PP15` | Lehti-Shiu et al., PP 2015 | 2 | 0 | 0 | 0 | 0 | — |
| `KUR17A` | Li et al., Nature 2017 | 2 | 0 | 0 | 0 | 0 | — |
| `MBE14` | Cheatle Jarvela et al., Mol Biol Evol 2014 | 2 | 0 | 0 | 0 | 0 | — |
| `GD13` | Soruco et al., Genes Dev 2013 | 1 | 0 | 0 | 0 | 0 | — |
| `NAR11` | De Masi et al., NAR 2011 | 1 | 0 | 0 | 0 | 0 | — |
| `PNAS12` | Busser et al., PNAS 2012 | 1 | 0 | 0 | 0 | 0 | — |
| `CB11` | Helfer et al., Curr Biol 2011 | 1 | 0 | 0 | 0 | 0 | — |
| `MMB08` | Pompeani et al., Mol Microbiol 2008 | 1 | 0 | 0 | 0 | 0 | — |
| `GD09` | Lesch et al., Genes Dev 2009 | 1 | 0 | 0 | 0 | 0 | — |
| `STI21B` | Stielow et al., Science Advances 2021 | 1 | 0 | 0 | 0 | 0 | — |
| `CR09` | Scharer et al., Cancer Res 2009 | 1 | 0 | 0 | 0 | 0 | — |
| `RAD13A` | Radke et al., PNAS 2013 | 1 | 0 | 0 | 0 | 0 | — |
| `MIZ19A` | Mizeracka et al. 2019 | 1 | 0 | 0 | 0 | 0 | — |

## What this settled

The brief asked for "any forkhead / ETS / bZIP family panels you find". They are:

- **ETS → `EMBO10`** (Wei et al., EMBO J 2010) — 20 of 22 proteins. Parsed.
- **forkhead → `PNAS13`** (Nakagawa et al., PNAS 2013) — parsed — and **`ROG18A`**
  (Rogers et al., Mol Cell 2019), 7 of 7 forkhead, not yet parsed.
- **bZIP → `MAR17A`** (Mariani et al., Cell Systems 2017), 8 bZIP alongside 8 forkhead,
  and `SCI09` (Badis et al., Science 2009), which spreads 104 mouse TFs across many
  families. Neither parsed yet.

**Noyes et al. 2008 is not in UniPROBE under any accession.** That study used bacterial
one-hybrid rather than PBM, so despite being the companion paper to Berger 2008 it belongs
to the Phase 3 B1H work, not to the UniPROBE panels.

## Obvious next candidates

| accession | why | proteins |
|---|---|---:|
| `SCI09` | broadest family spread in the database; Badis 2009 is the standard mouse PBM panel | 104 |
| `GR09` | Zhu et al. 2009, yeast TFs — adds a whole other clade of protein space | 89 |
| `MAR17A` | the clearest bZIP contribution, plus more forkhead | 61 |
| `SHO18A` | Shokri et al. 2019 | 45 |
| `ROG18A` | a second, independent forkhead panel — useful as a cross-source agreement check | 7 |

# UniPROBE accession survey

Every accession on <https://thebrain.bwh.harvard.edu/uniprobe/downloads.php>, with the
citation UniPROBE itself gives and a **heuristic** family count. Built in Phase 2 so that
the choice of panels is a decision on record rather than a guess.

The family columns count gene-name prefixes (`Fox*`, `Ets*/Elk*/Elf*`, `Atf*/Cebp*/Jun*`,
homeodomain stems). They are a **survey aid only** — the family actually written into
`dbd_family` comes from each detail page's Pfam `Domain` field at parse time.

Every accession exposes `<ACC>_contig8mers.zip`, the same contiguous-8-mer format used in
Phases 1-2. Downloads are not licence-gated (open item #6).

| accession | citation | proteins | forkhead | ETS | bZIP | homeodomain | outcome |
|---|---|---:|---:|---:|---:|---:|---|
| `Cell08` | Berger et al., Cell 2008 | 168 | 0 | 0 | 0 | 121 | **parsed** |
| `SCI09` | Badis et al., Science 2009 | 104 | 5 | 5 | 4 | 3 | **parsed** |
| `GR09` | Zhu et al., Genome Res 2009 | 89 | 0 | 0 | 1 | 0 | **parsed** |
| `MAR17A` | Mariani et al., Cell Systems 2017 | 61 | 8 | 1 | 8 | 3 | **parsed** |
| `SHO18A` | Shokri et al., Cell Reports 2019 | 45 | 0 | 0 | 2 | 0 | **parsed** |
| `BAR15A` | Barrera et al., Science 2016 | 41 | 1 | 0 | 0 | 20 | **parsed** |
| `GB11` | Gordan et al., Gen. Bio. 2011 | 27 | 0 | 0 | 0 | 0 | not usable |
| `EMBO10` | Wei et al., EMBO J 2010 | 22 | 0 | 20 | 0 | 0 | **parsed** |
| `Cell09` | Grove et al., Cell 2009 | 21 | 0 | 0 | 0 | 0 | **parsed** |
| `PNAS13` | Nakagawa et al., PNAS 2013 | 20 | 11 | 0 | 0 | 0 | **parsed** |
| `LAI20A` | Lai et al. 2020 | 20 | 3 | 0 | 1 | 3 | not usable |
| `Path10` | Campbell et al., PLoS Pathog 2010 | 19 | 0 | 0 | 0 | 0 | no admissible construct |
| `LIN14B` | Lindemose et al, Nucleic Acids Res. 2014 | 13 | 0 | 0 | 0 | 0 | **parsed** |
| `DEV12` | Busser et al., Development 2012 | 10 | 0 | 0 | 0 | 1 | **parsed** |
| `ROG18A` | Rogers et al., Mol Cell 2019 | 7 | 7 | 0 | 0 | 0 | **parsed** |
| `PO10` | Del Bianco et al., PLoS ONE 2010 | 6 | 0 | 0 | 0 | 0 | no admissible construct |
| `NBT06` | Berger et al., Nat Biotech 2006 | 5 | 0 | 0 | 0 | 0 | no admissible construct |
| `PNAS08` | De Silva et al., PNAS 2008 | 3 | 0 | 0 | 0 | 0 | **parsed** |
| `NAR10` | Alibes et al., NAR 2010 | 3 | 0 | 0 | 2 | 0 | not usable |
| `GD12` | Peterson et al., Genes Dev 2012 | 3 | 0 | 0 | 0 | 0 | no admissible construct |
| `LIU18B` | Liu et al., eLife 2018 | 2 | 0 | 0 | 0 | 0 | **parsed** |
| `LIU18A` | Liu et al., Cell 2018 | 2 | 0 | 0 | 0 | 0 | no admissible construct |
| `PP15` | Lehti-Shiu et al., PP 2015 | 2 | 0 | 0 | 0 | 0 | no admissible construct |
| `KUR17A` | Li et al., Nature 2017 | 2 | 0 | 0 | 0 | 0 | no admissible construct |
| `MBE14` | Cheatle Jarvela et al., Mol Biol Evol 2014 | 2 | 0 | 0 | 0 | 0 | **parsed** |
| `GD13` | Soruco et al., Genes Dev 2013 | 1 | 0 | 0 | 0 | 0 | no admissible construct |
| `NAR11` | De Masi et al., NAR 2011 | 1 | 0 | 0 | 0 | 0 | **parsed** |
| `PNAS12` | Busser et al., PNAS 2012 | 1 | 0 | 0 | 0 | 0 | no admissible construct |
| `CB11` | Helfer et al., Curr Biol 2011 | 1 | 0 | 0 | 0 | 0 | **parsed** |
| `MMB08` | Pompeani et al., Mol Microbiol 2008 | 1 | 0 | 0 | 0 | 0 | not usable |
| `GD09` | Lesch et al., Genes Dev 2009 | 1 | 0 | 0 | 0 | 0 | no admissible construct |
| `STI21B` | Stielow et al., Science Advances 2021 | 1 | 0 | 0 | 0 | 0 | no admissible construct |
| `CR09` | Scharer et al., Cancer Res 2009 | 1 | 0 | 0 | 0 | 0 | no admissible construct |
| `RAD13A` | Radke et al., PNAS 2013 | 1 | 0 | 0 | 0 | 0 | **parsed** |
| `MIZ19A` | Mizeracka et al. 2019 | 1 | 0 | 0 | 0 | 0 | not usable |

## Outcome legend

- **parsed** — contributed at least one admissible domain to the corpus.
- **no admissible construct** — retrieved and attempted, but every construct was rejected by
  the domain policy (`docs/DOMAIN_POLICY.md`), most often a C2H2 zinc-finger array or a
  multi-domain construct.
- **not usable** — could not be attempted at all: the server does not hold the advertised
  archive, the detail pages carry no sequence, or the 8-mer tables are truncated. Each case is
  recorded in `PROVENANCE.md`.

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
| `GB11` | Gordan 2011, bHLH — would most directly strengthen the smallest family (13 domains today) | 27 |
| `GR09` | Zhu et al. 2009, yeast TFs — adds a whole other clade of protein space | 89 |
| `MAR17A` | the clearest bZIP contribution, plus more forkhead | 61 |
| `SHO18A` | Shokri et al. 2019 | 45 |
| `ROG18A` | a second, independent forkhead panel — useful as a cross-source agreement check | 7 |

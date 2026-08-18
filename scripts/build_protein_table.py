#!/usr/bin/env python
"""Build data/interim/proteins/proteins.parquet from the admitted constructs.

    python scripts/build_protein_table.py

Re-runs the domain annotation over every source's constructs, keeps the admitted ones, and
fetches canonical UniProt sequences so embeddings can later be computed at the bare-domain,
padded-domain, construct or full-protein level. Fetched FASTA is cached under
data/external/uniprot/ so a re-run costs nothing.
"""

from __future__ import annotations

import time
import urllib.error
import urllib.request
from collections import namedtuple

import pandas as pd

from snp2prot import canonical, proteins, thresholds
from snp2prot.config import EXTERNAL_DIR, INTERIM_DIR, interim_table, raw_dir
from snp2prot.parsers import REGISTRY, _uniprobe, weirauch2014

#: Derived from the parser registry rather than hand-listed: a fixed list silently went
#: stale when the corpus grew from 9 sources to 18, so the protein table covered only half
#: of it while still reporting a plausible-looking total.
SOURCES = sorted(REGISTRY)
CACHE = EXTERNAL_DIR / "uniprot"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"

#: Accessions a deposit publishes that do not name the protein it assayed. Corrections only,
#: each with the evidence that settles it, never a guess (rule 1).
#:
#: `Q6P051` is a 305 aa TrEMBL entry, "SIX6 protein (Fragment)", from a cDNA clone; reviewed
#: `SIX6_HUMAN` is `O95475` at 246 aa. The stored `SIX6_REF` domain is a substring of **both**,
#: at offset 182 in the TrEMBL entry and 123 in the reviewed one — a 59-residue shift, which is
#: exactly the offset by which our `mut_positions` disagreed with the literature (`T22`).
ACCESSION_OVERRIDES = {"Q6P051": "O95475"}


#: One assayed construct, however its source happens to publish it.
Construct = namedtuple("Construct", "name sequence protein_id species")


def constructs_for(source: str) -> list[Construct]:
    """Every construct a source assayed, with whatever identity it publishes.

    Sources do not agree on where this lives. UniPROBE puts the clone insert on a per-gene
    HTML detail page; CIS-BP/Weirauch puts it in a supplementary spreadsheet keyed by plasmid.
    Deriving `SOURCES` from the parser registry is right, but it means a new source of a
    different shape lands here and must be handled rather than assumed away — before this
    existed, registering `weirauch2014` made this script raise `FileNotFoundError` looking for
    detail pages that were never going to exist.
    """
    if source == weirauch2014.SOURCE:
        clones = weirauch2014.load_constructs(raw_dir(source) / weirauch2014.CLONES_FILE)
        # Table S6 publishes no UniProt accession, so full-length mapping is unavailable for
        # this source until one is resolved from gene and species. Left empty, not guessed.
        return [Construct(c.gene, c.insert_aa, "", c.species) for c in clones.values()]

    pages = _uniprobe.load_details(raw_dir(source) / "details")
    return [
        Construct(name, seq, page.protein_id, page.species)
        for page in pages.values()
        for name, seq in page.inserts.items()
    ]


def fetch_uniprot(accession: str) -> tuple[str, bool, bool] | None:
    """`(sequence, reviewed, fragment)` for one accession, cached on disk.

    `reviewed` distinguishes Swiss-Prot from TrEMBL and `fragment` reads the header's own
    flag, because neither is a property a full-length sequence can be trusted without: 39 of
    the accessions the deposits publish are unreviewed and 10 say "(Fragment)".
    """
    if not accession or " " in accession:
        return None
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{accession}.fasta"
    if not path.exists():
        try:
            with urllib.request.urlopen(UNIPROT.format(accession), timeout=30) as fh:
                path.write_bytes(fh.read())
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            return None
        time.sleep(0.2)
    text = path.read_text()
    lines = text.splitlines()
    if not lines:
        return None
    header = lines[0]
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return (seq, header.startswith(">sp|"), "(Fragment)" in header) if seq else None


def locate(dbd: str, constructs: list[Construct]) -> tuple[Construct, int, int, str] | None:
    """Find the construct a stored domain came from, and where it sits inside it.

    A stored `dbd_seq` is usually a literal slice of its construct, and a substring search
    answers in microseconds. It is **not** a slice when the canonical sequence borrowed flank
    from a reference protein because the construct stopped short of the padding
    (`snp2prot.canonical`) — 111 domains of 1,334, which the substring test simply lost. Those
    are placed by alignment instead, and the domain's span is read off the alignment rather
    than assumed, so an indel between construct and canonical form cannot shift the window.
    """
    for c in constructs:
        if dbd in c.sequence:
            start = c.sequence.index(dbd) + 1
            return c, start, start + len(dbd) - 1, "substring"

    best: tuple[float, Construct, int, int] | None = None
    for c in constructs:
        # The stored domain plays the "construct" and the assayed construct the "reference":
        # `place` returns where the first argument sits inside the second. The envelope is the
        # whole domain, so every residue counts toward coverage.
        p = canonical.place(dbd, c.sequence, 1, len(dbd))
        if p is None or p.coverage < canonical.MIN_COVERAGE:
            continue
        if best is None or p.coverage > best[0]:
            best = (p.coverage, c, p.ref_start + 1, p.ref_end + 1)
    if best is None:
        return None
    return best[1], best[2], best[3], "aligned"


def main() -> None:
    """Build the protein table from the corpus, not by re-deriving domains.

    The parser is the authority on what `dbd_seq` is. BAR15A applies each cluster's
    *reference* domain boundaries to all of its alleles, so the cluster stays
    length-homogeneous and Hamming distance stays defined -- but a mutation can shift where
    the Pfam HMM aligns, so re-scanning a variant construct on its own yields a different
    envelope (PHOX2B: 15-71 for the reference, 18-71 for the variant). Re-deriving here
    therefore produced three domains the corpus does not contain, and missed three it does.

    So rows come from the interim tables, and each stored domain is located inside the
    construct it was cut from.
    """
    cfg = thresholds.load()["domain"]
    pad = int(cfg["padding_aa"])
    records: list[proteins.ProteinRecord] = []
    stats = {"rows": 0, "no_construct": 0, "uniprot_ok": 0, "mapped": 0,
             "substring": 0, "aligned": 0}  # fmt: skip

    for source in SOURCES:
        table = interim_table(source)
        if not table.exists():
            continue
        corpus = pd.read_parquet(
            table, columns=["dbd_seq", "dbd_family", "protein_id", "species"]
        ).drop_duplicates("dbd_seq")
        constructs = constructs_for(source)

        for _, row in corpus.iterrows():
            dbd = row["dbd_seq"]
            stats["rows"] += 1
            located = locate(dbd, constructs)
            if located is None:
                stats["no_construct"] += 1
                continue
            hit, start, end, how = located
            stats[how] += 1
            name, construct = hit.name, hit.sequence
            accession = ACCESSION_OVERRIDES.get(hit.protein_id, hit.protein_id)
            resolved = fetch_uniprot(accession) if accession else None
            full, reviewed, fragment = resolved if resolved else (None, None, None)
            if full:
                stats["uniprot_ok"] += 1
            span = proteins.locate_in_protein(construct, full or "", start, end)
            if span:
                stats["mapped"] += 1
            # The bare envelope is the stored slice minus the padding, clipped to it.
            bare_lo = min(start + pad, end)
            bare_hi = max(end - pad, start)
            records.append(
                proteins.ProteinRecord(
                    dbd_seq=dbd,
                    dbd_bare=construct[bare_lo - 1 : bare_hi],
                    pfam_family=row["dbd_family"],
                    pfam_start=bare_lo,
                    pfam_end=bare_hi,
                    dbd_start=start,
                    dbd_end=end,
                    construct_seq=construct,
                    full_seq=full,
                    dbd_start_protein=span[0] if span else None,
                    dbd_end_protein=span[1] if span else None,
                    protein_id=accession,
                    gene=name,
                    species=row["species"],
                    source_dataset=source,
                    placement=how,
                    uniprot_reviewed=reviewed,
                    uniprot_fragment=fragment,
                )
            )
        print(f"  {source}: {len(records)} records so far", flush=True)

    df = proteins.to_frame(records)
    out = INTERIM_DIR / "proteins"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "proteins.parquet"
    df.to_parquet(path, index=False)
    print(
        f"\n{stats['rows']} corpus domains, {stats['no_construct']} without a locatable "
        f"construct; UniProt resolved for {stats['uniprot_ok']}, domain located in the full "
        f"protein for {stats['mapped']}"
    )
    print(f"wrote {len(df):,} distinct domains -> {path}")


if __name__ == "__main__":
    main()

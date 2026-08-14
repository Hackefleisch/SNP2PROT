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

import pandas as pd

from snp2prot import proteins, thresholds
from snp2prot.config import EXTERNAL_DIR, INTERIM_DIR, interim_table, raw_dir
from snp2prot.parsers import REGISTRY, _uniprobe

#: Derived from the parser registry rather than hand-listed: a fixed list silently went
#: stale when the corpus grew from 9 sources to 18, so the protein table covered only half
#: of it while still reporting a plausible-looking total.
SOURCES = sorted(REGISTRY)
CACHE = EXTERNAL_DIR / "uniprot"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def fetch_uniprot(accession: str) -> str | None:
    """Canonical sequence for one accession, cached on disk. None if it does not resolve."""
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
    seq = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
    return seq or None


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
    stats = {"rows": 0, "no_construct": 0, "uniprot_ok": 0, "mapped": 0}

    for source in SOURCES:
        table = interim_table(source)
        if not table.exists():
            continue
        corpus = pd.read_parquet(
            table, columns=["dbd_seq", "dbd_family", "protein_id", "species"]
        ).drop_duplicates("dbd_seq")
        pages = _uniprobe.load_details(raw_dir(source) / "details")
        # construct sequence -> the page it came from, for locating a stored domain
        constructs = [
            (name, seq, page) for page in pages.values() for name, seq in page.inserts.items()
        ]

        for _, row in corpus.iterrows():
            dbd = row["dbd_seq"]
            stats["rows"] += 1
            hit = next(((n, c, pg) for n, c, pg in constructs if dbd in c), None)
            if hit is None:
                stats["no_construct"] += 1
                continue
            name, construct, page = hit
            start = construct.index(dbd) + 1
            end = start + len(dbd) - 1
            full = fetch_uniprot(page.protein_id)
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
                    protein_id=page.protein_id,
                    gene=name,
                    species=row["species"],
                    source_dataset=source,
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

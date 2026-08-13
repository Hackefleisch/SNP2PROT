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

from snp2prot import domains, proteins, thresholds
from snp2prot.config import EXTERNAL_DIR, INTERIM_DIR, raw_dir
from snp2prot.parsers import _uniprobe

SOURCES = ["BAR15A", "Cell08", "EMBO10", "PNAS13"]
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
    cfg = thresholds.load()["domain"]
    records: list[proteins.ProteinRecord] = []
    stats = {"admitted": 0, "rejected": 0, "uniprot_ok": 0, "mapped": 0}

    for source in SOURCES:
        pages = _uniprobe.load_details(raw_dir(source) / "details")
        constructs = {
            f"{g}|{allele}": seq
            for g, page in pages.items()
            for allele, seq in page.inserts.items()
        }
        hits = domains.scan(constructs)
        for key, construct in constructs.items():
            gene, _, allele = key.partition("|")
            call = domains.call_domain(construct, hits[key], cfg)
            if not call.ok:
                stats["rejected"] += 1
                continue
            stats["admitted"] += 1
            page = pages[gene]
            full = fetch_uniprot(page.protein_id)
            if full:
                stats["uniprot_ok"] += 1
            span = proteins.locate_in_protein(construct, full or "", call.start, call.end)
            if span:
                stats["mapped"] += 1
            bare_lo = min(h.start for h in call.hits)
            bare_hi = max(h.end for h in call.hits)
            records.append(
                proteins.ProteinRecord(
                    dbd_seq=call.sequence,
                    dbd_bare=construct[bare_lo - 1 : bare_hi],
                    pfam_family=call.family,
                    pfam_start=bare_lo,
                    pfam_end=bare_hi,
                    dbd_start=call.start,
                    dbd_end=call.end,
                    construct_seq=construct,
                    full_seq=full,
                    dbd_start_protein=span[0] if span else None,
                    dbd_end_protein=span[1] if span else None,
                    protein_id=page.protein_id,
                    gene=gene,
                    species=page.species,
                    source_dataset=source,
                )
            )
        print(f"  {source}: {len(records)} records so far")

    df = proteins.to_frame(records)
    out = INTERIM_DIR / "proteins"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "proteins.parquet"
    df.to_parquet(path, index=False)
    print(
        f"\n{stats['admitted']} admitted constructs, {stats['rejected']} rejected; "
        f"UniProt resolved for {stats['uniprot_ok']}, domain located in the full protein "
        f"for {stats['mapped']}"
    )
    print(f"wrote {len(df):,} distinct domains -> {path}")


if __name__ == "__main__":
    main()

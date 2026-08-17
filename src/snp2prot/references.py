"""Resolving a construct's reference protein sequence from whatever identifier it carries.

`snp2prot.canonical` needs a reference only where a construct fails to cover the canonical
window -- 357 of 1,477 admitted constructs. Every one of those carries an identifier, but not
the same kind: UniPROBE publishes a UniProt accession, while CIS-BP's Table S6 gives a `Gene ID`
drawn from whichever database that organism's genome project used. Across 131 species that is
Ensembl, FlyBase, Araport/TAIR, RefSeq, SGD, NCBI GeneID, dictyBase and a tail of
assembly-specific identifiers that resolve nowhere.

**Resolution is by cross-reference, never by gene name.** Searching UniProt for a gene symbol
plus an organism looks tempting and is wrong: `AT1G28420` resolves by cross-reference to a
1,705-residue protein and by gene-name search to a 317-residue one. Picking the wrong protein
is exactly the failure `PP15` was excluded to avoid (`docs/DECISIONS.md` §6), and a construct
that cannot be resolved is discarded rather than approximated.

Sequences are cached on disk, so a re-run costs nothing and the network is touched once.
"""

from __future__ import annotations

import re
import time
import urllib.error
import urllib.parse
import urllib.request

from snp2prot.config import EXTERNAL_DIR

CACHE = EXTERNAL_DIR / "uniprot"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"

#: Identifier shape -> the UniProt cross-reference database that indexes it. Ordered, because
#: the patterns are tested in sequence and the first match wins.
NAMESPACES: tuple[tuple[str, str], ...] = (
    (r"^ENS[A-Z]*[GTP]\d+", "ensembl"),
    (r"^FBgn\d+", "flybase"),
    (r"^AT\dG\d{5}", "araport"),
    (r"^(NP_|XP_|NM_|XM_)\w+", "refseq"),
    (r"^Y[A-P][LR]\d{3}[WC]", "sgd"),
    (r"^WBGene\d+", "wormbase"),
    (r"^DDB_G\d+", "dictybase"),
    (r"^\d+$", "geneid"),
)

#: A UniProt accession, which needs no cross-reference lookup at all. Some CIS-BP `Gene ID`
#: values are UniProt entry names in disguise (`Q0ZPQ8_NEMVE`), so the accession is taken from
#: the front rather than the whole string being rejected.
ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
)


def namespace(identifier: str) -> tuple[str | None, str]:
    """Which UniProt cross-reference indexes this identifier, and the id to look it up by.

    Returns `(None, identifier)` when the shape matches nothing known -- an assembly-specific
    identifier from a genome project UniProt does not cross-reference. Those are reported, not
    guessed at.
    """
    ident = (identifier or "").strip()
    if not ident:
        return None, ident
    if ACCESSION_RE.match(ident):
        return "uniprot", ident.split("_")[0]
    for pattern, db in NAMESPACES:
        if re.match(pattern, ident):
            # Araport ids arrive suffixed ("AT1G28420-TAIR-G"); the locus alone is the key.
            return db, ident.split("-")[0] if db == "araport" else ident
    return None, ident


def _cache_path(key: str):
    CACHE.mkdir(parents=True, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", key)
    return CACHE / f"{safe}.fasta"


def _get(url: str, timeout: int = 45) -> str | None:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            return fh.read().decode("utf8", errors="ignore")
    except (urllib.error.URLError, TimeoutError, OSError):
        return None


def _sequence(fasta: str | None) -> str | None:
    if not fasta or not fasta.startswith(">"):
        return None
    seq = "".join(line.strip() for line in fasta.splitlines()[1:] if not line.startswith(">"))
    return seq or None


def resolve(identifier: str, pause: float = 0.2) -> str | None:
    """The reference sequence for one identifier, or None if it does not resolve.

    Cached on disk under `data/external/uniprot/`, including negative results, so a rebuild
    does not re-ask the network about identifiers already known to be unresolvable.
    """
    db, ident = namespace(identifier)
    if db is None:
        return None
    path = _cache_path(f"{db}_{ident}")
    if path.exists():
        return _sequence(path.read_text()) if path.stat().st_size else None

    if db == "uniprot":
        text = _get(FASTA.format(ident))
    else:
        query = urllib.parse.urlencode(
            {"query": f"xref:{db}-{ident}", "format": "fasta", "size": "1"}
        )
        text = _get(f"{SEARCH}?{query}")
    time.sleep(pause)

    seq = _sequence(text)
    # An empty cache file records "asked, nothing there" so the miss is not repeated.
    path.write_text(text if seq else "")
    return seq

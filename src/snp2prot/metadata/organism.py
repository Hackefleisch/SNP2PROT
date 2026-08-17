"""Normalise the source organism to one convention across every deposit.

**The convention is the UniProt-style binomial**: `Genus species`, capitalised genus,
lower-case epithet, no strain, no abbreviation, no hybrid marker. Anything that is not an
organism is stored as the empty string rather than as a word that looks like one.

Two separate problems are handled here, both metadata-only — no label, cluster or domain
boundary depends on this column.

**Missing organisms** (`TODO.md` T5). UniPROBE renders each detail page from a template, and
on some pages the substitution did not run: the field arrives as the literal `$species`. Where
the construct carries an accession whose canonical FASTA is already cached under
`data/external/uniprot/`, the organism is read from that record's `OS=` field. Nothing here
fetches, so parsing stays offline and a cold cache degrades to an empty organism rather than
to a wrong one or to a network call inside a parser.

Guessing from the study's scope is deliberately not done. Every construct affected by the
template defect happens to come from one *C. elegans* panel, so guessing would be right today
and would be a rule that silently invents an organism the next time a mixed-species deposit
has the same defect.

**Inconsistent organisms** (`TODO.md` T5b). Deposits disagree about how to write a name that
is not in dispute. The corrections are enumerated rather than inferred: an alias table with a
reason per entry can be reviewed, whereas a rule that rewrites any string matching a pattern
will eventually rewrite a name that was correct. Every entry below was found by surveying the
distinct values actually present, not anticipated.

Deliberately **not** normalised:

- `Sarsia sp. Long Island Sound` — a genuinely unnamed species, not a formatting variant.
  Trimming it to a binomial would assert a species that has not been assigned.
- `Physcomitrella patens` — now *Physcomitrium patens* under current taxonomy. That is a
  revision of the name, not a disagreement about how to write it, and applying it here would
  put this module in the business of tracking taxonomy.
- `Acanthamoeba polyphaga mimivirus` — a virus, correctly named in three words.
"""

from __future__ import annotations

import re
from functools import cache

from snp2prot.config import EXTERNAL_DIR

#: Values that occupy the organism field without naming an organism, matched case-folded.
#: `$species` is UniPROBE's unsubstituted template variable; `Chimera` marks ROG18A's
#: engineered chimeras, which have no source organism by construction; `N/A` marks LIU18B's
#: reconstructed ancestors, likewise; `PBM CONSTRUCTS` is a stray heading from CIS-BP's Table
#: S6; `None Available` is how UniPROBE writes a missing value elsewhere on the same pages.
PLACEHOLDERS = frozenset(
    {
        "$species",
        "n/a",
        "na",
        "none",
        "none available",
        "unknown",
        "pbm constructs",
        "chimera",
    }
)

#: Deposit spelling -> the binomial, matched case-folded. One entry per variant actually
#: observed in the corpus, each with the reason it differs.
ALIASES = {
    # NAR11 abbreviates the genus on four constructs; every other deposit writes it out.
    "c. elegans": "Caenorhabditis elegans",
    # CIS-BP Table S6 transposes the 'h': the pea aphid is Acyrthosiphon, and the corpus
    # carries both spellings for what would otherwise be one organism.
    "acyrtosiphon pisum": "Acyrthosiphon pisum",
    # Hybrid marker. UniProt and NCBI both index the apple as Malus domestica.
    "malus x domestica": "Malus domestica",
}

#: Where `build_protein_table.py` caches canonical FASTA, one file per accession. A zero-byte
#: file is its negative cache: the accession did not resolve at UniProt.
CACHE = EXTERNAL_DIR / "uniprot"

#: `>sp|P46581|CND1_CAEEL ... OS=Caenorhabditis elegans OX=6239 GN=cnd-1 ...`
OS_FIELD = re.compile(r"\bOS=(.+?)\s+(?:OX|GN|PE|SV)=")


@cache
def from_cached_fasta(accession: str) -> str:
    """Organism named in a cached UniProt FASTA header. Empty if unavailable."""
    if not accession or accession.strip().casefold() in PLACEHOLDERS:
        return ""
    path = CACHE / f"{accession.strip()}.fasta"
    if not path.exists() or path.stat().st_size == 0:
        return ""
    header = path.read_text().split("\n", 1)[0]
    m = OS_FIELD.search(header)
    return m.group(1).strip() if m else ""


def resolve(raw: str, accession: str = "") -> str:
    """The organism to store for one construct, in the project's convention.

    A name that is already conventional is returned untouched. A known variant spelling is
    replaced by the binomial. A placeholder falls back to the cached UniProt record and, if
    that cannot answer, to the empty string — never to a guess.
    """
    value = re.sub(r"\s+", " ", (raw or "").strip())
    folded = value.casefold()
    if not value or folded in PLACEHOLDERS:
        return from_cached_fasta(accession)
    return ALIASES.get(folded, value)

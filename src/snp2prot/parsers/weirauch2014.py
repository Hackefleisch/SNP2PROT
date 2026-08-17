"""CIS-BP / Weirauch et al. 2014 — universal PBM over 1,032 cloned DNA-binding domains.

Not a UniPROBE accession, and not shaped like one. The measurements are a GEO SOFT family
record rather than a zip of per-gene folders, and the assayed sequence lives in a separate
supplementary spreadsheet rather than on a detail page. What the two sources share is the
assay, so everything assay-shaped comes from `_pbm` and only the file handling is here.

**The join is the whole reason this source is admissible.** Admission condition 1 needs the
sequence that was physically on the array, and Table S6 supplies it per plasmid — *"All
inserts were sequence verified in full"* (Weirauch et al. 2014, Extended Experimental
Procedures). Each GEO sample is titled `<plasmid>_<HK|ME>_8mer_<n>`, so the plasmid ID is the
key on both sides. It matches exhaustively: 1,032 plasmids in Table S6, 1,032 in GEO, none
unmatched either way.

Three properties of this source that a reader should not have to rediscover:

* **`GSE53348_RAW.tar` is not the data.** It holds probe-level intensities, 40,631 array spots
  per sample. The 8-mer E-scores exist only in the SOFT family record. See the acquisition
  note in `data/raw/weirauch2014/README.md`.
* **Two array designs, HK and ME, per plasmid.** They are replicates with different probe
  sequences, so they are binarized separately and combined — never averaged first.
* **Three construct architectures in one source**: 671 constructs carry 50 flanking residues,
  96 carry 15, and 265 carry none at all. A zero-flank construct yields a much shorter padded
  domain than a UniPROBE construct of the same protein, because our padding is clipped by the
  construct. That is recorded per construct, not silently absorbed.
"""

from __future__ import annotations

import gzip
import io
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from snp2prot import canonical, domains, references, schema, thresholds
from snp2prot.config import raw_dir
from snp2prot.metadata import organism
from snp2prot.parsers import _pbm

SOURCE = "weirauch2014"
SOFT_FILE = "GSE53348_family.soft.gz"
CLONES_FILE = "TabS6_DBD_clone_information.xlsx"
CLONES_SHEET = "Experimental constructs"
DBD_SOURCE = "pfam_hmmer_padded"

#: `!Sample_title = pTH1294_HK_8mer_593`. The plasmid ID is the first token and the array
#: design the second; both are needed, the first to attach a sequence and the second to know
#: these are two designs of one experiment rather than two experiments.
TITLE_RE = re.compile(r"^(?P<plasmid>pTH\d+)_(?P<array>HK|ME)_8mer_")


@dataclass(frozen=True)
class Construct:
    """One cloned insert, as Table S6 describes it."""

    plasmid: str
    gene: str
    species: str
    insert_aa: str
    #: Endogenous residues retained either side of the DBD: 50, 15, or 0.
    flanking_aa: int
    #: Table S6's `Gene ID`, in whichever namespace that organism's genome project used.
    #: Only consulted when the construct fails to cover the canonical window.
    gene_id: str = ""


@dataclass
class SkippedConstruct:
    plasmid: str
    reason: str


#: Table S6 columns this parser depends on. Named, never positional: reading a column by
#: index is what made the E-score wrong for three UniPROBE accessions (open item #36), and a
#: spreadsheet's column order is no more stable than a text file's.
CLONE_COLUMNS = ("Plasmid ID", "Gene Name", "Species", "Insert AA", "#Flanking AAs", "Gene ID")


def load_constructs(path: Path) -> dict[str, Construct]:
    """Table S6, keyed by plasmid ID."""
    df = pd.read_excel(path, sheet_name=CLONES_SHEET)
    absent = [c for c in CLONE_COLUMNS if c not in df.columns]
    if absent:
        raise ValueError(f"{path}: missing column(s) {absent}; have {list(df.columns)}")

    df = df[list(CLONE_COLUMNS)].dropna(subset=["Plasmid ID", "Insert AA"])
    out: dict[str, Construct] = {}
    for plasmid, gene, species, insert, flank, ident in df.itertuples(index=False, name=None):
        plasmid = str(plasmid).strip()
        if not plasmid.startswith("pTH"):
            continue
        out[plasmid] = Construct(
            plasmid=plasmid,
            gene=str(gene).strip(),
            species=str(species).strip(),
            insert_aa=str(insert).strip().upper(),
            flanking_aa=int(flank),
            gene_id="" if pd.isna(ident) else str(ident).strip(),
        )
    if not out:
        raise ValueError(f"{path}: no constructs with an insert sequence")
    return out


def iter_sample_tables(path: Path, wanted: set[str], require_complete: bool = True):
    """Stream the SOFT record, yielding `(plasmid, array, table)` for wanted plasmids.

    The file is ~1.7 GB compressed and holds 2,064 tables of 32,896 rows, so blocks for
    plasmids we are not going to use are consumed without being parsed. Each yielded table has
    `dna_seq` and `escore`; the E-score column is identified by `_pbm`, never by position.
    """
    with gzip.open(path, "rt", encoding="utf8", errors="ignore") as fh:
        plasmid = array = None
        block: list[str] | None = None
        for line in fh:
            if line.startswith("!Sample_title"):
                m = TITLE_RE.match(line.split("=", 1)[1].strip())
                plasmid, array = (m.group("plasmid"), m.group("array")) if m else (None, None)
            elif line.startswith("!sample_table_begin"):
                block = [] if (plasmid in wanted) else None
            elif line.startswith("!sample_table_end"):
                if block is not None and plasmid is not None:
                    yield (
                        plasmid,
                        array,
                        _read_block(block, f"{plasmid}_{array}", require_complete),
                    )
                block = None
            elif block is not None:
                block.append(line)


def _read_block(lines: list[str], label: str, require_complete: bool = True) -> pd.DataFrame:
    """One sample's data table. GEO writes a header, so it is used, then verified."""
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    escore = _pbm.escore_column(df, has_header=True, label=label)
    if require_complete:
        _pbm.check_complete(len(df), label)
    out = pd.DataFrame(
        {"dna_seq": df.iloc[:, 0].astype(str).str.strip().str.upper(), "escore": escore}
    )
    return out.dropna(subset=["escore"])


def _assign_wt_ids(groups: dict[str, dict]) -> dict[str, str]:
    """A cluster key per distinct domain sequence, unique and readable.

    Named for the gene, because that is what a reader recognises — but gene names are not
    unique here: 25 of them name a different sequence in a different organism (`FOXG1` in both
    human and medaka, `MATALPHA2` in three yeasts). A colliding name is therefore qualified
    with the plasmid, which is this source's own unique key, rather than silently fusing two
    orthologues into one cluster.
    """
    by_name: dict[str, list[str]] = defaultdict(list)
    for dbd, g in groups.items():
        by_name["/".join(sorted(g["genes"]))].append(dbd)
    out = {}
    for name, seqs in by_name.items():
        for dbd in seqs:
            plasmid = sorted(groups[dbd]["plasmids"])[0]
            out[dbd] = f"{SOURCE}:{name}" if len(seqs) == 1 else f"{SOURCE}:{name}_{plasmid}"
    return out


def parse_source(
    plasmids: list[str] | None = None,
    threshold_path: str | Path | None = None,
) -> tuple[pd.DataFrame, list[SkippedConstruct]]:
    """Parse the source. Returns the frame and the constructs that had to be skipped."""
    cfg = thresholds.for_assay("pbm", threshold_path)
    pos_cut, neg_cut = float(cfg["positive"]), float(cfg["negative"])
    if not cfg.get("per_experiment", False):
        raise ValueError("PBM sources binarize per experiment; thresholds.yaml disables it")

    root = raw_dir(SOURCE)
    constructs = load_constructs(root / CLONES_FILE)
    if plasmids:
        constructs = {k: v for k, v in constructs.items() if k in plasmids}

    # The construct that was on the array is what gets annotated, never the protein's own
    # domain annotation, which describes the full-length protein.
    domain_cfg = thresholds.load(threshold_path)["domain"]
    hits = domains.scan({p: c.insert_aa for p, c in constructs.items()})

    skipped: list[SkippedConstruct] = []
    calls: dict[str, domains.DomainCall] = {}
    canon: dict[str, str] = {}
    pad = canonical.padding(domain_cfg)
    for plasmid, c in constructs.items():
        call = domains.call_domain(c.insert_aa, hits[plasmid], domain_cfg)
        if not call.ok:
            skipped.append(SkippedConstruct(plasmid, f"rejected: {call.rejection} ({call.family})"))
            continue
        # The stored sequence is canonical, not the construct's own clipped window: 265 of
        # these constructs carry no flanking residues at all, so without this the same domain
        # is stored at one length here and another from UniPROBE.
        h = call.hits[0]
        ref = (
            None
            if canonical.covers(c.insert_aa, h.start, h.end, pad)
            else references.resolve(c.gene_id)
        )
        cc = canonical.canonicalise(c.insert_aa, h.start, h.end, reference=ref, pad=pad)
        if not cc.ok:
            skipped.append(SkippedConstruct(plasmid, f"not canonicalisable: {cc.rejection}"))
            continue
        calls[plasmid] = call
        canon[plasmid] = cc.sequence

    # Distinct plasmids resolving to one stored sequence are independent experiments on the
    # same protein sequence, so they are reconciled exactly like the HK/ME pair rather than
    # deduplicated -- keeping only the first would discard real data and hide disagreement.
    groups: dict[str, dict] = {}
    for plasmid, call in calls.items():
        c = constructs[plasmid]
        g = groups.setdefault(
            canon[plasmid], {"plasmids": [], "genes": set(), "call": call, "construct": c}
        )
        g["plasmids"].append(plasmid)
        g["genes"].add(c.gene)

    # Constructs sharing a canonical sequence but differing OUTSIDE it carry a signal the
    # model can never see. Bin the group rather than reconcile them as replicates.
    for dbd in [d for d, g in groups.items() if len(g["plasmids"]) > 1]:
        members = {p: constructs[p].insert_aa for p in groups[dbd]["plasmids"]}
        if canonical.conflicting_constructs(members):
            skipped.append(
                SkippedConstruct(
                    "/".join(sorted(members)),
                    "constructs share a canonical domain but differ outside it; discarded",
                )
            )
            del groups[dbd]

    wt_ids = _assign_wt_ids(groups)
    for g in groups.values():
        if len(g["plasmids"]) > 1:
            skipped.append(
                SkippedConstruct(
                    "/".join(sorted(g["plasmids"])),
                    f"{len(g['plasmids'])} plasmids share one DBD sequence; their experiments "
                    f"were reconciled as replicates",
                )
            )

    wanted = set(calls)
    tables: dict[str, list[tuple[str, pd.DataFrame]]] = defaultdict(list)
    for plasmid, array, table in iter_sample_tables(root / SOFT_FILE, wanted):
        tables[plasmid].append((f"{plasmid}_{array}", table))
    missing = wanted - set(tables)
    for plasmid in sorted(missing):
        skipped.append(SkippedConstruct(plasmid, "admitted, but no 8-mer table in the SOFT record"))

    # 884 domains x 32,896 8-mers is ~29 M rows, more than the rest of the corpus put
    # together, so memory is managed rather than assumed: each group's replicate tables are
    # released as soon as they are combined, and each frame is coerced to its final dtypes
    # (Int8 labels, Int32 lengths) before being appended rather than after concatenation.
    # Coercing per frame is equivalent -- every derived column is row-local.
    frames: list[pd.DataFrame] = []
    for dbd_seq, g in groups.items():
        replicates = [t for p in sorted(g["plasmids"]) for t in tables.get(p, [])]
        if not replicates:
            continue
        dna_seq, label, raw_score = _pbm.combine_replicates(replicates, pos_cut, neg_cut)
        for plasmid in g["plasmids"]:
            tables.pop(plasmid, None)
        c = g["construct"]
        frames.append(
            schema.coerce(
                _pbm.build_frame(
                    dna_seq=dna_seq,
                    label=label,
                    raw_score=raw_score,
                    dbd_seq=dbd_seq,
                    dbd_family=g["call"].family,
                    dbd_source=DBD_SOURCE,
                    wt_id=wt_ids[dbd_seq],
                    # No engineered variant series in this source: every construct is a distinct
                    # protein, so each is its own reference. Clustering near-identical orthologues
                    # is a separate decision (TODO.md T12), not something a parser may do.
                    n_mut_from_wt=0,
                    mut_positions="",
                    # Table S6 gives no UniProt accession; the protein table resolves one from the
                    # gene and species, and inventing one here would be a guess.
                    protein_id="",
                    species=organism.resolve(c.species),
                    pos_cut=pos_cut,
                    neg_cut=neg_cut,
                    source_dataset=SOURCE,
                    source_file=f"{SOURCE}/{SOFT_FILE}",
                )
            )
        )

    if not frames:
        raise ValueError(f"{SOURCE}: no experiments parsed")
    return pd.concat(frames, ignore_index=True), skipped


def parse() -> pd.DataFrame:
    """Parse the CIS-BP / Weirauch 2014 PBM set (Weirauch et al., Cell 2014; GEO GSE53348)."""
    return parse_source()[0]

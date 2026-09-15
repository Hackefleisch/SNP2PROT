"""Which TFs enter the GHT arm, and what protein sequence stands for each.

**The panel is decided by the admission policy, not by a family label.** `GHT_PLAN.md` §11 is
explicit about this and says why: the plan's own §4 counted the panel by string-matching the
`DBD(s)` column of the Codebook metadata — 139 benchmark TFs minus 92 labelled `C2H2 ZF` gives
47 — and that count is not what `snp2prot.domains` produces. It cannot be, because the label
records what the curators believe the protein contains while the policy records what a Pfam scan
finds *and* whether the found region is solely responsible, contiguous, and resolvable. The two
answers differ, and the policy's is the one every other source in this project was admitted by.

Measured: **39 of 139 admitted**, of which 6 carry a single `zf-C2H2` hit and are dropped by
`D10`, leaving **33**. The gap to 47 is not C2H2 — it is 11 constructs whose only Pfam hits are
families absent from `data/external/pfam/dbd_families.txt` (`FLYWCH`, `Myb_DNA-bind_4`,
`CGGBP1_N`, `Homeobox_KN`, …), 2 with no Pfam hit at all, and 2 genuinely two-domain constructs
(`MGA` T-box + HLH, `PAX7` PAX + homeodomain). See `reports/ght_panel.md`.

## Which construct stands for a TF

A TF can have several approved GHT-SELEX experiments run from different plasmids — full-length
(`FL`) and isolated-domain (`DBD`, `DBD.1`, …) constructs both appear, and the benchmark's
positives are pooled across replicates, so no single construct owns a peak file. The choice here
is **the longest full-length construct, falling back to the longest construct of any kind**, and
the reason is what `dbd_seq` is: the Pfam envelope *plus 10 residues each side*
(`docs/DOMAIN_POLICY.md`). A DBD-only construct can end inside that padding, so the same domain
would be stored shorter — a difference with nothing to do with the protein. The full-length
construct gives the envelope its padding, and which construct was used is recorded on every row.

**This is the construct-mismatch risk `GHT_PLAN.md` §12 names, made explicit rather than hidden:**
the model is fed the isolated padded domain while many of the peak files came from a full-length
protein. Nothing here can fix that; the column `insert_composition` is what lets a reader see it.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from snp2prot import domains
from snp2prot.ght import config

#: `GHT_PLAN.md` D10, decided 2026-09-15: a C2H2 zinc-finger array is 2-6 separate ~23-residue
#: folds on flexible linkers and fails admission condition 2. Six panel constructs pass the
#: policy's own array test only because HMMER's gathering threshold caught a single finger of an
#: 8-, 14- or 25-finger protein — the declared `DBD count` says otherwise, and the fingers the
#: scan missed still bind. Dropping them by family is the honest reading of D10.
DROPPED_FAMILY = "zf-C2H2"
DROPPED_BY_D10 = domains.Rejection("c2h2_array_d10")

APPROVED = "Approved"


@dataclass(frozen=True)
class Constructs:
    """The approved GHT-SELEX constructs of every benchmark TF, joined to their sequences."""

    frame: pd.DataFrame  # one row per (TF, plasmid) with an amino-acid sequence
    panel: list[str]  # the 139 benchmark TFs, in file order

    @property
    def n_with_sequence(self) -> int:
        return int(self.frame.tf.nunique())


def read_panel(path=None) -> list[str]:
    """`TFs_CHS+AFS.txt` — the benchmark's own TF list, in file order."""
    p = path or config.PANEL_FILE
    return [line.strip() for line in open(p) if line.strip()]


def approved_constructs(panel: list[str] | None = None) -> Constructs:
    """Every approved GHT-SELEX experiment of every panel TF, with its construct sequence.

    Three joins, all on the metadata's own keys: experiment -> plasmid -> insert. The insert
    join carries `TF` as well as `Insert Name` so that a name reused across TFs cannot silently
    attach one protein's sequence to another.
    """
    names = panel if panel is not None else read_panel()
    experiments = pd.read_excel(config.GHT_METADATA, "GHT_Experiment_Info")
    plasmids = pd.read_excel(config.TF_METADATA, "Plasmids")
    inserts = pd.read_excel(config.TF_METADATA, "Inserts")

    approved = experiments[experiments["Expert curation status"] == APPROVED]
    joined = (
        approved[approved.TF.isin(names)]
        .merge(plasmids[["Plasmid ID", "Insert Name"]], on="Plasmid ID", how="left")
        .merge(
            inserts[["Insert Name", "TF", "Insert composition", "Amino acid sequence"]],
            on=["Insert Name", "TF"],
            how="left",
        )
        .dropna(subset=["Amino acid sequence"])
    )
    frame = pd.DataFrame(
        {
            "tf": joined.TF.to_numpy(),
            "plasmid_id": joined["Plasmid ID"].to_numpy(),
            "insert_name": joined["Insert Name"].to_numpy(),
            "insert_composition": joined["Insert composition"].to_numpy(),
            "domain_label": joined["Domain"].to_numpy(),
            "construct_seq": joined["Amino acid sequence"].astype(str).str.strip().str.upper(),
        }
    ).drop_duplicates(["tf", "plasmid_id"])
    return Constructs(frame.reset_index(drop=True), names)


def choose_construct(constructs: Constructs) -> pd.DataFrame:
    """One row per TF: the longest full-length construct, else the longest of any kind.

    Full-length first because `dbd_seq` is the padded envelope and a DBD-only construct can end
    inside the padding — see the module docstring. Length breaks ties within a composition so the
    pick does not depend on row order.
    """
    frame = constructs.frame.assign(
        is_fl=lambda d: d.insert_composition.eq("FL"),
        length=lambda d: d.construct_seq.str.len(),
    )
    return (
        frame.sort_values(
            ["tf", "is_fl", "length", "plasmid_id"], ascending=[True, False, False, True]
        )
        .drop_duplicates("tf")
        .drop(columns=["is_fl"])
        .reset_index(drop=True)
    )


def declared_dbd_counts() -> pd.DataFrame:
    """The curators' own `DBD(s)` label and `DBD count`, carried through for the report.

    **Not used by the admission policy** — it is the thing §11 forbids deciding on. It is kept
    so a reader can see where the policy and the label disagree, which is the whole finding.
    """
    tfs = pd.read_excel(config.TF_METADATA, "TFs")
    return pd.DataFrame(
        {
            "tf": tfs.TF.to_numpy(),
            "declared_dbd": tfs["DBD(s)"].astype(str).to_numpy(),
            "declared_dbd_count": pd.to_numeric(tfs["DBD count"], errors="coerce").to_numpy(),
        }
    )


def build(
    panel: list[str] | None = None, drop_c2h2: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """`(admitted, rejected)` — the panel, through the same admission path as every other source.

    `drop_c2h2` applies `D10` on top of the policy. It is a separate flag from the policy's own
    `repeat_array` rejection because the two catch different things: the policy catches a
    construct where Pfam found more than one finger, `D10` catches the six where it found one
    finger of an array the curators counted 2 to 25 fingers in.
    """
    constructs = approved_constructs(panel)
    chosen = choose_construct(constructs)
    sequences = dict(zip(chosen.tf, chosen.construct_seq, strict=True))
    hits = domains.scan(sequences)
    declared = declared_dbd_counts().set_index("tf")

    rows = []
    for record in chosen.itertuples():
        call = domains.call_domain(record.construct_seq, hits[record.tf])
        rejection = call.rejection
        if call.ok and drop_c2h2 and call.family == DROPPED_FAMILY:
            rejection = DROPPED_BY_D10
        info = declared.loc[record.tf] if record.tf in declared.index else None
        rows.append(
            {
                "tf": record.tf,
                "plasmid_id": record.plasmid_id,
                "insert_name": record.insert_name,
                "insert_composition": record.insert_composition,
                "construct_seq": record.construct_seq,
                "construct_len": len(record.construct_seq),
                "dbd_seq": call.sequence or "",
                "dbd_family": call.family,
                "dbd_start": call.start,
                "dbd_end": call.end,
                "n_pfam_hits": len(call.hits),
                "pfam_hits": ";".join(f"{h.family}:{h.start}-{h.end}" for h in call.hits),
                "other_hits": ";".join(f"{h.family}:{h.start}-{h.end}" for h in call.other),
                "declared_dbd": "" if info is None else str(info.declared_dbd),
                "declared_dbd_count": (
                    float("nan") if info is None else float(info.declared_dbd_count)
                ),
                "rejection": str(rejection or ""),
            }
        )
    frame = pd.DataFrame(rows)
    admitted = frame[frame.rejection == ""].drop(columns="rejection").reset_index(drop=True)
    rejected = frame[frame.rejection != ""].reset_index(drop=True)
    return admitted, rejected


def load(path=None) -> pd.DataFrame:
    """The built panel. Sorted by `tf`, which is the order every downstream array is indexed by."""
    p = path or config.PANEL_TABLE
    if not p.exists():
        raise FileNotFoundError(f"no GHT panel at {p}\nrun scripts/build_ght_panel.py")
    return pd.read_parquet(p)

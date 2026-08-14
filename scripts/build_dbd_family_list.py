#!/usr/bin/env python
"""Generate the list of DNA-binding Pfam families used by the admission policy.

    python scripts/build_dbd_family_list.py

The whitelist decides which Pfam families take part in the admission policy, so a family
missing from it is rejected as `no_domain` — not because nothing was found, but because we
were not looking. That failure has now happened twice (open items #34 and T15), and both
times the cause was the same: the whitelist was derived from ONE database's vocabulary.

So it is built from every source database's own curation of what does the DNA binding:

  * **UniPROBE** publishes a `Domain` label per experiment on its browse page.
  * **CIS-BP / Weirauch 2014** publishes a `Pfam ID` per construct in Table S6.

Both are free text and both are stale in places, so **no label is trusted verbatim**. Each
token is resolved against Pfam-A's own `NAME` fields; where a label predates a Pfam rename it
is resolved through `MANUAL`, whose entries carry the accession they were checked against.

Tokens that resolve to nothing are printed, not silently dropped.
"""

from __future__ import annotations

import re

import pandas as pd

from snp2prot.config import EXTERNAL_DIR, PROJECT_ROOT, RAW_DIR

BROWSE = EXTERNAL_DIR / "uniprobe" / "browse.html"
CISBP = RAW_DIR / "weirauch2014" / "TabS6_DBD_clone_information.xlsx"
CISBP_SHEET = "Experimental constructs"
PFAM = EXTERNAL_DIR / "pfam" / "Pfam-A.hmm"
OUT = EXTERNAL_DIR / "pfam" / "dbd_families.txt"

#: Labels that name a fold or a concept rather than a Pfam family, mapped by hand. Each is a
#: judgement call and is listed explicitly so it can be audited.
MANUAL = {
    # UniPROBE labels predate several Pfam renames. Each was resolved by looking the Pfam
    # ACCESSION up in Pfam-A and reading its current NAME, not by guessing at the string.
    "homeobox": "Homeodomain",  # PF00046
    "fork_head": "Forkhead",  # PF00250
    "zn2cys6": "Zn_clus",  # PF00172
    "mads": "SRF-TF",  # PF00319
    "myb": "Myb_DNA-binding",  # PF00249
    "rfx": "RFX_DNA_binding",  # PF02257
    "e2f_tdp": "WHD_E2F_TDP",  # PF02319
    "nr": "zf-C4",  # PF00105
    "bright": "ARID",  # PF01388
    "ap-2": "TF_AP-2",  # PF03299 (mammalian TFAP2; plant "AP2" is separate)
    "pbx": "PBC",
    "zinc finger c2h2": "zf-C2H2",
    "c2h2 znf": "zf-C2H2",
    "gata type znf": "GATA",
    "zf-nhr/gata-type": "GATA",
    "hmg box": "HMG_box",
    "myb/sant-like": "Myb_DNA-binding",
    "sant": "Myb_DNA-binding",
    "tbox": "T-box",
    "znf_c4": "zf-C4",
    "nr_c4": "zf-C4",
    "bhlh": "HLH",
    "myc-type": "HLH",
    "basic": "HLH",
    "bzip": "bZIP_1",
    "brlz": "bZIP_1",
    "wh": "Forkhead",
    "hth apses-type": "KilA-N",
    # CIS-BP's Table S6 labels are older still. Each of these was resolved by scanning the
    # constructs CIS-BP gave the label to against the full Pfam-A library and reading off
    # which family actually hits them -- the sequence decides, not the string. The method
    # reproduces the three hand mappings above (Homeobox, Fork_head, E2F_TDP) exactly, which
    # is why it is trusted for the rest. Accession and description confirm each one.
    "cxc": "TCR",  # PF03638 "Tesmin/TSO1-like CXC domain"
    "duf260": "LOB",  # PF03195 "Lateral organ boundaries (LOB) domain"
    "duf573": "GeBP-like_DBD",  # PF04504 "...DBD domain"
    "ein3": "EIN3_DNA-bd",  # PF04873 -- the DNA-binding domain, NOT EIN3_N (PF27048)
    "rhd": "RHD_DNA_bind",  # PF00554 -- DNA-binding, NOT RHD_dimer (PF16179)
    "zf-dof": "Zn_ribbon_Dof",  # PF02701 "Dof domain, zinc finger"
}
#: Labels that are not DNA-binding domains at all, or are too vague to resolve.
IGNORE = {
    "",
    "distant similarity to homeodomain",
    "hth",
    "dwa",
    "nr lbd",
    "tig",
    "sam_pnt",
    "lim",
    "pwwp",
    "mtf2_c",
    "stb3",
    "vhr1",
    "gcr1_c",
    "ste",
}


def uniprobe_labels() -> list[str]:
    """UniPROBE's `Domain` column, one free-text label per experiment."""
    if not BROWSE.exists():
        print(f"note: UniPROBE browse page not cached ({BROWSE}); skipping that source")
        return []
    html = BROWSE.read_text(errors="ignore")
    return re.findall(r'value="[^"/]+/[^"/]+/[^"]*".*?</a></td>\s*<td>([^<]*)</td>', html, re.S)


def cisbp_labels() -> list[str]:
    """CIS-BP's `Pfam ID` column, one label per experimental construct.

    Multi-domain constructs are written as a comma-joined list ("Homeobox,Pou"), which the
    caller already splits -- the admission policy rejects those constructs anyway, but each
    of their domains is still a DNA-binding family and belongs on the whitelist.
    """
    if not CISBP.exists():
        print(f"note: CIS-BP Table S6 not present ({CISBP}); skipping that source")
        return []
    df = pd.read_excel(CISBP, sheet_name=CISBP_SHEET)
    return df["Pfam ID"].dropna().astype(str).tolist()


def pfam_names() -> set[str]:
    if not PFAM.exists():
        raise SystemExit(f"Pfam-A not found: {PFAM}")
    return {
        line.split(None, 1)[1].strip()
        for line in PFAM.read_text(errors="ignore").splitlines()
        if line.startswith("NAME")
    }


def main() -> None:
    sources = {"uniprobe": uniprobe_labels(), "cisbp": cisbp_labels()}
    if not any(sources.values()):
        raise SystemExit("no label source available; nothing to generate")

    names = pfam_names()
    lower = {n.lower(): n for n in names}
    resolved: dict[str, set[str]] = {}
    origin: dict[str, set[str]] = {}
    unresolved: dict[str, set[str]] = {}
    for source, rows in sources.items():
        for label in rows:
            for token in (t.strip() for t in label.split(",")):
                key = token.lower()
                if key in IGNORE:
                    continue
                family = MANUAL.get(key) or lower.get(key)
                if family is None:
                    unresolved.setdefault(source, set()).add(token)
                    continue
                resolved.setdefault(family, set()).add(token)
                origin.setdefault(family, set()).add(source)

    missing = sorted(f for f in resolved if f not in names)
    previous = set()
    if OUT.exists():
        previous = {
            ln.split("#")[0].strip()
            for ln in OUT.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        }

    lines = [
        "# DNA-binding Pfam families used by the domain admission policy.",
        "# Generated by scripts/build_dbd_family_list.py from the source databases' own",
        "# curation of what binds DNA -- UniPROBE's Domain labels and CIS-BP's Pfam ID column",
        "# -- each resolved against Pfam-A NAME fields. Do not hand-edit; re-run the script.",
        "#",
        "# The trailing comment records which database vouched for the family. A family listed",
        "# by only one of them is not more doubtful: they curate different organisms.",
        "",
    ]
    width = max(len(f) for f in resolved)
    lines += [f"{f:<{width}}  # {', '.join(sorted(origin[f]))}" for f in sorted(resolved)]
    OUT.write_text("\n".join(lines) + "\n")

    added = sorted(set(resolved) - previous)
    dropped = sorted(previous - set(resolved))
    print(f"{len(resolved)} DNA-binding families -> {OUT.relative_to(PROJECT_ROOT)}")
    for source, rows in sources.items():
        n = sum(1 for f, o in origin.items() if source in o)
        print(f"  {source:<9} {len(rows):>5} labels -> {n} families")
    if added:
        print(f"\nADDED {len(added)}: {added}")
    if dropped:
        print(f"REMOVED {len(dropped)}: {dropped}   <-- check this is intended")
    if missing:
        print(f"WARNING: {len(missing)} mapped names absent from Pfam-A: {missing}")
    for source, toks in unresolved.items():
        print(f"unresolved in {source} ({len(toks)}), not included: {sorted(toks)}")


if __name__ == "__main__":
    main()

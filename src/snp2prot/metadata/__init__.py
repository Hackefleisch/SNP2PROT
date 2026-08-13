"""Shared TF metadata used by every parser: DBD boundaries, families, canonical sequences.

Sources: CIS-BP (family + curated DBD), Pfam/InterPro + HMMER (HMM-derived DBD boundaries),
UniProt (canonical sequence, variant-notation resolution). Whichever route produced a given
domain is recorded per row in `dbd_source` — curated and HMM-derived boundaries are not
interchangeable and must stay distinguishable downstream.
"""

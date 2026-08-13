"""Featurization of the merged table into model inputs. Phase 6+.

Deliberately separate from `snp2prot.parsers`: padding and windowing are modeling
decisions and must never touch the stored table (brief §3, §5a). Two dataset variants
get built here and compared, not chosen between:

- `padded20`    — every site centred in a fixed 20 bp window, padded with a masked neutral
                  token, with `dna_context` passed through so the model can tell absent
                  evidence from an observed base.
- `common_core` — restricted to an 8-10 bp core shared by all sources, so `dna_len` cannot
                  leak assay identity at all.
"""

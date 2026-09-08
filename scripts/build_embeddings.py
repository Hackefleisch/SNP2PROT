#!/usr/bin/env python
"""Embed every domain with a protein language model and cache the pooled vectors.

    python scripts/build_embeddings.py --arm A1 [--batch-tokens 8192] [--device cuda]
    python scripts/build_embeddings.py --arm A1 --per-residue      # the un-pooled cache (T36)

Phase 7, build step 1 (`docs/ML_PLAN.md` §9). A one-off inference over the corpus's 1,338
canonical domains, mean-pooled to one vector per domain and written to
`data/processed/embeddings/<arm>.npz`. The language model never runs again — the training loop
reads this file (`ML_PLAN.md` §3.1).

**Weights are cached under `data/external/models/`**, not in the user's home directory, so that
a downloaded artifact the build depends on sits with the Pfam HMMs and the UniProt cache and can
carry a `PROVENANCE.md` row like everything else.

Order and identity: rows come out in `corpus.domains()` order and the domain list is written
into the file, so an embedding table built either side of a rebuild cannot be paired with a
mismatched matrix silently.

**`--per-residue` writes the un-pooled cache instead** — `<arm>_residues.npz`, every residue's
vector rather than their mean, about 600 MB. It is the same forward pass; the pooled run simply
averages the result and throws the rest away. It exists because that averaging was measured to
destroy the single-residue signal the project's central claim depends on
(`docs/ML_RESULTS.md` §9.3), and `scripts/check_unpooling.py` reads this file to find out whether
the signal survives before the mean is taken.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from snp2prot import corpus, embeddings
from snp2prot.config import (
    ESM_DBP_CHECKPOINT,
    MODEL_DIR,
    embedding_table,
    residue_embedding_table,
)


def load_model(arm: str, device: str):
    """The checkpoint for one arm, in eval mode on `device`, with its alphabet.

    Both sequence arms are the **same architecture** — ESM-2 650M, 33 layers, 1280-d — and
    differ only in the weights, which is what makes the `A1` -> `A4` delta isolate the
    pretraining corpus and nothing else (`ML_PLAN.md` §4.2). ESM-DBP publishes a bare
    `state_dict` rather than a loader, and it matches that architecture exactly: 572 tensors,
    none missing, none extra, no shape disagreement. So `A4` is `A1`'s architecture with the
    published weights loaded into it, and any future mismatch is an error rather than a
    `strict=False` shrug.
    """
    import esm
    import torch

    if arm not in embeddings.ARMS:
        raise SystemExit(f"unknown arm {arm!r}; have {sorted(embeddings.ARMS)}")

    # fair-esm downloads through torch.hub; point that at the project rather than ~/.cache.
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    torch.hub.set_dir(str(MODEL_DIR))
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()

    if arm == "A4":
        if not ESM_DBP_CHECKPOINT.exists():
            raise SystemExit(
                f"no ESM-DBP checkpoint at {ESM_DBP_CHECKPOINT}\nsee data/external/models/README.md"
            )
        state = torch.load(ESM_DBP_CHECKPOINT, map_location="cpu", weights_only=False)
        # Saved from a DataParallel wrapper, so every key carries a `module.` prefix.
        model.load_state_dict({strip_prefix(k): v for k, v in state.items()}, strict=True)

    return model.eval().to(device), alphabet, model.num_layers


def strip_prefix(key: str, prefix: str = "module.") -> str:
    return key[len(prefix) :] if key.startswith(prefix) else key


def pool(representation, lengths) -> np.ndarray:
    """Mean over each sequence's residues, dropping the BOS and EOS tokens.

    Those two carry sequence-level summary information of a kind every domain has equally, and
    averaging them in would dilute the per-residue signal pooling exists to capture
    (`snp2prot.embeddings`).
    """
    out = np.empty((representation.shape[0], representation.shape[2]), dtype=np.float32)
    for i, length in enumerate(lengths):
        out[i] = representation[i, 1 : length + 1].mean(0).float().cpu().numpy()
    return out


def residues_of(representation, lengths) -> list[np.ndarray]:
    """Every residue's vector, dropping BOS and EOS — `pool` without the mean.

    The same slice `pool` takes, so the two are the same numbers and
    `ResidueEmbeddings.pooled()` reproduces `DomainEmbeddings.vectors` exactly.
    """
    return [
        representation[i, 1 : length + 1].float().cpu().numpy().astype(np.float32)
        for i, length in enumerate(lengths)
    ]


def batches(sequences: list[str], batch_tokens: int):
    """Group sequences so that each batch is about `batch_tokens` residues, longest first.

    Longest-first keeps the padding waste in one place: a batch of similar lengths pads very
    little, and the largest batch is the one that decides peak memory, so it is met first
    rather than after an hour of work.
    """
    order = sorted(range(len(sequences)), key=lambda i: -len(sequences[i]))
    batch: list[int] = []
    for i in order:
        longest = max([len(sequences[j]) for j in batch] + [len(sequences[i])])
        if batch and longest * (len(batch) + 1) > batch_tokens:
            yield batch
            batch = []
        batch.append(i)
    if batch:
        yield batch


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="A1", choices=sorted(embeddings.ARMS))
    ap.add_argument("--batch-tokens", type=int, default=8192)
    ap.add_argument("--device", default=None, help="default: cuda when available")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--per-residue",
        action="store_true",
        help="write the un-pooled cache <arm>_residues.npz instead of the pooled table (T36)",
    )
    args = ap.parse_args()

    import torch

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    domains = corpus.domains()
    sequences = [str(s) for s in domains.dbd_seq]
    print(f"{len(sequences)} domains, {min(map(len, sequences))}-{max(map(len, sequences))} aa")
    print(f"model {embeddings.ARMS[args.arm]} on {device}")

    model, alphabet, layer = load_model(args.arm, device)
    converter = alphabet.get_batch_converter()

    vectors = np.empty((len(sequences), model.embed_dim), dtype=np.float32)
    # Batches come out longest-first, so the per-residue blocks are collected into a slot per
    # domain and concatenated in corpus order at the end rather than appended as they arrive.
    per_residue: list[np.ndarray | None] = [None] * len(sequences)
    started = time.time()
    done = 0
    with torch.no_grad():
        for group in batches(sequences, args.batch_tokens):
            _, _, tokens = converter([(str(i), sequences[i]) for i in group])
            out = model(tokens.to(device), repr_layers=[layer])["representations"][layer]
            lengths = [len(sequences[i]) for i in group]
            if args.per_residue:
                for i, block in zip(group, residues_of(out, lengths), strict=True):
                    per_residue[i] = block
            else:
                vectors[group] = pool(out, lengths)
            done += len(group)
            if done % 200 < len(group):
                print(f"  {done}/{len(sequences)}  ({time.time() - started:.0f}s)", flush=True)

    if args.per_residue:
        blocks = [b for b in per_residue if b is not None]
        if len(blocks) != len(sequences):
            raise SystemExit("some domains produced no representation; refusing to write")
        offsets = np.zeros(len(blocks) + 1, dtype=np.int64)
        np.cumsum([len(b) for b in blocks], out=offsets[1:])
        table = embeddings.ResidueEmbeddings(
            arm=args.arm,
            model=embeddings.ARMS[args.arm],
            domains=np.asarray(sequences, dtype=np.str_),
            offsets=offsets,
            vectors=np.concatenate(blocks, axis=0),
        )
        path = table.save(args.out or residue_embedding_table(args.arm))
        print(f"{table.width}-d over {offsets[-1]:,} residues in {len(blocks)} domains")
        print(f"wrote {path} ({path.stat().st_size / 1e6:.0f} MB) in {time.time() - started:.0f}s")
        return

    table = embeddings.DomainEmbeddings(
        arm=args.arm,
        model=embeddings.ARMS[args.arm],
        domains=np.asarray(sequences, dtype=np.str_),
        vectors=vectors,
    )
    path = table.save(args.out or embedding_table(args.arm))
    norms = np.linalg.norm(vectors, axis=1)
    print(f"{table.width}-d, norms {norms.min():.2f}-{norms.max():.2f}")
    print(f"wrote {path} ({path.stat().st_size / 1e6:.0f} MB) in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()

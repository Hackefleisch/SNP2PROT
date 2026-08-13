#!/usr/bin/env python
"""CLI entry point for training. Keeps all logic in `src/snp2prot`."""

from __future__ import annotations

import argparse

from snp2prot.config import load_config
from snp2prot.training.train import train
from snp2prot.utils.seed import set_seed


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a SNP2PROT model")
    parser.add_argument("--config", default="configs/default.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    set_seed(cfg.seed)
    train(cfg)


if __name__ == "__main__":
    main()

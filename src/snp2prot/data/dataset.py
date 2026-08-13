"""Dataset definitions. TODO."""

from __future__ import annotations


class SNPDataset:
    """Placeholder dataset. TODO: define inputs, targets, and splits."""

    def __init__(self, split: str = "train") -> None:
        self.split = split

    def __len__(self) -> int:
        raise NotImplementedError

    def __getitem__(self, idx: int):
        raise NotImplementedError

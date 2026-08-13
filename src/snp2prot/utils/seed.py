"""Reproducibility helpers."""

from __future__ import annotations

import os
import random


def set_seed(seed: int = 42, deterministic: bool = True) -> None:
    """Seed every RNG the project touches.

    Torch/numpy are imported lazily so this stays usable before they are deps.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)

    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass

    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass

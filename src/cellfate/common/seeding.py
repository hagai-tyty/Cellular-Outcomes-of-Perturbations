"""Global determinism. Seeds Python, NumPy and (if installed) PyTorch.

``torch`` is intentionally imported lazily so the foundation remains usable in a
lean ETL environment that has not installed the model extra.
"""

from __future__ import annotations

import os
import random

import numpy as np


def set_global_seed(seed: int, *, deterministic: bool = True) -> None:
    """Seed all RNGs. ``deterministic`` requests deterministic cuDNN kernels.

    Ensemble members must call this with *distinct* seeds (``base_seed + m``).
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch
    except ImportError:
        return  # lean environment: NumPy/Python seeding is sufficient

    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        # Opt-in; some ops have no deterministic implementation, so warn-only.
        try:
            torch.use_deterministic_algorithms(True, warn_only=True)
        except (AttributeError, RuntimeError):
            pass

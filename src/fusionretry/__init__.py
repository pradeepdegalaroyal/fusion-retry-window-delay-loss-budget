"""fusionretry: repeat-until-success fusion windows versus delay-line loss in a photonic fusion network."""

from .model import (
    FusionParams,
    best_window,
    erasure_probability,
    expected_attempts,
    lattice_blind_window,
    loss_budget_db,
)
from .lattice import failure_rate, winding_failure

__all__ = [
    "FusionParams",
    "best_window",
    "erasure_probability",
    "expected_attempts",
    "lattice_blind_window",
    "loss_budget_db",
    "failure_rate",
    "winding_failure",
]
__version__ = "0.1.0"

"""Threshold crossings and target-erasure inversion from Monte Carlo failure curves."""

from __future__ import annotations

import numpy as np


def crossing(e: np.ndarray, p_small: np.ndarray, p_large: np.ndarray) -> float:
    """Erasure at which the failure curve of the larger lattice crosses the smaller one.

    The larger lattice is below the smaller one under threshold and above it over
    threshold. Returns ``nan`` if no sign change of ``p_large - p_small`` exists.
    """
    diff = np.asarray(p_large) - np.asarray(p_small)
    found = float("nan")
    # Take the last negative-to-nonnegative change: Monte Carlo noise at very low
    # erasure (both curves near zero) must not be mistaken for the threshold.
    for i in range(len(e) - 1):
        if diff[i] < 0.0 <= diff[i + 1] and max(p_small[i + 1], p_large[i + 1]) > 0.02:
            t = -diff[i] / (diff[i + 1] - diff[i])
            found = float(e[i] + t * (e[i + 1] - e[i]))
    return found


def erasure_at_failure(e: np.ndarray, p: np.ndarray, target: float) -> float:
    """Erasure probability at which the (monotonised) failure curve reaches ``target``."""
    p = np.maximum.accumulate(np.asarray(p, dtype=float))
    if target <= p[0]:
        return float(e[0])
    if target >= p[-1]:
        return float(e[-1])
    return float(np.interp(target, p, e))

"""Analytic repeat-until-success (RUS) fusion model with delay-line loss.

One fusion consumes two photons. A fusion node may try up to ``R`` attempts
inside a synchronous window of ``R`` time slots. Attempt ``k`` (0-indexed) has
two independent heralded failure channels:

* the linear-optics fusion itself succeeds with probability ``p_success``;
* both photons must survive. Their transmission is
  ``10 ** (-(base_db + delay_db * (k + c * (R - 1))) / 10)`` where ``base_db``
  is the loss of a photon with no waiting (source, coupling, detector),
  ``delay_db`` is the delay-line loss per time slot, ``k`` slots are spent
  waiting for the retry, and the term ``c * (R - 1)`` is the idle loss that the
  synchronous clock imposes on every photon when the window is ``R`` slots long
  (``c = 1`` fully synchronous, ``c = 0`` retries do not stall the clock).

A fusion whose ``R`` attempts all fail leaves an *erased* fusion outcome.
Erasure is heralded, so the network decoder knows where it occurred.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class FusionParams:
    """Physical parameters of one fusion node.

    Attributes:
        p_success: Success probability of a single heralded fusion attempt.
        delay_db: Delay-line loss per time slot, in dB.
        window_coupling: ``c`` in [0, 1]; fraction of the window length that
            every photon of the network idles (synchronous-clock penalty).
        photons_per_attempt: Photons consumed by one attempt (cost model only).
    """

    p_success: float = 0.75
    delay_db: float = 0.1
    window_coupling: float = 1.0
    photons_per_attempt: int = 2


def attempt_success_probs(window: int, base_db: float, p: FusionParams) -> np.ndarray:
    """Probability that attempt ``k`` yields a usable fusion outcome, for k = 0..window-1."""
    if window < 1:
        raise ValueError("window must be >= 1")
    k = np.arange(window)
    slots = k + p.window_coupling * (window - 1)
    eta = 10.0 ** (-(base_db + p.delay_db * slots) / 10.0)
    return p.p_success * eta**2


def erasure_probability(window: int, base_db: float, p: FusionParams) -> float:
    """Probability that all ``window`` attempts fail, i.e. the fusion outcome is erased."""
    return float(np.prod(1.0 - attempt_success_probs(window, base_db, p)))


def expected_attempts(window: int, base_db: float, p: FusionParams) -> float:
    """Expected number of attempts actually used (<= window)."""
    q = attempt_success_probs(window, base_db, p)
    reach = np.concatenate(([1.0], np.cumprod(1.0 - q)[:-1]))
    return float(reach.sum())


def loss_budget_db(window: int, e_target: float, p: FusionParams, max_db: float = 30.0) -> float:
    """Largest base photon loss (dB) whose erasure probability stays <= ``e_target``.

    Returns ``nan`` if even a lossless photon (0 dB) exceeds the target.
    """
    if erasure_probability(window, 0.0, p) > e_target:
        return float("nan")
    if erasure_probability(window, max_db, p) <= e_target:
        return float(max_db)
    return float(brentq(lambda b: erasure_probability(window, b, p) - e_target, 0.0, max_db))


def best_window(e_target: float, p: FusionParams, r_max: int = 12) -> tuple[int, float, np.ndarray]:
    """Window length maximising the loss budget.

    Returns ``(R_star, budget_db, budgets)``; ``budgets[i]`` belongs to window ``i + 1``
    and is ``nan`` where the target is unreachable. ``R_star = 0`` if nothing is reachable.
    """
    budgets = np.array([loss_budget_db(r, e_target, p) for r in range(1, r_max + 1)])
    if np.all(np.isnan(budgets)):
        return 0, float("nan"), budgets
    idx = int(np.nanargmax(budgets))
    return idx + 1, float(budgets[idx]), budgets


def lattice_blind_window(base_db: float, p: FusionParams, r_max: int = 12) -> int:
    """Classical design rule: window minimising the independent-erasure probability.

    Needs no lattice or decoder knowledge. Used as the lattice-blind baseline.
    """
    e = [erasure_probability(r, base_db, p) for r in range(1, r_max + 1)]
    return int(np.argmin(e)) + 1

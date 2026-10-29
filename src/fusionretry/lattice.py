"""Erasure-correctability of a periodic cubic fusion network (wrapping-cycle criterion).

The network is an ``L x L x L`` periodic cubic lattice whose edges are fusion
outcomes. For a stabiliser code with a cubic-lattice cluster-state structure, a
set of erased edges is correctable (one logical sector) exactly when it contains
no non-contractible cycle. A union-find structure that tracks each node's
displacement to its root in the infinite cover detects such a cycle: an edge that
closes a loop with displacement different from the lattice step is a winding loop.

Erasures may be partly correlated: with fraction ``kappa`` of the marginal edge
erasure ``e`` coming from whole-site outages (all six edges of a node erased)
and the rest from independent edge erasures.
"""

from __future__ import annotations

import numpy as np

try:  # numba is optional; the pure-Python path is slow but identical.
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):  # type: ignore[misc]
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _find(parent, pot, x):
    root = x
    a0 = 0
    a1 = 0
    a2 = 0
    while parent[root] != root:
        a0 += pot[root, 0]
        a1 += pot[root, 1]
        a2 += pot[root, 2]
        root = parent[root]
    cur = x
    while parent[cur] != root and cur != root:
        nxt = parent[cur]
        p0 = pot[cur, 0]
        p1 = pot[cur, 1]
        p2 = pot[cur, 2]
        parent[cur] = root
        pot[cur, 0] = a0
        pot[cur, 1] = a1
        pot[cur, 2] = a2
        a0 -= p0
        a1 -= p1
        a2 -= p2
        cur = nxt
    return root


@njit(cache=True)
def _has_winding(L, erased, parent, pot):
    """True if the erased edge set (flat array of length 3 L^3) holds a winding cycle."""
    n = L * L * L
    for i in range(n):
        parent[i] = i
        pot[i, 0] = 0
        pot[i, 1] = 0
        pot[i, 2] = 0
    for eidx in range(3 * n):
        if not erased[eidx]:
            continue
        u = eidx // 3
        d = eidx % 3
        x = u // (L * L)
        y = (u // L) % L
        z = u % L
        if d == 0:
            x = (x + 1) % L
        elif d == 1:
            y = (y + 1) % L
        else:
            z = (z + 1) % L
        v = (x * L + y) * L + z
        ru = _find(parent, pot, u)
        rv = _find(parent, pot, v)
        if ru == rv:
            # displacement of v relative to u in the cover must equal the unit step
            for c in range(3):
                want = 1 if c == d else 0
                if pot[v, c] - pot[u, c] != want:
                    return True
        else:
            parent[rv] = ru
            for c in range(3):
                want = 1 if c == d else 0
                pot[rv, c] = pot[u, c] + want - pot[v, c]
    return False


@njit(cache=True)
def _failure_count(L, s_site, e_edge, trials, seed):
    np.random.seed(seed)
    n = L * L * L
    parent = np.empty(n, np.int64)
    pot = np.zeros((n, 3), np.int64)
    erased = np.zeros(3 * n, np.bool_)
    dead = np.zeros(n, np.bool_)
    fails = 0
    for _ in range(trials):
        for i in range(n):
            dead[i] = np.random.random() < s_site
        for eidx in range(3 * n):
            u = eidx // 3
            d = eidx % 3
            x = u // (L * L)
            y = (u // L) % L
            z = u % L
            if d == 0:
                x = (x + 1) % L
            elif d == 1:
                y = (y + 1) % L
            else:
                z = (z + 1) % L
            v = (x * L + y) * L + z
            erased[eidx] = dead[u] or dead[v] or (np.random.random() < e_edge)
        if _has_winding(L, erased, parent, pot):
            fails += 1
    return fails


def split_erasure(e: float, kappa: float) -> tuple[float, float]:
    """Split marginal edge erasure ``e`` into (site-outage prob, independent edge prob).

    ``kappa`` is the share of ``e`` produced by site outages; the marginal edge
    erasure probability is exactly ``e`` for every ``kappa``.
    """
    if not 0.0 <= e <= 1.0 or not 0.0 <= kappa <= 1.0:
        raise ValueError("e and kappa must lie in [0, 1]")
    site_part = kappa * e  # marginal edge erasure caused by outages = 1 - (1 - s)^2
    s = 1.0 - np.sqrt(1.0 - site_part)
    e_edge = 1.0 - (1.0 - e) / (1.0 - site_part) if site_part < 1.0 else 0.0
    return float(s), float(max(e_edge, 0.0))


def winding_failure(L: int, erased: np.ndarray) -> bool:
    """Deterministic check on an explicit flat boolean mask of ``3 * L**3`` erased edges.

    Edge ``3 * u + d`` joins node ``u = (x * L + y) * L + z`` to its +d neighbour.
    """
    n = L**3
    erased = np.asarray(erased, dtype=np.bool_)
    if erased.shape != (3 * n,):
        raise ValueError("erased must have length 3 * L**3")
    parent = np.empty(n, np.int64)
    pot = np.zeros((n, 3), np.int64)
    return bool(_has_winding(L, erased, parent, pot))


def failure_rate(L: int, e: float, trials: int, seed: int, kappa: float = 0.0) -> float:
    """Monte Carlo logical-failure probability at marginal edge erasure ``e``."""
    s, e_edge = split_erasure(e, kappa)
    return _failure_count(L, s, e_edge, trials, seed) / trials

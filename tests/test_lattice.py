import networkx as nx
import numpy as np

from fusionretry import failure_rate, winding_failure
from fusionretry.analysis import crossing, erasure_at_failure
from fusionretry.lattice import split_erasure


def _reference_winding(L, erased):
    """Independent check: BFS-lift every cluster and look for an inconsistent edge."""
    n = L**3
    g = nx.Graph()
    steps = {}
    for eidx in np.flatnonzero(erased):
        u, d = divmod(int(eidx), 3)
        x, y, z = u // (L * L), (u // L) % L, u % L
        c = [x, y, z]
        c[d] = (c[d] + 1) % L
        v = (c[0] * L + c[1]) * L + c[2]
        g.add_edge(u, v)
        steps.setdefault((u, v), []).append(d)
    pos = {}
    for comp in nx.connected_components(g):
        start = next(iter(comp))
        pos[start] = (0, 0, 0)
        stack = [start]
        while stack:
            a = stack.pop()
            for (u, v), ds in steps.items():
                if a not in (u, v):
                    continue
                for d in ds:
                    sign = 1 if a == u else -1
                    b = v if a == u else u
                    new = list(pos[a])
                    new[d] += sign
                    new = tuple(new)
                    if b not in pos:
                        pos[b] = new
                        stack.append(b)
                    elif pos[b] != new:
                        return True
    return False


def test_extremes():
    assert failure_rate(4, 0.0, 50, 0) == 0.0
    assert failure_rate(4, 1.0, 50, 0) == 1.0


def test_straight_noncontractible_line_fails_but_open_line_does_not():
    L = 3
    erased = np.zeros(3 * L**3, bool)
    for x in range(L):  # a loop wrapping around the x direction at y = z = 0
        erased[3 * ((x * L + 0) * L + 0) + 0] = True
    assert winding_failure(L, erased)
    erased[3 * ((1 * L + 0) * L + 0) + 0] = False  # break it
    assert not winding_failure(L, erased)


def test_contractible_square_is_harmless():
    L = 4
    erased = np.zeros(3 * L**3, bool)
    idx = lambda x, y, z, d: 3 * ((x * L + y) * L + z) + d
    for e in (idx(0, 0, 0, 0), idx(0, 0, 0, 1), idx(1, 0, 0, 1), idx(0, 1, 0, 0)):
        erased[e] = True
    assert not winding_failure(L, erased)


def test_matches_independent_bfs_reference():
    rng = np.random.default_rng(3)
    L = 3
    for p in (0.15, 0.3, 0.45):
        for _ in range(60):
            erased = rng.random(3 * L**3) < p
            assert winding_failure(L, erased) == _reference_winding(L, erased)


def test_failure_rate_increases_with_erasure():
    r = [failure_rate(4, e, 300, 5) for e in (0.1, 0.25, 0.4)]
    assert r[0] < r[1] < r[2]


def test_threshold_is_near_percolation_value():
    es = np.array([0.22, 0.25, 0.28])
    small = np.array([failure_rate(4, e, 800, 1) for e in es])
    large = np.array([failure_rate(8, e, 800, 1) for e in es])
    assert 0.22 <= crossing(es, small, large) <= 0.28


def test_split_erasure_keeps_marginal():
    for kappa in (0.0, 0.5, 1.0):
        s, e_edge = split_erasure(0.2, kappa)
        marginal = 1 - (1 - s) ** 2 * (1 - e_edge)
        assert abs(marginal - 0.2) < 1e-12


def test_erasure_at_failure_interpolates():
    e = np.array([0.1, 0.2, 0.3])
    p = np.array([0.0, 0.1, 0.9])
    assert abs(erasure_at_failure(e, p, 0.5) - 0.2 - 0.1 * 0.4 / 0.8) < 1e-12

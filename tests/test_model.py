import numpy as np
import pytest

from fusionretry import FusionParams, best_window, erasure_probability, expected_attempts, lattice_blind_window, loss_budget_db


def test_single_attempt_matches_closed_form():
    p = FusionParams(p_success=0.5, delay_db=0.3)
    e = erasure_probability(1, 3.0103, p)  # 3.0103 dB -> eta = 0.5
    assert e == pytest.approx(1 - 0.5 * 0.25, rel=1e-3)


def test_no_clock_penalty_makes_more_retries_always_help():
    p = FusionParams(p_success=0.5, delay_db=0.3, window_coupling=0.0)
    e = [erasure_probability(r, 1.0, p) for r in range(1, 10)]
    assert all(b < a for a, b in zip(e, e[1:]))


def test_clock_penalty_creates_finite_optimum():
    p = FusionParams(p_success=0.75, delay_db=0.3, window_coupling=1.0)
    r_star, budget, budgets = best_window(0.2488, p, r_max=12)
    assert 2 <= r_star < 12
    assert budget > 0.5
    assert np.isnan(budgets[-1])  # naive "retry forever" cannot reach the target


def test_budget_inverts_erasure():
    p = FusionParams(p_success=0.75, delay_db=0.05)
    b = loss_budget_db(4, 0.2, p)
    assert erasure_probability(4, b, p) == pytest.approx(0.2, abs=1e-6)


def test_unreachable_target_is_nan():
    assert np.isnan(loss_budget_db(1, 0.2, FusionParams(p_success=0.5)))


def test_expected_attempts_bounds():
    p = FusionParams(p_success=0.5, delay_db=0.1)
    n = expected_attempts(5, 0.5, p)
    assert 1.0 <= n <= 5.0


def test_lattice_blind_rule_matches_budget_optimum():
    """Erasure is monotone for any decoder, so the blind rule picks the same window."""
    p = FusionParams(p_success=0.75, delay_db=0.3)
    r_star, budget, _ = best_window(0.2488, p)
    assert lattice_blind_window(budget, p) == r_star

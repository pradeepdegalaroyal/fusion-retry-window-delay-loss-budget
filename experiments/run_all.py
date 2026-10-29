"""Reproduce every result and figure: ``python experiments/run_all.py`` (add ``--quick`` for a smoke run)."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from fusionretry import FusionParams, best_window, erasure_probability, failure_rate, lattice_blind_window
from fusionretry.analysis import crossing, erasure_at_failure
from fusionretry.model import loss_budget_db
from fusionretry import plotting
from summarise import headline

ROOT = Path(__file__).resolve().parents[1]
TARGET_PL = 0.05  # logical failure budget used to define the tolerable loss
PHOTONS = {0.5: 2, 0.75: 4}  # photons per attempt: bare pair; pair plus Bell-pair ancilla


def make_config(quick: bool) -> dict:
    if quick:
        return dict(
            sizes=(4, 6), kappas=(0.0, 1.0), e_grid=np.round(np.arange(0.14, 0.36, 0.03), 3),
            trials=120, seeds=2, l_ref=6, sm_lams=(0.05, 0.3), sm_windows=(1, 2, 4),
            sm_base=np.arange(0.0, 3.01, 0.5), sm_trials=100, r_max=12,
        )
    return dict(
        sizes=(4, 6, 8, 10), kappas=(0.0, 0.5, 1.0),
        e_grid=np.round(np.r_[np.arange(0.10, 0.20, 0.02), np.arange(0.20, 0.30, 0.01), np.arange(0.30, 0.37, 0.02)], 3),
        trials=1500, seeds=5, l_ref=8, sm_lams=(0.01, 0.05, 0.1, 0.3), sm_windows=(1, 2, 3, 4, 6, 8),
        sm_base=np.arange(0.0, 3.01, 0.25), sm_trials=800, r_max=30,
    )


def lattice_sweep(cfg: dict) -> pd.DataFrame:
    rows = []
    for kappa in cfg["kappas"]:
        for L in cfg["sizes"]:
            for e in cfg["e_grid"]:
                for s in range(cfg["seeds"]):
                    seed = 100003 * s + 101 * L + int(round(e * 1000)) + 7919 * int(round(kappa * 10)) + 1
                    rows.append(dict(kappa=kappa, L=L, e=float(e), seed=s,
                                     fail_rate=failure_rate(L, float(e), cfg["trials"], seed, kappa), trials=cfg["trials"]))
    return pd.DataFrame(rows)


def thresholds(pl: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    sizes = list(cfg["sizes"])
    pairs = list(zip(sizes[:-1], sizes[1:]))
    rows = []
    for kappa in cfg["kappas"]:
        for l1, l2 in pairs:
            cs = []
            for s in range(cfg["seeds"]):
                d = pl[(pl.kappa == kappa) & (pl.seed == s)]
                a = d[d.L == l1].sort_values("e")
                b = d[d.L == l2].sort_values("e")
                cs.append(crossing(a.e.to_numpy(), a.fail_rate.to_numpy(), b.fail_rate.to_numpy()))
            label = f"{l1}-{l2}"
            rows.append(dict(kappa=kappa, pair=label, crossing_mean=float(np.nanmean(cs)) if np.any(np.isfinite(cs)) else np.nan,
                             crossing_std=float(np.nanstd(cs)) if np.any(np.isfinite(cs)) else np.nan, n_valid=int(np.isfinite(cs).sum())))
            if (l1, l2) == pairs[-1]:
                rows.append({**rows[-1], "pair": "largest"})
    return pd.DataFrame(rows)


def targets(pl: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rows = []
    for kappa in cfg["kappas"]:
        for s in range(cfg["seeds"]):
            d = pl[(pl.kappa == kappa) & (pl.seed == s) & (pl.L == cfg["l_ref"])].sort_values("e")
            rows.append(dict(kappa=kappa, seed=s, e_target=erasure_at_failure(d.e.to_numpy(), d.fail_rate.to_numpy(), TARGET_PL)))
    return pd.DataFrame(rows)


def budget_grid(cfg: dict, tgt: pd.DataFrame, e_star: dict[float, float]) -> pd.DataFrame:
    lams = sorted(set(np.round(np.geomspace(0.01, 0.5, 10), 4)) | {0.01, 0.05, 0.1, 0.3, 0.5})
    rows = []
    combos = [(k, 1.0) for k in cfg["kappas"]] + [(0.0, c) for c in (0.0, 0.5)]
    for kappa, c in combos:
        for ps in (0.5, 0.75, 0.9):
            for lam in lams:
                p = FusionParams(ps, lam, c, PHOTONS.get(ps, 0))
                targets_ = [("threshold", -1, e_star[kappa])] + [
                    ("pl5", int(r.seed), float(r.e_target)) for r in tgt[tgt.kappa == kappa].itertuples()]
                for kind, seed, et in targets_:
                    for r in range(1, cfg["r_max"] + 1):
                        rows.append(dict(kappa=kappa, c=c, p_success=ps, delay_db=lam, kind=kind, seed=seed,
                                         e_target=et, window=r, budget_db=loss_budget_db(r, et, p),
                                         cost=r * PHOTONS.get(ps, np.nan)))
    return pd.DataFrame(rows)


def phase_grid(cfg: dict, e_target: float) -> pd.DataFrame:
    rows = []
    for ps in np.round(np.linspace(0.5, 0.95, 10), 3):
        for lam in np.geomspace(0.01, 0.5, 14):
            r, b, _ = best_window(e_target, FusionParams(float(ps), float(lam), 1.0), cfg["r_max"])
            rows.append(dict(p_success=float(ps), delay_db=float(lam), r_star=r, budget_db=b))
    return pd.DataFrame(rows)


def direct_small_multiples(cfg: dict) -> pd.DataFrame:
    rows = []
    for lam in cfg["sm_lams"]:
        p = FusionParams(0.75, lam, 1.0)
        for r in cfg["sm_windows"]:
            for base in cfg["sm_base"]:
                e = erasure_probability(r, float(base), p)
                for s in range(cfg["seeds"]):
                    seed = 7 + 31 * s + 1009 * r + int(1000 * lam) + int(100 * base)
                    rows.append(dict(delay_db=lam, window=r, base_db=float(base), e=e, seed=s,
                                     fail=failure_rate(cfg["l_ref"], e, cfg["sm_trials"], seed)))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="smoke run, well under a minute")
    args = ap.parse_args()
    cfg = make_config(args.quick)
    res, figs = ROOT / ("results_quick" if args.quick else "results"), ROOT / ("figures_quick" if args.quick else "figures")
    res.mkdir(exist_ok=True)
    figs.mkdir(exist_ok=True)
    t0 = time.time()

    pl = lattice_sweep(cfg)
    pl.to_csv(res / "lattice_failure.csv", index=False)
    thr = thresholds(pl, cfg)
    thr.to_csv(res / "lattice_thresholds.csv", index=False)
    tgt = targets(pl, cfg)
    tgt.to_csv(res / "lattice_targets.csv", index=False)
    print(f"[{time.time() - t0:5.0f}s] lattice sweep done ({len(pl)} points)")

    e_star = {}
    for k in cfg["kappas"]:
        v = thr[(thr.kappa == k) & (thr.pair == "largest")].crossing_mean.iloc[0]
        e_star[k] = float(v) if np.isfinite(v) else 0.2488
    bud = budget_grid(cfg, tgt, e_star)
    bud.to_csv(res / "budget_grid.csv", index=False)

    mean_target0 = float(tgt[tgt.kappa == 0.0].e_target.mean())
    phase = phase_grid(cfg, mean_target0)
    phase.to_csv(res / "phase_map.csv", index=False)
    sm = direct_small_multiples(cfg)
    sm.to_csv(res / "direct_failure_vs_loss.csv", index=False)
    print(f"[{time.time() - t0:5.0f}s] budgets, phase map, direct sweeps done")

    # interpolation check: direct Monte Carlo versus the tabulated failure curve
    curve = pl[(pl.kappa == 0.0) & (pl.L == cfg["l_ref"])].groupby("e").fail_rate.mean()
    dm = sm.groupby(["delay_db", "window", "base_db", "e"]).fail.mean().reset_index()
    interp = np.interp(dm.e, curve.index.to_numpy(), curve.to_numpy())
    mask = (dm.fail > 0.02) & (dm.fail < 0.98) & (dm.e >= curve.index.min()) & (dm.e <= curve.index.max())
    interp_mae = float(np.mean(np.abs(interp[mask] - dm.fail[mask]))) if mask.any() else float("nan")

    # per-cell optima (target: 5 percent logical failure at the reference lattice)
    pl5 = bud[(bud.kind == "pl5")]
    cells = []
    for (kappa, c, ps, lam, seed), g in pl5.groupby(["kappa", "c", "p_success", "delay_db", "seed"]):
        b = g.sort_values("window").budget_db.to_numpy()
        if np.all(np.isnan(b)):
            cells.append(dict(kappa=kappa, c=c, p_success=ps, delay_db=lam, seed=seed, r_star=0, budget_db=np.nan,
                              budget_r1=b[0], budget_rmax=b[-1]))
        else:
            i = int(np.nanargmax(b))
            cells.append(dict(kappa=kappa, c=c, p_success=ps, delay_db=lam, seed=seed, r_star=i + 1, budget_db=float(b[i]),
                              budget_r1=b[0], budget_rmax=b[-1]))
    cells = pd.DataFrame(cells)
    cells.to_csv(res / "optimal_window_cells.csv", index=False)
    summ = (cells.groupby(["kappa", "c", "p_success", "delay_db"])
            .agg(r_star_mean=("r_star", "mean"), r_star_std=("r_star", "std"), budget_mean=("budget_db", "mean"),
                 budget_std=("budget_db", "std"), r1_reachable=("budget_r1", lambda x: float(np.mean(~np.isnan(x)))),
                 rmax_reachable=("budget_rmax", lambda x: float(np.mean(~np.isnan(x))))).reset_index())
    summ.to_csv(res / "optimal_window_summary.csv", index=False)

    base = summ[(summ.kappa == 0.0) & (summ.c == 1.0)]
    # lattice-blind rule agreement
    agree, total = 0, 0
    for r in base.itertuples():
        if r.r_star_mean >= 1 and np.isfinite(r.budget_mean):
            p = FusionParams(r.p_success, r.delay_db, 1.0)
            rs, b, _ = best_window(mean_target0, p, cfg["r_max"])
            if rs:
                total += 1
                agree += int(lattice_blind_window(b, p, cfg["r_max"]) == rs)

    def cell(ps, lam, kappa=0.0, c=1.0):
        return summ[(summ.kappa == kappa) & (summ.c == c) & (summ.p_success == ps) & np.isclose(summ.delay_db, lam)].iloc[0]

    c0 = bud[(bud.kind == "pl5") & (bud.kappa == 0.0) & (bud.c == 0.0) & (bud.seed == 0)]
    monotone_c0 = bool(all(
        np.all(np.diff(g.sort_values("window").budget_db.dropna().to_numpy()) > 0) for _, g in c0.groupby(["p_success", "delay_db"])))
    rob_rows = []
    for ps in (0.5, 0.75):
        for k in cfg["kappas"]:
            g = cells[(cells.kappa == k) & (cells.c == 1.0) & (cells.p_success == ps) & np.isclose(cells.delay_db, 0.1)]
            rob_rows.append(dict(p_success=ps, kappa=k, budget_mean=float(g.budget_db.mean()), budget_std=float(g.budget_db.std()),
                                 r_star_mean=float(g.r_star.mean())))
    rob = pd.DataFrame(rob_rows)
    rob.to_csv(res / "correlated_outage_robustness.csv", index=False)

    summary = dict(
        target_logical_failure=TARGET_PL, reference_lattice=cfg["l_ref"], seeds=cfg["seeds"], trials_per_point=cfg["trials"],
        threshold_crossings={f"kappa={k}": {"largest_pair": thr[(thr.kappa == k) & (thr.pair == 'largest')].iloc[0][["crossing_mean", "crossing_std"]].to_dict()}
                             for k in cfg["kappas"]},
        e_target_mean_std={f"kappa={k}": [float(tgt[tgt.kappa == k].e_target.mean()), float(tgt[tgt.kappa == k].e_target.std())] for k in cfg["kappas"]},
        interpolation_mae_vs_direct_mc=interp_mae,
        lattice_blind_rule_agreement=[agree, total],
        c0_budget_monotone_in_window=monotone_c0,
        single_shot_reachable_fraction={str(ps): float(base[base.p_success == ps].r1_reachable.mean()) for ps in (0.5, 0.75, 0.9)},
        max_window_reachable_fraction={str(ps): float(base[base.p_success == ps].rmax_reachable.mean()) for ps in (0.5, 0.75, 0.9)},
        key_cells={f"ps={ps}_lam={lam}": cell(ps, lam)[["r_star_mean", "r_star_std", "budget_mean", "budget_std", "r1_reachable", "rmax_reachable"]].to_dict()
                   for ps in (0.5, 0.75, 0.9) for lam in (0.05, 0.1, 0.3)},
        runtime_seconds=None,
    )

    # figures
    plotting.fig_small_multiples(sm, TARGET_PL, cfg["l_ref"], figs / "fig1_logical_failure_small_multiples.png")
    pareto = bud[(bud.kind == "pl5") & (bud.kappa == 0.0) & (bud.c == 1.0) & bud.p_success.isin([0.5, 0.75])
                 & (np.isclose(bud.delay_db, 0.1) | np.isclose(bud.delay_db, 0.3))]
    pareto = pareto.groupby(["p_success", "delay_db", "window", "cost"]).budget_db.mean().reset_index()
    plotting.fig_pareto(pareto, figs / "fig2_cost_vs_loss_budget_pareto.png")
    plotting.fig_phase(phase, figs / "fig3_phase_map_window_and_budget.png")
    plotting.fig_threshold(pl, thr, figs / "fig4_threshold_scaling_correlated_outages.png")
    abl = bud[(bud.kind == "pl5") & (bud.kappa == 0.0)].groupby(["p_success", "delay_db", "c", "window"]).budget_db.mean().reset_index()
    plotting.fig_ablation(abl, rob, figs / "fig5_clock_ablation_and_robustness.png")

    summary["runtime_seconds"] = round(time.time() - t0, 1)
    (res / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    (res / "headline_numbers.json").write_text(json.dumps(headline(res), indent=2))
    print(f"[{time.time() - t0:5.0f}s] done; results in {res.name}/, figures in {figs.name}/")


if __name__ == "__main__":
    main()

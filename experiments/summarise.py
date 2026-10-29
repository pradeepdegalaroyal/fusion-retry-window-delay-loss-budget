"""Headline numbers for the README, computed only from the CSV files in a results folder."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _best(d: pd.DataFrame) -> dict:
    d = d.dropna(subset=["budget_db"])
    if d.empty:
        return {"window": 0, "budget_db": None}
    i = d.budget_db.idxmax()
    return {"window": int(d.loc[i, "window"]), "budget_db": round(float(d.loc[i, "budget_db"]), 3)}


def headline(res: Path) -> dict:
    b = pd.read_csv(res / "budget_grid.csv")
    pl5 = b[(b.kind == "pl5") & (b.kappa == 0.0)].groupby(["c", "p_success", "delay_db", "window"]).budget_db.mean().reset_index()
    thr = b[(b.kind == "threshold") & (b.kappa == 0.0) & (b.c == 1.0)]
    out: dict = {"clock_ablation": {}, "diminishing_returns": {}, "threshold_target": {}, "minimum_window": {}}
    for ps, lam in [(0.5, 0.1), (0.75, 0.1), (0.5, 0.3), (0.75, 0.3)]:
        key = f"ps={ps}_lam={lam}"
        sel = lambda c: pl5[(pl5.c == c) & (pl5.p_success == ps) & np.isclose(pl5.delay_db, lam)]
        out["clock_ablation"][key] = {f"c={c}": _best(sel(c)) for c in (0.0, 0.5, 1.0)}
        rmax = sel(1.0)
        out["clock_ablation"][key]["c=1_budget_at_R30"] = None if rmax[rmax.window == rmax.window.max()].budget_db.isna().all() else round(float(rmax[rmax.window == rmax.window.max()].budget_db.iloc[0]), 3)
        out["threshold_target"][key] = _best(thr[(thr.p_success == ps) & np.isclose(thr.delay_db, lam)].groupby("window").budget_db.mean().reset_index())
    for ps in (0.5, 0.75):
        d = pl5[(pl5.c == 1.0) & (pl5.p_success == ps) & np.isclose(pl5.delay_db, 0.1)].dropna(subset=["budget_db"])
        m = float(d.budget_db.max())
        r90 = int(d[d.budget_db >= 0.9 * m].window.min())
        out["diminishing_returns"][f"ps={ps}_lam=0.1"] = {
            "best": _best(d), "window_reaching_90pct_of_best": r90,
            "budget_at_that_window": round(float(d[d.window == r90].budget_db.iloc[0]), 3),
            "budget_at_R4": round(float(d[d.window == 4].budget_db.iloc[0]), 3),
            "budget_at_R8": round(float(d[d.window == 8].budget_db.iloc[0]), 3)}
        out["minimum_window"][f"ps={ps}_lam=0.1"] = {"window": int(d.window.min()), "budget_db": round(float(d.budget_db.iloc[0]), 3)}
    ph = pd.read_csv(res / "phase_map.csv")
    out["phase_map"] = {"fraction_unreachable": round(float((ph.r_star == 0).mean()), 3),
                        "max_budget_db": round(float(ph.budget_db.max()), 3), "cells": int(len(ph))}
    return out


if __name__ == "__main__":
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
    (folder / "headline_numbers.json").write_text(json.dumps(headline(folder), indent=2))
    print("wrote", folder / "headline_numbers.json")

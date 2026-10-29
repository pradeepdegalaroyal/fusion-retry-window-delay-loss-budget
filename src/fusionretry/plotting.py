"""Figures for the fusion-window study (earth-tone palette, 200 dpi)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

INK = "#2f2a24"
SAND = "#e8dcc2"
OCHRE = "#c8963e"
SIENNA = "#a0522d"
BROWN = "#6b4423"
TEAL = "#2a7f78"
EARTH_CMAP = LinearSegmentedColormap.from_list("earth", [SAND, OCHRE, SIENNA, BROWN, "#2f3e3a"])
TEAL_CMAP = LinearSegmentedColormap.from_list("teal_earth", [SAND, "#8fbfb0", TEAL, "#1c4f4b"])

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "#fbf8f1",
})


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_small_multiples(df: pd.DataFrame, target: float, lattice_l: int, out: Path) -> None:
    """Direct Monte Carlo logical failure versus base photon loss, one panel per delay loss."""
    lams = sorted(df["delay_db"].unique())
    windows = sorted(df["window"].unique())
    cols = [EARTH_CMAP(x) for x in np.linspace(0.12, 1.0, len(windows))]
    fig, axes = plt.subplots(1, len(lams), figsize=(3.0 * len(lams), 3.2), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, lam in zip(axes, lams):
        for r, c in zip(windows, cols):
            d = df[(df.delay_db == lam) & (df.window == r)].groupby("base_db")["fail"].agg(["mean", "std"]).reset_index()
            ax.plot(d.base_db, d["mean"], color=c, lw=1.6, label=f"R = {r}")
            ax.fill_between(d.base_db, d["mean"] - d["std"], d["mean"] + d["std"], color=c, alpha=0.18, lw=0)
        ax.axhline(target, color=TEAL, ls="--", lw=1)
        ax.set_title(f"delay loss {lam:g} dB per slot", fontsize=9)
        ax.set_xlabel("base photon loss (dB)")
    axes[0].set_ylabel(f"logical failure probability (L = {lattice_l})")
    axes[-1].legend(fontsize=7, frameon=False, title="window R", title_fontsize=7, loc="center right")
    fig.suptitle("Logical failure of the cubic fusion network, fusion success 0.75, synchronous clock", y=1.02)
    _save(fig, out)


def fig_pareto(df: pd.DataFrame, out: Path) -> None:
    """Hardware cost (photon-slots per fusion window) versus loss budget, Pareto front highlighted."""
    lams = sorted(df["delay_db"].unique())
    fig, axes = plt.subplots(1, len(lams), figsize=(4.4 * len(lams), 3.6), sharey=True)
    axes = np.atleast_1d(axes)
    colors = {0.5: TEAL, 0.75: OCHRE}
    for ax, lam in zip(axes, lams):
        d = df[(df.delay_db == lam) & df.budget_db.notna()]
        for ps, c in colors.items():
            s = d[d.p_success == ps]
            ax.scatter(s.cost, s.budget_db, s=26, color=c, label=f"fusion success {ps:g}", zorder=3)
            for _, row in s.iterrows():
                if int(row.window) in (1, 2, 3, 4, 6, 8, 12):
                    ax.annotate(f"R={int(row.window)}", (row.cost, row.budget_db), fontsize=6,
                                xytext=(3, -8), textcoords="offset points")
        pts = d.sort_values(["cost", "budget_db"], ascending=[True, False])
        front, best = [], -np.inf
        for _, row in pts.iterrows():
            if row.budget_db > best + 1e-12:
                front.append((row.cost, row.budget_db))
                best = row.budget_db
        if front:
            fx, fy = zip(*front)
            ax.step(fx, fy, where="post", color=BROWN, lw=1.5, label="Pareto front", zorder=2)
        ax.set_title(f"delay loss {lam:g} dB per slot", fontsize=9)
        ax.set_xlabel("cost: photons per attempt x window slots")
    axes[0].set_ylabel("tolerable base photon loss (dB)")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Cost versus loss budget at 5 percent logical failure, L = 8 (unreachable windows omitted)", y=1.03)
    _save(fig, out)


def fig_phase(df: pd.DataFrame, out: Path) -> None:
    """Filled contours over (delay loss, fusion success): optimal window and its loss budget."""
    lam = np.sort(df.delay_db.unique())
    ps = np.sort(df.p_success.unique())
    rstar = df.pivot(index="p_success", columns="delay_db", values="r_star").loc[ps, lam].to_numpy(float)
    bud = df.pivot(index="p_success", columns="delay_db", values="budget_db").loc[ps, lam].to_numpy(float)
    rstar = np.where(rstar == 0, np.nan, rstar)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True)
    levels = np.arange(0.5, np.nanmax(rstar) + 1.5, 1.0)
    cs = axes[0].contourf(lam, ps, rstar, levels=levels, cmap=EARTH_CMAP)
    cb = fig.colorbar(cs, ax=axes[0])
    cb.set_label("optimal window R* (slots)")
    axes[0].set_title("Optimal repeat-until-success window", fontsize=9)
    cs2 = axes[1].contourf(lam, ps, bud, levels=12, cmap=TEAL_CMAP)
    cb2 = fig.colorbar(cs2, ax=axes[1])
    cb2.set_label("tolerable base photon loss at R* (dB)")
    axes[1].contour(lam, ps, bud, levels=[0.5, 1.0, 2.0], colors=INK, linewidths=0.7)
    axes[1].set_title("Loss budget at the optimal window", fontsize=9)
    for ax in axes:
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("delay loss per slot (dB)")
    axes[0].set_ylabel("single-attempt fusion success probability")
    _save(fig, out)


def fig_threshold(pl: pd.DataFrame, thr: pd.DataFrame, out: Path) -> None:
    """Logical failure versus edge erasure for each lattice size, one panel per outage share kappa."""
    kappas = sorted(pl.kappa.unique())
    sizes = sorted(pl.L.unique())
    cols = [EARTH_CMAP(x) for x in np.linspace(0.2, 1.0, len(sizes))]
    fig, axes = plt.subplots(1, len(kappas), figsize=(3.5 * len(kappas), 3.3), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, kappas):
        for L, c in zip(sizes, cols):
            d = pl[(pl.kappa == k) & (pl.L == L)].groupby("e")["fail_rate"].agg(["mean", "std"]).reset_index()
            ax.errorbar(d.e, d["mean"], yerr=d["std"], color=c, lw=1.4, capsize=1.5, label=f"L = {L}")
        t = thr[(thr.kappa == k) & (thr.pair == "largest")]
        if len(t) and np.isfinite(t.crossing_mean.iloc[0]):
            x0 = t.crossing_mean.iloc[0]
            ax.axvline(x0, color=TEAL, ls="--", lw=1)
            ax.text(x0 + 0.004, 0.06, f"crossing {x0:.3f}", fontsize=7, color=TEAL)
        ax.set_title(f"site-outage share kappa = {k:g}", fontsize=9)
        ax.set_xlabel("marginal fusion erasure probability")
    axes[0].set_ylabel("logical failure probability")
    axes[0].legend(fontsize=7, frameon=False)
    _save(fig, out)


def fig_ablation(bud: pd.DataFrame, rob: pd.DataFrame, out: Path) -> None:
    """Clock-coupling ablation (two operating points) and correlated-outage robustness."""
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.5))
    cmap = {0.0: TEAL, 0.5: OCHRE, 1.0: BROWN}
    for ax, (ps, lam) in zip(axes[:2], [(0.5, 0.1), (0.75, 0.3)]):
        for c, col in cmap.items():
            d = bud[(bud.p_success == ps) & (np.isclose(bud.delay_db, lam)) & (bud.c == c)].sort_values("window")
            ax.plot(d.window, d.budget_db, "o-", ms=3.5, color=col, label=f"clock coupling c = {c:g}")
        ax.set_title(f"fusion success {ps:g}, delay {lam:g} dB per slot", fontsize=9)
        ax.set_xlabel("window R (slots)")
    axes[0].set_ylabel("tolerable base photon loss (dB)")
    axes[0].legend(fontsize=7, frameon=False)
    for ps, col in [(0.5, TEAL), (0.75, OCHRE)]:
        d = rob[rob.p_success == ps].sort_values("kappa")
        axes[2].errorbar(d.kappa, d.budget_mean, yerr=d.budget_std, fmt="o-", color=col, capsize=2,
                         label=f"fusion success {ps:g}")
    axes[2].set_title("Correlated site outages, delay 0.1 dB per slot", fontsize=9)
    axes[2].set_xlabel("site-outage share kappa")
    axes[2].set_ylabel("loss budget at best window (dB)")
    axes[2].legend(fontsize=7, frameon=False)
    _save(fig, out)

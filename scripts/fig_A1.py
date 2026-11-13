#!/usr/bin/env python3
"""Minimal draft of confirmatory Figure A1 (paired normal vs tumour HP load).

Official statistics: results/core/A1_confirmatory.json
Per-pair input:      results/core/korea_rna_per_sample.csv
(assembled and tested by scripts/core_a1.py / scripts/core_a1_stats.py)

Numbers drawn on the figure are verified against A1_confirmatory.json.
The only on-figure statistic is a compact one-sided Wilcoxon P value.
All other numerical details live in results/figures/draft/Fig_A1_caption.md.

Y-axis metric (as used by the Wilcoxon test):
    hp_log1p_rpm = ln(1 + HP RPM)
    HP RPM = Kraken2-confirmed taxid 210 clade reads
             / Salmon-mapped host fragments × 10^6

Single-threaded. Does not read or write results/prjeb9749.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
from scipy.stats import gaussian_kde, wilcoxon

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
B = Path("/REPLACE_WITH_PRIVATE_WORKDIR")
CORE = B / "results/core"
JSON_PATH = CORE / "A1_confirmatory.json"
DATA_PATH = CORE / "korea_rna_per_sample.csv"
OUTDIR = B / "results/figures/draft"
SEED = 20261030
N_BOOT = 10_000

# Okabe–Ito colour-blind-safe palette
COL_NGT = "#009E73"  # bluish green: N > T
COL_TGN = "#D55E00"  # vermillion: T > N
COL_TIE = "#999999"  # grey: tie
COL_EST = "#0072B2"  # blue: bootstrap estimate
COL_INK = "#222222"
COL_SPINE = "#333333"

# ---------------------------------------------------------------------------
# Fonts: static Source Sans Pro (editable TrueType in PDF/SVG)
# ---------------------------------------------------------------------------
FONT_DIR = Path("/usr/share/fonts/truetype/sand-box/custom/Source Sans Pro")
for fn in (
    "SourceSansPro-Regular.ttf",
    "SourceSansPro-Semibold.ttf",
    "SourceSansPro-Bold.ttf",
    "SourceSansPro-It.ttf",
    "SourceSansPro-BoldIt.ttf",
):
    font_manager.fontManager.addfont(str(FONT_DIR / fn))

mpl.rcParams.update(
    {
        "font.family": "Source Sans Pro",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.edgecolor": COL_SPINE,
        "axes.labelcolor": COL_INK,
        "xtick.color": COL_INK,
        "ytick.color": COL_INK,
        "text.color": COL_INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "pdf.compression": 6,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "none",
        "axes.unicode_minus": False,
    }
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def assert_match(name, computed, official, atol=0.0, rtol=0.0):
    if isinstance(computed, (float, np.floating)) or isinstance(
        official, (float, np.floating)
    ):
        ok = np.isclose(computed, official, atol=atol, rtol=rtol, equal_nan=True)
    else:
        ok = computed == official
    if not ok:
        raise SystemExit(
            f"MISMATCH vs A1_confirmatory.json: {name}: "
            f"computed={computed!r} official={official!r}"
        )


def half_violin(
    ax, values, x, side="left", width=0.22, color="#0072B2",
    alpha=0.22, lw=0.7, n=256, zorder=1,
):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3:
        return
    ymin, ymax = float(values.min()), float(values.max())
    pad = 0.04 * (ymax - ymin if ymax > ymin else 1.0)
    y = np.linspace(ymin - pad, ymax + pad, n)
    kde = gaussian_kde(values, bw_method="scott")
    v = kde(y)
    v = v / v.max() * width
    if side == "left":
        xx = x - v
        ax.fill_betweenx(
            y, xx, x, color=color, alpha=alpha, linewidth=0,
            zorder=zorder, clip_on=False,
        )
        ax.plot(xx, y, color=color, lw=lw, alpha=0.85,
                zorder=zorder + 1, solid_capstyle="round")
    else:
        xx = x + v
        ax.fill_betweenx(
            y, x, xx, color=color, alpha=alpha, linewidth=0,
            zorder=zorder, clip_on=False,
        )
        ax.plot(xx, y, color=color, lw=lw, alpha=0.85,
                zorder=zorder + 1, solid_capstyle="round")


def iqr_median_marker(ax, x, values, color=COL_INK, zorder=6):
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    ax.plot([x, x], [q1, q3], color=color, lw=1.5, solid_capstyle="round", zorder=zorder)
    ax.scatter(
        [x], [med],
        s=26, marker="o",
        facecolor="white", edgecolor=color, linewidth=1.15,
        zorder=zorder + 1,
    )
    return float(med)


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="both", which="both", length=3.0, width=0.6, pad=2.5)
    ax.tick_params(axis="x", which="minor", length=0)


def density_beeswarm(y, x_centre=0.0, max_spread=0.16, seed=SEED + 1):
    y = np.asarray(y, dtype=float)
    n = y.size
    if n == 0:
        return np.full(0, x_centre)
    kde = gaussian_kde(y, bw_method="scott")
    dens = kde(y)
    dens = dens / dens.max()
    order = np.argsort(y, kind="mergesort")
    pos = np.empty(n, dtype=int)
    pos[order] = np.arange(n)
    side = np.where(pos % 2 == 0, 1.0, -1.0)
    spread = 0.016 + max_spread * dens
    rng = np.random.default_rng(seed)
    return x_centre + side * spread + rng.uniform(-0.005, 0.005, size=n)


def letter_at(fig, ax, letter, dx=-0.12, dy=0.04):
    bbox = ax.get_position()
    fig.text(
        bbox.x0 + dx * bbox.width,
        bbox.y1 + dy * bbox.height,
        letter,
        fontsize=12, fontweight="bold", color=COL_INK,
        ha="left", va="bottom",
    )


# ---------------------------------------------------------------------------
# Load data and verify against the official JSON
# ---------------------------------------------------------------------------
def load_and_verify():
    official = json.loads(JSON_PATH.read_text())
    D = pd.read_csv(DATA_PATH)

    assert_match("n_runs", len(D), official["n_runs"])
    assert_match("n_patients", int(D.pid.nunique()), official["n_patients"])
    types = D.type.value_counts().to_dict()
    assert_match("types.T", int(types.get("T", 0)), official["types"]["T"])
    assert_match("types.N", int(types.get("N", 0)), official["types"]["N"])
    assert_match("n_hp_ge10", int((D.hp_k2_confirmed >= 10).sum()), official["n_hp_ge10"])
    assert_match("n_hp_ge1000", int((D.hp_k2_confirmed >= 1000).sum()), official["n_hp_ge1000"])
    assert_match("max_hp", int(D.hp_k2_confirmed.max()), official["max_hp"])

    w = D.pivot_table(index="pid", columns="type", values="hp_log1p_rpm")
    n_single = int(w.isna().any(axis=1).sum())
    assert_match(
        "patients_with_only_one_tissue",
        n_single,
        official["patients_with_only_one_tissue"],
    )
    w = w.dropna()
    assert_match("pairs", len(w), official["pairs"])

    n_gt = int((w.N > w["T"]).sum())
    t_gt = int((w["T"] > w.N).sum())
    ties = int((w.N == w["T"]).sum())
    assert_match("N_gt_T", n_gt, official["N_gt_T"])
    assert_match("T_gt_N", t_gt, official["T_gt_N"])
    assert_match("ties", ties, official["ties"])

    pw = wilcoxon(w.N, w["T"], alternative="greater")
    pt = wilcoxon(w.N, w["T"])
    assert_match(
        "wilcoxon_P_one_sided_N_gt_T",
        float(pw.pvalue),
        official["wilcoxon_P_one_sided_N_gt_T"],
    )
    assert_match("wilcoxon_P_two_sided", float(pt.pvalue), official["wilcoxon_P_two_sided"])

    dd = w.N - w["T"]
    d2 = dd[dd != 0]
    n = len(d2)
    rk = d2.abs().rank()
    rp = rk[d2 > 0].sum()
    tot = n * (n + 1) / 2
    r_rb = float((rp - (tot - rp)) / tot)
    assert_match("rank_biserial_r_N_vs_T", r_rb, official["rank_biserial_r_N_vs_T"])

    W = D.pivot_table(index="pid", columns="type", values="hp_rpm_host").loc[w.index]
    assert_match("median_rpm_N", float(W.N.median()), official["median_rpm_N"])
    assert_match("median_rpm_T", float(W["T"].median()), official["median_rpm_T"])

    tie_idx = w.index[w.N == w["T"]]
    ties_both_zero = int(((W.loc[tie_idx, "N"] == 0) & (W.loc[tie_idx, "T"] == 0)).sum())

    pairs = pd.DataFrame(
        {
            "pid": w.index.astype(str),
            "N": w.N.to_numpy(),
            "T": w["T"].to_numpy(),
            "rpm_N": W.N.to_numpy(),
            "rpm_T": W["T"].to_numpy(),
        }
    )
    pairs["delta_T_minus_N"] = pairs["T"] - pairs["N"]
    pairs["direction"] = np.where(
        pairs["N"] > pairs["T"],
        "N>T",
        np.where(pairs["T"] > pairs["N"], "T>N", "tie"),
    )

    stats = dict(
        official=official,
        n_pairs=len(pairs),
        n_gt=n_gt,
        t_gt=t_gt,
        ties=ties,
        ties_both_zero=ties_both_zero,
        p_one=float(pw.pvalue),
        p_two=float(pt.pvalue),
        r_rb=r_rb,
        r_ci=official["rank_biserial_95CI_bootstrap"],
        median_rpm_N=float(W.N.median()),
        median_rpm_T=float(W["T"].median()),
        median_log1p_N=float(w.N.median()),
        median_log1p_T=float(w["T"].median()),
        n_single=n_single,
    )
    return pairs, stats


def bootstrap_delta(values, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    draws = rng.choice(values, size=(n_boot, values.size), replace=True)
    boot_median = np.median(draws, axis=1)
    boot_mean = draws.mean(axis=1)
    return dict(
        boot_median=boot_median,
        median=float(np.median(values)),
        median_ci=np.percentile(boot_median, [2.5, 97.5]),
        boot_mean=boot_mean,
        mean=float(np.mean(values)),
        mean_ci=np.percentile(boot_mean, [2.5, 97.5]),
    )


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def draw_figure(pairs, stats, boot):
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Compact double-column layout, no title band.
    fig = plt.figure(figsize=(6.90, 4.55), facecolor="white")
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[3.35, 0.52],
        width_ratios=[1.25, 1.00],
        left=0.10, right=0.985,
        top=0.88, bottom=0.10,
        wspace=0.28, hspace=0.55,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    colour_of = {"N>T": COL_NGT, "T>N": COL_TGN, "tie": COL_TIE}
    alpha_of = {"N>T": 0.32, "T>N": 0.55, "tie": 0.45}
    order = ["tie", "N>T", "T>N"]  # minority T>N drawn last

    # ======================================================================
    # Panel a — paired slope plot
    # ======================================================================
    x0, x1 = 0.0, 1.0
    half_violin(
        ax_a, pairs["N"].values, x0, side="left", width=0.20,
        color=COL_EST, alpha=0.18, lw=0.6, zorder=1,
    )
    half_violin(
        ax_a, pairs["T"].values, x1, side="right", width=0.20,
        color=COL_EST, alpha=0.18, lw=0.6, zorder=1,
    )

    for direc in order:
        sub = pairs[pairs["direction"] == direc]
        segs = np.stack(
            [
                np.column_stack([np.full(len(sub), x0), sub["N"].to_numpy()]),
                np.column_stack([np.full(len(sub), x1), sub["T"].to_numpy()]),
            ],
            axis=1,
        )
        ax_a.add_collection(
            LineCollection(
                segs, colors=colour_of[direc], linewidths=0.70,
                alpha=alpha_of[direc], zorder=2, capstyle="round",
            )
        )
        pt_kw = dict(
            s=10 if direc == "tie" else 9,
            color=colour_of[direc],
            alpha=0.80 if direc != "tie" else 0.70,
            edgecolors="white" if direc == "tie" else "none",
            linewidths=0.25 if direc == "tie" else 0,
            zorder=3,
        )
        ax_a.scatter(np.full(len(sub), x0), sub["N"], **pt_kw)
        ax_a.scatter(np.full(len(sub), x1), sub["T"], **pt_kw)

    med_n = iqr_median_marker(ax_a, x0 - 0.07, pairs["N"].values)
    med_t = iqr_median_marker(ax_a, x1 + 0.07, pairs["T"].values)
    ax_a.plot(
        [x0 - 0.07, x1 + 0.07], [med_n, med_t],
        color=COL_INK, lw=1.35, zorder=5, solid_capstyle="round",
    )

    ymax = max(pairs["N"].max(), pairs["T"].max())
    # room for the single P-value bracket above the data
    y_br = ymax + 0.55
    h_br = 0.22
    ax_a.plot(
        [x0, x0, x1, x1],
        [y_br, y_br + h_br, y_br + h_br, y_br],
        color=COL_INK, lw=0.7, clip_on=False, solid_capstyle="butt",
        solid_joinstyle="miter",
    )
    ax_a.text(
        0.5, y_br + h_br + 0.10,
        r"$P = 1.4 \times 10^{-18}$",
        ha="center", va="bottom", fontsize=7.5, color=COL_INK, clip_on=False,
    )

    ax_a.set_xlim(-0.48, 1.48)
    ax_a.set_ylim(-0.20, y_br + h_br + 0.85)
    ax_a.set_xticks([x0, x1])
    ax_a.set_xticklabels(["Adjacent\nnormal", "Tumour"])
    ax_a.set_ylabel(r"HP load,  $\ln(1+\mathrm{RPM})$")
    ax_a.yaxis.set_major_locator(MultipleLocator(2))
    ax_a.yaxis.set_minor_locator(MultipleLocator(1))
    ax_a.tick_params(axis="y", which="minor", length=1.8, width=0.5)
    style_ax(ax_a)

    # ======================================================================
    # Panel b — Gardner–Altman estimation plot (tumour − normal)
    # ======================================================================
    deltas = pairs["delta_T_minus_N"].to_numpy()
    x_obs, x_est = 0.0, 1.15
    x_jit = density_beeswarm(deltas, x_centre=x_obs, max_spread=0.13)

    for direc in order:
        mask = pairs["direction"].to_numpy() == direc
        ax_b.scatter(
            x_jit[mask],
            pairs.loc[mask, "delta_T_minus_N"],
            s=11 if direc == "tie" else 10,
            color=colour_of[direc],
            alpha=0.78 if direc != "tie" else 0.80,
            edgecolors="white" if direc == "tie" else "none",
            linewidths=0.25 if direc == "tie" else 0,
            zorder=4 if direc == "tie" else 3,
        )

    half_violin(
        ax_b, boot["boot_median"], x_est, side="right",
        width=0.32, color=COL_EST, alpha=0.35, lw=0.8, zorder=2,
    )
    med, (lo, hi) = boot["median"], boot["median_ci"]
    cap = 0.07
    ax_b.plot([x_est, x_est], [lo, hi], color=COL_EST, lw=2.0,
              solid_capstyle="round", zorder=5)
    ax_b.plot([x_est - cap, x_est + cap], [lo, lo], color=COL_EST, lw=1.4, zorder=5)
    ax_b.plot([x_est - cap, x_est + cap], [hi, hi], color=COL_EST, lw=1.4, zorder=5)
    ax_b.scatter(
        [x_est], [med],
        s=38, marker="o", color=COL_EST, edgecolors="white", linewidths=0.7, zorder=6,
    )

    ax_b.axhline(0.0, color=COL_SPINE, lw=0.6, ls=(0, (3, 2)), zorder=1)

    ax_b.set_xlim(-0.40, 1.68)
    dmin, dmax = float(deltas.min()), float(deltas.max())
    ax_b.set_ylim(dmin - 0.40, dmax + 0.45)
    ax_b.set_xticks([x_obs, x_est])
    ax_b.set_xticklabels(["Observed\npairwise Δ", "Bootstrap\nmedian"])
    ax_b.set_ylabel(r"Δ HP load  (tumour − normal)")
    ax_b.yaxis.set_major_locator(MultipleLocator(2))
    ax_b.yaxis.set_minor_locator(MultipleLocator(1))
    ax_b.tick_params(axis="y", which="minor", length=1.8, width=0.5)
    style_ax(ax_b)

    # ======================================================================
    # Panel c — single slim stacked bar of direction counts
    # ======================================================================
    n_gt, t_gt, ties = stats["n_gt"], stats["t_gt"], stats["ties"]
    total = n_gt + t_gt + ties
    bar_h = 0.46
    left = 0.0
    for val, col in ((n_gt, COL_NGT), (t_gt, COL_TGN), (ties, COL_TIE)):
        ax_c.barh(
            0.0, val, left=left, height=bar_h,
            color=col, edgecolor="white", linewidth=0.4, align="center",
        )
        if val >= 8:
            ax_c.text(
                left + val / 2.0, 0.0,
                f"{val}",
                ha="center", va="center",
                fontsize=7.5, color="white", fontweight="semibold",
            )
        left += val
    ax_c.text(
        total + 3, 0.0, f"$n$ = {total}",
        ha="left", va="center", fontsize=8, color=COL_INK,
    )

    ax_c.set_yticks([])
    ax_c.set_xlabel("Paired patients")
    ax_c.set_xlim(0, total + 28)
    ax_c.set_ylim(-0.55, 0.55)
    ax_c.xaxis.set_major_locator(MultipleLocator(20))
    style_ax(ax_c)
    ax_c.spines["left"].set_visible(False)
    ax_c.tick_params(axis="y", length=0)

    # Shared colour legend, above the panels (outside the violins)
    handles = [
        Patch(facecolor=COL_NGT, edgecolor="none", label="N > T"),
        Patch(facecolor=COL_TGN, edgecolor="none", label="T > N"),
        Patch(facecolor=COL_TIE, edgecolor="none", label="Tie"),
        plt.Line2D(
            [0], [0], color=COL_INK, lw=1.4, marker="o", markersize=4.5,
            markerfacecolor="white", markeredgecolor=COL_INK, markeredgewidth=1.1,
            label="Median + IQR",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.10, 0.995),
        ncol=4, frameon=False,
        handlelength=1.15, handleheight=0.7, handletextpad=0.4,
        columnspacing=1.15, fontsize=7.5, borderaxespad=0.0,
    )

    letter_at(fig, ax_a, "a", dx=-0.13, dy=0.03)
    letter_at(fig, ax_b, "b", dx=-0.15, dy=0.03)
    letter_at(fig, ax_c, "c", dx=-0.055, dy=0.28)

    pdf = OUTDIR / "Fig_A1_draft.pdf"
    svg = OUTDIR / "Fig_A1_draft.svg"
    png = OUTDIR / "Fig_A1_draft.png"
    fig.savefig(pdf)
    fig.savefig(svg)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return pdf, svg, png


def main():
    pairs, stats = load_and_verify()
    boot = bootstrap_delta(pairs["delta_T_minus_N"].to_numpy())
    pdf, svg, png = draw_figure(pairs, stats, boot)
    print("OK: official A1_confirmatory.json fields matched exactly.")
    print(f"pairs: {stats['n_pairs']}  N>T: {stats['n_gt']}  T>N: {stats['t_gt']}  ties: {stats['ties']}")
    print(f"P_one_sided: {stats['p_one']}")
    print(f"r: {stats['r_rb']}  CI: {stats['r_ci']}")
    print(f"median_rpm_N/T: {stats['median_rpm_N']} / {stats['median_rpm_T']}")
    print(f"median_delta: {boot['median']}  CI: {list(boot['median_ci'])}")
    print(f"mean_delta: {boot['mean']}  CI: {list(boot['mean_ci'])}")
    print(f"wrote: {pdf}")
    print(f"wrote: {svg}")
    print(f"wrote: {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

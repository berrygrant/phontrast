#!/usr/bin/env python3
"""Reproduce the final LabPhon 2026 poster analysis figures.

The script uses the processed inputs committed under analysis/labphon_2026/data
and writes regenerated tables/figures under analysis/labphon_2026/outputs.
It intentionally excludes MMO from the overlap-estimator robustness comparison;
the faithful PB52 MMO analysis is handled separately as a talking-point analysis.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, t
from shapely.geometry import Polygon
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import ConvexHull


ANALYSIS_DIR = Path(__file__).resolve().parent
DATA = ANALYSIS_DIR / "data"
OUT = ANALYSIS_DIR / "outputs"

MFCC_TOKENS = DATA / "rerun_mfcc_sampled_tokens_20260624.csv"
MFCC_METRICS = DATA / "rerun_pairwise_metrics_mfcc13_20260624.csv"
PB52_METRICS = DATA / "rerun_pb52_pairwise_20260624.csv"
STANLEY_SUMMARY = DATA / "stanley_jsd_sample_size_summary.csv"
CONSONANT_SUMMARY = DATA / "final_consonant_jsd_generalization.csv"
PB52_MMO_COMPARISON = DATA / "pb52_E_I_jsd_mmo_comparison.csv"
PB52_MMO_DRAWS = DATA / "pb52_E_I_mmo_ba_draws.csv"
PB52_MMO_SUMMARY = DATA / "pb52_E_I_mmo_ba_summary.csv"
PB52_MMO_DIAGNOSTICS = DATA / "pb52_E_I_mmo_diagnostics.csv"
PB52_MMO_MODEL_DATA = DATA / "pb52_E_I_mmo_model_data.csv"

PAIR_VALUES_OUT = OUT / "audited_overlap_estimator_values_no_mmo_20260625.csv"
PERF_OUT = OUT / "audited_overlap_estimator_metric_performance_no_mmo_20260625.csv"
AUDIT_OUT = OUT / "audited_overlap_estimator_method_notes_no_mmo_20260625.csv"
SUMMARY_OUT = OUT / "poster_result_summary_20260625.csv"

MAIN_FIG = OUT / "poster_main_metric_comparison_20260624.png"
MAIN_FIG_PDF = OUT / "poster_main_metric_comparison_20260624.pdf"
MAIN_FIG_SVG = OUT / "poster_main_metric_comparison_20260624.svg"
ROBUST_FIG = OUT / "poster_overlap_estimator_audit_no_mmo_20260625.png"
ROBUST_FIG_PDF = OUT / "poster_overlap_estimator_audit_no_mmo_20260625.pdf"
ROBUST_FIG_SVG = OUT / "poster_overlap_estimator_audit_no_mmo_20260625.svg"
STANLEY_FIG = OUT / "stanley_jsd_sample_size_robustness.png"
STANLEY_FIG_PDF = OUT / "stanley_jsd_sample_size_robustness.pdf"
STANLEY_FIG_SVG = OUT / "stanley_jsd_sample_size_robustness.svg"
CONSONANT_FIG = OUT / "final_consonant_jsd_generalization.png"
CONSONANT_FIG_PDF = OUT / "final_consonant_jsd_generalization.pdf"
CONSONANT_FIG_SVG = OUT / "final_consonant_jsd_generalization.svg"
PB52_MMO_COMPARISON_FIG = OUT / "poster_pb52_jsd_mmo_comparison_20260625.png"
PB52_MMO_COMPARISON_FIG_PDF = OUT / "poster_pb52_jsd_mmo_comparison_20260625.pdf"
PB52_MMO_COMPARISON_FIG_SVG = OUT / "poster_pb52_jsd_mmo_comparison_20260625.svg"
PB52_MMO_FIG = OUT / "poster_pb52_mmo_20260625.png"
PB52_MMO_FIG_PDF = OUT / "poster_pb52_mmo_20260625.pdf"
PB52_MMO_FIG_SVG = OUT / "poster_pb52_mmo_20260625.svg"

FEATURES = [f"mfcc{i}" for i in range(1, 14)]
GLOBAL_SEED = 2026

COLORS = {
    "JSD": "#D62828",
    "Pillai": "#168C84",
    "Bhattacharyya distance": "#405260",
    "Bhattacharyya affinity": "#7A5C9E",
}


def stable_seed(*parts: object) -> int:
    text = "|".join(str(p) for p in parts)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return (int(digest[:12], 16) + GLOBAL_SEED) % (2**32 - 1)


def standardize_pair(X0: np.ndarray, X1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    X = np.vstack([X0, X1])
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    Xz = (X - mu) / sd
    return Xz[: len(X0)], Xz[len(X0) :]


def balanced_pair(tokens: pd.DataFrame, v1: str, v2: str) -> tuple[np.ndarray, np.ndarray]:
    d0 = tokens.loc[tokens["vowel"] == v1, FEATURES].to_numpy(float)
    d1 = tokens.loc[tokens["vowel"] == v2, FEATURES].to_numpy(float)
    n = min(len(d0), len(d1))
    rng = np.random.default_rng(stable_seed(v1, v2, "balance"))
    if len(d0) > n:
        d0 = d0[rng.choice(len(d0), n, replace=False)]
    if len(d1) > n:
        d1 = d1[rng.choice(len(d1), n, replace=False)]
    d0, d1 = standardize_pair(d0, d1)
    X = np.vstack([d0, d1])
    y = np.r_[np.zeros(len(d0), dtype=int), np.ones(len(d1), dtype=int)]
    return X, y


def normalize_log_density(lp: np.ndarray) -> np.ndarray:
    return np.exp(lp - logsumexp(lp))


def knn_density_overlap(X: np.ndarray, y: np.ndarray, k: int = 25) -> float:
    """Shared-mass overlap from kNN density estimates evaluated on pooled points."""
    X0 = X[y == 0]
    X1 = X[y == 1]
    if len(X0) <= k + 2 or len(X1) <= k + 2:
        return np.nan
    d = X.shape[1]

    def log_density(train: np.ndarray, eval_points: np.ndarray) -> np.ndarray:
        nn = NearestNeighbors(n_neighbors=k + 1)
        nn.fit(train)
        dist = nn.kneighbors(eval_points, return_distance=True)[0][:, k]
        dist = np.maximum(dist, np.finfo(float).eps)
        # Unit-ball constants cancel after normalization on the shared support.
        return math.log(k) - math.log(len(train)) - d * np.log(dist)

    p = normalize_log_density(log_density(X0, X))
    q = normalize_log_density(log_density(X1, X))
    return float(np.minimum(p, q).sum())


def mvn_monte_carlo_overlap(
    X: np.ndarray,
    y: np.ndarray,
    n_mc: int = 10000,
    ridge: float = 1e-3,
) -> float:
    """Posterior-overlap integral under one full-covariance Gaussian per category."""
    X0 = X[y == 0]
    X1 = X[y == 1]
    d = X.shape[1]
    if len(X0) <= d + 2 or len(X1) <= d + 2:
        return np.nan
    mu0 = X0.mean(axis=0)
    mu1 = X1.mean(axis=0)
    cov0 = np.cov(X0, rowvar=False) + np.eye(d) * ridge
    cov1 = np.cov(X1, rowvar=False) + np.eye(d) * ridge
    rng = np.random.default_rng(stable_seed(len(X0), len(X1), X0[0, 0], X1[0, 0], "mvn10000"))
    n0 = n_mc // 2
    n1 = n_mc - n0
    try:
        sample = np.vstack(
            [
                rng.multivariate_normal(mu0, cov0, size=n0),
                rng.multivariate_normal(mu1, cov1, size=n1),
            ]
        )
        f0 = multivariate_normal.pdf(sample, mean=mu0, cov=cov0, allow_singular=False)
        f1 = multivariate_normal.pdf(sample, mean=mu1, cov=cov1, allow_singular=False)
    except Exception:
        return np.nan
    denom = f0 + f1
    ok = np.isfinite(denom) & (denom > 0) & np.isfinite(f0) & np.isfinite(f1)
    if ok.sum() < 10:
        return np.nan
    posterior_1 = f1[ok] / denom[ok]
    return float(np.mean(2.0 * np.minimum(posterior_1, 1.0 - posterior_1)))


def pca2_grid_overlap(X: np.ndarray, y: np.ndarray, bins: int = 35) -> float:
    """Shared histogram mass after projecting the pair to two PCA dimensions."""
    if len(y) < 10:
        return np.nan
    Z = PCA(n_components=2, random_state=GLOBAL_SEED).fit_transform(X)
    ranges = []
    for j in range(2):
        lo = float(Z[:, j].min())
        hi = float(Z[:, j].max())
        width = hi - lo
        ranges.append((lo - 0.02 * width - 1e-6, hi + 0.02 * width + 1e-6))
    h0, _, _ = np.histogram2d(Z[y == 0, 0], Z[y == 0, 1], bins=bins, range=ranges)
    h1, _, _ = np.histogram2d(Z[y == 1, 0], Z[y == 1, 1], bins=bins, range=ranges)
    p = h0.ravel().astype(float)
    q = h1.ravel().astype(float)
    if p.sum() <= 0 or q.sum() <= 0:
        return np.nan
    p /= p.sum()
    q /= q.sum()
    return float(np.minimum(p, q).sum())


def pca2_convex_hull_overlap(X: np.ndarray, y: np.ndarray) -> float:
    """Convex-hull support intersection after projecting to two PCA dimensions."""
    if len(y) < 10:
        return np.nan
    Z = PCA(n_components=2, random_state=GLOBAL_SEED).fit_transform(X)
    try:
        z0 = Z[y == 0]
        z1 = Z[y == 1]
        h0 = ConvexHull(z0)
        h1 = ConvexHull(z1)
        p0 = Polygon(z0[h0.vertices])
        p1 = Polygon(z1[h1.vertices])
        if not p0.is_valid or not p1.is_valid:
            return np.nan
        denom = min(p0.area, p1.area)
        if denom <= 0:
            return np.nan
        return float(p0.intersection(p1).area / denom)
    except Exception:
        return np.nan


def compute_overlap_values() -> pd.DataFrame:
    tokens = pd.read_csv(MFCC_TOKENS)
    metrics = pd.read_csv(MFCC_METRICS)
    if "bhatt_affinity" not in metrics.columns:
        metrics["bhatt_affinity"] = np.exp(-metrics["bhatt_dist"])

    rows = []
    for i, row in metrics.iterrows():
        v1, v2 = row["v1"], row["v2"]
        X, y = balanced_pair(tokens, v1, v2)
        rows.append(
            {
                "v1": v1,
                "v2": v2,
                "n_balanced_per_vowel": int(min((y == 0).sum(), (y == 1).sum())),
                "overlap_kde": row["overlap_kde"],
                "overlap_knn_density_k10": knn_density_overlap(X, y, k=10),
                "overlap_knn_density_k25": knn_density_overlap(X, y, k=25),
                "overlap_knn_density_k50": knn_density_overlap(X, y, k=50),
                "overlap_mvn_mc": mvn_monte_carlo_overlap(X, y),
                "overlap_pca2_grid": pca2_grid_overlap(X, y),
                "overlap_pca2_hull": pca2_convex_hull_overlap(X, y),
            }
        )
        if (i + 1) % 50 == 0:
            print(f"computed alternative overlap for {i + 1}/{len(metrics)} pairs", flush=True)

    alt = pd.DataFrame(rows)
    pair_df = metrics.merge(alt, on=["v1", "v2", "overlap_kde"], how="left")
    pair_df.to_csv(PAIR_VALUES_OUT, index=False)
    return pair_df


def correlation_table(pair_df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = {
        "JSD": "jsd",
        "Pillai": "pillai",
        "Bhattacharyya distance": "bhatt_dist",
        "Bhattacharyya affinity": "bhatt_affinity",
    }
    overlap_cols = {
        "KDE density overlap (13D)": "overlap_kde",
        "kNN density overlap k=10 (13D)": "overlap_knn_density_k10",
        "kNN density overlap k=25 (13D)": "overlap_knn_density_k25",
        "kNN density overlap k=50 (13D)": "overlap_knn_density_k50",
        "MVN Monte Carlo overlap (13D)": "overlap_mvn_mc",
        "PCA-2D grid overlap": "overlap_pca2_grid",
        "PCA-2D convex hull support": "overlap_pca2_hull",
    }
    rows = []
    for overlap_name, overlap_col in overlap_cols.items():
        for metric_name, metric_col in metric_cols.items():
            x = pair_df[metric_col]
            y = pair_df[overlap_col]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 3:
                x_ok = x[ok].to_numpy(float)
                y_ok = y[ok].to_numpy(float)
                pearson = float(pd.Series(x_ok).corr(pd.Series(y_ok), method="pearson"))
                spearman = float(pd.Series(x_ok).corr(pd.Series(y_ok), method="spearman"))
                boot_low, boot_high = bootstrap_abs_pearson_ci(
                    x_ok,
                    y_ok,
                    seed=stable_seed(overlap_name, metric_name, "abs_pearson_boot"),
                )
            else:
                pearson = np.nan
                spearman = np.nan
                boot_low = np.nan
                boot_high = np.nan
            rows.append(
                {
                    "overlap_estimator": overlap_name,
                    "metric": metric_name,
                    "n_pairs": int(ok.sum()),
                    "pearson_r": pearson,
                    "abs_pearson_r": abs(pearson) if math.isfinite(pearson) else np.nan,
                    "abs_pearson_boot_low": boot_low,
                    "abs_pearson_boot_high": boot_high,
                    "spearman_rho": spearman,
                    "abs_spearman_rho": abs(spearman) if math.isfinite(spearman) else np.nan,
                }
            )
    perf = pd.DataFrame(rows)
    perf.to_csv(PERF_OUT, index=False)
    return perf


def bootstrap_abs_pearson_ci(
    x: np.ndarray,
    y: np.ndarray,
    seed: int,
    n_boot: int = 3000,
) -> tuple[float, float]:
    """Bootstrap 95% interval for |Pearson r| across pairwise contrasts."""
    n = len(x)
    if n < 4:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        xb = x[idx]
        yb = y[idx]
        if np.std(xb) <= 0 or np.std(yb) <= 0:
            vals[i] = np.nan
        else:
            vals[i] = abs(float(np.corrcoef(xb, yb)[0, 1]))
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    return float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def write_method_notes() -> None:
    notes = pd.DataFrame(
        [
            {
                "overlap_estimator": "KDE density overlap (13D)",
                "definition": "Shared mass from the same KDE framework used for the main overlap benchmark.",
                "status": "primary but circular with KDE-based JSD",
            },
            {
                "overlap_estimator": "kNN density overlap k=10/25/50 (13D)",
                "definition": "Shared mass from k-nearest-neighbor class-conditional density estimates on pooled points.",
                "status": "independent density-family sensitivity check",
            },
            {
                "overlap_estimator": "MVN Monte Carlo overlap (13D)",
                "definition": "Posterior overlap under one full-covariance multivariate normal per category; 10,000 Monte Carlo draws.",
                "status": "parametric Gaussian reference; close to Bhattacharyya assumptions",
            },
            {
                "overlap_estimator": "PCA-2D grid overlap",
                "definition": "Shared histogram mass after pairwise projection to two principal components.",
                "status": "low-dimensional projection diagnostic",
            },
            {
                "overlap_estimator": "PCA-2D convex hull support",
                "definition": "Intersection area of PCA-2D convex-hull supports divided by smaller hull area.",
                "status": "support-overlap diagnostic, not a density-overlap estimate",
            },
        ]
    )
    notes.to_csv(AUDIT_OUT, index=False)


def _fit_line_with_ci(ax, x: np.ndarray, y: np.ndarray, color: str) -> None:
    if len(x) < 2:
        return
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
    y_hat = slope * xs + intercept
    ax.plot(xs, y_hat, color=color, linewidth=2.4)
    n = len(x)
    if n <= 2:
        return
    y_fit = slope * x + intercept
    resid = y - y_fit
    x_mean = x.mean()
    sxx = np.sum((x - x_mean) ** 2)
    if sxx <= 0:
        return
    mse = np.sum(resid ** 2) / (n - 2)
    se_mean = np.sqrt(mse * (1.0 / n + ((xs - x_mean) ** 2) / sxx))
    crit = t.ppf(0.975, df=n - 2)
    ax.fill_between(
        xs,
        y_hat - crit * se_mean,
        y_hat + crit * se_mean,
        color=color,
        alpha=0.16,
        linewidth=0,
    )


def make_main_metric_figure() -> None:
    mfcc = pd.read_csv(MFCC_METRICS)
    pb52 = pd.read_csv(PB52_METRICS)

    datasets = [
        ("PB52 F1/F2", pb52, "overlap", {"Bhatt. distance": "bhatt", "JSD": "jsd", "Pillai": "pillai"}),
        (
            "SBCSAE 13-MFCC",
            mfcc,
            "overlap_kde",
            {"Bhatt. distance": "bhatt_dist", "JSD": "jsd", "Pillai": "pillai"},
        ),
    ]
    display_to_color = {
        "Bhatt. distance": COLORS["Bhattacharyya distance"],
        "JSD": COLORS["JSD"],
        "Pillai": COLORS["Pillai"],
    }

    fig, axes = plt.subplots(2, 3, figsize=(16.8, 11.8), dpi=220)
    fig.patch.set_facecolor("#F6F9FC")

    for row_i, (dataset_label, df, overlap_col, specs) in enumerate(datasets):
        for col_i, (label, metric_col) in enumerate(specs.items()):
            ax = axes[row_i, col_i]
            ax.set_facecolor("white")
            x = df[metric_col].to_numpy(float)
            y = df[overlap_col].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            x = x[ok]
            y = y[ok]
            color = display_to_color[label]
            ax.scatter(x, y, s=23 if row_i == 0 else 18, alpha=0.66, color=color, edgecolors="none")
            _fit_line_with_ci(ax, x, y, color)
            r = pd.Series(x).corr(pd.Series(y), method="pearson")
            rho = pd.Series(x).corr(pd.Series(y), method="spearman")
            ax.text(
                0.04,
                0.95,
                f"r = {r:.2f}\n\u03c1 = {rho:.2f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=14,
                color="#2E3A4A",
                bbox=dict(boxstyle="round,pad=0.28", fc="#F6F9FC", ec="#C9D6E2"),
            )
            if row_i == 0:
                ax.set_title(label, fontsize=18, fontweight="bold", color="#12213A", pad=10)
            if col_i == 0:
                ax.set_ylabel(f"{dataset_label}\nKDE overlap", fontsize=17, fontweight="bold", color="#12213A")
            else:
                ax.set_ylabel("")
            ax.set_xlabel("Metric value", fontsize=15, fontweight="bold", color="#12213A")
            ax.grid(color="#DFE7F0", linewidth=0.9, alpha=0.75)
            ax.tick_params(labelsize=12, colors="#12213A")
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=1.2)
    fig.savefig(MAIN_FIG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(MAIN_FIG_PDF, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(MAIN_FIG_SVG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def make_robustness_figure(perf: pd.DataFrame) -> None:
    display_order = [
        "KDE density overlap (13D)",
        "kNN density overlap k=25 (13D)",
        "MVN Monte Carlo overlap (13D)",
        "PCA-2D grid overlap",
        "PCA-2D convex hull support",
    ]
    xlabels = ["KDE\n13D", "kNN density\n13D", "MVN MC\n13D", "PCA grid\n2D", "PCA hull\n2D"]
    metrics = ["JSD", "Pillai", "Bhattacharyya distance"]
    label_map = {"Bhattacharyya distance": "Bhatt. distance"}
    markers = {"JSD": "o", "Pillai": "s", "Bhattacharyya distance": "^"}
    x_offsets = {"JSD": -0.18, "Pillai": 0.00, "Bhattacharyya distance": 0.18}

    fig, ax = plt.subplots(figsize=(17.2, 7.8), dpi=220)
    fig.patch.set_facecolor("#F6F9FC")
    ax.set_facecolor("white")
    xs = np.arange(len(display_order))
    for metric in metrics:
        sub = (
            perf[(perf["metric"] == metric) & (perf["overlap_estimator"].isin(display_order))]
            .set_index("overlap_estimator")
            .reindex(display_order)
        )
        ys = sub["abs_pearson_r"].to_numpy(float)
        lows = sub["abs_pearson_boot_low"].to_numpy(float)
        highs = sub["abs_pearson_boot_high"].to_numpy(float)
        yerr = np.vstack([np.maximum(0, ys - lows), np.maximum(0, highs - ys)])
        color = COLORS[metric]
        metric_xs = xs + x_offsets[metric]
        ax.errorbar(
            metric_xs,
            ys,
            yerr=yerr,
            marker=markers[metric],
            linewidth=3.2,
            markersize=9.5,
            color=color,
            ecolor=color,
            elinewidth=1.7,
            capsize=4.5,
            capthick=1.7,
            label=label_map.get(metric, metric),
        )
        for xi, yi in zip(metric_xs, ys):
            if np.isfinite(yi):
                va = "bottom"
                offset = 0.024
                if metric == "Bhattacharyya distance":
                    va = "top"
                    offset = -0.034
                ax.text(xi, yi + offset, f"{yi:.2f}", ha="center", va=va, fontsize=11.8, fontweight="bold", color=color)

    ax.set_ylim(0.58, 1.04)
    ax.set_xlim(-0.35, len(display_order) - 0.65)
    ax.set_xticks(xs)
    ax.set_xticklabels(xlabels, fontsize=13)
    ax.set_ylabel("|Pearson r| with reference overlap", fontsize=18, fontweight="bold", color="#12213A")
    ax.set_xlabel("Reference overlap estimator", fontsize=18, fontweight="bold", color="#12213A")
    ax.grid(axis="y", color="#D8E2EE", linewidth=1.0)
    ax.tick_params(axis="y", labelsize=13, colors="#12213A")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02), ncol=3, frameon=False, fontsize=14)
    ax.text(
        0.995,
        0.02,
        "KDE shares density-estimation machinery with JSD; parametric references favor Gaussian/Bhattacharyya assumptions",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11.5,
        color="#526173",
    )
    fig.tight_layout(pad=1.0)
    fig.savefig(ROBUST_FIG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(ROBUST_FIG_PDF, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(ROBUST_FIG_SVG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def make_sample_size_figure() -> None:
    df = pd.read_csv(STANLEY_SUMMARY)
    fig, ax = plt.subplots(figsize=(10.8, 5.8), dpi=220)
    fig.patch.set_facecolor("#F6F9FC")
    ax.set_facecolor("white")

    x = df["n_per_category"].to_numpy(float)
    mean = df["mean"].to_numpy(float)
    low = df["p05"].to_numpy(float)
    high = df["p95"].to_numpy(float)
    ax.fill_between(x, low, high, color=COLORS["JSD"], alpha=0.17, linewidth=0, label="95% simulation interval")
    ax.plot(x, mean, color=COLORS["JSD"], linewidth=3.0, label="Mean JSD under complete overlap")
    ax.axvline(25, color="#526173", linewidth=1.8, linestyle="--")
    n25 = df.loc[df["n_per_category"] == 25].iloc[0]
    ax.scatter([25], [n25["mean"]], color=COLORS["JSD"], s=70, zorder=3)
    ax.text(
        25.8,
        n25["mean"] + 0.035,
        f"n=25 mean={n25['mean']:.2f}\n95% sim int [{n25['p05']:.2f}, {n25['p95']:.2f}]",
        ha="left",
        va="bottom",
        fontsize=11.5,
        color="#12213A",
        bbox={"boxstyle": "round,pad=0.30", "facecolor": "white", "edgecolor": "#C9D6E2"},
    )
    ax.set_xlim(5, 100)
    ax.set_ylim(0, max(0.62, float(high.max()) + 0.04))
    ax.set_xlabel("Tokens per category", fontsize=15, fontweight="bold", color="#12213A")
    ax.set_ylabel("JSD estimate under complete overlap", fontsize=15, fontweight="bold", color="#12213A")
    ax.grid(color="#DFE7F0", linewidth=0.9)
    ax.tick_params(labelsize=12, colors="#12213A")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=11.5)
    fig.tight_layout()
    fig.savefig(STANLEY_FIG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(STANLEY_FIG_PDF, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(STANLEY_FIG_SVG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def make_consonant_figure() -> None:
    df = pd.read_csv(CONSONANT_SUMMARY)
    df["display"] = df["domain"] + "\n" + df["contrast"]
    fig, ax = plt.subplots(figsize=(11.8, 5.8), dpi=220)
    fig.patch.set_facecolor("#F6F9FC")
    ax.set_facecolor("white")

    colors = {"Stop voicing": "#2D6A8E", "Fricative place": "#D62828"}
    xs = np.arange(len(df))
    for i, row in df.iterrows():
        color = colors.get(row["domain"], "#405260")
        y = float(row["jsd_point"])
        low = float(row["jsd_interval_low"])
        high = float(row["jsd_interval_high"])
        ax.errorbar(
            i,
            y,
            yerr=[[max(0, y - low)], [max(0, high - y)]],
            marker="o",
            markersize=9,
            linewidth=2.4,
            capsize=6,
            color=color,
            ecolor=color,
        )
        ax.text(i, min(1.03, high + 0.04), f"{y:.2f}", ha="center", va="bottom", fontsize=11.5, fontweight="bold", color=color)

    ax.set_xticks(xs)
    ax.set_xticklabels(df["display"], fontsize=11.5)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("JSD separation", fontsize=15, fontweight="bold", color="#12213A")
    ax.set_xlabel("Exploratory SBCSAE contrast", fontsize=15, fontweight="bold", color="#12213A")
    ax.grid(axis="y", color="#DFE7F0", linewidth=0.9)
    ax.tick_params(axis="y", labelsize=12, colors="#12213A")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(
        0.01,
        0.98,
        "VOT stop voicing remains low; fricative place is high in the specified spectral feature space",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11.5,
        color="#526173",
    )
    fig.tight_layout()
    fig.savefig(CONSONANT_FIG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(CONSONANT_FIG_PDF, facecolor=fig.get_facecolor(), bbox_inches="tight")
    fig.savefig(CONSONANT_FIG_SVG, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def _draw_mmo_comparison_panel(ax, df: pd.DataFrame, title: str, ylabel: str, ylim: tuple[float, float]) -> None:
    panel_colors = {
        "JSD": "#D7263D",
        "1 - MMO BA": "#2D6A8E",
        "KDE overlap": "#D7263D",
        "MMO BA": "#2D6A8E",
    }
    ax.set_facecolor("white")
    xs = np.arange(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        color = panel_colors[row["measure"]]
        y = float(row["point"])
        ax.scatter(i, y, s=170, color=color, zorder=3)
        if np.isfinite(row["low"]) and np.isfinite(row["high"]):
            ax.errorbar(
                i,
                y,
                yerr=[[y - row["low"]], [row["high"] - y]],
                color=color,
                linewidth=2.8,
                capsize=7,
                capthick=2.2,
                zorder=2,
            )
        label = f"{y:.3f}" if row["measure"] == "KDE overlap" else f"{y:.2f}"
        label_y = row["high"] + 0.025 if np.isfinite(row["high"]) else y + 0.055
        ax.text(i, label_y, label, ha="center", va="bottom", fontsize=13.5, fontweight="bold", color=color)
    ax.set_xticks(xs)
    ax.set_xticklabels(df["measure"], fontsize=12.5, fontweight="bold")
    ax.set_xlim(-0.55, len(df) - 0.45)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=14, fontweight="bold", color="#12213A")
    ax.grid(axis="y", color="#DFE7F0", linewidth=1.0)
    ax.tick_params(axis="y", labelsize=11.5, colors="#12213A")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold", color="#12213A", pad=10)


def make_pb52_jsd_mmo_comparison_figure() -> None:
    df = pd.read_csv(PB52_MMO_COMPARISON)
    sep = df[df["estimand"] == "separation"].copy()
    overlap = df[df["estimand"] == "overlap"].copy()
    sep["order"] = sep["measure"].map({"JSD": 0, "1 - MMO BA": 1})
    overlap["order"] = overlap["measure"].map({"KDE overlap": 0, "MMO BA": 1})
    sep = sep.sort_values("order")
    overlap = overlap.sort_values("order")

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.7), dpi=220, gridspec_kw={"wspace": 0.34})
    fig.patch.set_facecolor("#F6F9FC")
    _draw_mmo_comparison_panel(axes[0], sep, "Same PB52 /E/-/I/ contrast", "Separation scale", (0, 1.02))
    axes[0].text(
        0.5,
        0.08,
        "JSD: 0 = identical, 1 = maximally distinct\nMMO shown as 1 - Bhattacharyya affinity",
        transform=axes[0].transAxes,
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="#526173",
    )
    _draw_mmo_comparison_panel(axes[1], overlap, "Overlap estimates", "Overlap scale", (0, 1.02))
    fig.text(
        0.5,
        0.02,
        "Intervals: JSD uses a 95% nonparametric bootstrap CI; MMO uses a 95% posterior credible interval.",
        ha="center",
        va="bottom",
        fontsize=11.3,
        color="#526173",
    )
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.22, top=0.88)
    fig.savefig(PB52_MMO_COMPARISON_FIG, facecolor=fig.get_facecolor())
    fig.savefig(PB52_MMO_COMPARISON_FIG_PDF, facecolor=fig.get_facecolor())
    fig.savefig(PB52_MMO_COMPARISON_FIG_SVG, facecolor=fig.get_facecolor())
    plt.close(fig)


def covariance_ellipse(x: np.ndarray, y: np.ndarray, level: float = 0.68) -> Ellipse:
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    radius = np.sqrt(-2 * np.log(1 - level))
    width, height = 2 * radius * np.sqrt(vals)
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    return Ellipse((np.mean(x), np.mean(y)), width=width, height=height, angle=angle)


def make_pb52_mmo_figure() -> None:
    df = pd.read_csv(PB52_MMO_MODEL_DATA)
    ba = pd.read_csv(PB52_MMO_DRAWS)["bhattacharyya_affinity"].to_numpy(float)
    summary = pd.read_csv(PB52_MMO_SUMMARY).iloc[0]
    diag = pd.read_csv(PB52_MMO_DIAGNOSTICS).iloc[0]
    vowel_colors = {"/ɛ/": "#D7263D", "/ɪ/": "#15928A", "/E/": "#D7263D", "/I/": "#15928A"}

    fig, axes = plt.subplots(1, 2, figsize=(15.4, 6.2), dpi=220, gridspec_kw={"width_ratios": [1.12, 1.0], "wspace": 0.28})
    fig.patch.set_facecolor("#F6F9FC")
    ax = axes[0]
    ax.set_facecolor("white")
    for label, sub in df.groupby("vowel_ipa"):
        color = vowel_colors.get(label, "#405260")
        ax.scatter(sub["f2_z"], sub["f1_z"], s=26, alpha=0.62, color=color, edgecolor="white", linewidth=0.3, label=label)
        ell = covariance_ellipse(sub["f2_z"].to_numpy(), sub["f1_z"].to_numpy())
        ell.set_facecolor("none")
        ell.set_edgecolor(color)
        ell.set_linewidth(2.4)
        ax.add_patch(ell)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.grid(color="#DFE7F0", linewidth=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("F2, speaker-normalized", fontsize=15, fontweight="bold", color="#12213A")
    ax.set_ylabel("F1, speaker-normalized", fontsize=15, fontweight="bold", color="#12213A")
    ax.tick_params(labelsize=11.5, colors="#12213A")
    ax.legend(title="Vowel", frameon=False, fontsize=13, title_fontsize=13, loc="lower left")
    ax.text(0.02, 0.98, "PB52 /E/-/I/ tokens", transform=ax.transAxes, va="top", ha="left", fontsize=16, fontweight="bold", color="#12213A")

    ax = axes[1]
    ax.set_facecolor("white")
    ax.hist(ba, bins=42, color="#2D6A8E", edgecolor="white", linewidth=0.4)
    ax.axvline(summary["ba_median"], color="#D7263D", linewidth=3.2)
    ax.axvspan(summary["ba_q025"], summary["ba_q975"], color="#D7263D", alpha=0.12)
    ax.set_xlim(0, 1)
    ax.grid(color="#DFE7F0", linewidth=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Modelled overlap (Bhattacharyya affinity)", fontsize=15, fontweight="bold", color="#12213A")
    ax.set_ylabel("Posterior draws", fontsize=15, fontweight="bold", color="#12213A")
    ax.tick_params(labelsize=11.5, colors="#12213A")
    stat_text = (
        f"Median = {summary['ba_median']:.2f}\n"
        f"95% CrI [{summary['ba_q025']:.2f}, {summary['ba_q975']:.2f}]\n"
        f"Rhat max = {diag['max_rhat']:.3f}; divergences = {int(diag['divergences'])}"
    )
    ax.text(
        0.98,
        0.83,
        stat_text,
        transform=ax.transAxes,
        va="top",
        ha="right",
        fontsize=13.2,
        color="#12213A",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#B8C7D9"},
    )
    ax.text(0.02, 0.98, "MMO posterior overlap", transform=ax.transAxes, va="top", ha="left", fontsize=16, fontweight="bold", color="#12213A")
    fig.text(
        0.5,
        0.02,
        "Bayesian multivariate mixed-effects model: F1/F2 ~ vowel + talker type + repetition + (1 + vowel | speaker).",
        ha="center",
        va="bottom",
        fontsize=11.5,
        color="#526173",
    )
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.15, top=0.95)
    fig.savefig(PB52_MMO_FIG, facecolor=fig.get_facecolor())
    fig.savefig(PB52_MMO_FIG_PDF, facecolor=fig.get_facecolor())
    fig.savefig(PB52_MMO_FIG_SVG, facecolor=fig.get_facecolor())
    plt.close(fig)


def write_result_summary(perf: pd.DataFrame) -> None:
    pb52 = pd.read_csv(PB52_METRICS)
    mfcc = pd.read_csv(MFCC_METRICS)
    sample = pd.read_csv(STANLEY_SUMMARY)
    consonants = pd.read_csv(CONSONANT_SUMMARY)
    mmo = pd.read_csv(PB52_MMO_COMPARISON)

    rows: list[dict[str, object]] = []
    metric_specs = [
        ("PB52 F1/F2", pb52, "overlap", {"JSD": "jsd", "Pillai": "pillai", "Bhattacharyya distance": "bhatt"}),
        ("SBCSAE 13-MFCC", mfcc, "overlap_kde", {"JSD": "jsd", "Pillai": "pillai", "Bhattacharyya distance": "bhatt_dist"}),
    ]
    for dataset, df, overlap_col, specs in metric_specs:
        for metric, col in specs.items():
            x = df[col].to_numpy(float)
            y = df[overlap_col].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            rows.append(
                {
                    "section": "main_metric_comparison",
                    "dataset": dataset,
                    "measure": metric,
                    "statistic": "Pearson r with KDE overlap",
                    "value": float(pd.Series(x[ok]).corr(pd.Series(y[ok]), method="pearson")),
                    "n": int(ok.sum()),
                }
            )
            rows.append(
                {
                    "section": "main_metric_comparison",
                    "dataset": dataset,
                    "measure": metric,
                    "statistic": "Spearman rho with KDE overlap",
                    "value": float(pd.Series(x[ok]).corr(pd.Series(y[ok]), method="spearman")),
                    "n": int(ok.sum()),
                }
            )

    for _, row in perf.iterrows():
        if row["metric"] in {"JSD", "Pillai", "Bhattacharyya distance"}:
            rows.append(
                {
                    "section": "overlap_estimator_robustness",
                    "dataset": "SBCSAE 13-MFCC",
                    "measure": row["metric"],
                    "statistic": f"abs Pearson r with {row['overlap_estimator']}",
                    "value": row["abs_pearson_r"],
                    "n": row["n_pairs"],
                }
            )

    for n in [10, 20, 25]:
        row = sample.loc[sample["n_per_category"] == n].iloc[0]
        rows.append(
            {
                "section": "sample_size_simulation",
                "dataset": "simulated complete overlap",
                "measure": "JSD",
                "statistic": f"n={n} mean JSD",
                "value": row["mean"],
                "n": row["reps"],
            }
        )
        rows.append(
            {
                "section": "sample_size_simulation",
                "dataset": "simulated complete overlap",
                "measure": "JSD",
                "statistic": f"n={n} 95% simulation interval",
                "value": f"[{row['p05']:.6f}, {row['p95']:.6f}]",
                "n": row["reps"],
            }
        )

    for _, row in consonants.iterrows():
        rows.append(
            {
                "section": "contrast_generalization",
                "dataset": "SBCSAE exploratory consonants",
                "measure": row["domain"],
                "statistic": row["contrast"],
                "value": row["jsd_point"],
                "n": row["n_label"],
            }
        )

    for _, row in mmo.iterrows():
        rows.append(
            {
                "section": "pb52_mmo_talking_point",
                "dataset": "PB52 /E/-/I/",
                "measure": row["measure"],
                "statistic": row["interval"],
                "value": row["point"],
                "n": row["n_tokens"],
            }
        )

    pd.DataFrame(rows).to_csv(SUMMARY_OUT, index=False)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    pair_df = compute_overlap_values()
    perf = correlation_table(pair_df)
    write_method_notes()
    make_main_metric_figure()
    make_robustness_figure(perf)
    make_sample_size_figure()
    make_consonant_figure()
    make_pb52_jsd_mmo_comparison_figure()
    make_pb52_mmo_figure()
    write_result_summary(perf)
    print(PAIR_VALUES_OUT)
    print(PERF_OUT)
    print(AUDIT_OUT)
    print(SUMMARY_OUT)
    print(MAIN_FIG)
    print(ROBUST_FIG)
    print(STANLEY_FIG)
    print(CONSONANT_FIG)
    print(PB52_MMO_COMPARISON_FIG)
    print(PB52_MMO_FIG)


if __name__ == "__main__":
    main()

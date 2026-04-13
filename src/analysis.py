"""
Connectivity-tax benchmark — analysis and plotting.

Produces the five main plots:
1. Metric overview vs alpha
2. Connectivity tax ratio curves
3. Heat maps (2Q depth: unitary vs dynamic)
4. Phase diagram (feasibility regions)
5. N-dependence comparison
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

from .metrics import CircuitMetrics


# ---------------------------------------------------------------------------
# Data wrangling
# ---------------------------------------------------------------------------

def metrics_to_dataframe(metrics: list[CircuitMetrics]) -> pd.DataFrame:
    """Convert list of CircuitMetrics to a pandas DataFrame."""
    return pd.DataFrame([m.to_dict() for m in metrics])


def compute_tax_ratios(df: pd.DataFrame) -> pd.DataFrame:
    """Add connectivity-tax ratio columns to the DataFrame.

    Ratios are computed relative to the 'ideal' variant for each
    (n_qubits, alpha, n_steps, theta) group.
    """
    group_cols = ["n_qubits", "alpha", "n_steps", "theta"]
    merged = df.copy()

    ideal = df[df["variant"] == "ideal"].set_index(group_cols)
    for col in ["two_q_depth", "two_q_count", "total_depth", "cx_count"]:
        ratio_col = f"{col}_tax"
        ideal_map = ideal[col].to_dict()
        merged[ratio_col] = merged.apply(
            lambda row: (
                row[col] / ideal_map.get(
                    tuple(row[c] for c in group_cols), np.nan
                )
                if ideal_map.get(tuple(row[c] for c in group_cols), 0) > 0
                else np.nan
            ),
            axis=1,
        )
    return merged


# ---------------------------------------------------------------------------
# Plot 1: Metric overview vs alpha
# ---------------------------------------------------------------------------

def plot_metrics_vs_alpha(
    df: pd.DataFrame,
    n_qubits: int = 6,
    n_steps: int = 1,
    theta: float | None = None,
    metrics: list[str] | None = None,
    figsize: tuple = (14, 10),
):
    """Line plots of circuit metrics vs alpha for all three variants.

    Parameters
    ----------
    df : DataFrame with columns from CircuitMetrics + tax ratios.
    n_qubits : filter to this system size.
    n_steps : filter to this Trotter step count.
    theta : filter to this theta value (if None, uses the first available).
    metrics : list of column names to plot.  Default: 4 key metrics.
    """
    if metrics is None:
        metrics = ["two_q_depth", "two_q_count", "swap_count", "total_depth"]

    sub = df[
        (df["n_qubits"] == n_qubits)
        & (df["n_steps"] == n_steps)
    ]
    if theta is not None:
        sub = sub[sub["theta"] == theta]
    else:
        theta_val = sub["theta"].iloc[0] if len(sub) > 0 else 0
        sub = sub[sub["theta"] == theta_val]

    fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True)
    axes = axes.ravel()

    variant_styles = {
        "ideal":   {"color": "tab:green",  "marker": "o", "ls": "-"},
        "unitary": {"color": "tab:red",    "marker": "s", "ls": "--"},
        "dynamic": {"color": "tab:blue",   "marker": "^", "ls": "-."},
    }

    for ax, metric in zip(axes, metrics):
        for variant, style in variant_styles.items():
            data = sub[sub["variant"] == variant].sort_values("alpha")
            if len(data) == 0:
                continue
            ax.plot(data["alpha"], data[metric], label=variant, **style, markersize=7)
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    for ax in axes[-2:]:
        ax.set_xlabel(r"$\alpha$")

    fig.suptitle(
        f"Circuit metrics vs power-law exponent  (N={n_qubits}, steps={n_steps})",
        fontsize=13,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 2: Connectivity tax ratios
# ---------------------------------------------------------------------------

def plot_tax_ratios(
    df: pd.DataFrame,
    n_qubits: int = 6,
    n_steps: int = 1,
    theta: float | None = None,
    metric: str = "two_q_depth",
    figsize: tuple = (8, 5),
):
    """Plot tax ratio (variant / ideal) vs alpha."""
    tax_col = f"{metric}_tax"
    if tax_col not in df.columns:
        df = compute_tax_ratios(df)

    sub = df[(df["n_qubits"] == n_qubits) & (df["n_steps"] == n_steps)]
    if theta is not None:
        sub = sub[sub["theta"] == theta]
    else:
        theta_val = sub["theta"].iloc[0] if len(sub) > 0 else 0
        sub = sub[sub["theta"] == theta_val]

    fig, ax = plt.subplots(figsize=figsize)

    for variant, color in [("unitary", "tab:red"), ("dynamic", "tab:blue")]:
        data = sub[sub["variant"] == variant].sort_values("alpha")
        if len(data) == 0:
            continue
        ax.plot(data["alpha"], data[tax_col], "o-", color=color, label=variant, markersize=7)

    ax.axhline(1.0, color="tab:green", ls="--", alpha=0.7, label="ideal (ratio=1)")
    ax.set_xlabel(r"$\alpha$", fontsize=12)
    ax.set_ylabel(f"{metric.replace('_', ' ').title()} — Tax Ratio", fontsize=12)
    ax.set_title(f"Connectivity Tax  (N={n_qubits}, steps={n_steps})", fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 3: Heat maps (alpha vs Trotter steps)
# ---------------------------------------------------------------------------

def plot_heatmaps(
    df: pd.DataFrame,
    n_qubits: int = 6,
    theta: float | None = None,
    metric: str = "two_q_depth",
    t2_us: float = 200.0,
    figsize: tuple = (14, 5),
):
    """Side-by-side heat maps: unitary (left) vs dynamic (right).

    Color = metric value.  Axes: alpha (y) vs Trotter steps (x).
    Optional contour at estimated duration = T2.
    """
    sub = df[df["n_qubits"] == n_qubits]
    if theta is not None:
        sub = sub[sub["theta"] == theta]
    else:
        theta_val = sub["theta"].iloc[0] if len(sub) > 0 else 0
        sub = sub[sub["theta"] == theta_val]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for ax, variant, title in [
        (ax1, "unitary", "Heavy-hex Unitary"),
        (ax2, "dynamic", "Heavy-hex Dynamic"),
    ]:
        data = sub[sub["variant"] == variant]
        if len(data) == 0:
            ax.set_title(f"{title} — no data")
            continue

        pivot = data.pivot_table(
            index="alpha", columns="n_steps", values=metric, aggfunc="mean"
        )
        alphas = pivot.index.values
        steps = pivot.columns.values
        Z = pivot.values

        im = ax.imshow(
            Z,
            aspect="auto",
            origin="lower",
            extent=[steps.min() - 0.5, steps.max() + 0.5,
                    alphas.min(), alphas.max()],
            cmap="YlOrRd",
        )
        ax.set_xlabel("Trotter Steps")
        ax.set_title(title)
        plt.colorbar(im, ax=ax, label=metric.replace("_", " ").title())

    ax1.set_ylabel(r"$\alpha$")
    fig.suptitle(
        f"{metric.replace('_', ' ').title()}  (N={n_qubits})",
        fontsize=13,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 4: Phase diagram (feasibility regions)
# ---------------------------------------------------------------------------

def plot_phase_diagram(
    df: pd.DataFrame,
    n_qubits: int = 6,
    theta: float | None = None,
    metric: str = "two_q_depth",
    threshold: float | None = None,
    figsize: tuple = (8, 6),
):
    """Phase diagram: green / yellow / red regions.

    - Green:  both unitary and dynamic below threshold
    - Yellow: unitary above, dynamic below threshold
    - Red:    both above threshold

    If threshold is None, uses a heuristic: 3x the ideal metric value
    as a rough "hardware-feasible" cutoff.
    """
    sub = df[df["n_qubits"] == n_qubits]
    if theta is not None:
        sub = sub[sub["theta"] == theta]
    else:
        theta_val = sub["theta"].iloc[0] if len(sub) > 0 else 0
        sub = sub[sub["theta"] == theta_val]

    # Determine threshold from ideal if not given
    if threshold is None:
        ideal_max = sub[sub["variant"] == "ideal"][metric].max()
        threshold = max(3 * ideal_max, 10)

    alphas = np.sort(sub["alpha"].unique())
    steps = np.sort(sub["n_steps"].unique())

    phase = np.full((len(alphas), len(steps)), np.nan)

    for ai, a in enumerate(alphas):
        for si, s in enumerate(steps):
            u_val = sub[
                (sub["variant"] == "unitary") & (sub["alpha"] == a) & (sub["n_steps"] == s)
            ][metric]
            d_val = sub[
                (sub["variant"] == "dynamic") & (sub["alpha"] == a) & (sub["n_steps"] == s)
            ][metric]

            u = u_val.values[0] if len(u_val) > 0 else np.nan
            d = d_val.values[0] if len(d_val) > 0 else np.nan

            if np.isnan(u) or np.isnan(d):
                phase[ai, si] = np.nan
            elif u <= threshold and d <= threshold:
                phase[ai, si] = 0  # green
            elif (u > threshold and d <= threshold) or (u <= threshold and d > threshold):
                phase[ai, si] = 1  # yellow
            else:
                phase[ai, si] = 2  # red

    cmap = mcolors.ListedColormap(["#4CAF50", "#FFC107", "#F44336"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        phase,
        aspect="auto",
        origin="lower",
        extent=[steps.min() - 0.5, steps.max() + 0.5,
                alphas.min(), alphas.max()],
        cmap=cmap,
        norm=norm,
    )

    legend_elements = [
        Patch(facecolor="#4CAF50", label="Both feasible"),
        Patch(facecolor="#FFC107", label="Only dynamic feasible"),
        Patch(facecolor="#F44336", label="Neither feasible (need all-to-all)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9)
    ax.set_xlabel("Trotter Steps", fontsize=12)
    ax.set_ylabel(r"$\alpha$", fontsize=12)
    ax.set_title(
        f"Feasibility Phase Diagram  (N={n_qubits}, threshold={threshold})",
        fontsize=13,
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot 5: N-dependence panels
# ---------------------------------------------------------------------------

def plot_n_dependence(
    df: pd.DataFrame,
    n_values: list[int] = [4, 6, 8],
    theta: float | None = None,
    metric: str = "two_q_depth",
    threshold: float | None = None,
    figsize: tuple = (18, 5),
):
    """Side-by-side phase diagrams for different N."""
    fig, axes = plt.subplots(1, len(n_values), figsize=figsize, sharey=True)
    if len(n_values) == 1:
        axes = [axes]

    for ax, n in zip(axes, n_values):
        # Re-use phase diagram logic inline
        sub = df[df["n_qubits"] == n]
        if theta is not None:
            sub = sub[sub["theta"] == theta]
        else:
            theta_val = sub["theta"].iloc[0] if len(sub) > 0 else 0
            sub = sub[sub["theta"] == theta_val]

        if threshold is None:
            ideal_max = sub[sub["variant"] == "ideal"][metric].max()
            thr = max(3 * ideal_max, 10)
        else:
            thr = threshold

        alphas = np.sort(sub["alpha"].unique())
        steps = np.sort(sub["n_steps"].unique())

        phase = np.full((len(alphas), len(steps)), np.nan)
        for ai, a in enumerate(alphas):
            for si, s in enumerate(steps):
                u = sub[(sub["variant"] == "unitary") & (sub["alpha"] == a) & (sub["n_steps"] == s)][metric]
                d = sub[(sub["variant"] == "dynamic") & (sub["alpha"] == a) & (sub["n_steps"] == s)][metric]
                uv = u.values[0] if len(u) > 0 else np.nan
                dv = d.values[0] if len(d) > 0 else np.nan
                if np.isnan(uv) or np.isnan(dv):
                    continue
                elif uv <= thr and dv <= thr:
                    phase[ai, si] = 0
                elif (uv > thr and dv <= thr) or (uv <= thr and dv > thr):
                    phase[ai, si] = 1
                else:
                    phase[ai, si] = 2

        cmap = mcolors.ListedColormap(["#4CAF50", "#FFC107", "#F44336"])
        norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

        if len(alphas) > 0 and len(steps) > 0:
            ax.imshow(
                phase, aspect="auto", origin="lower",
                extent=[steps.min() - 0.5, steps.max() + 0.5,
                        alphas.min(), alphas.max()],
                cmap=cmap, norm=norm,
            )
        ax.set_title(f"N = {n}")
        ax.set_xlabel("Trotter Steps")

    axes[0].set_ylabel(r"$\alpha$")

    legend_elements = [
        Patch(facecolor="#4CAF50", label="Both feasible"),
        Patch(facecolor="#FFC107", label="Only dynamic feasible"),
        Patch(facecolor="#F44336", label="Neither (need all-to-all)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3, fontsize=10)
    fig.suptitle(f"Phase Diagram vs System Size  (metric: {metric})", fontsize=13)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    return fig

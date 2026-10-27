"""
================================================================================
Script : generate_publication_figures.py
Purpose: Automatically generate all 12 standalone publication-ready vector
         figures (.pdf and .svg) required by manuscript.tex directly from the
         outputs of wa.R (Cam-CAN Watershed SEM Pipeline, N = 513).
Usage  : python3 generate_publication_figures.py
         python3 generate_publication_figures.py --data_dir output --output_dir output_figures
================================================================================
"""

import argparse
import os
import sys

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patches

# ------------------------------------------------------------------------------
# 1. Unified Publication Palette & Typography System
# ------------------------------------------------------------------------------
PALETTE = {
    "NAVY_DARK": "#0f172a",  # Titles, primary labels, dark text
    "SLATE_DARK": "#334155",  # Spines, structural borders, arrows
    "SLATE_MUTED": "#64748b",  # Reference lines, non-significant items
    "SLATE_LIGHT": "#cbd5e1",  # Grid lines, legend/container borders
    "BG_SUBTLE": "#f8fafc",  # Subtle container background
    "BLUE_PRIMARY": "#1d4ed8",  # Primary series (ATR, Forward Model, CFA, Simple RT)
    "BLUE_DARK": "#1e3a8a",  # Dark blue edges & text
    "BLUE_FILL": "#dbeafe",  # Soft blue box highlight
    "TEAL_SECONDARY": "#0f766e",  # Secondary series (TLI, Bimodal RT, Adjusted R2)
    "AMBER_ACCENT": "#d97706",  # Accent series (Forceps Minor, Choice RT)
    "AMBER_DARK": "#b45309",  # Dark amber text/edges
    "AMBER_FILL": "#fef3c7",  # Soft amber box highlight
    "GREEN_WIN": "#15803d",  # Winning Model B / Positive paths / Final cohort
    "GREEN_DARK": "#14532d",  # Dark green text/edges
    "GREEN_FILL": "#f0fdf4",  # Soft green inclusion zone shading
    "RED_NEG": "#b91c1c",  # Negative/suppressor paths, Inverted model, Thresholds
    "RED_LIGHT": "#fca5a5",  # Soft red border
    "RED_FILL": "#fef2f2",  # Soft red callout background
}

TYPO = {
    "W_SINGLE": 7.0,  # Uniform width (in) for all single-column figures
    "W_DOUBLE": 14.5,  # Uniform width (in) for double-column Fig 1
    "TITLE": 10.6,  # Uniform title size across single-col plots
    "AXIS": 9.4,  # Uniform axis label size
    "TICK": 8.6,  # Uniform tick label size
    "ANNOT": 7.9,  # Uniform data annotation size
    "LEGEND": 8.0,  # Uniform legend font size
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": TYPO["TICK"],
        "axes.titlesize": TYPO["TITLE"],
        "axes.titleweight": "bold",
        "axes.titlecolor": PALETTE["NAVY_DARK"],
        "axes.labelsize": TYPO["AXIS"],
        "axes.labelweight": "bold",
        "axes.labelcolor": PALETTE["NAVY_DARK"],
        "xtick.labelsize": TYPO["TICK"],
        "ytick.labelsize": TYPO["TICK"],
        "xtick.color": PALETTE["NAVY_DARK"],
        "ytick.color": PALETTE["NAVY_DARK"],
        "legend.fontsize": TYPO["LEGEND"],
        "figure.titlesize": TYPO["TITLE"],
        "figure.titleweight": "bold",
        "figure.dpi": 300,
        "mathtext.fontset": "dejavusans",
        "axes.edgecolor": PALETTE["SLATE_DARK"],
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": PALETTE["SLATE_LIGHT"],
        "grid.linestyle": ":",
        "grid.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
    }
)


def apply_title_and_top_legend(ax, title, ncol=2, n_rows=1, handles=None, labels=None):
    """Place a non-overlapping legend bar cleanly between the axis top and title."""
    pad_pt = 32 if n_rows == 1 else 50
    ax.set_title(
        title,
        loc="left",
        fontsize=TYPO["TITLE"],
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
        pad=pad_pt,
    )
    kwargs = {
        "loc": "lower center",
        "bbox_to_anchor": (0.5, 1.02),
        "ncol": ncol,
        "frameon": True,
        "fontsize": TYPO["LEGEND"],
        "facecolor": "#ffffff",
        "edgecolor": PALETTE["SLATE_LIGHT"],
        "framealpha": 0.98,
        "borderpad": 0.45,
        "columnspacing": 1.1,
        "handletextpad": 0.5,
    }
    if handles is not None and labels is not None:
        ax.legend(handles, labels, **kwargs)
    elif handles is not None:
        ax.legend(handles=handles, **kwargs)
    else:
        ax.legend(**kwargs)


def save_figure(fig, out_dir, name):
    """Export both PDF and SVG vector formats."""
    pdf_path = os.path.join(out_dir, f"{name}.pdf")
    svg_path = os.path.join(out_dir, f"{name}.svg")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", pad_inches=0.12)
    fig.savefig(svg_path, format="svg", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(f"  [SAVED] {pdf_path} | {svg_path}")


# ------------------------------------------------------------------------------
# 2. Figure 01: Watershed Path Architecture Diagram (Model B, Double-Column)
# ------------------------------------------------------------------------------
def plot_fig01(std_df, paths_df, r2_unadj_df, output_dir):
    std_b = std_df[std_df["model"] == "Model_B"].copy()
    paths_b = paths_df[paths_df["model"] == "Model_B"].copy()
    r2_b = (
        r2_unadj_df[r2_unadj_df["model"] == "Model_B"].set_index("var")["r2"].to_dict()
    )

    fig, ax = plt.subplots(figsize=(TYPO["W_DOUBLE"], 8.6), dpi=300)
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(-0.12, 1.07)
    ax.axis("off")

    ax.text(
        0.00,
        1.045,
        "Empirical Watershed SEM Architecture (Model B: Free-Water-Eliminated Tissue FA, N = 513, df = 60)",
        ha="left",
        va="center",
        fontsize=12.0,
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
    )

    layer_specs = [
        (
            0.00,
            -0.02,
            0.28,
            1.03,
            PALETTE["BG_SUBTLE"],
            PALETTE["SLATE_LIGHT"],
            "TIER 3: UPSTREAM WHITE MATTER\n10 Bilateral JHU Tracts (FWE Tissue FA)",
            PALETTE["NAVY_DARK"],
        ),
        (
            0.36,
            -0.02,
            0.25,
            1.03,
            PALETTE["BG_SUBTLE"],
            PALETTE["SLATE_LIGHT"],
            "TIER 2: INTERMEDIATE COGNITIVE\nProcessing Speed Primitives",
            PALETTE["NAVY_DARK"],
        ),
        (
            0.70,
            -0.02,
            0.30,
            1.03,
            "#eff6ff",
            "#93c5fd",
            "TIER 1: DOWNSTREAM PHENOTYPE\nFluid Intelligence (Cattell CFA)",
            PALETTE["BLUE_DARK"],
        ),
    ]
    for lx, ly, lw, lh, fc, ec, title, tc in layer_specs:
        ax.add_patch(
            patches.FancyBboxPatch(
                (lx, ly),
                lw,
                lh,
                boxstyle="round,pad=0.008,rounding_size=0.018",
                facecolor=fc,
                edgecolor=ec,
                linewidth=1.2,
                linestyle="--",
            )
        )
        ax.text(
            lx + lw / 2,
            ly + lh - 0.025,
            title,
            ha="center",
            va="top",
            fontsize=9.0,
            fontweight="bold",
            color=tc,
            linespacing=1.3,
        )

    tracts = [
        ("ATR", "Anterior Thalamic Rad. (ATR)", 0.865),
        ("CING", "Cingulate Cingulum (CING)", 0.770),
        ("FMin", "Forceps Minor — Genu (FMin)", 0.675),
        ("IFOF", "Inf. Fronto-Occipital (IFOF)", 0.580),
        ("CINGHipp", "Hippocampal Cingulum (CINGHipp)", 0.485),
        ("UNC", "Uncinate Fasciculus (UNC)", 0.390),
        ("SLF", "Sup. Longitudinal Fasc. (SLF)", 0.295),
        ("CST", "Corticospinal Tract (CST)", 0.200),
        ("FMaj", "Forceps Major — Splenium (FMaj)", 0.105),
        ("ILF", "Inf. Longitudinal Fasc. (ILF)", 0.010),
    ]

    speed_meta = [
        ("SRTspeed", "Simple RT Speed", 0.81),
        ("SRTcons", "Simple RT Consistency", 0.66),
        ("CRTspeed", "Choice RT Speed", 0.51),
        ("CRTcons", "Choice RT Consistency", 0.36),
        ("AVspeed", "Audio-Visual Speed", 0.21),
        ("AVcons", "Audio-Visual Consistency", 0.06),
    ]

    sig_tracts = set(paths_b[paths_b["pvalue"] < 0.05]["tract"].unique())
    tract_coords = {}
    for code, full_label, y in tracts:
        is_sig = code in sig_tracts
        fc = PALETTE["BLUE_FILL"] if is_sig else "#ffffff"
        ec = PALETTE["BLUE_PRIMARY"] if is_sig else PALETTE["SLATE_MUTED"]
        rect = patches.FancyBboxPatch(
            (0.018, y - 0.028),
            0.244,
            0.056,
            boxstyle="round,pad=0.005,rounding_size=0.010",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.5 if is_sig else 0.9,
        )
        ax.add_patch(rect)
        ax.text(
            0.140,
            y,
            full_label,
            ha="center",
            va="center",
            fontsize=8.1,
            fontweight="bold" if is_sig else "normal",
            color=PALETTE["NAVY_DARK"],
        )
        tract_coords[code] = (0.267, y)

    speed_coords = {}
    for code, name_lbl, y in speed_meta:
        is_crt = "CRT" in code
        r2_pct = r2_b.get(code, 0.0) * 100.0
        label = f"{name_lbl}\n(R² = {r2_pct:.1f}%)"
        fc = PALETTE["AMBER_FILL"] if is_crt else "#ffffff"
        ec = PALETTE["AMBER_ACCENT"] if is_crt else PALETTE["SLATE_MUTED"]
        rect = patches.FancyBboxPatch(
            (0.380, y - 0.048),
            0.210,
            0.096,
            boxstyle="round,pad=0.005,rounding_size=0.014",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.5 if is_crt else 1.0,
        )
        ax.add_patch(rect)
        ax.text(
            0.485,
            y,
            label,
            ha="center",
            va="center",
            fontsize=8.2,
            fontweight="bold" if is_crt else "normal",
            color=PALETTE["NAVY_DARK"],
            linespacing=1.25,
        )
        speed_coords[code] = (0.375, y, 0.595, y)

    fi_cx, fi_cy = 0.780, 0.450
    fi_rx, fi_ry = 0.062, 0.125
    fi_r2_pct = r2_b.get("FI", 0.5777) * 100.0
    fi_ellipse = patches.Ellipse(
        (fi_cx, fi_cy),
        width=2 * fi_rx,
        height=2 * fi_ry,
        facecolor=PALETTE["BLUE_PRIMARY"],
        edgecolor=PALETTE["BLUE_DARK"],
        linewidth=2.0,
        zorder=4,
    )
    ax.add_patch(fi_ellipse)
    ax.text(
        fi_cx,
        fi_cy,
        f"Latent FI\n(Fluid Intelligence)\n\nR² = {fi_r2_pct:.1f}%",
        ha="center",
        va="center",
        fontsize=8.3,
        fontweight="bold",
        color="#ffffff",
        linespacing=1.2,
        zorder=5,
    )

    cfa_b = std_b[(std_b["lhs"] == "FI") & (std_b["op"] == "=~")].set_index("rhs")
    indicators = [
        ("C1", "Series", 0.66),
        ("C2", "Classif.", 0.52),
        ("C3", "Matrices", 0.38),
        ("C4", "Topology", 0.24),
    ]
    for code, short_lbl, y in indicators:
        lam = cfa_b.loc[code, "est.std"]
        rect = patches.FancyBboxPatch(
            (0.892, y - 0.042),
            0.090,
            0.084,
            boxstyle="round,pad=0.005,rounding_size=0.010",
            facecolor="#ffffff",
            edgecolor=PALETTE["SLATE_DARK"],
            linewidth=1.1,
        )
        ax.add_patch(rect)
        ax.text(
            0.937,
            y,
            f"{code}: {short_lbl}\n(λ = {lam:.2f})",
            ha="center",
            va="center",
            fontsize=7.8,
            fontweight="bold",
            color=PALETTE["NAVY_DARK"],
            linespacing=1.2,
        )
        dy_norm = (y - fi_cy) / 0.28
        y_start = fi_cy + dy_norm * (fi_ry * 0.65)
        x_start = fi_cx + fi_rx * np.sqrt(
            max(0.0, 1.0 - ((y_start - fi_cy) / fi_ry) ** 2)
        )
        ax.annotate(
            "",
            xy=(0.887, y),
            xytext=(x_start, y_start),
            arrowprops={
                "arrowstyle": "-|>,head_length=0.5,head_width=0.28",
                "color": PALETTE["SLATE_DARK"],
                "lw": 1.3,
            },
        )

    tier2_b = std_b[(std_b["lhs"] == "FI") & (std_b["op"] == "~")].set_index("rhs")
    label_pos_t = {
        "SRTspeed": 0.50,
        "SRTcons": 0.46,
        "CRTspeed": 0.45,
        "CRTcons": 0.45,
        "AVspeed": 0.42,
        "AVcons": 0.54,
    }
    for scode, (sx, sy, ex, ey) in speed_coords.items():
        beta = tier2_b.loc[scode, "est.std"]
        pval = tier2_b.loc[scode, "pvalue"]
        sig_star = (
            "***"
            if pval < 0.001
            else ("**" if pval < 0.01 else ("*" if pval < 0.05 else ""))
        )
        lbl_txt = f"{beta:+.2f}{sig_star}" if pval < 0.05 else "ns"
        label_t = label_pos_t[scode]
        color = (
            PALETTE["SLATE_MUTED"]
            if pval >= 0.05
            else (PALETTE["RED_NEG"] if beta < 0 else PALETTE["GREEN_WIN"])
        )
        lw = 1.0 if pval >= 0.05 else 1.2 + abs(beta) * 4.0
        ls = "--" if pval >= 0.05 else "-"
        y_targ = fi_cy + (ey - fi_cy) * 0.26
        x_targ = (
            fi_cx
            - fi_rx * np.sqrt(max(0.0, 1.0 - ((y_targ - fi_cy) / fi_ry) ** 2))
            - 0.003
        )
        ax.annotate(
            "",
            xy=(x_targ, y_targ),
            xytext=(ex, ey),
            arrowprops={
                "arrowstyle": "-|>,head_length=0.55,head_width=0.30",
                "color": color,
                "lw": lw,
                "linestyle": ls,
                "shrinkA": 2,
                "shrinkB": 2,
            },
        )
        if pval < 0.05:
            mx = ex * (1.0 - label_t) + x_targ * label_t
            my = ey * (1.0 - label_t) + y_targ * label_t
            ax.text(
                mx,
                my,
                lbl_txt,
                fontsize=7.6,
                fontweight="bold",
                color=color,
                ha="center",
                va="center",
                bbox={
                    "boxstyle": "round,pad=0.14",
                    "fc": "#ffffff",
                    "ec": "none",
                    "alpha": 0.92,
                },
                zorder=6,
            )

    sig_paths_b = paths_b[paths_b["pvalue"] < 0.05].copy()
    for scode, grp in sig_paths_b.groupby("lhs"):
        n_g = len(grp)
        offsets = np.linspace(0.015, -0.015, n_g) if n_g > 1 else [0.0]
        for idx, (_, row) in enumerate(grp.iterrows()):
            tcode = row["tract"]
            beta = row["est.std"]
            p_holm = row["p_holm"]
            tx, ty = tract_coords[tcode]
            sx, sy, _, _ = speed_coords[scode]
            color = (
                PALETTE["RED_NEG"]
                if beta < 0
                else (
                    PALETTE["BLUE_DARK"] if p_holm < 0.05 else PALETTE["BLUE_PRIMARY"]
                )
            )
            lw = (1.6 if p_holm < 0.05 else 1.1) + abs(beta) * 3.8
            rad = 0.025 if ty > sy else -0.025
            ax.annotate(
                "",
                xy=(sx, sy + offsets[idx]),
                xytext=(tx, ty),
                arrowprops={
                    "arrowstyle": "-|>,head_length=0.5,head_width=0.28",
                    "color": color,
                    "lw": lw,
                    "shrinkA": 2,
                    "shrinkB": 2,
                    "connectionstyle": f"arc3,rad={rad}",
                },
            )

    leg_handles = [
        patches.Patch(
            facecolor=PALETTE["BLUE_FILL"],
            edgecolor=PALETTE["BLUE_PRIMARY"],
            lw=1.4,
            label="Significant Upstream Tract (p < 0.05)",
        ),
        patches.Patch(
            facecolor=PALETTE["AMBER_FILL"],
            edgecolor=PALETTE["AMBER_ACCENT"],
            lw=1.4,
            label="Executive Choice RT Primitive",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["BLUE_PRIMARY"],
            lw=2.4,
            label="Tract → Speed (Positive β; Holm-sig ATR in thick navy)",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["GREEN_WIN"],
            lw=2.5,
            label="Speed → FI (Positive β, p < 0.05)",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["RED_NEG"],
            lw=2.2,
            label="Negative / Suppressor Path (p < 0.05)",
        ),
    ]
    ax.legend(
        handles=leg_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.105),
        ncol=5,
        frameon=True,
        fontsize=TYPO["LEGEND"],
        facecolor="#ffffff",
        edgecolor=PALETTE["SLATE_LIGHT"],
        borderpad=0.55,
    )
    save_figure(fig, output_dir, "fig01_watershed_path_diagram")


# ------------------------------------------------------------------------------
# 3. Figure 02: Fluid Intelligence Measurement Model (CFA)
# ------------------------------------------------------------------------------
def plot_fig02(std_df, cfa_fit_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 4.4), dpi=300)
    loadings = (
        std_df[
            (std_df["model"] == "Model_B")
            & (std_df["lhs"] == "FI")
            & (std_df["op"] == "=~")
        ]
        .copy()
        .reset_index(drop=True)
    )
    label_map = {
        "C1": "Series Completion (C1)",
        "C2": "Classification (C2)",
        "C3": "Matrices (C3)",
        "C4": "Conditions / Topology (C4)",
    }
    loadings["label"] = loadings["rhs"].map(label_map)
    y_pos = np.arange(len(loadings))

    ax.axvline(
        0.70,
        color=PALETTE["SLATE_MUTED"],
        linestyle=":",
        linewidth=1.2,
        label="Strong Loading Benchmark (λ = 0.70)",
    )
    ax.errorbar(
        loadings["est.std"],
        y_pos,
        xerr=[
            loadings["est.std"] - loadings["ci.lower"],
            loadings["ci.upper"] - loadings["est.std"],
        ],
        fmt="o",
        color=PALETTE["BLUE_PRIMARY"],
        ecolor=PALETTE["BLUE_PRIMARY"],
        elinewidth=2.2,
        capsize=5,
        capthick=1.5,
        markersize=7.5,
        zorder=3,
        label="Factor Loading (λ ± 95% CI)",
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(loadings["label"], fontweight="bold")
    ax.set_ylim(len(loadings) - 0.35, -0.65)
    ax.set_xlabel("Standardized Factor Loading (λ ± 95% CI)")
    ax.set_xlim(0.60, 1.05)
    ax.grid(True, axis="x", alpha=0.7)

    for i, row in loadings.iterrows():
        r2_ind = (row["est.std"] ** 2) * 100.0
        ax.text(
            row["ci.upper"] + 0.010,
            i,
            f"λ = {row['est.std']:.3f} [{row['ci.lower']:.3f}, {row['ci.upper']:.3f}]\n(z = {row['z']:.2f}, R² = {r2_ind:.1f}%)",
            ha="left",
            va="center",
            fontsize=TYPO["ANNOT"],
            color=PALETTE["NAVY_DARK"],
            fontweight="bold",
            linespacing=1.22,
        )

    c_row = cfa_fit_df.iloc[0]
    fit_txt = (
        f"Standalone 1-Factor CFA Fit (N = {int(c_row['n'])}, df = {int(c_row['df'])}): χ² = {c_row['chisq.scaled']:.2f} (p = {c_row['pvalue.scaled']:.3f})\n"
        f"Robust CFI = {c_row['cfi.robust']:.4f}, Robust RMSEA = {c_row['rmsea.robust']:.4f}, SRMR = {c_row['srmr']:.4f}"
    )
    ax.text(
        0.5,
        -0.24,
        fit_txt,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.6,
        color=PALETTE["SLATE_DARK"],
        bbox={
            "boxstyle": "round,pad=0.3",
            "fc": PALETTE["BG_SUBTLE"],
            "ec": PALETTE["SLATE_LIGHT"],
        },
    )

    apply_title_and_top_legend(
        ax,
        "Tier 1: Fluid Intelligence Reflective Measurement Model (CFA)",
        ncol=2,
        n_rows=1,
    )
    save_figure(fig, output_dir, "fig02_fi_measurement_cfa")


# ------------------------------------------------------------------------------
# 4. Figure 03: Speed Primitives to Latent FI Paths
# ------------------------------------------------------------------------------
def plot_fig03(std_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 4.7), dpi=300)
    sp_paths = (
        std_df[
            (std_df["model"] == "Model_B")
            & (std_df["lhs"] == "FI")
            & (std_df["op"] == "~")
        ]
        .copy()
        .reset_index(drop=True)
    )
    label_map = {
        "CRTspeed": "Choice RT Speed (CRTspeed)",
        "AVcons": "Audio-Visual Consistency (AVcons)",
        "CRTcons": "Choice RT Consistency (CRTcons)",
        "SRTcons": "Simple RT Consistency (SRTcons)",
        "SRTspeed": "Simple RT Speed (SRTspeed)",
        "AVspeed": "Audio-Visual Speed (AVspeed)",
    }
    order = ["CRTspeed", "AVcons", "CRTcons", "SRTcons", "SRTspeed", "AVspeed"]
    sp_paths = sp_paths.set_index("rhs").loc[order].reset_index()
    sp_paths["label"] = sp_paths["rhs"].map(label_map)

    y_sp = np.arange(len(sp_paths))
    colors = [
        PALETTE["GREEN_WIN"]
        if b > 0 and p < 0.05
        else PALETTE["RED_NEG"]
        if b < 0 and p < 0.05
        else PALETTE["SLATE_MUTED"]
        for b, p in zip(sp_paths["est.std"], sp_paths["pvalue"])
    ]

    ax.axvline(0, color=PALETTE["SLATE_DARK"], linestyle="--", linewidth=1.1, zorder=1)
    for i in range(len(sp_paths)):
        ax.errorbar(
            sp_paths["est.std"].iloc[i],
            y_sp[i],
            xerr=[
                [sp_paths["est.std"].iloc[i] - sp_paths["ci.lower"].iloc[i]],
                [sp_paths["ci.upper"].iloc[i] - sp_paths["est.std"].iloc[i]],
            ],
            fmt="o",
            color=colors[i],
            ecolor=colors[i],
            elinewidth=2.2,
            capsize=5,
            capthick=1.5,
            markersize=7.5,
            zorder=3,
        )

    ax.set_yticks(y_sp)
    ax.set_yticklabels(sp_paths["label"], fontweight="bold")
    ax.set_ylim(len(sp_paths) - 0.4, -0.6)
    ax.set_xlabel("Standardized Path Coefficient (β ± 95% CI)")
    ax.set_xlim(-0.42, 1.06)
    ax.grid(True, axis="x", alpha=0.7)

    for i, row in sp_paths.iterrows():
        sig_str = "p < 0.0001" if row["pvalue"] < 0.0001 else f"p = {row['pvalue']:.4f}"
        x_anchor = max(row["ci.upper"], 0.01) + 0.020
        ax.text(
            x_anchor,
            i,
            f"β = {row['est.std']:+.3f} ({sig_str})",
            va="center",
            ha="left",
            fontsize=TYPO["ANNOT"],
            color=colors[i],
            fontweight="bold" if row["pvalue"] < 0.05 else "normal",
            bbox={
                "boxstyle": "round,pad=0.15",
                "fc": "#ffffff",
                "ec": "none",
                "alpha": 0.88,
            },
        )

    leg_handles = [
        mlines.Line2D(
            [],
            [],
            color=PALETTE["GREEN_WIN"],
            marker="o",
            lw=2.0,
            markersize=6,
            label="Positive Path (p < 0.05)",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["RED_NEG"],
            marker="o",
            lw=2.0,
            markersize=6,
            label="Suppressor Path (p < 0.05)",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["SLATE_MUTED"],
            marker="o",
            lw=2.0,
            markersize=6,
            label="Non-Significant (ns)",
        ),
    ]
    apply_title_and_top_legend(
        ax,
        "Tier 2: Processing Speed Paths to Latent FI (R² = 57.77%)",
        ncol=3,
        n_rows=1,
        handles=leg_handles,
    )
    save_figure(fig, output_dir, "fig03_speed_to_fi_forest")


# ------------------------------------------------------------------------------
# 5. Figure 04: Incremental Fit Indices (CFI & TLI Across Modalities)
# ------------------------------------------------------------------------------
def plot_fig04(comp_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 4.8), dpi=300)
    comp_sorted = comp_df.sort_values("Rank").reset_index(drop=True)
    short_labels = {
        "Model B (FWE Tissue FA)": "Model B\n(FWE Tissue FA)",
        "Model D (Free-Water FW)": "Model D\n(Free-Water FW)",
        "Model A (Conventional FA)": "Model A\n(Conventional FA)",
        "Model C (WMTI AWF)": "Model C\n(WMTI AWF)",
    }
    models_order = [short_labels.get(m, m) for m in comp_sorted["model"].tolist()]
    x = np.arange(len(models_order))
    width = 0.34

    cfi = comp_sorted["cfi.robust"].values
    tli = comp_sorted["tli.robust"].values

    rects1 = ax.bar(
        x - width / 2,
        cfi,
        width,
        label="Robust CFI",
        color=PALETTE["BLUE_PRIMARY"],
        edgecolor=PALETTE["BLUE_DARK"],
        linewidth=0.8,
    )
    rects2 = ax.bar(
        x + width / 2,
        tli,
        width,
        label="Robust TLI",
        color=PALETTE["TEAL_SECONDARY"],
        edgecolor="#115e59",
        linewidth=0.8,
    )

    ax.axhline(
        0.95,
        color=PALETTE["RED_NEG"],
        linestyle="--",
        linewidth=1.2,
        label="Acceptable Threshold (0.95)",
    )
    ax.axhline(
        0.97,
        color=PALETTE["SLATE_MUTED"],
        linestyle=":",
        linewidth=1.1,
        label="Close-Fit Threshold (0.97)",
    )
    ax.set_ylim(0.940, 1.014)
    ax.set_xticks(x)
    ax.set_xticklabels(models_order, fontweight="bold")
    ax.set_ylabel("Incremental Fit Index")
    ax.grid(axis="y", alpha=0.7)

    for r in rects1:
        h = r.get_height()
        ax.text(
            r.get_x() + r.get_width() / 2,
            h + 0.0014,
            f"{h:.4f}",
            ha="center",
            va="bottom",
            fontsize=TYPO["ANNOT"],
            fontweight="bold",
            color=PALETTE["BLUE_DARK"],
        )
    for r in rects2:
        h = r.get_height()
        ax.text(
            r.get_x() + r.get_width() / 2,
            h + 0.0014,
            f"{h:.4f}",
            ha="center",
            va="bottom",
            fontsize=TYPO["ANNOT"],
            fontweight="bold",
            color=PALETTE["TEAL_SECONDARY"],
        )

    apply_title_and_top_legend(
        ax,
        "Global Incremental Fit (Robust CFI & TLI Across Modalities, df = 60)",
        ncol=2,
        n_rows=2,
    )
    save_figure(fig, output_dir, "fig04_modality_cfi_tli")


# ------------------------------------------------------------------------------
# 6. Figure 05: Approximation Error (RMSEA 90% CI Across Modalities)
# ------------------------------------------------------------------------------
def plot_fig05(comp_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 4.6), dpi=300)
    comp_sorted = comp_df.sort_values("Rank").reset_index(drop=True)
    short_labels = {
        "Model B (FWE Tissue FA)": "Model B\n(FWE Tissue FA)",
        "Model D (Free-Water FW)": "Model D\n(Free-Water FW)",
        "Model A (Conventional FA)": "Model A\n(Conventional FA)",
        "Model C (WMTI AWF)": "Model C\n(WMTI AWF)",
    }
    models_order = [short_labels.get(m, m) for m in comp_sorted["model"].tolist()]
    y = np.arange(len(models_order))

    rmsea = comp_sorted["rmsea.robust"].values
    rmsea_low = comp_sorted["rmsea.ci.lower.robust"].values
    rmsea_upp = comp_sorted["rmsea.ci.upper.robust"].values
    srmr = comp_sorted["srmr"].values

    ax.axvline(
        0.05,
        color=PALETTE["RED_NEG"],
        linestyle="--",
        linewidth=1.2,
        label="Close Fit Threshold (≤ 0.05)",
    )
    for i in range(len(models_order)):
        is_winner = i == 0
        color = PALETTE["GREEN_WIN"] if is_winner else PALETTE["BLUE_PRIMARY"]
        ax.errorbar(
            rmsea[i],
            y[i],
            xerr=[[rmsea[i] - rmsea_low[i]], [rmsea_upp[i] - rmsea[i]]],
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=2.2,
            capsize=5,
            capthick=1.5,
            markersize=7.5,
            label="Rank 1: Model B (FWE Tissue FA)"
            if is_winner
            else ("Candidate Modality" if i == 1 else None),
        )
        fw = "bold" if is_winner else "normal"
        ax.text(
            rmsea_upp[i] + 0.0018,
            y[i],
            f"RMSEA = {rmsea[i]:.4f} [{rmsea_low[i]:.4f}, {rmsea_upp[i]:.4f}]\n(SRMR = {srmr[i]:.4f})",
            va="center",
            ha="left",
            fontsize=TYPO["ANNOT"],
            color=PALETTE["GREEN_DARK"] if is_winner else PALETTE["NAVY_DARK"],
            fontweight=fw,
            linespacing=1.2,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(models_order, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Robust RMSEA (90% Confidence Interval)")
    ax.set_xlim(0.012, 0.088)
    ax.set_ylim(len(models_order) - 0.4, -0.6)
    ax.grid(axis="x", alpha=0.7)

    apply_title_and_top_legend(
        ax,
        "Approximation & Residual Error Across Modalities (df = 60)",
        ncol=3,
        n_rows=1,
    )
    save_figure(fig, output_dir, "fig05_modality_rmsea_ci")


# ------------------------------------------------------------------------------
# 7. Figure 06: Directional Hierarchy Inversion Fit Degradation (Models B & A)
# ------------------------------------------------------------------------------
def plot_fig06(dir_df, output_dir):
    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(TYPO["W_SINGLE"], 4.8),
        dpi=300,
        gridspec_kw={"width_ratios": [1.0, 1.35], "wspace": 0.36},
    )

    d_map = dir_df.set_index("model")
    b_fwd = d_map.loc["Model_B_forward"]
    b_inv = d_map.loc["Model_B_inverted_hierarchy"]
    a_fwd = d_map.loc["Model_A_forward"]
    a_inv = d_map.loc["Model_A_inverted_hierarchy"]

    # Panel A: Robust CFI across Model B and Model A
    x1 = np.arange(2)
    w = 0.34
    b1 = ax1.bar(
        x1 - w / 2,
        [b_fwd["cfi.robust"], a_fwd["cfi.robust"]],
        w,
        color=PALETTE["BLUE_PRIMARY"],
        edgecolor=PALETTE["BLUE_DARK"],
        linewidth=0.8,
    )
    b2 = ax1.bar(
        x1 + w / 2,
        [b_inv["cfi.robust"], a_inv["cfi.robust"]],
        w,
        color=PALETTE["RED_NEG"],
        edgecolor="#7f1d1d",
        linewidth=0.8,
    )
    ax1.axhline(0.97, color=PALETTE["SLATE_MUTED"], linestyle=":", linewidth=1.1)
    ax1.set_ylim(0.92, 1.018)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(["Model B\n(FWE FA)", "Model A\n(Conv. FA)"], fontweight="bold")
    ax1.set_ylabel("Robust CFI (Higher = Better)")
    ax1.grid(axis="y", alpha=0.7)
    for bars, col in [(b1, PALETTE["BLUE_DARK"]), (b2, PALETTE["RED_NEG"])]:
        for bar in bars:
            h = bar.get_height()
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.002,
                f"{h:.3f}",
                ha="center",
                va="bottom",
                fontsize=7.4,
                fontweight="bold",
                color=col,
            )

    # Panel B: Robust RMSEA & SRMR for Winning Model B (with Model A annotated)
    x2 = np.arange(2)
    b3 = ax2.bar(
        x2 - w / 2,
        [b_fwd["rmsea.robust"], b_fwd["srmr"]],
        w,
        color=PALETTE["BLUE_PRIMARY"],
        edgecolor=PALETTE["BLUE_DARK"],
        linewidth=0.8,
        label="Forward (WM → PS → FI, df = 60)",
    )
    b4 = ax2.bar(
        x2 + w / 2,
        [b_inv["rmsea.robust"], b_inv["srmr"]],
        w,
        color=PALETTE["RED_NEG"],
        edgecolor="#7f1d1d",
        linewidth=0.8,
        label="Inverted (PS → WM → FI, df = 56)",
    )
    ax2.axhline(
        0.05,
        color=PALETTE["RED_NEG"],
        linestyle="--",
        linewidth=1.1,
        label="Close Error Threshold (≤ 0.05)",
    )
    ax2.plot(
        [],
        [],
        color=PALETTE["SLATE_MUTED"],
        linestyle=":",
        linewidth=1.1,
        label="Close CFI Threshold (≥ 0.97)",
    )
    ax2.set_ylim(0.0, 0.098)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(
        ["Model B RMSEA\n(Lower = Better)", "Model B SRMR\n(Lower = Better)"],
        fontweight="bold",
    )
    ax2.set_ylabel("Approximation / Residual Error")
    ax2.grid(axis="y", alpha=0.7)
    for bars, col in [(b3, PALETTE["BLUE_DARK"]), (b4, PALETTE["RED_NEG"])]:
        for bar in bars:
            h = bar.get_height()
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.002,
                f"{h:.4f}",
                ha="center",
                va="bottom",
                fontsize=TYPO["ANNOT"],
                fontweight="bold",
                color=col,
            )

    fig.suptitle(
        "Hierarchy Inversion Test: Fit Degradation (df = 60 vs. df = 56)",
        x=0.01,
        y=1.12,
        ha="left",
        fontsize=TYPO["TITLE"],
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
    )
    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.07),
        ncol=2,
        frameon=True,
        fontsize=TYPO["LEGEND"],
        facecolor="#ffffff",
        edgecolor=PALETTE["SLATE_LIGHT"],
        framealpha=0.98,
        borderpad=0.45,
    )
    save_figure(fig, output_dir, "fig06_inversion_fit_degradation")


# ------------------------------------------------------------------------------
# 8. Figure 07: Likelihood Evidence & Delta AIC Decision (Models B & A)
# ------------------------------------------------------------------------------
def plot_fig07(dir_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 4.7), dpi=300)
    d_map = dir_df.set_index("model")

    b_delta = d_map.loc["Model_B_inverted_hierarchy", "delta_aic_vs_forward"]
    a_delta = d_map.loc["Model_A_inverted_hierarchy", "delta_aic_vs_forward"]
    b_rev = d_map.loc["Model_B_reversed_WM_PS_arrows", "delta_aic_vs_forward"]

    specs = [
        "Model B Inverted\n(PS → WM → FI, df = 56)",
        "Model A Inverted\n(PS → WM → FI, df = 56)",
        "Reversed Lower Block\n(WM ← PS, df = 60, Markov Eq.)",
    ]
    deltas = [b_delta, a_delta, b_rev]
    colors = [
        PALETTE["RED_NEG"],
        PALETTE["AMBER_ACCENT"],
        PALETTE["SLATE_MUTED"],
    ]

    y = np.arange(len(specs))
    bars = ax.barh(
        y,
        deltas,
        color=colors,
        height=0.52,
        edgecolor=PALETTE["SLATE_DARK"],
        linewidth=0.9,
    )
    ax.axvline(
        10.0,
        color=PALETTE["RED_NEG"],
        linestyle="--",
        linewidth=1.2,
        label="Decisive Support Threshold (ΔAIC > +10)",
    )
    ax.axvline(0.0, color=PALETTE["SLATE_DARK"], linestyle="-", linewidth=1.0)

    ax.set_yticks(y)
    ax.set_yticklabels(specs, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel(
        "AIC Penalty Relative to Forward Watershed (ΔAIC = AIC_alt − AIC_forward)"
    )
    ax.set_xlim(-8, 248)
    ax.grid(axis="x", alpha=0.7)

    chi_info = [
        f"ΔAIC = +{b_delta:.2f}  (χ² = 241.17 vs. 92.27)",
        f"ΔAIC = +{a_delta:.2f}  (χ² = 205.52 vs. 97.56)",
        "ΔAIC = 0.00  (Identical Joint Likelihood)",
    ]
    for i, b in enumerate(bars):
        w = b.get_width()
        ax.text(
            w + 4.0,
            b.get_y() + b.get_height() / 2,
            chi_info[i],
            va="center",
            ha="left",
            fontsize=TYPO["ANNOT"],
            fontweight="bold",
            color=colors[i] if i < 2 else PALETTE["SLATE_DARK"],
        )

    leg_handles = [
        patches.Patch(
            facecolor=PALETTE["RED_NEG"],
            edgecolor=PALETTE["SLATE_DARK"],
            label="Model B Inverted (+160.26 AIC)",
        ),
        patches.Patch(
            facecolor=PALETTE["AMBER_ACCENT"],
            edgecolor=PALETTE["SLATE_DARK"],
            label="Model A Inverted (+118.48 AIC)",
        ),
        mlines.Line2D(
            [],
            [],
            color=PALETTE["RED_NEG"],
            linestyle="--",
            lw=1.3,
            label="Decisive Threshold (ΔAIC = +10)",
        ),
    ]
    apply_title_and_top_legend(
        ax,
        "Directional Hierarchy Comparison: Joint AIC Penalties (N = 513)",
        ncol=3,
        n_rows=1,
        handles=leg_handles,
    )
    save_figure(fig, output_dir, "fig07_inversion_delta_aic")


# ------------------------------------------------------------------------------
# 9. Figure 08: Tract-to-Speed Functional Specificity Matrix (Heatmap)
# ------------------------------------------------------------------------------
def plot_fig08(paths_df, r2_unadj_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 6.7), dpi=300)
    paths_b = paths_df[paths_df["model"] == "Model_B"].copy()
    r2_b = (
        r2_unadj_df[r2_unadj_df["model"] == "Model_B"].set_index("var")["r2"].to_dict()
    )

    tracts_order = [
        "ATR",
        "FMin",
        "CING",
        "IFOF",
        "CINGHipp",
        "UNC",
        "SLF",
        "CST",
        "FMaj",
        "ILF",
    ]
    tract_names_map = {
        "ATR": "Anterior Thalamic Rad. (ATR)",
        "FMin": "Forceps Minor — Genu (FMin)",
        "CING": "Cingulate Cingulum (CING)",
        "IFOF": "Inf. Fronto-Occipital Fasc. (IFOF)",
        "CINGHipp": "Hippocampal Cingulum (CINGHipp)",
        "UNC": "Uncinate Fasciculus (UNC)",
        "SLF": "Sup. Longitudinal Fasc. (SLF)",
        "CST": "Corticospinal Tract (CST)",
        "FMaj": "Forceps Major — Splenium (FMaj)",
        "ILF": "Inf. Longitudinal Fasc. (ILF)",
    }
    speed_order = [
        "CRTspeed",
        "CRTcons",
        "SRTspeed",
        "SRTcons",
        "AVcons",
        "AVspeed",
    ]
    speed_names = {
        "CRTspeed": "Choice RT\nSpeed",
        "CRTcons": "Choice RT\nConsistency",
        "SRTspeed": "Simple RT\nSpeed",
        "SRTcons": "Simple RT\nConsistency",
        "AVcons": "Audio-Visual\nConsistency",
        "AVspeed": "Audio-Visual\nSpeed",
    }
    speed_labels = [
        f"{speed_names[s]}\n(R²={r2_b.get(s, 0.0) * 100:.1f}%)" for s in speed_order
    ]

    beta_mat = np.zeros((len(tracts_order), len(speed_order)))
    pval_mat = np.ones((len(tracts_order), len(speed_order)))
    holm_mat = np.ones((len(tracts_order), len(speed_order)))

    for i, t in enumerate(tracts_order):
        for j, s in enumerate(speed_order):
            row = paths_b[(paths_b["tract"] == t) & (paths_b["lhs"] == s)]
            if len(row) > 0:
                beta_mat[i, j] = row["est.std"].values[0]
                pval_mat[i, j] = row["pvalue"].values[0]
                holm_mat[i, j] = row["p_holm"].values[0]

    vlag_cmap = sns.color_palette("vlag", as_cmap=True)
    im = ax.imshow(beta_mat, cmap=vlag_cmap, vmin=-0.38, vmax=0.38, aspect="auto")

    ax.set_xticks(np.arange(len(speed_order) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(tracts_order) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="#ffffff", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_xticks(np.arange(len(speed_order)))
    ax.set_xticklabels(speed_labels, fontsize=8.0, fontweight="bold")
    ax.set_yticks(np.arange(len(tracts_order)))
    ax.set_yticklabels(
        [tract_names_map[t] for t in tracts_order],
        fontsize=8.4,
        fontweight="bold",
    )
    ax.set_title(
        "Upstream FWE Tissue FA → Processing Speed Beta Matrix (β)",
        loc="left",
        fontsize=TYPO["TITLE"],
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
        pad=14,
    )

    for i in range(len(tracts_order)):
        for j in range(len(speed_order)):
            val = beta_mat[i, j]
            pval = pval_mat[i, j]
            pholm = holm_mat[i, j]
            sig = (
                "***"
                if pval < 0.001
                else "**"
                if pval < 0.01
                else "*"
                if pval < 0.05
                else "†"
                if pval < 0.10
                else ""
            )
            text_color = "#ffffff" if abs(val) > 0.20 else PALETTE["NAVY_DARK"]
            fw = "bold" if pval < 0.05 else "normal"
            label_cell = f"{val:+.3f}{sig}"
            if pholm < 0.05:
                ax.add_patch(
                    patches.Rectangle(
                        (j - 0.48, i - 0.48),
                        0.96,
                        0.96,
                        fill=False,
                        edgecolor=PALETTE["NAVY_DARK"],
                        linewidth=2.2,
                        zorder=4,
                    )
                )
            ax.text(
                j,
                i,
                label_cell,
                ha="center",
                va="center",
                fontsize=TYPO["ANNOT"],
                color=text_color,
                fontweight=fw,
            )

    cbar = fig.colorbar(
        im,
        ax=ax,
        orientation="horizontal",
        pad=0.14,
        fraction=0.046,
        aspect=26,
        shrink=0.92,
    )
    cbar.set_label(
        "Standardized Beta (β)  [*p < 0.05, **p < 0.01, ***p < 0.001, †p < 0.10; dark border = p_Holm < 0.05]",
        fontsize=7.6,
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
    )
    save_figure(fig, output_dir, "fig08_tract_speed_beta_matrix")


# ------------------------------------------------------------------------------
# 10. Figure 09: Anterior System Specificity (ATR vs. FMin in Model A vs. Model B)
# ------------------------------------------------------------------------------
def plot_fig09(paths_df, wald_df, boot_df, output_dir):
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(TYPO["W_SINGLE"], 5.2),
        dpi=300,
        sharey=True,
        gridspec_kw={"wspace": 0.10},
    )
    speed_order = [
        "CRTspeed",
        "CRTcons",
        "SRTspeed",
        "SRTcons",
        "AVcons",
        "AVspeed",
    ]
    speed_labels = [
        "Choice RT Speed",
        "Choice RT Consistency",
        "Simple RT Speed",
        "Simple RT Consistency",
        "Audio-Visual Consistency",
        "Audio-Visual Speed",
    ]
    y_f = np.arange(len(speed_order))
    offset = 0.18

    boot_unadj = boot_df[
        (boot_df["specification"] == "unadjusted")
        & (boot_df["contrast"] == "change_in_contrast_B_vs_A")
    ].set_index("outcome")
    wald_map = wald_df.set_index(["model", "outcome"])

    for ax, mod_key, panel_title in [
        (ax_a, "Model_A", "Model A: Conventional FA"),
        (ax_b, "Model_B", "Model B: FWE Tissue FA"),
    ]:
        sub = paths_df[paths_df["model"] == mod_key].copy()
        atr_rows = sub[sub["tract"] == "ATR"].set_index("lhs").loc[speed_order]
        fmin_rows = sub[sub["tract"] == "FMin"].set_index("lhs").loc[speed_order]

        ax.axvline(0, color=PALETTE["SLATE_DARK"], linestyle="--", linewidth=1.0)
        ax.errorbar(
            atr_rows["est.std"],
            y_f - offset,
            xerr=[
                atr_rows["est.std"] - atr_rows["ci.lower"],
                atr_rows["ci.upper"] - atr_rows["est.std"],
            ],
            fmt="o",
            color=PALETTE["BLUE_PRIMARY"],
            ecolor=PALETTE["BLUE_PRIMARY"],
            elinewidth=2.0,
            capsize=3.5,
            markersize=6.0,
            label="Anterior Thalamic Radiations (ATR)",
        )
        ax.errorbar(
            fmin_rows["est.std"],
            y_f + offset,
            xerr=[
                fmin_rows["est.std"] - fmin_rows["ci.lower"],
                fmin_rows["ci.upper"] - fmin_rows["est.std"],
            ],
            fmt="s",
            color=PALETTE["AMBER_ACCENT"],
            ecolor=PALETTE["AMBER_ACCENT"],
            elinewidth=2.0,
            capsize=3.5,
            markersize=6.0,
            label="Forceps Minor — Genu (FMin)",
        )

        for i, scode in enumerate(speed_order):
            r_atr = atr_rows.loc[scode]
            r_fmin = fmin_rows.loc[scode]
            sig_a = (
                "***"
                if r_atr["pvalue"] < 0.001
                else (
                    "**"
                    if r_atr["pvalue"] < 0.01
                    else ("*" if r_atr["pvalue"] < 0.05 else "")
                )
            )
            sig_f = (
                "***"
                if r_fmin["pvalue"] < 0.001
                else (
                    "**"
                    if r_fmin["pvalue"] < 0.01
                    else ("*" if r_fmin["pvalue"] < 0.05 else "")
                )
            )
            ax.text(
                r_atr["ci.upper"] + 0.012,
                y_f[i] - offset,
                f"{r_atr['est.std']:+.2f}{sig_a}",
                va="center",
                ha="left",
                fontsize=7.3,
                fontweight="bold" if r_atr["pvalue"] < 0.05 else "normal",
                color=PALETTE["BLUE_DARK"],
            )
            ax.text(
                r_fmin["ci.upper"] + 0.012,
                y_f[i] + offset,
                f"{r_fmin['est.std']:+.2f}{sig_f}",
                va="center",
                ha="left",
                fontsize=7.3,
                fontweight="bold" if r_fmin["pvalue"] < 0.05 else "normal",
                color=PALETTE["AMBER_DARK"],
            )

        ax.set_title(
            panel_title,
            fontsize=9.3,
            fontweight="bold",
            color=PALETTE["NAVY_DARK"],
            pad=8,
        )
        ax.set_xlim(-0.22, 0.68)
        ax.set_ylim(len(speed_order) - 0.35, -0.65)
        ax.grid(axis="x", alpha=0.7)
        ax.set_xlabel("Standardized β (±95% CI)", fontsize=8.6)

    ax_a.set_yticks(y_f)
    ax_a.set_yticklabels(speed_labels, fontweight="bold")
    ax_a.invert_yaxis()

    w_crt = wald_map.loc[("Model_B", "CRTspeed")]
    b_crt = boot_unadj.loc["CRTspeed"]
    b_cons = boot_unadj.loc["CRTcons"]
    annot_txt = (
        f"Cross-Modality Shift (B vs. A, 5,000 Bootstraps):\n"
        f"• CRT Speed: ΔΔβ = {b_crt['estimate']:+.3f} [{b_crt['ci_lo']:+.3f}, {b_crt['ci_hi']:+.3f}], "
        f"Model B Wald χ² = {w_crt['wald_chisq']:.2f} (p = {w_crt['p']:.4f})\n"
        f"• CRT Cons:  ΔΔβ = {b_cons['estimate']:+.3f} [{b_cons['ci_lo']:+.3f}, {b_cons['ci_hi']:+.3f}]"
    )
    fig.text(
        0.5,
        -0.06,
        annot_txt,
        ha="center",
        va="top",
        fontsize=7.5,
        color=PALETTE["NAVY_DARK"],
        linespacing=1.3,
        bbox={
            "boxstyle": "round,pad=0.35",
            "fc": PALETTE["BG_SUBTLE"],
            "ec": PALETTE["SLATE_LIGHT"],
        },
    )

    fig.suptitle(
        "Anterior White Matter Specificity Across Modalities: ATR vs. Forceps Minor",
        x=0.01,
        y=1.10,
        ha="left",
        fontsize=TYPO["TITLE"],
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
    )
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.05),
        ncol=2,
        frameon=True,
        fontsize=TYPO["LEGEND"],
        facecolor="#ffffff",
        edgecolor=PALETTE["SLATE_LIGHT"],
    )
    save_figure(fig, output_dir, "fig09_anterior_system_specificity")


# ------------------------------------------------------------------------------
# 11. Figure 10: Structural Variance Funneling (Unadjusted vs. Adjusted R²)
# ------------------------------------------------------------------------------
def plot_fig10(r2_unadj_df, r2_adj_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 5.0), dpi=300)
    u_map = (
        r2_unadj_df[r2_unadj_df["model"] == "Model_B"].set_index("var")["r2"].to_dict()
    )
    a_map = (
        r2_adj_df[r2_adj_df["model"] == "adj_age_sex_vent_Model_B"]
        .set_index("var")["r2"]
        .to_dict()
    )

    items = [
        ("FI", "Downstream Latent Fluid Intelligence (FI)"),
        ("CRTspeed", "Choice RT Speed (CRTspeed)"),
        ("CRTcons", "Choice RT Consistency (CRTcons)"),
        ("SRTspeed", "Simple RT Speed (SRTspeed)"),
        ("SRTcons", "Simple RT Consistency (SRTcons)"),
        ("AVcons", "Audio-Visual Consistency (AVcons)"),
        ("AVspeed", "Audio-Visual Speed (AVspeed)"),
    ]
    codes = [x[0] for x in items]
    labels = [x[1] for x in items]
    unadj_vals = [u_map.get(c, 0.0) * 100.0 for c in codes]
    adj_vals = [a_map.get(c, 0.0) * 100.0 for c in codes]

    y = np.arange(len(items))
    h = 0.36
    b1 = ax.barh(
        y - h / 2,
        unadj_vals,
        height=h,
        color=PALETTE["BLUE_PRIMARY"],
        edgecolor=PALETTE["BLUE_DARK"],
        linewidth=0.8,
        label="Unadjusted Lifespan Model B (N = 513)",
    )
    b2 = ax.barh(
        y + h / 2,
        adj_vals,
        height=h,
        color=PALETTE["TEAL_SECONDARY"],
        edgecolor="#115e59",
        linewidth=0.8,
        label="Age + Age² + Sex + Ventricle/eTIV Adjusted",
    )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Variance Explained (R² %)")
    ax.set_xlim(0, 72)
    ax.grid(axis="x", alpha=0.7)

    for bar in b1:
        w = bar.get_width()
        ax.text(
            w + 1.0,
            bar.get_y() + bar.get_height() / 2,
            f"{w:.2f}%",
            va="center",
            ha="left",
            fontsize=7.6,
            fontweight="bold",
            color=PALETTE["BLUE_DARK"],
        )
    for bar in b2:
        w = bar.get_width()
        ax.text(
            w + 1.0,
            bar.get_y() + bar.get_height() / 2,
            f"{w:.2f}%",
            va="center",
            ha="left",
            fontsize=7.6,
            fontweight="bold",
            color=PALETTE["TEAL_SECONDARY"],
        )

    apply_title_and_top_legend(
        ax,
        "Watershed Variance Funneling: Lifespan vs. Adjusted R² (%)",
        ncol=2,
        n_rows=1,
    )
    save_figure(fig, output_dir, "fig10_variance_funneling_r2")


# ------------------------------------------------------------------------------
# 12. Figure 11: CONSORT-Style Cohort Derivation & Attrition Audit Flowchart
# ------------------------------------------------------------------------------
def plot_fig11(demo_df, sex_df, att_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 9.0), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(
        "CONSORT-Style Cohort Derivation and Harmonization Flowchart",
        loc="left",
        fontsize=TYPO["TITLE"],
        fontweight="bold",
        color=PALETTE["NAVY_DARK"],
        pad=12,
    )

    d_row = demo_df.iloc[0]
    sex_counts = sex_df.set_index("Var1")["Freq"].to_dict()
    n_fem = int(sex_counts.get("F", 253))
    n_mal = int(sex_counts.get("M", 260))

    att_map = att_df.set_index("group")
    ex_all_age = att_map.loc["excluded_all", "mean_age"]
    ex_all_p = att_map.loc["excluded_all", "wilcox_p_vs_retained"]
    ex_tr_age = att_map.loc["excluded_tract_outlier", "mean_age"]
    ex_dwi_age = att_map.loc["excluded_dwi_incomplete", "mean_age"]

    stages = [
        (
            0.905,
            0.050,
            "1. Initial Cam-CAN Stage 2 Intake Cohort",
            "N = 653 participants with demographic & structural MRI records",
            PALETTE["BLUE_DARK"],
            "#ffffff",
            PALETTE["NAVY_DARK"],
            "-",
        ),
        (
            0.740,
            0.055,
            "2. Cross-Modality Sequential Inner Merge",
            "Demographics + Cattell FI + Simple/Choice/AV RT + JHU DWI + FreeSurfer aseg/eTIV\nN = 550 merged cross-modal participants (561 features)",
            PALETTE["BLUE_PRIMARY"],
            "#ffffff",
            PALETTE["BLUE_DARK"],
            "-",
        ),
        (
            0.565,
            0.053,
            "3. Behavioral Quality Control Screening",
            "Valid task completion flags across Simple, Choice, and Audio-Visual RT paradigms\nN = 547 behaviorally valid participants (3 excluded; 99.45% retained)",
            PALETTE["BLUE_PRIMARY"],
            "#ffffff",
            PALETTE["BLUE_DARK"],
            "-",
        ),
        (
            0.350,
            0.084,
            "4. Multi-Compartment Fitting & Outlier Attrition Audit",
            f"• 4-SD outlier screen (|z| > 4.0): 25 unique participants flagged (4 behavioral, 22 tract; 1 overlap)\n"
            f"• Multi-compartment DWI convergence audit: 9 participants excluded (FWE / WMTI boundary failures)\n"
            f"• Attrition Age Audit (Wilcoxon vs. retained M = {d_row['age_mean']:.2f} yrs): Excluded n = 34 older\n"
            f"  (M = {ex_all_age:.2f} yrs, p = {ex_all_p:.2e}; tract outliers M = {ex_tr_age:.2f} yrs; DWI incomplete M = {ex_dwi_age:.2f} yrs)",
            PALETTE["BG_SUBTLE"],
            PALETTE["NAVY_DARK"],
            PALETTE["SLATE_DARK"],
            "--",
        ),
        (
            0.095,
            0.068,
            "5. Final Harmonized Complete-Case SEM Cohort",
            f"Complete 4-Modality DWI + FreeSurfer Ventricles/eTIV + Valid Behavior + Zero 4-SD Outliers\n"
            f"N = {int(d_row['n'])} participants ({n_fem} females, {n_mal} males; Age: {d_row['age_mean']:.2f} ± {d_row['age_sd']:.2f} yrs, Range: {d_row['age_min']:.1f}–{d_row['age_max']:.1f})\n"
            f"Strict sample parity across Models A–D, Covariate-Adjusted SEMs, and Sensitivity Regimes",
            PALETTE["GREEN_WIN"],
            "#ffffff",
            PALETTE["GREEN_DARK"],
            "-",
        ),
    ]

    for cy, hh, title, desc, bg, fg, ec, ls in stages:
        box = patches.FancyBboxPatch(
            (0.03, cy - hh),
            0.94,
            2 * hh,
            boxstyle="round,pad=0.008,rounding_size=0.018",
            facecolor=bg,
            edgecolor=ec,
            linewidth=1.5,
            linestyle=ls,
        )
        ax.add_patch(box)
        title_y = cy + hh * 0.54 if hh > 0.06 else cy + hh * 0.44
        desc_y = cy - hh * 0.20 if hh > 0.06 else cy - hh * 0.24
        ax.text(
            0.50,
            title_y,
            title,
            ha="center",
            va="center",
            fontsize=9.1,
            fontweight="bold",
            color=fg,
        )
        ax.text(
            0.50,
            desc_y,
            desc,
            ha="center",
            va="center",
            fontsize=6.9 if hh > 0.07 else 7.8,
            color=fg if fg == PALETTE["NAVY_DARK"] else "#f1f5f9",
            linespacing=1.32,
        )

    pad = 0.010
    for i in range(len(stages) - 1):
        y_start = stages[i][0] - stages[i][1] - pad
        y_end = stages[i + 1][0] + stages[i + 1][1] + pad
        ax.annotate(
            "",
            xy=(0.50, y_end),
            xytext=(0.50, y_start),
            arrowprops={
                "arrowstyle": "-|>,head_length=0.6,head_width=0.35",
                "color": PALETTE["SLATE_DARK"],
                "lw": 1.8,
            },
        )
    save_figure(fig, output_dir, "fig11_cohort_attrition_flowchart")


# ------------------------------------------------------------------------------
# 13. Figure 12: Standardized Variable Extremes Relative to 4-SD Rule
# ------------------------------------------------------------------------------
def plot_fig12(prep_df, output_dir):
    fig, ax = plt.subplots(figsize=(TYPO["W_SINGLE"], 7.2), dpi=300)
    feature_labels = [
        ("C1", "Cattell C1: Series Completion"),
        ("C2", "Cattell C2: Classification"),
        ("C3", "Cattell C3: Matrices"),
        ("C4", "Cattell C4: Conditions / Topology"),
        ("SRTspeed", "Simple RT Speed"),
        ("SRTcons", "Simple RT Consistency"),
        ("CRTspeed", "Choice RT Speed"),
        ("CRTcons", "Choice RT Consistency"),
        ("AVspeed", "Audio-Visual Cued Speed"),
        ("AVcons", "Audio-Visual Cued Consistency"),
        ("FWE_FAtissue_FMin", "FAtissue: Forceps Minor (FMin)"),
        ("FWE_FAtissue_FMaj", "FAtissue: Forceps Major (FMaj)"),
        ("FWE_FAtissue_CST", "FAtissue: Corticospinal Tract (CST)"),
        ("FWE_FAtissue_ATR", "FAtissue: Anterior Thalamic (ATR)"),
        ("FWE_FAtissue_SLF", "FAtissue: Sup. Longitudinal (SLF)"),
        ("FWE_FAtissue_IFOF", "FAtissue: Inf. Fronto-Occipital (IFOF)"),
        ("FWE_FAtissue_ILF", "FAtissue: Inf. Longitudinal (ILF)"),
        ("FWE_FAtissue_UNC", "FAtissue: Uncinate Fasciculus (UNC)"),
        ("FWE_FAtissue_CING", "FAtissue: Cingulate Cingulum (CING)"),
        ("FWE_FAtissue_CINGHipp", "FAtissue: Hippocampal Cingulum"),
    ]
    codes = [x[0] for x in feature_labels]
    labels = [x[1] for x in feature_labels]
    prep_sub = prep_df.set_index("Feature").loc[codes].reset_index()
    y_p = np.arange(len(prep_sub))

    ax.axvspan(
        -4.0,
        4.0,
        color=PALETTE["GREEN_FILL"],
        alpha=0.85,
        zorder=0,
        label="Included Harmonized Zone (|z| ≤ 4.0 SD)",
    )
    ax.axvline(
        -4.0,
        color=PALETTE["RED_NEG"],
        linestyle="--",
        linewidth=1.3,
        label="±4.0 SD Outlier Screening Threshold",
    )
    ax.axvline(+4.0, color=PALETTE["RED_NEG"], linestyle="--", linewidth=1.3)
    ax.axvline(0.0, color=PALETTE["SLATE_MUTED"], linestyle=":", linewidth=1.0)

    for i, row in prep_sub.iterrows():
        ax.hlines(
            y_p[i],
            row["Min"],
            row["Max"],
            color=PALETTE["BLUE_PRIMARY"],
            linewidth=2.3,
            zorder=2,
        )
        ax.plot(
            row["Min"],
            y_p[i],
            "|",
            color=PALETTE["RED_NEG"] if row["Min"] < -4.0 else PALETTE["BLUE_DARK"],
            markersize=8,
            markeredgewidth=2.0,
            zorder=3,
        )
        ax.plot(
            row["Max"],
            y_p[i],
            "|",
            color=PALETTE["RED_NEG"] if row["Max"] > 4.0 else PALETTE["BLUE_DARK"],
            markersize=8,
            markeredgewidth=2.0,
            zorder=3,
        )
        ax.plot(
            row["Mean"],
            y_p[i],
            "o",
            color=PALETTE["NAVY_DARK"],
            markersize=4.0,
            zorder=4,
        )

    ax.set_yticks(y_p)
    ax.set_yticklabels(labels, fontsize=8.4, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Pre-Filter Standardized Z-Score Range [Min to Max] (N = 547)")
    ax.set_xlim(-14.8, 7.2)
    ax.set_ylim(len(prep_sub) - 0.4, -0.6)
    ax.grid(axis="x", alpha=0.7)

    apply_title_and_top_legend(
        ax,
        "Pre-Filter Manifest Z-Score Extremes Relative to ±4.0 SD Rule",
        ncol=2,
        n_rows=1,
    )
    save_figure(fig, output_dir, "fig12_outlier_4sd_range_profile")


# ------------------------------------------------------------------------------
# 14. Main Orchestrator
# ------------------------------------------------------------------------------
def resolve_data_dir(user_dir):
    """Automatically locate the directory containing wa.R CSV outputs."""
    for candidate in [user_dir, "output", "."]:
        if os.path.exists(os.path.join(candidate, "watershed_model_comparison.csv")):
            return candidate
    return user_dir


def main():
    parser = argparse.ArgumentParser(
        description="Generate all 12 standalone vector figures (.pdf and .svg) from wa.R outputs."
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        default="output",
        help="Directory containing wa.R output CSVs (default: auto-detects 'output' or '.')",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output_figures",
        help="Target directory for exported PDF and SVG figures (default: output_figures)",
    )
    args = parser.parse_args()

    data_dir = resolve_data_dir(args.data_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"[SETUP] Reading CSVs from: {data_dir}")
    print(f"[SETUP] Saving figures to: {args.output_dir}")

    req_files = {
        "comp": "watershed_model_comparison.csv",
        "dir": "directionality_model_comparison.csv",
        "std": "standardized_solutions_all_models.csv",
        "paths": "tract_to_speed_paths_unadjusted.csv",
        "wald": "wald_ATR_vs_FMin_unadjusted.csv",
        "boot": "bootstrap_ATR_FMin_contrasts.csv",
        "r2_u": "rsquared_unadjusted.csv",
        "r2_a": "rsquared_adjusted.csv",
        "prep": "preprocessing_summary.csv",
        "cfa": "cfa_fi_fit.csv",
        "demo": "cohort_demographics.csv",
        "sex": "cohort_sex_counts.csv",
        "att": "excluded_vs_retained_age.csv",
    }

    dfs = {}
    for key, fname in req_files.items():
        fpath = os.path.join(data_dir, fname)
        if not os.path.exists(fpath):
            print(f"[ERROR] Required file missing: {fpath}", file=sys.stderr)
            sys.exit(1)
        dfs[key] = pd.read_csv(fpath)

    print("[DATA] Generating 12 standalone vector figures (.pdf + .svg)...")
    plot_fig01(dfs["std"], dfs["paths"], dfs["r2_u"], args.output_dir)
    plot_fig02(dfs["std"], dfs["cfa"], args.output_dir)
    plot_fig03(dfs["std"], args.output_dir)
    plot_fig04(dfs["comp"], args.output_dir)
    plot_fig05(dfs["comp"], args.output_dir)
    plot_fig06(dfs["dir"], args.output_dir)
    plot_fig07(dfs["dir"], args.output_dir)
    plot_fig08(dfs["paths"], dfs["r2_u"], args.output_dir)
    plot_fig09(dfs["paths"], dfs["wald"], dfs["boot"], args.output_dir)
    plot_fig10(dfs["r2_u"], dfs["r2_a"], args.output_dir)
    plot_fig11(dfs["demo"], dfs["sex"], dfs["att"], args.output_dir)
    plot_fig12(dfs["prep"], args.output_dir)

    print("\n[SUCCESS] All 12 figures exported in PDF and SVG formats.")


if __name__ == "__main__":
    main()

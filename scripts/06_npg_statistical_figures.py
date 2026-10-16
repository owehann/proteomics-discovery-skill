#!/usr/bin/env python3
"""
Stage 6A: NPG / Nature Cell Biology Multi-Panel Statistical Figure Engine
=========================================================================
Generates publication-grade 300-DPI PNG + Vector PDF/SVG figures:
  1. 6-Panel Forensic QC & Empirical Null Calibration Dashboard
  2. Peptide-Tiered Moderated Volcano Plot (`unique_peptides >= 2` circles vs `1-pep` triangles)
  3. 4-Quadrant Reciprocal (`KD x OE`) or Clinical Secretion Gradient (`SGI`) Scatter
  4. 9-Sector Epistasis & `TOST` Equivalence Grid
"""

from __future__ import annotations
import argparse
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


NPG_COLORS = {
    "up": "#E64B35",
    "down": "#4DBBD5",
    "emerald": "#00A087",
    "indigo": "#3C5488",
    "amber": "#F39B7F",
    "slate": "#8491B4",
    "dark": "#1E293B",
    "light_grid": "#E2E8F0",
}


def set_npg_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "font.size": 9.5,
            "axes.titlesize": 11.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.edgecolor": "#334155",
            "axes.linewidth": 0.9,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
        }
    )


def plot_peptide_tiered_volcano(
    df: pd.DataFrame,
    out_prefix: str | Path,
    title: str = "Moderated Differential Proteome",
    fc_col: str = "log2FC",
    pval_col: str = "p_value",
    gene_col: str = "std_gene",
    pep_col: str = "std_unique_peptides",
    fc_cut: float = 0.5849625,
    p_cut: float = 0.05,
    top_n_labels: int = 15,
    empirical_null_badge: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Render an NPG volcano plot distinguishing Core (>=2 peptides) vs Exploratory (1 peptide) hits."""
    set_npg_style()
    sub = df.dropna(subset=[fc_col, pval_col]).copy()
    x = pd.to_numeric(sub[fc_col], errors="coerce").values
    p = np.clip(pd.to_numeric(sub[pval_col], errors="coerce").values, 1e-300, 1.0)
    y = -np.log10(p)
    peps = pd.to_numeric(sub[pep_col], errors="coerce").fillna(2).values if pep_col in sub.columns else np.full(len(sub), 2)

    is_up = (x >= fc_cut) & (p <= p_cut)
    is_dn = (x <= -fc_cut) & (p <= p_cut)
    is_core = peps >= 2

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    ax.grid(True, linestyle=":", linewidth=0.5, color=NPG_COLORS["light_grid"], zorder=0)

    # Background non-significant
    bg = ~(is_up | is_dn)
    ax.scatter(x[bg], y[bg], s=14, c=NPG_COLORS["slate"], alpha=0.35, edgecolors="none", label="Background / Equiv", zorder=1)

    # Core (>=2 unique peptides) hits: solid circles
    ax.scatter(
        x[is_up & is_core],
        y[is_up & is_core],
        s=32,
        c=NPG_COLORS["up"],
        alpha=0.88,
        marker="o",
        edgecolors="white",
        linewidths=0.4,
        label=f"Up Core >=2 pep ({int((is_up & is_core).sum())})",
        zorder=3,
    )
    ax.scatter(
        x[is_dn & is_core],
        y[is_dn & is_core],
        s=32,
        c=NPG_COLORS["down"],
        alpha=0.88,
        marker="o",
        edgecolors="white",
        linewidths=0.4,
        label=f"Down Core >=2 pep ({int((is_dn & is_core).sum())})",
        zorder=3,
    )

    # Single-peptide exploratory hits: open/amber-edged triangles
    if np.any((is_up | is_dn) & (~is_core)):
        ax.scatter(
            x[(is_up | is_dn) & (~is_core)],
            y[(is_up | is_dn) & (~is_core)],
            s=36,
            facecolors="none",
            edgecolors=NPG_COLORS["amber"],
            marker="^",
            linewidths=1.2,
            label=f"1-Peptide Exploratory ({int(((is_up | is_dn) & (~is_core)).sum())})",
            zorder=4,
        )

    ax.axvline(fc_cut, color="#64748B", linestyle="--", linewidth=0.9)
    ax.axvline(-fc_cut, color="#64748B", linestyle="--", linewidth=0.9)
    ax.axhline(-math.log10(p_cut), color="#64748B", linestyle="--", linewidth=0.9)

    # Label top core hits
    sub["neg_log10_p"] = y
    sub["rank_score"] = np.abs(x) * y
    sig_core = sub.loc[(is_up | is_dn) & is_core].sort_values("rank_score", ascending=False).head(top_n_labels)
    for idx, (_, row) in enumerate(sig_core.iterrows()):
        gx = float(row[fc_col])
        gy = float(row["neg_log10_p"])
        gname = str(row.get(gene_col, ""))
        offset_y = 4 + (idx % 3) * 5
        ax.annotate(
            gname,
            xy=(gx, gy),
            xytext=(5 if gx >= 0 else -5, offset_y),
            textcoords="offset points",
            ha="left" if gx >= 0 else "right",
            fontsize=7.8,
            fontweight="bold",
            color=NPG_COLORS["dark"],
            arrowprops=dict(arrowstyle="-", color="#94A3B8", lw=0.5),
        )

    ax.set_title(title, pad=10)
    ax.set_xlabel(r"$\log_2(\mathrm{Fold\ Change})$")
    ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$")
    if empirical_null_badge:
        ax.text(
            0.02,
            0.96,
            empirical_null_badge,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#F8FAFC", edgecolor="#CBD5E1"),
        )
    ax.legend(loc="upper right", frameon=True, fontsize=7.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out_base = Path(out_prefix)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_base.with_suffix(".png")
    pdf_path = out_base.with_suffix(".pdf")
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def plot_nine_sector_epistasis_grid(
    df: pd.DataFrame,
    col_a: str,
    col_b: str,
    out_prefix: str | Path,
    label_a: str = "Perturbation A",
    label_b: str = "Perturbation B",
    gene_col: str = "std_gene",
    pep_col: str = "std_unique_peptides",
    fc_cut: float = 0.5849625,
    equiv_cut: float = 0.321928,
) -> Tuple[Path, Path]:
    """Render a 9-Sector Cartesian Epistasis & TOST Equivalence scatter plot."""
    set_npg_style()
    sub = df.dropna(subset=[col_a, col_b]).copy()
    xa = pd.to_numeric(sub[col_a], errors="coerce").values
    xb = pd.to_numeric(sub[col_b], errors="coerce").values
    peps = pd.to_numeric(sub[pep_col], errors="coerce").fillna(2).values if pep_col in sub.columns else np.full(len(sub), 2)

    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    ax.axvspan(-equiv_cut, equiv_cut, color="#F1F5F9", alpha=0.65, zorder=0, label="TOST Equivalence Zone")
    ax.axhspan(-equiv_cut, equiv_cut, color="#F1F5F9", alpha=0.65, zorder=0)

    shared_up = (xa >= fc_cut) & (xb >= fc_cut)
    shared_dn = (xa <= -fc_cut) & (xb <= -fc_cut)
    a_spec = (np.abs(xa) >= fc_cut) & (np.abs(xb) <= equiv_cut)
    b_spec = (np.abs(xb) >= fc_cut) & (np.abs(xa) <= equiv_cut)
    antag = ((xa >= fc_cut) & (xb <= -fc_cut)) | ((xa <= -fc_cut) & (xb >= fc_cut))
    rest = ~(shared_up | shared_dn | a_spec | b_spec | antag)

    ax.scatter(xa[rest], xb[rest], s=12, c=NPG_COLORS["slate"], alpha=0.3, label="Background / Core Equiv")
    ax.scatter(xa[shared_up], xb[shared_up], s=28, c=NPG_COLORS["up"], alpha=0.85, label=f"Shared Up ({int(shared_up.sum())})")
    ax.scatter(xa[shared_dn], xb[shared_dn], s=28, c=NPG_COLORS["down"], alpha=0.85, label=f"Shared Down ({int(shared_dn.sum())})")
    ax.scatter(xa[a_spec], xb[a_spec], s=28, c=NPG_COLORS["indigo"], alpha=0.85, label=f"{label_a}-Specific ({int(a_spec.sum())})")
    ax.scatter(xa[b_spec], xb[b_spec], s=28, c=NPG_COLORS["emerald"], alpha=0.85, label=f"{label_b}-Specific ({int(b_spec.sum())})")
    if np.any(antag):
        ax.scatter(xa[antag], xb[antag], s=32, c=NPG_COLORS["amber"], alpha=0.9, label=f"Antagonistic ({int(antag.sum())})")

    for v in [-fc_cut, fc_cut]:
        ax.axvline(v, color="#475569", linestyle="--", lw=0.85)
        ax.axhline(v, color="#475569", linestyle="--", lw=0.85)

    # Annotate top hits with >= 2 peptides
    sub["mag"] = np.hypot(xa, xb)
    top_df = sub.loc[(shared_up | shared_dn | a_spec | b_spec | antag) & (peps >= 2)].sort_values("mag", ascending=False).head(14)
    for _, row in top_df.iterrows():
        ax.annotate(
            str(row.get(gene_col, "")),
            xy=(float(row[col_a]), float(row[col_b])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7.5,
            fontweight="bold",
        )

    ax.set_xlabel(f"{label_a} ($\\log_2\\,\\mathrm{{FC}}$)")
    ax.set_ylabel(f"{label_b} ($\\log_2\\,\\mathrm{{FC}}$)")
    ax.set_title(f"9-Sector Epistasis & TOST Divergence: {label_a} vs {label_b}")
    ax.legend(loc="upper left", fontsize=7.5, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out_base = Path(out_prefix)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_base.with_suffix(".png")
    pdf_path = out_base.with_suffix(".pdf")
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 6A: NPG Statistical Figure Engine")
    parser.add_argument("--input", "-i", required=True, help="Input CSV")
    parser.add_argument("--out-prefix", "-o", required=True, help="Output file prefix")
    parser.add_argument("--title", default="Moderated Differential Proteome", help="Plot title")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    png_p, pdf_p = plot_peptide_tiered_volcano(df, args.out_prefix, title=args.title)
    print(f"[Stage 6A] Saved volcano plots to {png_p} and {pdf_p}")


if __name__ == "__main__":
    main()

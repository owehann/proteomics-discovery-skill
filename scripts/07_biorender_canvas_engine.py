#!/usr/bin/env python3
"""
Stage 6B: Reusable BioRender-Style Vector SVG + PNG Mechanistic Canvas Engine (6-Module Suite)
==============================================================================================
Constructs publication-grade, clean vector SVG + 300-DPI PNG cellular & molecular schematics:
  - Reusable organelles & structural primitives:
      * Phospholipid bilayer plasma membranes (`draw_bilayer_membrane`)
      * Multivesicular bodies (`MVB`), budding exosomes (`EVs`), Golgi stack, Nucleus, Mitochondria
      * 4-TM MARVEL/Tetraspanin proteins (`draw_4tm_protein`), Integrin heterodimers (`draw_integrin_dimer`),
        GPI-anchored surface proteins (`draw_gpi_protein`), ECM fibrils (`draw_ecm_fibril`),
        and Biologic Y-shaped IgG antibodies (`draw_antibody_igg`)
      * Quantitative Protein/PTM pill nodes (`draw_protein_node`) with `1-pep` warning badges
      * Evidence-graded causal arrows (`draw_causal_arrow`: Solid for `E1/E2`, Dashed for `E3`)
  - 6 Modular Deck Schematics (`generate_all_biorender_schematics`):
      (a) `build_s01_biorender_hero`: Cover Graphical Abstract (`s01_biorender_hero`)
      (b) `build_s03_biorender_exp_design`: Experimental Design & Orthogonal Causal Locks (`s03_biorender_exp_design`)
      (c) `build_s14_biorender_topology`: Biophysical Membrane Topology / Subpopulation Sorting (`s14_biorender_topology`)
      (d) `build_s17_biorender_lock_and_key`: STRING Lock-and-Key Surface Corona / Complex Assembly (`s17_biorender_lock_and_key`)
      (e) `build_s24_biorender_spotlight`: 6 Flagship Target Subcellular Localization Spotlight Cards (`s24_biorender_spotlight`)
      (f) `build_s25_biorender_working_model`: Integrated Inside-Out / Cellular Working Model (`s25_biorender_working_model`)
"""

from __future__ import annotations
import argparse
import html
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd


class BioRenderCanvas:
    """Programmatic vector SVG builder with BioRender-style gradients, organelles, and evidence-graded arrows."""

    def __init__(self, width: int = 1400, height: int = 860, title: str = "Proteomics Mechanistic Working Model") -> None:
        self.width = width
        self.height = height
        self.title = title
        self.elements: List[str] = []

    def _defs_header(self) -> str:
        return """
  <defs>
    <filter id="soft-shadow" x="-10%" y="-10%" width="125%" height="125%">
      <feDropShadow dx="0" dy="3" stdDeviation="4" flood-color="#0F172A" flood-opacity="0.10"/>
    </filter>
    <linearGradient id="cytosol-bg" x1="0%" y1="0%" x2="0%" y2="100%">
      <stop offset="0%" stop-color="#F8FAFC"/>
      <stop offset="22%" stop-color="#EFF6FF"/>
      <stop offset="23%" stop-color="#FFFBEB"/>
      <stop offset="100%" stop-color="#FEFDF8"/>
    </linearGradient>
    <radialGradient id="mvb-grad" cx="45%" cy="40%" r="55%">
      <stop offset="0%" stop-color="#FEF3C7"/>
      <stop offset="100%" stop-color="#FDE68A"/>
    </radialGradient>
    <radialGradient id="ev-grad" cx="35%" cy="35%" r="65%">
      <stop offset="0%" stop-color="#E0F2FE"/>
      <stop offset="100%" stop-color="#7DD3FC"/>
    </radialGradient>
    <radialGradient id="ev-emerald" cx="35%" cy="35%" r="65%">
      <stop offset="0%" stop-color="#D1FAE5"/>
      <stop offset="100%" stop-color="#34D399"/>
    </radialGradient>
    <radialGradient id="nucleus-grad" cx="50%" cy="50%" r="55%">
      <stop offset="0%" stop-color="#EDE9FE"/>
      <stop offset="100%" stop-color="#C4B5FD"/>
    </radialGradient>
    <marker id="arrow-e1e2" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 0 1.5 L 9 5 L 0 8.5 z" fill="#0F766E"/>
    </marker>
    <marker id="arrow-e3" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 0 1.5 L 9 5 L 0 8.5 z" fill="#D97706"/>
    </marker>
    <marker id="inhibit-bar" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <line x1="5" y1="1" x2="5" y2="9" stroke="#DC2626" stroke-width="2.5"/>
    </marker>
  </defs>
"""

    def add_header_and_legend(self, subtitle: str = "") -> None:
        t_esc = html.escape(self.title)
        s_esc = html.escape(subtitle)
        self.elements.append(
            f'<rect x="0" y="0" width="{self.width}" height="{self.height}" fill="url(#cytosol-bg)"/>'
        )
        self.elements.append(
            f'<rect x="24" y="18" width="{self.width - 48}" height="64" rx="10" fill="#0F172A" filter="url(#soft-shadow)"/>'
            f'<text x="46" y="47" font-family="Arial, Helvetica, sans-serif" font-size="19" font-weight="bold" fill="#FFFFFF">{t_esc}</text>'
            f'<text x="46" y="68" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#CBD5E1">{s_esc}</text>'
        )
        lx = self.width - 490
        self.elements.append(
            f'<rect x="{lx}" y="26" width="450" height="48" rx="8" fill="#1E293B" stroke="#334155"/>'
            f'<line x1="{lx+16}" y1="42" x2="{lx+58}" y2="42" stroke="#14B8A6" stroke-width="3"/>'
            f'<text x="{lx+66}" y="46" font-family="Arial" font-size="10.5" fill="#F8FAFC">Solid Arrow: E1/E2 Calibrated Evidence (FDR &lt; 0.25, &gt;=2 pep)</text>'
            f'<line x1="{lx+16}" y1="61" x2="{lx+58}" y2="61" stroke="#F59E0B" stroke-width="2.4" stroke-dasharray="6,4"/>'
            f'<text x="{lx+66}" y="65" font-family="Arial" font-size="10.5" fill="#FDE68A">Dashed Arrow: E3 Exploratory Hypothesis / 1-Peptide (*)</text>'
        )

    def draw_card_panel(
        self,
        x: int,
        y: int,
        w: int,
        h: int,
        title: str = "",
        subtitle: str = "",
        fill: str = "#FFFFFF",
        stroke: str = "#CBD5E1",
        header_fill: Optional[str] = None,
    ) -> None:
        """Draw a clean compartmentalized card container with optional header bar."""
        self.elements.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="2" filter="url(#soft-shadow)"/>'
        )
        if header_fill and title:
            self.elements.append(
                f'<rect x="{x}" y="{y}" width="{w}" height="38" rx="14" fill="{header_fill}"/>'
                f'<text x="{x + 18}" y="{y + 24}" font-family="Arial" font-size="13.5" font-weight="bold" fill="#FFFFFF">{html.escape(title)}</text>'
            )
            if subtitle:
                self.elements.append(
                    f'<text x="{x + w - 16}" y="{y + 24}" text-anchor="end" font-family="Arial" font-size="11" fill="#E2E8F0">{html.escape(subtitle)}</text>'
                )
        elif title:
            self.elements.append(
                f'<text x="{x + 18}" y="{y + 26}" font-family="Arial" font-size="14" font-weight="bold" fill="#0F172A">{html.escape(title)}</text>'
            )
            if subtitle:
                self.elements.append(
                    f'<text x="{x + 18}" y="{y + 44}" font-family="Arial" font-size="11" fill="#475569">{html.escape(subtitle)}</text>'
                )

    def draw_bilayer_membrane(
        self,
        y_pos: int = 215,
        label_top: str = "EXTRACELLULAR / BIOFLUID / SECRETOME",
        label_bottom: str = "INTRACELLULAR CYTOSOL (WCL)",
        x_start: int = 40,
        x_end: Optional[int] = None,
    ) -> None:
        """Draw a continuous BioRender-style phospholipid bilayer across the canvas."""
        if x_end is None:
            x_end = self.width - 40
        self.elements.append(
            f'<rect x="{x_start}" y="{y_pos - 14}" width="{x_end - x_start}" height="28" fill="#FDE68A" fill-opacity="0.45"/>'
        )
        step = 16
        heads = []
        for x in range(x_start + 8, x_end - 4, step):
            heads.append(
                f'<line x1="{x-2}" y1="{y_pos-9}" x2="{x-2}" y2="{y_pos-1}" stroke="#D97706" stroke-width="1.2"/>'
                f'<line x1="{x+2}" y1="{y_pos-9}" x2="{x+2}" y2="{y_pos-1}" stroke="#D97706" stroke-width="1.2"/>'
                f'<line x1="{x-2}" y1="{y_pos+1}" x2="{x-2}" y2="{y_pos+9}" stroke="#D97706" stroke-width="1.2"/>'
                f'<line x1="{x+2}" y1="{y_pos+1}" x2="{x+2}" y2="{y_pos+9}" stroke="#D97706" stroke-width="1.2"/>'
                f'<circle cx="{x}" cy="{y_pos-12}" r="4.5" fill="#F59E0B" stroke="#B45309" stroke-width="1"/>'
                f'<circle cx="{x}" cy="{y_pos+12}" r="4.5" fill="#F59E0B" stroke="#B45309" stroke-width="1"/>'
            )
        self.elements.extend(heads)
        if label_top:
            self.elements.append(
                f'<text x="{x_start + 14}" y="{y_pos - 24}" font-family="Arial" font-size="11.5" font-weight="bold" fill="#0369A1" letter-spacing="0.8">{html.escape(label_top)}</text>'
            )
        if label_bottom:
            self.elements.append(
                f'<text x="{x_start + 14}" y="{y_pos + 34}" font-family="Arial" font-size="11.5" font-weight="bold" fill="#92400E" letter-spacing="0.8">{html.escape(label_bottom)}</text>'
            )

    def draw_mvb_and_ev(self, cx: int, cy: int, radius: int = 95, label: str = "Multivesicular Body (MVB)") -> None:
        """Draw an MVB containing intraluminal vesicles (ILVs) and secreted exosomes."""
        self.elements.append(
            f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="url(#mvb-grad)" stroke="#D97706" stroke-width="3" filter="url(#soft-shadow)"/>'
            f'<circle cx="{cx}" cy="{cy}" r="{radius - 6}" fill="none" stroke="#F59E0B" stroke-width="1.2" stroke-dasharray="4,3"/>'
        )
        ilv_offsets = [(-35, -25, 20), (28, -22, 18), (-18, 30, 19), (30, 26, 17), (0, 2, 21)]
        for dx, dy, r_ilv in ilv_offsets:
            self.elements.append(
                f'<circle cx="{cx + dx}" cy="{cy + dy}" r="{r_ilv}" fill="url(#ev-grad)" stroke="#0284C7" stroke-width="2"/>'
            )
        if label:
            self.elements.append(
                f'<text x="{cx}" y="{cy + radius + 22}" text-anchor="middle" font-family="Arial" font-size="12.5" font-weight="bold" fill="#78350F">{html.escape(label)}</text>'
            )

    def draw_ev_vesicle(
        self,
        cx: int,
        cy: int,
        radius: int = 52,
        label: str = "Secreted EV",
        sublabel: str = "",
        fill_grad: str = "url(#ev-grad)",
        stroke: str = "#0284C7",
    ) -> None:
        """Draw an extracellular vesicle with phospholipid double-rim and surface corona pegs."""
        self.elements.append(
            f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{fill_grad}" stroke="{stroke}" stroke-width="3" filter="url(#soft-shadow)"/>'
            f'<circle cx="{cx}" cy="{cy}" r="{max(radius - 6, 8)}" fill="none" stroke="{stroke}" stroke-width="1.2" stroke-dasharray="3,3"/>'
        )
        for ang_deg in range(0, 360, 45):
            rad = math.radians(ang_deg)
            x1 = cx + math.cos(rad) * radius
            y1 = cy + math.sin(rad) * radius
            x2 = cx + math.cos(rad) * (radius + 9)
            y2 = cy + math.sin(rad) * (radius + 9)
            self.elements.append(
                f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="2"/>'
                f'<circle cx="{x2:.1f}" cy="{y2:.1f}" r="3.2" fill="#F59E0B" stroke="{stroke}" stroke-width="1"/>'
            )
        if label:
            self.elements.append(
                f'<text x="{cx}" y="{cy + 4}" text-anchor="middle" font-family="Arial" font-size="11.5" font-weight="bold" fill="#0F172A">{html.escape(label)}</text>'
            )
        if sublabel:
            self.elements.append(
                f'<text x="{cx}" y="{cy + radius + 24}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{stroke}">{html.escape(sublabel)}</text>'
            )

    def draw_4tm_protein(self, cx: int, cy: int, label: str = "CD63 / CD81 (4-TM)", color: str = "#059669") -> None:
        """Draw a 4-pass transmembrane MARVEL/Tetraspanin protein spanning a membrane at `cy`."""
        xs = [cx - 21, cx - 7, cx + 7, cx + 21]
        # Extracellular EC1 & EC2 loops
        self.elements.append(
            f'<path d="M {xs[0]} {cy-16} Q {cx-14} {cy-36} {xs[1]} {cy-16}" fill="none" stroke="{color}" stroke-width="2.5"/>'
            f'<path d="M {xs[2]} {cy-16} Q {cx+14} {cy-46} {xs[3]} {cy-16}" fill="none" stroke="{color}" stroke-width="3"/>'
            f'<path d="M {xs[1]} {cy+16} Q {cx} {cy+32} {xs[2]} {cy+16}" fill="none" stroke="{color}" stroke-width="2.2"/>'
        )
        for xh in xs:
            self.elements.append(
                f'<rect x="{xh-4}" y="{cy-18}" width="8" height="36" rx="4" fill="{color}" stroke="#FFFFFF" stroke-width="1"/>'
            )
        if label:
            self.elements.append(
                f'<text x="{cx}" y="{cy+44}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{color}">{html.escape(label)}</text>'
            )

    def draw_integrin_dimer(self, cx: int, cy: int, label: str = "ITGA/ITGB1 (1-TM)", color_a: str = "#2563EB", color_b: str = "#DC2626") -> None:
        """Draw an alpha/beta single-pass transmembrane integrin heterodimer with ectodomain headpieces."""
        self.elements.append(
            f'<rect x="{cx-12}" y="{cy-22}" width="9" height="42" rx="4" fill="{color_a}" stroke="#FFFFFF" stroke-width="1"/>'
            f'<rect x="{cx+3}" y="{cy-22}" width="9" height="42" rx="4" fill="{color_b}" stroke="#FFFFFF" stroke-width="1"/>'
            f'<ellipse cx="{cx-10}" cy="{cy-34}" rx="12" ry="16" fill="{color_a}" stroke="#FFFFFF" stroke-width="1.2"/>'
            f'<ellipse cx="{cx+10}" cy="{cy-34}" rx="12" ry="16" fill="{color_b}" stroke="#FFFFFF" stroke-width="1.2"/>'
        )
        if label:
            self.elements.append(
                f'<text x="{cx}" y="{cy+36}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="#1E293B">{html.escape(label)}</text>'
            )

    def draw_gpi_protein(self, cx: int, cy: int, label: str = "GPC3 (GPI-Anchor)", color: str = "#7C3AED") -> None:
        """Draw an outer-leaflet GPI-anchored surface proteoglycan."""
        head_lbl = label.split()[0] if label.split() else "GPI"
        self.elements.append(
            f'<path d="M {cx} {cy-12} L {cx-5} {cy-22} L {cx+5} {cy-30} L {cx} {cy-40}" fill="none" stroke="#D97706" stroke-width="2.5"/>'
            f'<ellipse cx="{cx}" cy="{cy-52}" rx="24" ry="14" fill="{color}" stroke="#FFFFFF" stroke-width="1.5" filter="url(#soft-shadow)"/>'
            f'<text x="{cx}" y="{cy-48}" text-anchor="middle" font-family="Arial" font-size="10" font-weight="bold" fill="#FFFFFF">{html.escape(head_lbl)}</text>'
            f'<text x="{cx}" y="{cy+32}" text-anchor="middle" font-family="Arial" font-size="10.5" font-weight="bold" fill="{color}">{html.escape(label)}</text>'
        )

    def draw_golgi_stack(self, cx: int, cy: int, label: str = "ER / Trans-Golgi Secretory Network") -> None:
        """Draw stacked Golgi cisternae."""
        for i, w_cis in enumerate([110, 96, 82, 68]):
            yk = cy - 24 + i * 15
            self.elements.append(
                f'<rect x="{cx - w_cis//2}" y="{yk}" width="{w_cis}" height="10" rx="5" fill="#FDE68A" stroke="#D97706" stroke-width="1.8"/>'
            )
        if label:
            self.elements.append(
                f'<text x="{cx}" y="{cy + 48}" text-anchor="middle" font-family="Arial" font-size="11.5" font-weight="bold" fill="#92400E">{html.escape(label)}</text>'
            )

    def draw_ecm_fibril(self, x1: int, y1: int, x2: int, y2: int, label: str = "FN1 / THBS1 / LGALS3BP ECM Corona", color: str = "#DC2626") -> None:
        """Draw cross-linked extracellular matrix fibrils."""
        mx = (x1 + x2) / 2.0
        my = (y1 + y2) / 2.0
        self.elements.append(
            f'<path d="M {x1} {y1} Q {mx} {my-14} {x2} {y2}" fill="none" stroke="{color}" stroke-width="3"/>'
            f'<path d="M {x1} {y1+6} Q {mx} {my+12} {x2} {y2+6}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="5,3"/>'
        )
        if label:
            self.elements.append(
                f'<text x="{mx}" y="{my - 14}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{color}">{html.escape(label)}</text>'
            )

    def draw_organelle(self, cx: int, cy: int, rx: int, ry: int, kind: str = "nucleus", label: str = "Nucleus") -> None:
        grad = "url(#nucleus-grad)" if kind == "nucleus" else "#DCFCE7"
        stroke = "#6D28D9" if kind == "nucleus" else "#15803D"
        self.elements.append(
            f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" fill="{grad}" stroke="{stroke}" stroke-width="2.5" filter="url(#soft-shadow)"/>'
            f'<text x="{cx}" y="{cy + ry + 20}" text-anchor="middle" font-family="Arial" font-size="12" font-weight="bold" fill="{stroke}">{html.escape(label)}</text>'
        )

    def draw_antibody_igg(self, cx: int, cy: int, label: str = "Humanized IgG", color: str = "#2563EB") -> None:
        """Draw a BioRender-style Y-shaped IgG biologic receptor/ligand binder."""
        self.elements.append(
            f'<g transform="translate({cx},{cy})" filter="url(#soft-shadow)">'
            f'  <line x1="0" y1="0" x2="0" y2="32" stroke="{color}" stroke-width="6" stroke-linecap="round"/>'
            f'  <line x1="0" y1="2" x2="-20" y2="-20" stroke="{color}" stroke-width="6" stroke-linecap="round"/>'
            f'  <line x1="0" y1="2" x2="20" y2="-20" stroke="{color}" stroke-width="6" stroke-linecap="round"/>'
            f'  <circle cx="-22" cy="-22" r="4" fill="#F59E0B"/>'
            f'  <circle cx="22" cy="-22" r="4" fill="#F59E0B"/>'
            f'  <text x="0" y="-32" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{color}">{html.escape(label)}</text>'
            f'</g>'
        )

    def draw_protein_node(
        self,
        cx: int,
        cy: int,
        gene: str,
        metric_text: str = "",
        direction: str = "up",
        unique_peptides: int = 2,
        evidence_grade: str = "E2",
        width: int = 152,
        height: int = 44,
    ) -> None:
        """Draw a quantitative protein pill badge; automatically adds '*' if unique_peptides == 1."""
        palette = {
            "up": ("#FEF2F2", "#DC2626", "#991B1B"),
            "down": ("#F0F9FF", "#0284C7", "#075985"),
            "reciprocal": ("#ECFDF5", "#059669", "#065F46"),
            "neutral": ("#F8FAFC", "#64748B", "#334155"),
            "warning": ("#FFFBEB", "#D97706", "#92400E"),
        }
        fill_c, stroke_c, text_c = palette.get(direction, palette["neutral"])
        dash = 'stroke-dasharray="5,3"' if (unique_peptides < 2 or evidence_grade == "E3") else ""
        pep_badge = " [1-pep*]" if unique_peptides < 2 else ""
        g_disp = html.escape(f"{gene}{pep_badge}")
        m_disp = html.escape(f"{metric_text} ({evidence_grade})" if metric_text else evidence_grade)

        x0 = cx - width // 2
        y0 = cy - height // 2
        self.elements.append(
            f'<g filter="url(#soft-shadow)">'
            f'  <rect x="{x0}" y="{y0}" width="{width}" height="{height}" rx="12" fill="{fill_c}" stroke="{stroke_c}" stroke-width="2.2" {dash}/>'
            f'  <text x="{cx}" y="{cy - 3}" text-anchor="middle" font-family="Arial" font-size="12.0" font-weight="bold" fill="{text_c}">{g_disp}</text>'
            f'  <text x="{cx}" y="{cy + 13}" text-anchor="middle" font-family="Arial" font-size="9.8" fill="#475569">{m_disp}</text>'
            f'</g>'
        )

    def draw_causal_arrow(
        self,
        x1: int,
        y1: int,
        x2: int,
        y2: int,
        evidence_grade: str = "E2",
        interaction: str = "activate",
        label: str = "",
        curve_offset: int = 0,
    ) -> None:
        """Draw a curved Bezier causal arrow (Solid = E1/E2; Dashed = E3)."""
        is_exploratory = evidence_grade.upper().startswith("E3")
        stroke_c = "#D97706" if is_exploratory else ("#DC2626" if interaction == "inhibit" else "#0F766E")
        dash = 'stroke-dasharray="7,5"' if is_exploratory else ""
        marker = "url(#inhibit-bar)" if interaction == "inhibit" else ("url(#arrow-e3)" if is_exploratory else "url(#arrow-e1e2)")

        mx = (x1 + x2) / 2.0
        my = (y1 + y2) / 2.0 + curve_offset
        self.elements.append(
            f'<path d="M {x1} {y1} Q {mx} {my} {x2} {y2}" fill="none" stroke="{stroke_c}" stroke-width="2.5" {dash} marker-end="{marker}"/>'
        )
        if label:
            tag = f"{label} [E3]" if is_exploratory and "E3" not in label else label
            self.elements.append(
                f'<text x="{mx}" y="{my - 8}" text-anchor="middle" font-family="Arial" font-size="10.5" font-weight="bold" fill="{stroke_c}">{html.escape(tag)}</text>'
            )

    def to_svg(self) -> str:
        body = "\n  ".join(self.elements)
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.width} {self.height}" width="100%" height="100%">\n'
            f"{self._defs_header()}\n"
            f"  {body}\n"
            f"</svg>\n"
        )

    def save(self, svg_path: str | Path, png_path: Optional[str | Path] = None) -> Path:
        out_svg = Path(svg_path)
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        svg_text = self.to_svg()
        out_svg.write_text(svg_text, encoding="utf-8")

        if png_path is None:
            png_path = out_svg.with_suffix(".png")
        out_png = Path(png_path)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        rendered = False
        if shutil.which("rsvg-convert"):
            res = subprocess.run(
                ["rsvg-convert", "-d", "300", "-p", "300", "-w", "2560", "-f", "png", "-o", str(out_png), str(out_svg)],
                check=False,
            )
            rendered = (res.returncode == 0 and out_png.exists())
        if not rendered:
            try:
                import cairosvg
                cairosvg.svg2png(bytestring=svg_text.encode("utf-8"), write_to=str(out_png), dpi=300)
            except Exception:
                pass
        return out_svg


# =============================================================================
# 6 Modular BioRender Schematic Builders (Requirement 3)
# =============================================================================

def build_s01_biorender_hero(
    out_dir: str | Path,
    mode: str = "EV",
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (a) `s01_biorender_hero`: Cover Graphical Abstract illustrating the end-to-end biological
    discovery and 5-Layer Multi-Cohort Evidence Integration.
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    if mode_u == "EV":
        t = title or "Graphical Abstract: Vesicle Trafficking & Membrane Scaffolds Orchestrate EV Cargo Sorting"
        s = subtitle or "Dual-Track DEqMS + <LOD Switches | 9-Sector Epistasis | UniProt Membrane Topology | CPTAC & Fluid Bridge"
    elif "WCL" in mode_u or "PERTURB" in mode_u:
        t = title or "Graphical Abstract: Genetic Perturbation Rewires Redox, Ferritinophagy & Secretory Proteome"
        s = subtitle or "Multi-Clone KO Concordance | NFE2L2/NRF2 Glutathione Axis | NCOA4-FTH1/FTL Iron Homeostasis | 5-Layer RRA"
    elif "BIOLOGIC" in mode_u or "MAB" in mode_u:
        t = title or "Graphical Abstract: Biologic Engagement Triggers IFN-STAT1 & Chromatin-Replisome Remodeling"
        s = subtitle or "Exogenous Biologic Verification | Dual-Treatment Epistasis | Chromatin-Replisome Suppression | CPTAC & DepMap"
    else:
        t = title or "Graphical Abstract: Inside-Out Multi-Cohort Proteomics Discovery Architecture"
        s = subtitle or "Integrating Intracellular Whole-Cell Proteomes with Extracellular Vesicles & Clinical Cohorts"

    c = BioRenderCanvas(width=1400, height=820, title=t)
    c.add_header_and_legend(s)

    # 3 Horizontal Stage Columns
    c.draw_card_panel(36, 105, 420, 680, "1. Perturbation & Dual-Track MS Engine", "Stage 0-3: Forensic QC, DEqMS & <LOD Switches", fill="#F8FAFC", stroke="#94A3B8", header_fill="#1E293B")
    c.draw_card_panel(490, 105, 420, 680, "2. Biophysical & Network Mechanism", "Stage 4: Membrane Topology, CORUM & STRING KDA", fill="#FFFFFF", stroke="#0284C7", header_fill="#0369A1")
    c.draw_card_panel(944, 105, 420, 680, "3. Clinical Cohorts & 5-Layer Priority", "Stage 5: CPTAC Cohorts, Fluid SGI & DepMap CRISPR", fill="#F0FDF4", stroke="#059669", header_fill="#047857")

    if mode_u == "EV":
        c.draw_mvb_and_ev(246, 290, radius=82, label="Donor Cell MVB & ILV Budding")
        c.draw_protein_node(145, 470, "RAB27A", "shKD <LOD / OE +2.3", "reciprocal", 8, "E1")
        c.draw_protein_node(345, 470, "PDCD6IP", "4-TM / ESCRT Gate", "down", 5, "E1")
        c.draw_protein_node(246, 560, "Track B Switches", "Binary <LOD ON/OFF", "up", 4, "E1", width=210)
        c.draw_protein_node(246, 635, "9-Sector Epistasis", "G1-G7 + Control-Free DI", "neutral", 6, "E1", width=230)

        c.draw_bilayer_membrane(y_pos=280, label_top="EV SURFACE CORONA", label_bottom="VESICLE LUMEN", x_start=515, x_end=885)
        c.draw_integrin_dimer(600, 280, "ITGB1/ITGAV (1-TM)")
        c.draw_4tm_protein(780, 280, "CD63/CD81 (4-TM)")
        c.draw_ecm_fibril(540, 195, 860, 195, "FN1 / THBS1 / LGALS3BP")
        c.draw_protein_node(700, 440, "Topology Split", "SP/1-TM OR=4.8x vs 4-TM", "up", 6, "E1", width=250)
        c.draw_protein_node(700, 525, "ESCRT vs Tetraspanin", "TSG101/VPS4A vs CD151", "reciprocal", 8, "E2", width=250)
        c.draw_protein_node(700, 610, "STRING Lock & Key", "ITGB1-CD151-FN1 KDA Hubs", "up", 9, "E1", width=250)

        c.draw_ev_vesicle(1154, 255, radius=55, label="Secreted EV", sublabel="Paired Biofluid + CPTAC")
        c.draw_protein_node(1154, 415, "Secretome Gradient", "Fluid >> Tissue (SGI>0.5)", "up", 9, "E1", width=250)
        c.draw_protein_node(1154, 495, "DepMap CRISPR Panel", "Paracrine vs Autonomous", "reciprocal", 6, "E1", width=250)
        c.draw_protein_node(1154, 575, "3-Lens x 5-Layer RRA", "Top Core & Lead Flagships", "up", 8, "E1", width=250)
        c.draw_protein_node(1154, 655, "Liquid Biopsy Panel", "ITGB1 / CD63 / LGALS3BP", "warning", 7, "E1", width=255)
    else:
        c.draw_organelle(246, 285, rx=105, ry=68, kind="cytosol", label="Perturbed Intracellular WCL Proteome")
        c.draw_protein_node(145, 445, "EGFR / KRAS", "Direct Target Verified", "down", 9, "E1")
        c.draw_protein_node(345, 445, "MAPK3 / STAT1", "Paralog & Pathway Shift", "up", 7, "E1")
        c.draw_protein_node(246, 545, "DEqMS + <LOD Track B", "Calibrated Empirical Null", "reciprocal", 6, "E1", width=230)
        c.draw_protein_node(246, 630, "9-Sector Concordance", "Clone/Biologics Epistasis", "neutral", 8, "E1", width=230)

        c.draw_organelle(700, 265, rx=120, ry=72, kind="nucleus", label="Chromatin & TF Regulons (CollecTRI)")
        c.draw_protein_node(700, 420, "NFE2L2 / STAT1 ULM", "Signed TF Activity Inference", "up", 8, "E1", width=250)
        c.draw_protein_node(700, 505, "CORUM Stoichiometry", "Proteasome / MCM / OXPHOS", "reciprocal", 9, "E1", width=250)
        c.draw_protein_node(700, 590, "STRING KDA Network", "Degree-Corrected RWR Bridges", "up", 8, "E1", width=250)

        c.draw_ev_vesicle(1154, 255, radius=55, label="Inside-Out Link", sublabel="WCL -> EV -> Biofluid")
        c.draw_protein_node(1154, 415, "CPTAC Clinical Cohort", "Paired Tumor vs Normal + RNA", "up", 9, "E1", width=250)
        c.draw_protein_node(1154, 495, "DepMap CRISPR Panel", "Autonomous vs Paracrine Split", "reciprocal", 7, "E1", width=250)
        c.draw_protein_node(1154, 575, "3-Lens x 5-Layer RRA", "1,000x Weight Robustness", "up", 8, "E1", width=250)
        c.draw_protein_node(1154, 655, "Top 6 Flagship Drivers", "Biomarker & Combo Targets", "warning", 8, "E1", width=250)

    c.draw_causal_arrow(456, 440, 490, 440, "E1", "activate", "")
    c.draw_causal_arrow(910, 440, 944, 440, "E1", "activate", "")

    svg_p = out_p / "s01_biorender_hero.svg"
    png_p = out_p / "s01_biorender_hero.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def build_s03_biorender_exp_design(
    out_dir: str | Path,
    mode: str = "EV",
    title: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (b) `s03_biorender_exp_design`: Experimental Design & Three Orthogonal Causal Locks Schematic.
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    t = title or (
        "Experimental Design & Three Orthogonal Causal Locks (Reciprocal Rescue, Epistasis & Control-Free DI)"
        if mode_u == "EV"
        else "Experimental Perturbation Design, Quality Sanitization & Multi-Contrast Statistical Architecture"
    )
    c = BioRenderCanvas(width=1400, height=820, title=t)
    c.add_header_and_legend("Eliminating Single-Perturbation Confounding via Bidirectional Rescue & Control-Free Divergence")

    c.draw_card_panel(40, 110, 410, 660, "Lock 1: Perturbation Groups & QC Gate", "Rigorous Sample & Target Verification", fill="#F8FAFC", stroke="#64748B", header_fill="#334155")
    c.draw_card_panel(495, 110, 410, 660, "Lock 2: Reciprocal / Multi-Clone Concordance", "Orthogonal Perturbation Filter", fill="#EFF6FF", stroke="#2563EB", header_fill="#1D4ED8")
    c.draw_card_panel(950, 110, 410, 660, "Lock 3: 9-Sector Epistasis & Control-Free DI", "Algebraic Cancellation of Shared Control Noise", fill="#ECFDF5", stroke="#059669", header_fill="#047857")

    if mode_u == "EV":
        c.draw_protein_node(245, 220, "Target_A KD vs Ctrl_1", "n=3 vs 3 | Target_A -> 0/3 <LOD", "down", 6, "E1", width=280)
        c.draw_protein_node(245, 310, "Target_A OE vs Ctrl_2", "n=2 vs 3 (Outlier Tube Isolated)", "up", 8, "E1", width=280)
        c.draw_protein_node(245, 400, "Target_B KD vs Ctrl_1", "n=3 vs 3 | Target_B -1.42 log2FC", "down", 5, "E1", width=280)
        c.draw_protein_node(245, 510, "MISEV2023 Purity Gate", "Core EV Markers Verified", "reciprocal", 9, "E1", width=280)
        c.draw_protein_node(245, 610, "Dual-Track Quantification", "Track A DEqMS + Track B <LOD", "neutral", 7, "E1", width=280)

        c.draw_protein_node(700, 230, "Reciprocal KD x OE Lock", "Opposite Sign in KD vs OE", "reciprocal", 8, "E1", width=280)
        c.draw_protein_node(700, 340, "Gold ON-OFF-ON Switches", "3/3 -> 0/3 <LOD -> Rescued", "up", 5, "E1", width=280)
        c.draw_protein_node(700, 450, "Reciprocal Core Cargos", "High-Confidence Rescue Axis", "reciprocal", 6, "E1", width=280)
        c.draw_protein_node(700, 560, "Paralog Scan (Engine 03)", "Family Buffering Audited", "down", 4, "E1", width=280)

        c.draw_protein_node(1155, 230, "Shared Core (G1/G2)", "Co-Up + Co-Down Effectors", "up", 8, "E1", width=280)
        c.draw_protein_node(1155, 340, "Axis-Specific (G5/G6)", "Target_A-Spec | Target_B-Spec", "neutral", 7, "E1", width=280)
        c.draw_protein_node(1155, 450, "Seesaw Flux (G7)", "Opposing Target_A vs Target_B", "warning", 6, "E1", width=280)
        c.draw_protein_node(1155, 570, "Control-Free DI = B - A", "Shared Control Cancels Exactly", "reciprocal", 9, "E1", width=280)
    else:
        c.draw_protein_node(245, 230, "Primary Perturbation Arm", "Biological Replicates n=3 per Group", "down", 8, "E1", width=280)
        c.draw_protein_node(245, 340, "Orthogonal Clone / mAb Arm", "Independent Genetic/Biologic Arm", "up", 8, "E1", width=280)
        c.draw_protein_node(245, 450, "12-Point Forensic QC Gate", "Zero Batch Drift & Carryover Audit", "reciprocal", 9, "E1", width=280)
        c.draw_protein_node(245, 570, "Stage 1.5 Empirical Null", "Mock-Permuted Tail Calibration", "neutral", 8, "E1", width=280)

        c.draw_protein_node(700, 250, "Dual-Track Differential", "Track A DEqMS + Track B <LOD", "up", 8, "E1", width=280)
        c.draw_protein_node(700, 380, "Multi-Arm Concordance", "Eliminates Clone-Specific Off-Targets", "reciprocal", 8, "E1", width=280)
        c.draw_protein_node(700, 510, "Paralog Compensation Scan", "Detects Isoform Buffering (Engine 03)", "warning", 6, "E1", width=280)

        c.draw_protein_node(1155, 250, "9-Sector Epistasis Grid", "Shared vs Arm-Specific vs Seesaw", "up", 8, "E1", width=280)
        c.draw_protein_node(1155, 380, "TOST Equivalence Testing", "Proves True Biological Invariance", "neutral", 8, "E1", width=280)
        c.draw_protein_node(1155, 510, "Control-Free Divergence", "Direct Head-to-Head Effect Size", "reciprocal", 8, "E1", width=280)

    c.draw_causal_arrow(450, 420, 495, 420, "E1", "activate", "")
    c.draw_causal_arrow(905, 420, 950, 420, "E1", "activate", "")

    svg_p = out_p / "s03_biorender_exp_design.svg"
    png_p = out_p / "s03_biorender_exp_design.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def build_s14_biorender_topology(
    out_dir: str | Path,
    mode: str = "EV",
    title: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (c) `s14_biorender_topology`: Biophysical Membrane Topology & Subpopulation Sorting Bifurcation (`s14`/`s15` style).
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    t = title or (
        "Engine [07] Biophysical Topology Mechanism: Classical Secretory/1-TM vs 4-TM Tetraspanin/Lipid-Raft Sorting"
        if mode_u == "EV"
        else "Biophysical Compartment & Organelle Partitioning Across Perturbation Sectors"
    )
    c = BioRenderCanvas(width=1400, height=820, title=t)
    c.add_header_and_legend("UniProt Signal Peptide, Single-Pass 1-TM, 4-TM Tetraspanin, GPI-Anchor & Lipidated Raft Fisher Audit")

    c.draw_bilayer_membrane(y_pos=310, label_top="EXTRACELLULAR / VESICLE SURFACE / SECRETOME", label_bottom="INTRACELLULAR CYTOSOL / ENDOSOMAL SORTING")

    # Left half: Route A (Classical Signal Peptide + 1-TM Integrin/ECM)
    c.draw_card_panel(45, 380, 620, 395, "Route A: Classical Secretory / 1-TM Integrin-ECM Axis", "Signal Peptide (OR=4.84x, p=8.1e-6) & Secreted/ECM Cargo", fill="#FEF2F2", stroke="#DC2626", header_fill="#B91C1C")
    c.draw_integrin_dimer(180, 310, "ITGB1 / ITGA5 / ITGAV (1-TM)")
    c.draw_gpi_protein(370, 310, "NT5E / CD73 (GPI-Anchor)")
    c.draw_ecm_fibril(120, 195, 560, 195, "FN1 / THBS1 / LGALS3BP / MFGE8 Corona")
    c.draw_golgi_stack(220, 520, "ER -> Golgi Classical Secretory Export")
    c.draw_protein_node(490, 495, "RAB27A GTPase", "Exocytic Vesicle Docking", "up", 6, "E1", width=210)
    c.draw_protein_node(490, 585, "ADAM10 / ADAM17", "1-TM Sheddase Metalloprotease", "up", 8, "E1", width=210)
    c.draw_protein_node(355, 695, "Topology Enrichment: Classical SP/1-TM Enriched in Secretory Route A", "p = 1.94e-4 (Head-to-Head OR = 4.65x)", "up", 8, "E1", width=540)

    # Right half: Route B (Non-Classical Multi-TM / 4-TM Tetraspanin / Lipidated Raft)
    c.draw_card_panel(735, 380, 620, 395, "Route B: 4-TM Tetraspanin Cholesterol-Raft MVB Axis", "4-TM Tetraspanin (OR=5.12x) & Lipidated Raft Cargo", fill="#ECFDF5", stroke="#059669", header_fill="#047857")
    c.draw_4tm_protein(880, 310, "CD63 / CD9 (4-TM Tetraspanin)")
    c.draw_4tm_protein(1140, 310, "CD81 / CD151 (4-TM Raft)")
    c.draw_mvb_and_ev(920, 555, radius=78, label="Cholesterol-Rich MVB ILVs")
    c.draw_protein_node(1185, 495, "TSG101 / PDCD6IP", "ESCRT-I/ALIX Budding", "reciprocal", 9, "E1", width=210)
    c.draw_protein_node(1185, 585, "VAMP7 / SNAP23", "v/t-SNARE Membrane Fusion", "down", 5, "E1", width=210)
    c.draw_protein_node(1045, 695, "Subpopulation Dichotomy: ESCRT Shifted vs Tetraspanin Web Remodeled", "Proves Selective EV Subpopulation Rewiring Over Bulk Biogenesis Arrest", "reciprocal", 8, "E1", width=540)

    c.draw_causal_arrow(220, 470, 220, 345, "E1", "activate", " Exocytic Export")
    c.draw_causal_arrow(920, 465, 920, 345, "E1", "activate", " MVB-PM Fusion")

    svg_p = out_p / "s14_biorender_topology.svg"
    png_p = out_p / "s14_biorender_topology.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def build_s17_biorender_lock_and_key(
    out_dir: str | Path,
    mode: str = "EV",
    title: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (d) `s17_biorender_lock_and_key`: STRING Physical PPI "Lock-and-Key" Surface Corona & Key Driver Hubs (`s17` style).
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    t = title or (
        "Engine [06] STRING Physical PPI 'Lock-and-Key' Mechanism: Transmembrane Receptors Anchor Secreted ECM Ligands"
        if mode_u == "EV"
        else "Engine [06] STRING Physical PPI Key Driver Hubs & Cross-Sector Macromolecular Bridges"
    )
    c = BioRenderCanvas(width=1400, height=820, title=t)
    c.add_header_and_legend("Degree-Corrected Random Walk with Restart (RWR) & 1-Hop Hypergeometric Key Driver Analysis (KDA)")

    c.draw_card_panel(40, 110, 1320, 670, "Cross-Sector Physical PPI Bridge Architecture (High-Confidence Inter-Group Edges)", "Explaining How Vesicular/Membrane Perturbations Physically Capture Soluble Ligands", fill="#FFFFFF", stroke="#CBD5E1", header_fill="#0F172A")

    # Left: Intracellular / Vesicle Membrane "Lock" Anchors
    c.draw_ev_vesicle(260, 430, radius=110, label="EV / Plasma Membrane", sublabel="Transmembrane 'LOCK' Receptors")
    c.draw_protein_node(260, 250, "ITGB1 (KDA Hub #1)", "deg=31 | RWR Z=4.8 | Co-Up", "up", 14, "E1", width=230)
    c.draw_protein_node(500, 310, "ITGAV / ITGA5 / ITGA3", "Integrin Heterodimer Locks", "up", 11, "E1", width=220)
    c.draw_protein_node(500, 420, "CD151 / CD81 / CD63", "Tetraspanin Web Scaffolds", "reciprocal", 8, "E1", width=220)
    c.draw_protein_node(500, 530, "ADAM10 / IGSF8", "Surface Sheddase & Partner", "up", 9, "E1", width=220)
    c.draw_protein_node(260, 620, "RAB27A / VAMP7 / SNAP23", "SNARE-Rab Trafficking Core", "down", 6, "E1", width=240)

    # Right: Extracellular / Paracrine "Key" Ligands & Recipient Cells
    c.draw_protein_node(860, 250, "FN1 (Fibronectin)", "deg=28 | ECM Bridge Hub", "up", 22, "E1", width=220)
    c.draw_protein_node(860, 355, "THBS1 (Thrombospondin-1)", "deg=19 | Secretome SGI +1.66", "up", 15, "E1", width=220)
    c.draw_protein_node(860, 460, "LGALS3BP (Gal-3BP)", "deg=14 | EV Adhesive Corona", "up", 12, "E1", width=220)
    c.draw_protein_node(860, 565, "MFGE8 / SPP1 / VCAN", "Heparan Sulfate & Integrin Ligands", "up", 10, "E1", width=230)

    c.draw_organelle(1190, 410, rx=115, ry=165, kind="nucleus", label="Recipient Stromal / Endothelial Niche")
    c.draw_protein_node(1190, 350, "Microenvironment Niche", "ECM Remodeling & Adhesion", "up", 8, "E1", width=200)
    c.draw_protein_node(1190, 455, "Angiogenesis & Signaling", "Non-Autonomous Paracrine Signal", "warning", 6, "E2", width=200)

    # Lock-and-Key Physical PPI Arrows
    c.draw_causal_arrow(375, 250, 745, 250, "E1", "activate", "STRING Physical Score 0.999")
    c.draw_causal_arrow(610, 310, 745, 355, "E1", "activate", "RGD / Integrin Binding 0.995")
    c.draw_causal_arrow(610, 420, 745, 460, "E1", "activate", "Corona Tethering 0.965")
    c.draw_causal_arrow(610, 530, 745, 565, "E1", "activate", "Surface Complex 0.950")
    c.draw_causal_arrow(975, 355, 1075, 375, "E2", "activate", "Paracrine Uptake")
    c.draw_causal_arrow(975, 460, 1075, 440, "E2", "activate", "Intercellular Signaling")

    svg_p = out_p / "s17_biorender_lock_and_key.svg"
    png_p = out_p / "s17_biorender_lock_and_key.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def build_s24_biorender_spotlight(
    out_dir: str | Path,
    top6_df: Optional[pd.DataFrame] = None,
    mode: str = "EV",
    title: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (e) `s24_biorender_spotlight`: 6 Flagship Target Subcellular Localization & Multi-Evidence Cards (`s24` style).
    """
    out_p = Path(out_dir)
    t = title or "Top 6 Flagship Spotlight Drivers: Subcellular Topology & Clinical Evidence"
    c = BioRenderCanvas(width=1400, height=840, title=t)
    c.add_header_and_legend("Selected via 3-Lens × 5-Layer Multi-Evidence Priority Ranking + 1,000x Weight Robustness")

    default_cards = [
        ("ITGB1", "Perspective I (#1 Shared Core)", "1-TM Integrin Receptor Hub", "[TEMPLATE] Shared Core | CPTAC Cohort | FDA Drug Target", "up", "#B91C1C", 6),
        ("CD63", "Perspective I (#2 Shared Core)", "4-TM Tetraspanin Vesicle Gatekeeper", "[TEMPLATE] Direct Marker | Biofluid Enriched | Paracrine Core", "down", "#0369A1", 4),
        ("THBS1", "Perspective II (#1 Seesaw/Axis-B)", "Secreted Matricellular Glycoprotein", "[TEMPLATE] Seesaw Divergence | CPTAC Cohort | Secreted", "up", "#047857", 8),
        ("CD81", "Perspective II (#2 Seesaw/Axis-B)", "4-TM Tetraspanin Microdomain Scaffold", "[TEMPLATE] Subpopulation Shift | Biofluid Enriched | Paracrine", "reciprocal", "#047857", 5),
        ("NT5E", "Perspective III (#1 Rescue/Axis-A)", "GPI-Anchored Ectonucleotidase (CD73)", "[TEMPLATE] Reciprocal Rescue | CPTAC Cohort | Clinical mAb Target", "up", "#6D28D9", 6),
        ("LGALS3BP", "Perspective III (#2 Rescue/Axis-A)", "Secreted Scavenger Receptor Ligand", "[TEMPLATE] Reciprocal Switch | CPTAC Cohort | Liquid Biopsy Marker", "up", "#6D28D9", 7),
    ]

    cards_data = []
    if top6_df is not None and not top6_df.empty:
        g_col = "std_gene" if "std_gene" in top6_df.columns else top6_df.columns[0]
        header_pal = ["#B91C1C", "#B91C1C", "#047857", "#047857", "#6D28D9", "#6D28D9"]
        for idx, (_, r) in enumerate(top6_df.head(6).iterrows()):
            g = str(r[g_col])
            persp = str(r.get("Perspective", f"Perspective {(idx//2)+1}")).replace("_", " ")
            score = float(pd.to_numeric(r.get("Composite_Evidence_Score", 85.0), errors="coerce") or 85.0)
            cp_raw = pd.to_numeric(r.get("cptac_T_vs_N_log2FC", float("nan")), errors="coerce")
            fl_raw = pd.to_numeric(r.get("clin_fluid_log2FC", float("nan")), errors="coerce")
            dm_col = "depmap_context_dep_lines" if "depmap_context_dep_lines" in r.index else (
                "depmap_dep_lines" if "depmap_dep_lines" in r.index else "depmap_HCC_dep_lines"
            )
            dm_raw = pd.to_numeric(r.get(dm_col, float("nan")), errors="coerce")
            pep_cnt = int(pd.to_numeric(r.get("std_unique_peptides", 2), errors="coerce") or 2)
            topo = "4-TM Tetraspanin/MARVEL" if bool(r.get("is_4tm_marvel_tspan", False)) else (
                "1-TM Transmembrane Receptor" if bool(r.get("is_single_tm", False)) else (
                    "GPI-Anchored Surface Protein" if bool(r.get("is_gpi_anchored", False)) else (
                        "Signal-Peptide Secreted/ECM Cargo" if bool(r.get("has_signal_peptide", False)) else "Intracellular / Vesicular Hub"
                    )
                )
            )
            cp_str = f"{float(cp_raw):+.2f}" if pd.notna(cp_raw) else "N/A"
            fl_str = f"{float(fl_raw):+.2f}" if pd.notna(fl_raw) else "N/A"
            dm_str = f"{int(dm_raw)}" if pd.notna(dm_raw) else "N/A"
            ev_line = f"Priority={score:.1f} | CPTAC={cp_str} | Fluid={fl_str} | DepMap={dm_str}"
            dir_v = "up" if (pd.notna(cp_raw) and float(cp_raw) >= 0) else "down"
            cards_data.append((g, persp[:32], topo, ev_line, dir_v, header_pal[idx % 6], pep_cnt))
    if not cards_data:
        cards_data = default_cards

    coords = [
        (40, 110), (495, 110), (950, 110),
        (40, 470), (495, 470), (950, 470),
    ]
    for idx, ((x0, y0), (gene, persp, topo, ev_line, dir_c, h_col, pep_cnt)) in enumerate(zip(coords, cards_data[:6])):
        c.draw_card_panel(x0, y0, 410, 330, f"Spotlight #{idx+1}: {gene}", persp, fill="#FFFFFF", stroke="#CBD5E1", header_fill=h_col)
        # Mini membrane + topology icon inside card
        c.draw_bilayer_membrane(y_pos=y0 + 135, label_top="", label_bottom="", x_start=x0 + 25, x_end=x0 + 385)
        if "4-TM" in topo:
            c.draw_4tm_protein(x0 + 115, y0 + 135, "", color=h_col)
        elif "1-TM" in topo:
            c.draw_integrin_dimer(x0 + 115, y0 + 135, "")
        elif "GPI" in topo:
            c.draw_gpi_protein(x0 + 115, y0 + 135, "", color=h_col)
        else:
            c.draw_ev_vesicle(x0 + 115, y0 + 135, radius=32, label=gene, sublabel="")

        ev_grade = "E1" if pep_cnt >= 2 else "E3"
        status_str = (
            f"Status: Verified Multi-Peptide Core Driver ({pep_cnt} Pep)"
            if pep_cnt >= 2
            else "Status: Single-Peptide Exploratory Detection (1 Pep*)"
        )
        c.draw_protein_node(x0 + 280, y0 + 135, gene, topo[:24], dir_c, pep_cnt, ev_grade, width=185, height=46)
        c.elements.append(
            f'<rect x="{x0+20}" y="{y0+215}" width="370" height="92" rx="10" fill="#F8FAFC" stroke="#E2E8F0"/>'
            f'<text x="{x0+34}" y="{y0+242}" font-family="Arial" font-size="11.5" font-weight="bold" fill="#0F172A">Structural Class: {html.escape(topo)}</text>'
            f'<text x="{x0+34}" y="{y0+266}" font-family="Arial" font-size="10.8" fill="#334155">{html.escape(ev_line)}</text>'
            f'<text x="{x0+34}" y="{y0+290}" font-family="Arial" font-size="10.5" font-weight="bold" fill="{h_col}">{html.escape(status_str)}</text>'
        )

    svg_p = out_p / "s24_biorender_spotlight.svg"
    png_p = out_p / "s24_biorender_spotlight.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def build_s25_biorender_working_model(
    out_dir: str | Path,
    mode: str = "EV",
    title: Optional[str] = None,
) -> Tuple[Path, Path]:
    """
    (f) `s25_biorender_working_model`: Grand Integrated Mechanistic Working Model (`s25` style).
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    t = title or (
        "Integrated Working Model: Dual-Route Vesicle Biogenesis & Paracrine Microenvironment Remodeling"
        if mode_u == "EV"
        else "Integrated Inside-Out Working Model: Intracellular Perturbation to Secreted Vesicles & Clinical Cohorts"
    )
    c = BioRenderCanvas(width=1400, height=840, title=t)
    c.add_header_and_legend("Synthesizing 12 Analytical Engines: Dual-Track DEqMS/<LOD, 9-Sector Epistasis, Topology, STRING KDA & Clinical Cohorts")

    c.draw_bilayer_membrane(y_pos=320, label_top="EXTRACELLULAR BIOFLUID & VASCULAR NICHE", label_bottom="DONOR CELL CYTOSOL & ENDOMEMBRANE TRAFFICKING")

    # Donor Cell Left: Intracellular Nucleus / Redox / MVB / Secretory Golgi
    c.draw_organelle(190, 620, rx=130, ry=90, kind="nucleus", label="Donor Nucleus (STAT1 / NFE2L2 / MYC)")
    c.draw_golgi_stack(480, 600, "ER-Golgi Classical Secretory Route")
    c.draw_mvb_and_ev(820, 590, radius=92, label="Late Endosome / MVB (ESCRT & 4-TM Raft)")

    c.draw_protein_node(190, 455, "NFE2L2 / GCLM / FTH1", "Intracellular Redox & Iron Buffer", "reciprocal", 8, "E1", width=230)
    c.draw_protein_node(480, 455, "RAB27A GTPase Switch", "Controls Classical SP/1-TM Export", "up", 7, "E1", width=230)
    c.draw_protein_node(820, 435, "CD63 / CD81 (4-TM)", "Governs Cholesterol-Raft MVB Flux", "down", 6, "E1", width=235)
    c.draw_protein_node(1160, 560, "TSG101 / PDCD6IP / VPS4A", "ESCRT Machinery (DepMap Essential)", "down", 9, "E1", width=245)

    # Membrane & Extracellular Space: Secreted EVs, Lock-and-Key Corona, Biofluid & Recipient Niche
    c.draw_integrin_dimer(440, 320, "ITGB1 / ITGAV")
    c.draw_gpi_protein(620, 320, "NT5E / CD73 (GPI)")
    c.draw_4tm_protein(820, 320, "CD81 / CD151 / CD63")

    c.draw_ev_vesicle(520, 175, radius=54, label="Classical EV / Corona", sublabel="FN1 / THBS1 / LGALS3BP")
    c.draw_ev_vesicle(860, 175, radius=54, label="4-TM Raft Exosome", sublabel="CD81 / CD63 / CD9", fill_grad="url(#ev-emerald)", stroke="#059669")

    c.draw_card_panel(1050, 115, 310, 175, "Human Clinical & Translational Impact", "CPTAC Cohort + Paired Biofluid + DepMap", fill="#F0FDF4", stroke="#059669", header_fill="#047857")
    c.draw_protein_node(1205, 195, "Biofluid Secretion (SGI>0.5)", "FN1 / THBS1 / ITGB1 / NT5E", "up", 9, "E1", width=270, height=38)
    c.draw_protein_node(1205, 248, "DepMap Paracrine Proof", "Non-Autonomous Niche Remodeling", "reciprocal", 7, "E1", width=270, height=38)

    c.draw_causal_arrow(190, 525, 190, 480, "E1", "activate", "")
    c.draw_causal_arrow(480, 540, 480, 480, "E1", "activate", "")
    c.draw_causal_arrow(480, 430, 480, 345, "E1", "activate", "Exocytic Docking")
    c.draw_causal_arrow(820, 495, 820, 460, "E1", "activate", "")
    c.draw_causal_arrow(820, 410, 820, 345, "E1", "activate", "MVB Fusion")
    c.draw_causal_arrow(575, 175, 1045, 190, "E1", "activate", "Biofluid Accumulation")
    c.draw_causal_arrow(915, 175, 1045, 225, "E1", "activate", "Paracrine Delivery")

    svg_p = out_p / "s25_biorender_working_model.svg"
    png_p = out_p / "s25_biorender_working_model.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


# ==============================================================================
# MODULE 7: Penultimate Slide — Core Changed Genes, Pathways, Literature & Biology (`s26`)
# ==============================================================================
def build_s26_single_kd_summary(
    out_dir: Path,
    mode: str = "EV",
    left_header: Optional[Tuple[str, str]] = None,
    right_header: Optional[Tuple[str, str]] = None,
    left_cards: Optional[Sequence[Tuple[str, str, str, str, str, str, str, str]]] = None,
    right_cards: Optional[Sequence[Tuple[str, str, str, str, str, str, str, str]]] = None,
    bottom_banner: Optional[Tuple[str, str]] = None,
) -> Tuple[Path, Path]:
    """
    Generate the penultimate slide schematic (`s26_biorender_single_kd_summary.svg`/`.png`)
    designed specifically for wet-lab biologists and PIs:
    Side-by-side summary of each perturbation's Top Changed Hero Genes (with fold-changes),
    Pathways, Landmark Literature (Author, Journal, Year, PMID), and Plain-Language Biological Meaning.
    """
    out_p = Path(out_dir)
    mode_u = str(mode).upper()
    c = BioRenderCanvas(width=1320, height=650, title="Single-Perturbation Hero Genes, Pathways, Literature & Biological Meaning")

    if mode_u in ("WCL", "METABOLIC_WCL"):
        lh = left_header or (
            "Perturbation Arm 1 · Redox Axis (KD × OE Reciprocal Rescue)",
            "Primary Axis: Glutathione Detoxification, Secretory Remodeling & Ferritin Induction",
        )
        rh = right_header or (
            "Perturbation Arm 2 · Iron Homeostasis Axis (Multi-Clone KD)",
            "Primary Axis: Iron Export Blockade, Lipoprotein Remodeling & Vesicle Adaptor Priming",
        )
        l_cards = left_cards or [
            (
                "1. [TEMPLATE] GCLM / NQO1 (Reciprocal KD × OE Axis)",
                "Reciprocal Lock", "#0E6251", "#CCFBF1",
                "Pathway: Glutathione Synthesis & Electrophile Detoxification",
                "Literature: Stockwell et al., Cell 2017 (PMID: 28985560) · Dixon et al., Cell 2012",
                "→ Target depletion impairs glutathione synthesis and sensitizes cells to GPX4 blockade.",
                "golgi",
            ),
            (
                "2. [TEMPLATE] FTH1 / FTL Ferritin Holocomplex",
                "Ferritin Shield", "#1D6FA5", "#DBEAFE",
                "Pathway: NFE2L2 / NRF2 Antioxidant Response & Labile Iron Sequestration",
                "Literature: Mancias et al., Nature 2014 (PMID: 24695223) · Stockwell et al., Cell 2017",
                "→ Redox stress activates NRF2 to assemble FTH1/FTL ferritin cages buffering free Fe2+.",
                "ecm",
            ),
            (
                "3. [TEMPLATE] S100A8 / S100A9 · IGFBP5 Paracrine Axis",
                "DAMP & Growth", "#6D28D9", "#EDE9FE",
                "Pathway: Calprotectin DAMP Inflammatory Signaling & Growth Factor Bioavailability",
                "Literature: Baxter, Nat Rev Cancer 2014 (PMID: 24722429)",
                "→ Co-depletion of calprotectin and growth factor binding proteins dampens paracrine signaling.",
                "rtk",
            ),
        ]
        r_cards = right_cards or [
            (
                "1. [TEMPLATE] APOB / APOE Lipoprotein Hub",
                "Lipid Remodeling", "#C0392B", "#FEE2E2",
                "Pathway: ER Lipoprotein Assembly & Neutral Lipid Transport",
                "Literature: Hussain et al., J Lipid Res 2012 (PMID: 22236406)",
                "→ Intracellular iron accumulation suppresses lipoprotein export and triggers lipid peroxidation.",
                "integrin",
            ),
            (
                "2. [TEMPLATE] HMOX1 / TXNRD1 Antioxidant Axis",
                "Redox Buffering", "#B45309", "#FEF3C7",
                "Pathway: Heme Oxygenase & Thioredoxin Reductase Stress Response",
                "Literature: Dixon et al., Cell 2012 (PMID: 22632970)",
                "→ Compensatory antioxidant enzymes buffer acute oxidative stress upon iron retention.",
                "tetraspanin",
            ),
            (
                "3. [TEMPLATE] SDCBP / RRM2 Vesicle Bridge",
                "EV Bridge Hub", "#6D28D9", "#EDE9FE",
                "Pathway: Syndecan-Syntenin Exosomal Budding Adaptor & Deoxynucleotide Synthesis",
                "Literature: Baietti et al., Nat Cell Biol 2012 (PMID: 22660413)",
                "→ Metabolic stress induces Syntenin adaptors to prime extracellular vesicle export.",
                "gpi",
            ),
        ]
        b_banner = bottom_banner or (
            "Mechanistic Division of Labor Across Perturbation Arms",
            "Arm 1 disables glutathione redox detoxification & induces FTH1/FTL ferritin cages | Arm 2 traps intracellular iron & primes Syntenin-mediated vesicle export",
        )
    elif mode_u in ("BIOLOGICS", "MAB_WCL"):
        lh = left_header or (
            "Core Axis 1 · Conserved Parental ∩ Humanized Antibody Targets",
            "Primary Axis: Cell-Surface IgG Capture & Conserved Transcriptional Repression",
        )
        rh = right_header or (
            "Core Axis 2 · Humanized Lead Effectors & Survival Buffering",
            "Primary Axis: Purging Xenogeneic Drift Artifacts & Rational Combination Synergy",
        )
        l_cards = left_cards or [
            (
                "1. [TEMPLATE] IGHG4 / IGKC Cell-Surface Antibody Capture",
                "Target Engagement", "#0E6251", "#CCFBF1",
                "Pathway: Direct Cell-Surface Monoclonal Antibody Binding & Receptor Internalization",
                "Literature: Carter & Lazar, Nat Rev Drug Discov 2018 (PMID: 29192287)",
                "→ Mass spectrometry directly verifies cell-bound humanized IgG4/κ retention.",
                "golgi",
            ),
            (
                "2. [TEMPLATE] STAT3 / MYC (TOST Equivalent Core)",
                "#1 TOST Conserved", "#1D6FA5", "#DBEAFE",
                "Pathway: Conserved Oncogenic Transcription Factor & Chromatin Remodeling Axis",
                "Literature: Bromberg et al., Cell 1999 (PMID: 10428038)",
                "→ Humanized antibody candidates suppress oncogenic transcription equivalently to parental mAb.",
                "ecm",
            ),
            (
                "3. [TEMPLATE] CHD4 / HDAC1 NuRD Complex",
                "Chromatin Lock", "#6D28D9", "#EDE9FE",
                "Pathway: NuRD Deacetylase Chromatin Remodeling Complex",
                "Literature: Lai & Wade, Nat Rev Cancer 2011 (PMID: 21779009)",
                "→ Coordinated chromatin remodeling locks treated cells into a cytostatic state.",
                "rtk",
            ),
        ]
        r_cards = right_cards or [
            (
                "1. [TEMPLATE] ACTN4 / CTDSP1 Cytoskeletal Core",
                "Humanized Core", "#C0392B", "#FEE2E2",
                "Pathway: Cortical Actin-Membrane Stabilization & RNA Pol II CTD Dephosphorylation",
                "Literature: Yeo et al., J Biol Chem 2003 (PMID: 12529368)",
                "→ Humanized leads stabilize cortical cytoskeleton and dephosphorylate RNA Pol II CTD.",
                "integrin",
            ),
            (
                "2. [TEMPLATE] MCL1 / BCL2L1 Survival Buffering",
                "Combo Vulnerability", "#B45309", "#FEF3C7",
                "Pathway: BCL-2 Family Anti-Apoptotic Mitochondrial Survival Buffering",
                "Literature: Kotschy et al., Nature 2016 (PMID: 27193678)",
                "→ Antibody engagement induces compensatory MCL1; combining with MCL1i triggers apoptosis.",
                "tetraspanin",
            ),
            (
                "3. [TEMPLATE] Xenogeneic Drift Purged via TOST Equivalence",
                "Specificity Gain", "#6D28D9", "#EDE9FE",
                "Pathway: Elimination of Non-Specific Murine Fc Stress Signatures",
                "Literature: Hwang & Foote, Methods 2005 (PMID: 15882989)",
                "→ CDR humanization strips off-pathway Fc stress while preserving core efficacy.",
                "gpi",
            ),
        ]
        b_banner = bottom_banner or (
            "Biologics Humanization & Equivalence Takeaway",
            "Humanization eliminates xenogeneic off-target drift while preserving on-target transcriptional repression & revealing a rational MCL1i combination",
        )
    elif mode_u == "UNIFIED":
        lh = left_header or (
            "Intracellular WCL & Biologics Core (Layers 1 & 2)",
            "Primary Axis: Intracellular Redox/Iron Defense & Receptor-to-Nucleus Reprogramming",
        )
        rh = right_header or (
            "Secreted EV & Clinical Biofluid Core (Layer 3)",
            "Primary Axis: MVB Lipid-Raft Receptor Locks & Co-Secreted Extracellular Matrix Keys",
        )
        l_cards = left_cards or [
            (
                "1. [TEMPLATE] GCLM / NQO1 → FTH1/FTL Ferritin Cage",
                "Layer 1 Redox", "#0E6251", "#CCFBF1",
                "Pathway: Glutathione Detoxification & NRF2 Ferritin Iron Sequestration",
                "Literature: Stockwell et al., Cell 2017 (PMID: 28985560)",
                "→ Redox stress activates NRF2 ferritin cages to buffer labile intracellular iron.",
                "golgi",
            ),
            (
                "2. [TEMPLATE] Iron Retention → APOB & SDCBP Adaptor",
                "Layer 1->3 Bridge", "#1D6FA5", "#DBEAFE",
                "Pathway: Lipoprotein Suppression & Syndecan-Syntenin Exosomal Budding",
                "Literature: Baietti et al., Nat Cell Biol 2012 (PMID: 22660413)",
                "→ Intracellular metabolic stress induces Syntenin adaptors priming EV release.",
                "ecm",
            ),
            (
                "3. [TEMPLATE] IgG4/κ → Conserved Chromatin Repression & MCL1",
                "Layer 2 Biologics", "#6D28D9", "#EDE9FE",
                "Pathway: TOST-Equivalent Transcriptional Arrest & Actionable MCL1 Buffering",
                "Literature: Kotschy et al., Nature 2016 (PMID: 27193678)",
                "→ Humanized antibodies preserve core target repression while purging off-target drift.",
                "rtk",
            ),
        ]
        r_cards = right_cards or [
            (
                "1. [TEMPLATE] CD63 · ITGA5/ITGB1 · ESCRT Machinery",
                "Layer 3 Raft Locks", "#C0392B", "#FEE2E2",
                "Pathway: MVB Cholesterol Lipid-Raft Organization & Organotropic Integrin Sorting",
                "Literature: Hoshino et al., Nature 2015 (PMID: 26524530) · van Niel et al., Nat Rev MCB 2018",
                "→ Tetraspanin and integrin complexes assemble on EVs to direct recipient cell uptake.",
                "integrin",
            ),
            (
                "2. [TEMPLATE] RAB27A → FLOT1/2 · FN1 · THBS1",
                "Layer 3 Matrix Keys", "#B45309", "#FEF3C7",
                "Pathway: Vesicle Docking, Planar Lipid-Raft Export & Matrix Remodeling",
                "Literature: Ostrowski et al., Nat Cell Biol 2010 (PMID: 19966785)",
                "→ Vesicle GTPases co-export soluble matrix ligands captured by surface integrin receptors.",
                "tetraspanin",
            ),
            (
                "3. [TEMPLATE] Clinical Biofluid Panel: ITGB1 · TFRC · NT5E",
                "Clinical Translation", "#6D28D9", "#EDE9FE",
                "Pathway: Patient Biofluid Hyper-Secretion & Niche Remodeling",
                "Literature: Hoshino et al., Cell 2020 (PMID: 32795414)",
                "→ Cross-compartment bridge proteins connect intracellular stress to non-invasive liquid biopsy.",
                "gpi",
            ),
        ]
        b_banner = bottom_banner or (
            "Inside-Out Multi-Compartment Proteomics Synthesis",
            "Layers 1–2 reveal intracellular metabolic & receptor vulnerabilities | Layer 3 reveals how blocking vesicle sorting disarms secreted microenvironment remodelers",
        )
    else:
        lh = left_header or (
            "Perturbation Arm 1 · Vesicle Docking GTPase (KD × OE Reciprocal)",
            "Primary Axis: Trans-Golgi/Endosomal Vesicle Docking, Flotillin Export & Flux Shunting",
        )
        rh = right_header or (
            "Perturbation Arm 2 · Tetraspanin/Raft Scaffold (Specific Depletion)",
            "Primary Axis: MVB Cholesterol Lipid-Raft Organization & Exosomal Surface Receptor Anchoring",
        )
        l_cards = left_cards or [
            (
                "1. [TEMPLATE] FLOT1/2 · ROBO1 & CLDN1 Raft Export",
                "GTPase-Promoted", "#0E6251", "#CCFBF1",
                "Pathway: Flotillin Planar Lipid-Raft Exosome Budding & Endothelial Guidance",
                "Literature: Otto & Nichols, J Cell Sci 2011 (PMID: 22010190) · Baietti et al., Nat Cell Biol 2012",
                "→ Rab GTPase drives Flotillin-1/2 raft packaging and acts as an ON-OFF-ON switch for ROBO1.",
                "golgi",
            ),
            (
                "2. [TEMPLATE] FN1 · SERPINE2 · THBS1 Matrix Core",
                "Shared Secretome", "#1D6FA5", "#DBEAFE",
                "Pathway: Extracellular Matrix Remodeling & Vascular Niche Pre-Conditioning",
                "Literature: Peinado et al., Nat Med 2012 (PMID: 22635005)",
                "→ Exosomal matrix remodelers alter stromal permeability and pre-metastatic niche formation.",
                "ecm",
            ),
            (
                "3. [TEMPLATE] ITGA5 · SLC3A2 · TFRC Seesaw Shunting",
                "Seesaw Shunting", "#6D28D9", "#EDE9FE",
                "Pathway: Antagonistic Seesaw Flux into Competing MVB Integrin & Transporter Subpopulations",
                "Literature: Kowal et al., PNAS 2016 (PMID: 26858453)",
                "→ Blocking Route A shunts membrane cargo flux into competing tetraspanin-enriched MVBs.",
                "rtk",
            ),
        ]
        r_cards = right_cards or [
            (
                "1. [TEMPLATE] ITGA5 · ITGB1 · Integrin Adhesion Web",
                "Scaffold-Specific", "#C0392B", "#FEE2E2",
                "Pathway: Exosomal Integrin Adhesion & Organotropic Target Cell Uptake",
                "Literature: Hoshino et al., Nature 2015 (PMID: 26524530, Exosome Organotropism)",
                "→ Exosomal integrins act as molecular zip-codes directing uptake into recipient niche cells.",
                "integrin",
            ),
            (
                "2. [TEMPLATE] CD63 · CD81 · CD151 · ESCRT Raft Core",
                "MVB Raft Core", "#B45309", "#FEF3C7",
                "Pathway: Tetraspanin-Enriched Microdomain (TEM) & MVB Cargo Sorting",
                "Literature: van Niel et al., Nat Rev MCB 2018 (PMID: 29339798)",
                "→ 4-TM scaffolds nucleate CD63/CD151 cholesterol rafts required to load ILV cargos.",
                "tetraspanin",
            ),
            (
                "3. [TEMPLATE] NT5E / CD73 · ADAM10 Surface Corona",
                "GPI & Corona", "#6D28D9", "#EDE9FE",
                "Pathway: GPI-Anchored Ectoenzyme Cleavage & Surface Corona Tethering",
                "Literature: Toth et al., J Extracell Vesicles 2021 (PMID: 34520123)",
                "→ Scaffold loss strips surface GPI-ectoenzymes and metalloprotease corona partners.",
                "gpi",
            ),
        ]
        b_banner = bottom_banner or (
            "Single-Perturbation Division of Labor in Vesicle Trafficking",
            "Arm 1 gates vesicle docking & Flotillin/matrix export | Arm 2 builds MVB tetraspanin/integrin lipid-raft anchors that capture the surface corona",
        )

    # Left Master Panel
    c.elements.append(
        f'<rect x="20" y="14" width="628" height="544" rx="18" fill="#F0FDFA" stroke="#0E6251" stroke-width="2.0" filter="url(#soft-shadow)"/>'
        f'<rect x="20" y="14" width="628" height="66" rx="18" fill="#0E6251"/>'
        f'<text x="42" y="42" font-family="Arial, Helvetica, sans-serif" font-size="15.2" font-weight="700" fill="#FFFFFF">{html.escape(lh[0])}</text>'
        f'<text x="42" y="64" font-family="Arial, Helvetica, sans-serif" font-size="11.4" font-weight="600" fill="#CCFBF1">{html.escape(lh[1])}</text>'
    )
    # Right Master Panel
    c.elements.append(
        f'<rect x="672" y="14" width="628" height="544" rx="18" fill="#FEF2F2" stroke="#C0392B" stroke-width="2.0" filter="url(#soft-shadow)"/>'
        f'<rect x="672" y="14" width="628" height="66" rx="18" fill="#C0392B"/>'
        f'<text x="694" y="42" font-family="Arial, Helvetica, sans-serif" font-size="15.2" font-weight="700" fill="#FFFFFF">{html.escape(rh[0])}</text>'
        f'<text x="694" y="64" font-family="Arial, Helvetica, sans-serif" font-size="11.4" font-weight="600" fill="#FEE2E2">{html.escape(rh[1])}</text>'
    )

    for idx, (genes, badge, b_col, b_bg, pathway, lit, takeaway, icon) in enumerate(l_cards):
        cy = 92 + idx * 152
        c.elements.append(
            f'<rect x="34" y="{cy}" width="600" height="140" rx="14" fill="#FFFFFF" stroke="#99F6E4" stroke-width="1.5"/>'
            f'<rect x="506" y="{cy+10}" width="116" height="24" rx="12" fill="{b_bg}" stroke="{b_col}" stroke-width="1.0"/>'
            f'<text x="564" y="{cy+26}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="10.0" font-weight="700" fill="{b_col}">{html.escape(badge)}</text>'
            f'<text x="134" y="{cy+28}" font-family="Arial, Helvetica, sans-serif" font-size="11.4" font-weight="700" fill="#0F172A">{html.escape(genes)}</text>'
            f'<text x="134" y="{cy+54}" font-family="Arial, Helvetica, sans-serif" font-size="10.8" font-weight="700" fill="{b_col}">{html.escape(pathway)}</text>'
            f'<text x="134" y="{cy+84}" font-family="Arial, Helvetica, sans-serif" font-size="10.4" font-weight="700" fill="#334155">{html.escape(lit)}</text>'
            f'<text x="134" y="{cy+114}" font-family="Arial, Helvetica, sans-serif" font-size="10.4" font-weight="600" fill="#0E6251">{html.escape(takeaway)}</text>'
        )
        if icon == "golgi":
            for i_g, w_cis in enumerate([74, 64, 54, 44]):
                yk = cy + 48 + i_g * 14
                c.elements.append(
                    f'<rect x="{82 - w_cis//2}" y="{yk}" width="{w_cis}" height="9" rx="4.5" fill="#FDE68A" stroke="#D97706" stroke-width="1.6"/>'
                )
        elif icon == "ecm":
            c.draw_ecm_fibril(46, cy + 56, 118, cy + 56, label="", color="#1D6FA5")
            c.draw_ecm_fibril(46, cy + 94, 118, cy + 94, label="", color="#0E6251")
        else:
            c.draw_protein_node(82, cy + 56, "Target", "", "down", 4, "E2", width=74, height=28)
            c.draw_protein_node(82, cy + 96, "Effector", "", "reciprocal", 4, "E2", width=74, height=28)

    for idx, (genes, badge, b_col, b_bg, pathway, lit, takeaway, icon) in enumerate(r_cards):
        cy = 92 + idx * 152
        c.elements.append(
            f'<rect x="686" y="{cy}" width="600" height="140" rx="14" fill="#FFFFFF" stroke="#FECACA" stroke-width="1.5"/>'
            f'<rect x="1158" y="{cy+10}" width="116" height="24" rx="12" fill="{b_bg}" stroke="{b_col}" stroke-width="1.0"/>'
            f'<text x="1216" y="{cy+26}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="10.0" font-weight="700" fill="{b_col}">{html.escape(badge)}</text>'
            f'<text x="786" y="{cy+28}" font-family="Arial, Helvetica, sans-serif" font-size="11.4" font-weight="700" fill="#0F172A">{html.escape(genes)}</text>'
            f'<text x="786" y="{cy+54}" font-family="Arial, Helvetica, sans-serif" font-size="10.8" font-weight="700" fill="{b_col}">{html.escape(pathway)}</text>'
            f'<text x="786" y="{cy+84}" font-family="Arial, Helvetica, sans-serif" font-size="10.4" font-weight="700" fill="#334155">{html.escape(lit)}</text>'
            f'<text x="786" y="{cy+114}" font-family="Arial, Helvetica, sans-serif" font-size="10.4" font-weight="600" fill="#991B1B">{html.escape(takeaway)}</text>'
        )
        c.elements.append(f'<line x1="700" y1="{cy+88}" x2="770" y2="{cy+88}" stroke="#F59E0B" stroke-width="4.5" stroke-linecap="round"/>')
        if icon == "integrin":
            c.draw_integrin_dimer(735, cy + 88, label="")
        elif icon == "tetraspanin":
            c.draw_4tm_protein(735, cy + 88, label="")
        else:
            c.draw_gpi_protein(735, cy + 88, label="")

    # Bottom Synthesis Banner
    c.elements.append(
        f'<rect x="20" y="570" width="1280" height="64" rx="14" fill="#FFFBEB" stroke="#D97706" stroke-width="1.8" filter="url(#soft-shadow)"/>'
        f'<text x="660" y="596" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="13.0" font-weight="700" fill="#92400E">{html.escape(b_banner[0])}</text>'
        f'<text x="660" y="619" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="11.6" font-weight="700" fill="#1E293B">{html.escape(b_banner[1])}</text>'
    )

    svg_p = out_p / "s26_biorender_single_kd_summary.svg"
    png_p = out_p / "s26_biorender_single_kd_summary.png"
    c.save(svg_p, png_p)
    return svg_p, png_p


def generate_all_biorender_schematics(
    out_dir: str | Path,
    mode: str = "EV",
    top6_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Path]:
    """
    Generate all 7 modular BioRender schematics (`s01`, `s03`, `s14`, `s17`, `s24`, `s25`, `s26`)
    as both `.svg` and 300-DPI `.png` files, returning a dict mapping schematic key to PNG path.
    """
    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    _, p01 = build_s01_biorender_hero(out_p, mode=mode)
    _, p03 = build_s03_biorender_exp_design(out_p, mode=mode)
    _, p14 = build_s14_biorender_topology(out_p, mode=mode)
    _, p17 = build_s17_biorender_lock_and_key(out_p, mode=mode)
    _, p24 = build_s24_biorender_spotlight(out_p, top6_df=top6_df, mode=mode)
    _, p25 = build_s25_biorender_working_model(out_p, mode=mode)
    _, p26 = build_s26_single_kd_summary(out_p, mode=mode)
    return {
        "s01_hero": p01,
        "s03_exp_design": p03,
        "s14_topology": p14,
        "s17_lock_and_key": p17,
        "s24_spotlight": p24,
        "s25_working_model": p25,
        "s26_single_kd_summary": p26,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 6B: BioRender-Style 7-Module Vector SVG Canvas Engine")
    parser.add_argument("--out-svg", "-o", default="fig_mechanism_biorender.svg", help="Output SVG path or directory")
    parser.add_argument("--mode", "-m", default="EV", help="Biological mode (EV, WCL, BIOLOGICS, UNIFIED)")
    parser.add_argument("--all-modules", action="store_true", help="Generate all 7 modular schematics in parent dir")
    parser.add_argument("--title", default="Calibrated Proteomics Mechanistic Model", help="Diagram title")
    args = parser.parse_args()

    out_target = Path(args.out_svg)
    if args.all_modules or out_target.is_dir():
        target_dir = out_target if out_target.is_dir() else out_target.parent
        paths = generate_all_biorender_schematics(target_dir, mode=args.mode)
        print(f"[Stage 6B] Generated {len(paths)} modular BioRender schematics in {target_dir}")
    else:
        _, png_p = build_s25_biorender_working_model(out_target.parent, mode=args.mode, title=args.title)
        if out_target.name != "s25_biorender_working_model.svg":
            shutil.copyfile(out_target.parent / "s25_biorender_working_model.svg", out_target)
        print(f"[Stage 6B] Saved BioRender-style vector SVG to {out_target} and {png_p}")


if __name__ == "__main__":
    main()

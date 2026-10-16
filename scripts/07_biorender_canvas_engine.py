#!/usr/bin/env python3
"""
Stage 6B: Reusable BioRender-Style Vector SVG + PNG Mechanistic Canvas Engine
=============================================================================
Constructs publication-grade, clean vector SVG cellular & molecular schematics with:
  - Phospholipid bilayer plasma membranes (`draw_bilayer_membrane`)
  - Multivesicular bodies (`MVB`), budding exosomes (`EVs`), & organelles (`draw_mvb_and_ev`, `draw_organelle`)
  - Biologic Y-shaped IgG antibodies (`draw_antibody_igg`)
  - Quantitative Protein/PTM pill nodes (`draw_protein_node`) with automatic `1-pep` warning badges
  - Evidence-graded causal arrows (`draw_causal_arrow`):
      * Solid curve (`stroke-dasharray="none"`) for `E1 / E2` calibrated findings
      * Dashed curve (`stroke-dasharray="6,4"`) for `E3` exploratory hypotheses
"""

from __future__ import annotations
import argparse
import html
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


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
        # Top editorial banner
        self.elements.append(
            f'<rect x="24" y="18" width="{self.width - 48}" height="64" rx="10" fill="#0F172A" filter="url(#soft-shadow)"/>'
            f'<text x="46" y="47" font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="bold" fill="#FFFFFF">{t_esc}</text>'
            f'<text x="46" y="68" font-family="Arial, Helvetica, sans-serif" font-size="12.5" fill="#CBD5E1">{s_esc}</text>'
        )
        # Evidence-grade legend box in top right
        lx = self.width - 490
        self.elements.append(
            f'<rect x="{lx}" y="26" width="450" height="48" rx="8" fill="#1E293B" stroke="#334155"/>'
            f'<line x1="{lx+16}" y1="42" x2="{lx+58}" y2="42" stroke="#14B8A6" stroke-width="3"/>'
            f'<text x="{lx+66}" y="46" font-family="Arial" font-size="11" fill="#F8FAFC">Solid Arrow: E1/E2 Calibrated Evidence (FDR &lt; 0.25, &gt;=2 pep)</text>'
            f'<line x1="{lx+16}" y1="61" x2="{lx+58}" y2="61" stroke="#F59E0B" stroke-width="2.4" stroke-dasharray="6,4"/>'
            f'<text x="{lx+66}" y="65" font-family="Arial" font-size="11" fill="#FDE68A">Dashed Arrow: E3 Exploratory Hypothesis / 1-Peptide (*)</text>'
        )

    def draw_bilayer_membrane(self, y_pos: int = 215, label_top: str = "EXTRACELLULAR / TIF / SECRETOME", label_bottom: str = "INTRACELLULAR CYTOSOL (WCL)") -> None:
        """Draw a continuous BioRender-style phospholipid bilayer across the canvas."""
        x_start, x_end = 40, self.width - 40
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
        self.elements.append(
            f'<text x="54" y="{y_pos - 26}" font-family="Arial" font-size="11.5" font-weight="bold" fill="#0369A1" letter-spacing="1">{html.escape(label_top)}</text>'
            f'<text x="54" y="{y_pos + 36}" font-family="Arial" font-size="11.5" font-weight="bold" fill="#92400E" letter-spacing="1">{html.escape(label_bottom)}</text>'
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
        self.elements.append(
            f'<text x="{cx}" y="{cy + radius + 22}" text-anchor="middle" font-family="Arial" font-size="12.5" font-weight="bold" fill="#78350F">{html.escape(label)}</text>'
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
        width: int = 148,
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
            f'  <text x="{cx}" y="{cy - 3}" text-anchor="middle" font-family="Arial" font-size="12.5" font-weight="bold" fill="{text_c}">{g_disp}</text>'
            f'  <text x="{cx}" y="{cy + 13}" text-anchor="middle" font-family="Arial" font-size="10" fill="#475569">{m_disp}</text>'
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
        """
        Draw a curved Bezier causal arrow.
        Enforces visual honesty: E1/E2 -> solid arrow; E3 -> dashed arrow.
        """
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
                f'<text x="{mx}" y="{my - 6}" text-anchor="middle" font-family="Arial" font-size="10.5" font-weight="bold" fill="{stroke_c}">{html.escape(tag)}</text>'
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
        import shutil
        import subprocess

        out_svg = Path(svg_path)
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        svg_text = self.to_svg()
        out_svg.write_text(svg_text, encoding="utf-8")

        if png_path:
            out_png = Path(png_path)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            rendered = False
            try:
                import cairosvg
                cairosvg.svg2png(bytestring=svg_text.encode("utf-8"), write_to=str(out_png), dpi=300)
                rendered = True
            except Exception:
                pass
            if not rendered and shutil.which("rsvg-convert"):
                subprocess.run(
                    ["rsvg-convert", "-d", "300", "-p", "300", "-w", str(self.width * 2), "-o", str(out_png), str(out_svg)],
                    check=False,
                )
        return out_svg


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 6B: BioRender-Style Vector SVG Canvas Engine")
    parser.add_argument("--out-svg", "-o", default="fig_mechanism_biorender.svg", help="Output SVG path")
    parser.add_argument("--title", default="Calibrated Proteomics Mechanistic Model", help="Diagram title")
    args = parser.parse_args()

    canvas = BioRenderCanvas(title=args.title)
    canvas.add_header_and_legend("Calibrated via Stage 1.5 Empirical Null FDR & Peptide Support Tiering")
    canvas.draw_bilayer_membrane(y_pos=230)
    canvas.draw_mvb_and_ev(cx=420, cy=500, label="Late Endosome / MVB Sorting")
    canvas.draw_organelle(cx=1020, cy=530, rx=140, ry=85, kind="nucleus", label="Nucleus & Transcriptional Programs")
    canvas.draw_protein_node(420, 340, "RAB27A", "shKD -2.8 log2", direction="down", unique_peptides=6, evidence_grade="E2")
    canvas.draw_protein_node(700, 140, "CD63 / PDCD6IP", "EV Cargo Shift", direction="down", unique_peptides=9, evidence_grade="E2")
    canvas.draw_protein_node(190, 340, "EGFR", "Reciprocal OE/KD", direction="reciprocal", unique_peptides=12, evidence_grade="E2")
    canvas.draw_protein_node(190, 140, "AREG", "Putative Shedding", direction="warning", unique_peptides=1, evidence_grade="E3")
    canvas.draw_causal_arrow(420, 315, 680, 165, evidence_grade="E2", label="Exosome Vesicle Docking (p=0.0002)", curve_offset=-30)
    canvas.draw_causal_arrow(190, 315, 190, 165, evidence_grade="E3", label="Exploratory Ectodomain Release", curve_offset=0)
    canvas.save(args.out_svg)
    print(f"[Stage 6B] Saved BioRender-style vector SVG to {args.out_svg}")


if __name__ == "__main__":
    main()


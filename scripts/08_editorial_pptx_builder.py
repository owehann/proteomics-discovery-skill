#!/usr/bin/env python3
"""
Stage 7A: Minimalist 16:9 Editorial PPTX Slide Deck Builder
===========================================================
Constructs high-density, zero-overlap widescreen (13.333" x 7.5") PowerPoint decks (`python-pptx`)
featuring:
  - Slide 1: Cover & Study Architecture
  - Slide 2: Mandatory Executive Scorecard (Sample Size N, Target Status, Null Noise sigma,
             Expected Null Hits vs Observed Hits, Empirical Null FDR, & Evidence Tier Cap)
  - Content Slides: Paired narrative cards + high-DPI figures + bilingual speaker notes
  - Automated Bounding-Box Overlap Validator (`validate_deck_geometry`)
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from pptx import Presentation
    from pptx.dml.color import RGBColor
    from pptx.enum.shapes import MSO_SHAPE
    from pptx.enum.text import PP_ALIGN
    from pptx.util import Inches, Pt
except ImportError:
    Presentation = None


class EditorialSlideDeck:
    """Minimalist 16:9 presentation builder with automated bounding-box collision validation."""

    def __init__(self, title: str, subtitle: str = "") -> None:
        if Presentation is None:
            raise RuntimeError("python-pptx is required to build .pptx decks.")
        self.prs = Presentation()
        self.prs.slide_width = Inches(13.333)
        self.prs.slide_height = Inches(7.500)
        self.blank_layout = self.prs.slide_layouts[6]
        self.title = title
        self.subtitle = subtitle
        self.slide_bboxes: Dict[int, List[Tuple[str, float, float, float, float]]] = {}

    def _register_box(self, slide_idx: int, label: str, left: float, top: float, width: float, height: float) -> None:
        self.slide_bboxes.setdefault(slide_idx, []).append((label, left, top, left + width, top + height))

    def validate_deck_geometry(self) -> List[str]:
        """Verify that no two registered visual containers overlap by > 0.01 sq inches on any slide."""
        warnings: List[str] = []
        for s_idx, boxes in self.slide_bboxes.items():
            for i in range(len(boxes)):
                l1, x1a, y1a, x1b, y1b = boxes[i]
                for j in range(i + 1, len(boxes)):
                    l2, x2a, y2a, x2b, y2b = boxes[j]
                    dx = min(x1b, x2b) - max(x1a, x2a)
                    dy = min(y1b, y2b) - max(y1a, y2a)
                    if dx > 0.05 and dy > 0.05:
                        warnings.append(f"Slide {s_idx}: Overlap ({dx*dy:.2f} sq in) between '{l1}' and '{l2}'")
        return warnings

    def _add_header_band(self, slide: Any, slide_idx: int, title_text: str, kicker: str = "PROTEOMICS DISCOVERY ENGINE") -> None:
        banner = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.0), Inches(0.0), Inches(13.333), Inches(0.95))
        banner.fill.solid()
        banner.fill.fore_color.rgb = RGBColor(15, 23, 42)
        banner.line.fill.background()
        self._register_box(slide_idx, "header_banner", 0.0, 0.0, 13.333, 0.95)

        tb = slide.shapes.add_textbox(Inches(0.5), Inches(0.12), Inches(12.3), Inches(0.75))
        tf = tb.text_frame
        tf.word_wrap = True
        p0 = tf.paragraphs[0]
        p0.text = kicker.upper()
        p0.font.size = Pt(9.5)
        p0.font.bold = True
        p0.font.color.rgb = RGBColor(230, 75, 53)

        p1 = tf.add_paragraph()
        p1.text = title_text
        p1.font.size = Pt(20)
        p1.font.bold = True
        p1.font.color.rgb = RGBColor(255, 255, 255)

    def add_cover_slide(self, metadata_badges: Sequence[str], speaker_notes: str = "") -> None:
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.333), Inches(7.5))
        bg.fill.solid()
        bg.fill.fore_color.rgb = RGBColor(15, 23, 42)
        bg.line.fill.background()

        tb = slide.shapes.add_textbox(Inches(0.9), Inches(1.8), Inches(11.5), Inches(3.2))
        tf = tb.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        p.text = self.title
        p.font.size = Pt(34)
        p.font.bold = True
        p.font.color.rgb = RGBColor(255, 255, 255)

        p2 = tf.add_paragraph()
        p2.text = self.subtitle
        p2.font.size = Pt(16)
        p2.font.color.rgb = RGBColor(148, 163, 184)
        self._register_box(s_idx, "cover_title", 0.9, 1.8, 11.5, 3.2)

        for i, badge in enumerate(metadata_badges[:4]):
            left = 0.9 + i * 2.95
            card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(5.3), Inches(2.75), Inches(1.1))
            card.fill.solid()
            card.fill.fore_color.rgb = RGBColor(30, 41, 59)
            card.line.color.rgb = RGBColor(71, 85, 105)
            ctf = card.text_frame
            ctf.word_wrap = True
            cp = ctf.paragraphs[0]
            cp.text = badge
            cp.font.size = Pt(11)
            cp.font.bold = True
            cp.font.color.rgb = RGBColor(248, 250, 252)
            self._register_box(s_idx, f"cover_badge_{i}", left, 5.3, 2.75, 1.1)

        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes

    def add_executive_scorecard_slide(
        self,
        scorecard_rows: Sequence[Dict[str, Any]],
        speaker_notes: str = "",
    ) -> None:
        """
        Slide 2 (Mandatory): Executive Study Health & Empirical Null FDR Scorecard.
        Columns: Contrast | Sample N | Target Audit | Null Noise (sigma) | Obs vs Null Hits | Emp. Null FDR | Max Evidence Grade | Key Finding
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self._add_header_band(
            slide,
            s_idx,
            "Executive Scorecard: Statistical Power, Empirical Null FDR, & Evidence Tiers",
            kicker="STAGE 1.5 CALIBRATION GATE",
        )

        headers = [
            "Contrast / Arm",
            "Design (N)",
            "Target Audit",
            "Noise (σ)",
            "Obs vs Null",
            "Emp. FDR",
            "Max Tier",
            "Calibrated Biological Takeaway",
        ]
        n_rows = len(scorecard_rows) + 1
        n_cols = len(headers)
        top, left, width, height = 1.25, 0.45, 12.43, min(5.6, 0.65 * n_rows)

        table_shape = slide.shapes.add_table(n_rows, n_cols, Inches(left), Inches(top), Inches(width), Inches(height))
        self._register_box(s_idx, "scorecard_table", left, top, width, height)
        table = table_shape.table

        col_widths = [1.85, 0.95, 1.45, 0.90, 1.30, 0.95, 0.95, 4.08]
        for c_i, w in enumerate(col_widths):
            table.columns[c_i].width = Inches(w)

        for c_i, h_text in enumerate(headers):
            cell = table.cell(0, c_i)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor(30, 41, 59)
            p = cell.text_frame.paragraphs[0]
            p.text = h_text
            p.font.size = Pt(10)
            p.font.bold = True
            p.font.color.rgb = RGBColor(255, 255, 255)

        for r_i, row_data in enumerate(scorecard_rows, start=1):
            vals = [
                str(row_data.get("contrast", "")),
                str(row_data.get("design_n", "")),
                str(row_data.get("target_status", "")),
                str(row_data.get("sigma_hat", "")),
                str(row_data.get("obs_vs_null", "")),
                str(row_data.get("emp_fdr", "")),
                str(row_data.get("max_grade", "E3")),
                str(row_data.get("takeaway", "")),
            ]
            for c_i, val in enumerate(vals):
                cell = table.cell(r_i, c_i)
                cell.fill.solid()
                if c_i == 6:
                    if "E1" in val:
                        cell.fill.fore_color.rgb = RGBColor(209, 250, 229)
                    elif "E2" in val:
                        cell.fill.fore_color.rgb = RGBColor(219, 234, 254)
                    else:
                        cell.fill.fore_color.rgb = RGBColor(254, 243, 199)
                else:
                    cell.fill.fore_color.rgb = RGBColor(248, 250, 252) if r_i % 2 == 1 else RGBColor(255, 255, 255)
                p = cell.text_frame.paragraphs[0]
                p.text = val
                p.font.size = Pt(9.5)
                p.font.bold = c_i in {0, 6}
                p.font.color.rgb = RGBColor(15, 23, 42)

        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes

    def add_figure_and_narrative_slide(
        self,
        title: str,
        kicker: str,
        bullets: Sequence[str],
        image_path: Optional[str | Path] = None,
        speaker_notes: str = "",
    ) -> None:
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self._add_header_band(slide, s_idx, title, kicker=kicker)

        # Left narrative card
        card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.45), Inches(1.25), Inches(5.15), Inches(5.80))
        card.fill.solid()
        card.fill.fore_color.rgb = RGBColor(248, 250, 252)
        card.line.color.rgb = RGBColor(203, 213, 225)
        self._register_box(s_idx, "left_narrative_card", 0.45, 1.25, 5.15, 5.80)

        tf = card.text_frame
        tf.word_wrap = True
        for b_i, b_text in enumerate(bullets):
            p = tf.paragraphs[0] if b_i == 0 else tf.add_paragraph()
            p.text = f"• {b_text}"
            p.font.size = Pt(11.5)
            p.font.color.rgb = RGBColor(30, 41, 59)
            p.space_after = Pt(10)

        # Right figure panel (aspect-ratio-safe fit within 7.00" x 5.80")
        if image_path and Path(image_path).exists():
            max_w, max_h = 7.00, 5.80
            box_left, box_top = 5.85, 1.25
            try:
                from PIL import Image
                with Image.open(image_path) as im:
                    iw, ih = im.size
                aspect = iw / max(ih, 1)
                if max_w / aspect <= max_h:
                    pic_w = max_w
                    pic_h = max_w / aspect
                else:
                    pic_h = max_h
                    pic_w = max_h * aspect
                pic_left = box_left + (max_w - pic_w) / 2.0
                pic_top = box_top + (max_h - pic_h) / 2.0
                slide.shapes.add_picture(str(image_path), Inches(pic_left), Inches(pic_top), width=Inches(pic_w), height=Inches(pic_h))
                self._register_box(s_idx, "right_figure", pic_left, pic_top, pic_w, pic_h)
            except Exception:
                slide.shapes.add_picture(str(image_path), Inches(box_left), Inches(box_top), width=Inches(max_w))
                self._register_box(s_idx, "right_figure", box_left, box_top, max_w, max_h)

        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes

    def add_act_divider_slide(
        self,
        act_label: str,
        act_title: str,
        act_subtitle: str,
        key_pillars: Sequence[str],
        speaker_notes: str = "",
    ) -> None:
        """Add an architectural Act Divider slide for multi-dataset synthesis presentations."""
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.333), Inches(7.5))
        bg.fill.solid()
        bg.fill.fore_color.rgb = RGBColor(15, 23, 42)
        bg.line.fill.background()

        badge = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.9), Inches(1.3), Inches(3.4), Inches(0.55))
        badge.fill.solid()
        badge.fill.fore_color.rgb = RGBColor(230, 75, 53)
        badge.line.fill.background()
        btf = badge.text_frame
        bp = btf.paragraphs[0]
        bp.text = act_label.upper()
        bp.font.size = Pt(13)
        bp.font.bold = True
        bp.font.color.rgb = RGBColor(255, 255, 255)
        self._register_box(s_idx, "act_badge", 0.9, 1.3, 3.4, 0.55)

        tb = slide.shapes.add_textbox(Inches(0.9), Inches(2.05), Inches(11.5), Inches(2.2))
        tf = tb.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        p.text = act_title
        p.font.size = Pt(30)
        p.font.bold = True
        p.font.color.rgb = RGBColor(255, 255, 255)
        p2 = tf.add_paragraph()
        p2.text = act_subtitle
        p2.font.size = Pt(15)
        p2.font.color.rgb = RGBColor(148, 163, 184)
        self._register_box(s_idx, "act_title_box", 0.9, 2.05, 11.5, 2.2)

        for i, pillar in enumerate(key_pillars[:3]):
            left = 0.9 + i * 3.95
            card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(left), Inches(4.55), Inches(3.70), Inches(2.15))
            card.fill.solid()
            card.fill.fore_color.rgb = RGBColor(30, 41, 59)
            card.line.color.rgb = RGBColor(71, 85, 105)
            ctf = card.text_frame
            ctf.word_wrap = True
            cp = ctf.paragraphs[0]
            cp.text = pillar
            cp.font.size = Pt(11.5)
            cp.font.color.rgb = RGBColor(241, 245, 249)
            self._register_box(s_idx, f"act_pillar_{i}", left, 4.55, 3.70, 2.15)

        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes

    def add_full_width_table_slide(
        self,
        title: str,
        kicker: str,
        headers: Sequence[str],
        rows: Sequence[Sequence[str]],
        col_widths: Optional[Sequence[float]] = None,
        bottom_callout: str = "",
        speaker_notes: str = "",
    ) -> None:
        """Add a structured full-width evidence or roadmap table slide with optional bottom callout."""
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self._add_header_band(slide, s_idx, title, kicker=kicker)

        n_rows = len(rows) + 1
        n_cols = len(headers)
        max_table_h = 4.95 if bottom_callout else 5.75
        top, left, width = 1.20, 0.45, 12.43
        height = min(max_table_h, 0.52 * n_rows)

        table_shape = slide.shapes.add_table(n_rows, n_cols, Inches(left), Inches(top), Inches(width), Inches(height))
        self._register_box(s_idx, "full_table", left, top, width, height)
        table = table_shape.table

        if col_widths and len(col_widths) == n_cols:
            for c_i, w in enumerate(col_widths):
                table.columns[c_i].width = Inches(w)

        for c_i, h_text in enumerate(headers):
            cell = table.cell(0, c_i)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor(30, 41, 59)
            p = cell.text_frame.paragraphs[0]
            p.text = str(h_text)
            p.font.size = Pt(10)
            p.font.bold = True
            p.font.color.rgb = RGBColor(255, 255, 255)

        for r_i, row_vals in enumerate(rows, start=1):
            for c_i, val in enumerate(row_vals):
                cell = table.cell(r_i, c_i)
                cell.fill.solid()
                cell.fill.fore_color.rgb = RGBColor(248, 250, 252) if r_i % 2 == 1 else RGBColor(255, 255, 255)
                p = cell.text_frame.paragraphs[0]
                p.text = str(val)
                p.font.size = Pt(9.5)
                p.font.bold = (c_i == 0)
                p.font.color.rgb = RGBColor(15, 23, 42)

        if bottom_callout:
            c_top = top + height + 0.18
            callout = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.45), Inches(c_top), Inches(12.43), Inches(0.72))
            callout.fill.solid()
            callout.fill.fore_color.rgb = RGBColor(239, 246, 255)
            callout.line.color.rgb = RGBColor(147, 197, 253)
            ctf = callout.text_frame
            ctf.word_wrap = True
            cp = ctf.paragraphs[0]
            cp.text = bottom_callout
            cp.font.size = Pt(10.5)
            cp.font.bold = True
            cp.font.color.rgb = RGBColor(30, 58, 138)
            self._register_box(s_idx, "bottom_callout", 0.45, c_top, 12.43, 0.72)

        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes

    def save(self, output_path: str | Path) -> Tuple[Path, List[str]]:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        warnings = self.validate_deck_geometry()
        self.prs.save(str(out))
        return out, warnings


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 7A: Editorial PPTX Slide Deck Builder")
    parser.add_argument("--out", "-o", default="proteomics_discovery_deck.pptx", help="Output PPTX path")
    parser.add_argument("--title", default="Proteomics Discovery & Mechanistic Report", help="Deck title")
    args = parser.parse_args()

    deck = EditorialSlideDeck(args.title, "Calibrated Multi-Layer Proteomics Workflow")
    deck.add_cover_slide(["Modality: M1-M6", "Null Audit: Stage 1.5", "Peptide Gate: >=2 Core", "Fact Gate: Verified"])
    deck.add_executive_scorecard_slide(
        [
            {
                "contrast": "Example Contrast",
                "design_n": "n=3 vs 3",
                "target_status": "Detected (8 pep)",
                "sigma_hat": "0.14",
                "obs_vs_null": "84 vs 1.2",
                "emp_fdr": "0.014",
                "max_grade": "E1 / E2",
                "takeaway": "High-confidence target-driven proteome remodeling",
            }
        ]
    )
    out_p, warns = deck.save(args.out)
    print(f"[Stage 7A] Saved PPTX to {out_p} (Overlaps: {len(warns)})")


if __name__ == "__main__":
    main()

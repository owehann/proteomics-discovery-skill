#!/usr/bin/env python3
"""
Step 7A: Editorial Minimalist 16:9 Group-Meeting & Self-Review PPTX Builder
===========================================================================
Constructs zero-overlap, publication-grade widescreen (13.333" x 7.5") PowerPoint decks (`python-pptx`)
governed by strict visual design constraints and collision-proof text flow:

1. Collision-Proof Unified `TextFrame` Architecture:
   - Every multi-line component (`stat()` KPI callout, Cover badge, Layout D mechanism card,
     Layout E roadmap row) places its bold headline (`Paragraph 0`) and muted caption (`Paragraph 1`)
     inside ONE unified `TextFrame` (`add_textbox`), so PowerPoint natively flows the caption
     below the headline even if the headline wraps to 2 lines — making vertical text collisions
     between a headline and its caption physically impossible.
   - Automatic single-line font scaling (`_auto_scale_headline_pt`) dynamically adjusts `stat()`
     and `action_title` font sizes based on string character width and container width `w` so
     long gene lists or metrics stay on a single line without mid-word truncation.

2. Airy Editorial Header (No Heavy Dark Top Banner):
   - Pure white (#FFFFFF) or warm paper (#FAFAF7) background.
   - Muted uppercase step/section overline at y=0.24" (10.5pt Arial Bold, #8A8984).
   - Single-line concise Action Title at y=0.50" (auto-scaled 16.5-23pt Arial Bold, #1F1F1D).
   - Ultra-thin hairline divider (0.012" height, #D9D8D3) at y=1.02".
   - Minimalist bottom footer at y=7.14" (study title left, slide counter right).

3. 6 Alternating Layout Archetypes (Layouts A-F, max 2 consecutive repeats):
   - Layout A: Full-Bleed Visual Canvas (Cover Graphical Abstract, Multi-Panel Spotlights, BioRender Model)
   - Layout B: Wide Figure Top (12.23" x 4.62") + 3 Horizontal `stat()` Callouts Bottom
   - Layout C: Large Figure Left (8.45" x 5.82") + 4 Vertical `stat()` Callouts Right
   - Layout D: Dual Side-by-Side Split (Quantitative Data Left + BioRender Mechanism / Native Cards Right)
   - Layout E: 100% Native Editable Vector Workflow / 5-Step Numbered Wet-Lab Validation Roadmap
   - Layout F: 100% Native Editable Evidence Tier & Clinical Tractability Matrix Table

4. Direct Bilingual Notes (No Ceremonial Greetings):
   - Embeds `【中文深度解析与汇报要点】` and `[English Detailed Walkthrough]` directly into `notes_slide`
     and exports a companion `*_Notes_Bilingual.md` without ceremonial openings ("各位老师同学好").
"""

from __future__ import annotations
import argparse
import math
import sys
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
    RGBColor = None
    MSO_SHAPE = None
    PP_ALIGN = None
    Inches = None
    Pt = None


# ── Editorial Color Tokens (NPG / Warm-Minimalist Academic Palette) ──────────
C_BG_WHITE = RGBColor(255, 255, 255) if RGBColor else None       # #FFFFFF
C_BG_PAPER = RGBColor(250, 250, 247) if RGBColor else None       # #FAFAF7
C_CARD_BG = RGBColor(246, 245, 240) if RGBColor else None        # #F6F5F0
C_INK = RGBColor(31, 31, 29) if RGBColor else None               # #1F1F1D
C_MUTED = RGBColor(90, 89, 85) if RGBColor else None             # #5A5955
C_LIGHT = RGBColor(138, 137, 132) if RGBColor else None          # #8A8984
C_RULE = RGBColor(217, 216, 211) if RGBColor else None           # #D9D8D3

# Semantic Accent Tokens
C_CORAL = RGBColor(200, 62, 43) if RGBColor else None            # #C83E2B (Upregulated / Activated)
C_SLATE_BLUE = RGBColor(42, 98, 143) if RGBColor else None       # #2A628F (Downregulated / Depleted)
C_EMERALD = RGBColor(42, 123, 92) if RGBColor else None          # #2A7B5C (QC Pass / Rescue / E1)
C_AMBER = RGBColor(217, 119, 6) if RGBColor else None            # #D97706 (PTM / Shedding / E3)
C_PLUM = RGBColor(109, 40, 217) if RGBColor else None            # #6D28D9 (Regulator / Complex)

SEMANTIC_COLOR_MAP = {
    "coral": C_CORAL,
    "red": C_CORAL,
    "up": C_CORAL,
    "blue": C_SLATE_BLUE,
    "down": C_SLATE_BLUE,
    "teal": C_EMERALD,
    "emerald": C_EMERALD,
    "green": C_EMERALD,
    "amber": C_AMBER,
    "gold": C_AMBER,
    "plum": C_PLUM,
    "purple": C_PLUM,
    "ink": C_INK,
}

FORBIDDEN_GREETINGS = (
    "各位老师同学好",
    "各位老师、同学好",
    "各位老师好",
    "大家好",
    "hello everyone",
    "good morning everyone",
    "good afternoon everyone",
)

# Cryptic standalone bioinformatic jargon phrases banned from uncontextualized slide titles/headlines
BANNED_OPAQUE_JARGON_HEADLINES = (
    "kolde beta order-statistic",
    "dirichlet stability",
    "archetype a1",
    "archetype a2",
    "archetype a3",
    "archetype a4",
    "archetype a5",
    "ulm z-score",
)

JARGON_TO_PLAIN_BIOLOGY_MAP = {
    "Kolde Beta Order-Statistic": "5-Evidence Priority Score",
    "1,000x Dirichlet Stability": "1,000x Weight Robustness",
    "Stage 1.5 Shared-Control Null": "False-Positive Noise Audit",
    "Degree-Corrected RWR": "Physical Hub & Bridge Scan",
    "CollecTRI ULM": "Master Regulator Activity",
    "DEqMS Variance Shrinkage": "Multi-Peptide Noise Filter",
    "TOST Equivalence": "Equivalent within ±δ",
}


def polish_jargon_to_biology(text: str) -> str:
    """Replace opaque statistical jargon with biologist-friendly plain English terms."""
    out = str(text)
    for k, v in JARGON_TO_PLAIN_BIOLOGY_MAP.items():
        if k in out:
            out = out.replace(k, v)
    return out


def _effective_char_units(text: str) -> float:
    """Estimate relative visual width of a string in Arial Bold (uppercase/CJK/digits are wider)."""
    total = 0.0
    for ch in str(text):
        if ord(ch) > 127 and ch not in "−–—±×·":
            total += 1.65
        elif ch.isupper() or ch in "MWT%@+#/":
            total += 1.15
        elif ch in "il1.,:;|! ":
            total += 0.55
        else:
            total += 0.95
    return max(total, 1.0)


def _auto_scale_headline_pt(text: str, box_width_in: float, max_pt: float = 22.0, min_pt: float = 13.0) -> float:
    """
    Compute the largest font size in [min_pt, max_pt] that allows `text` to fit on a SINGLE line
    inside `box_width_in` inches without wrapping.
    """
    usable_w = max(box_width_in - 0.10, 0.8)
    eff_len = _effective_char_units(text)
    # 1 unit in Arial Bold at `pt` points is approx `0.58 * pt / 72.0` inches
    fit_pt = (usable_w * 72.0) / (eff_len * 0.58)
    return round(min(max_pt, max(min_pt, fit_pt)), 1)


def _estimate_rendered_height_in(paragraphs: Sequence[Tuple[str, float]], box_width_in: float, extra_pad_in: float = 0.08) -> float:
    """Estimate actual vertical height occupied by wrapped paragraphs inside a PowerPoint TextFrame."""
    usable_w = max(box_width_in - 0.08, 0.6)
    total_h = extra_pad_in
    for text, pt_size in paragraphs:
        eff_w_in = _effective_char_units(text) * (0.56 * pt_size / 72.0)
        n_lines = max(1, int(math.ceil(eff_w_in / usable_w)))
        line_h_in = (pt_size * 1.28) / 72.0
        total_h += n_lines * line_h_in
    return total_h


def format_bilingual_notes(zh_points: Sequence[str] | str, en_points: Sequence[str] | str) -> str:
    """
    Format structured bilingual slide notes for both solo reading and group-meeting reference.
    Strictly omits ceremonial greetings ('各位老师同学好' / 'Hello everyone').
    """
    if isinstance(zh_points, str):
        zh_body = zh_points.strip()
    else:
        zh_body = "\n".join(f"- {pt.strip()}" for pt in zh_points if str(pt).strip())

    if isinstance(en_points, str):
        en_body = en_points.strip()
    else:
        en_body = "\n".join(f"- {pt.strip()}" for pt in en_points if str(pt).strip())

    return (
        "【中文深度解析与看图指南 (Why · What to Look At · Core Genes & Pathways · Next Step)】\n"
        f"{zh_body}\n\n"
        "[English Visual Reading Guide & Biological Walkthrough]\n"
        f"{en_body}"
    )


class EditorialSlideDeck:
    """
    Minimalist 16:9 scientific presentation builder enforcing:
      - Airy pure-white header (no heavy dark top bar)
      - Single-TextFrame unboxed `stat()` KPI callouts (100% immune to headline/caption overlap)
      - Automatic single-line font scaling for long gene names or metric strings
      - 6 alternating layout archetypes (Layouts A-F)
      - Real wrapped-height bounding-box collision & design constraint verification
    """

    def __init__(
        self,
        title: str,
        subtitle: str = "",
        total_expected_slides: int = 22,
        enforce_comprehensive: bool = False,
        min_comprehensive_slides: int = 20,
    ) -> None:
        if Presentation is None:
            raise RuntimeError("python-pptx is required to build .pptx decks.")
        self.prs = Presentation()
        self.prs.slide_width = Inches(13.333)
        self.prs.slide_height = Inches(7.500)
        self.blank_layout = self.prs.slide_layouts[6]
        self.title = title
        self.subtitle = subtitle
        self.total_expected_slides = total_expected_slides
        self.enforce_comprehensive = enforce_comprehensive
        self.min_comprehensive_slides = min_comprehensive_slides
        self.slide_bboxes: Dict[int, List[Tuple[str, float, float, float, float]]] = {}
        self.slide_layouts_used: List[str] = []
        self.slide_titles: List[str] = []
        self.slide_stat_labels: List[Tuple[int, str, str]] = []
        self.slide_images_used: List[Tuple[int, str, str]] = []
        self.bilingual_notes_log: List[Dict[str, Any]] = []

    def _resolve_color(self, color_or_key: Any) -> Any:
        if isinstance(color_or_key, str):
            return SEMANTIC_COLOR_MAP.get(color_or_key.lower(), C_SLATE_BLUE)
        return color_or_key if color_or_key is not None else C_SLATE_BLUE

    def _register_box(self, slide_idx: int, label: str, left: float, top: float, width: float, height: float) -> None:
        self.slide_bboxes.setdefault(slide_idx, []).append((label, left, top, left + width, top + height))

    def _set_bg(self, slide: Any, color: Any = None) -> None:
        bg = slide.background
        fill = bg.fill
        fill.solid()
        fill.fore_color.rgb = color if color is not None else C_BG_WHITE

    def _add_textbox(
        self,
        slide: Any,
        slide_idx: int,
        label: str,
        x: float,
        y: float,
        w: float,
        h: float,
        text: str,
        size: float = 13,
        bold: bool = False,
        color: Any = None,
        align: Any = None,
        register: bool = True,
    ) -> Any:
        box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        if register:
            real_h = max(h, _estimate_rendered_height_in([(str(text), size)], w, extra_pad_in=0.04))
            self._register_box(slide_idx, label, x, y, w, real_h)
        tf = box.text_frame
        tf.word_wrap = True
        tf.margin_left = Inches(0.02)
        tf.margin_right = Inches(0.02)
        tf.margin_top = Inches(0.01)
        tf.margin_bottom = Inches(0.01)
        p = tf.paragraphs[0]
        p.text = str(text)
        p.font.name = "Arial"
        p.font.size = Pt(size)
        p.font.bold = bold
        p.font.color.rgb = self._resolve_color(color) if color is not None else C_INK
        p.alignment = align if align is not None else PP_ALIGN.LEFT
        return box

    def _add_stacked_textbox(
        self,
        slide: Any,
        slide_idx: int,
        label: str,
        x: float,
        y: float,
        w: float,
        h: float,
        headline_text: str,
        caption_text: str,
        headline_pt: float = 20.0,
        caption_pt: float = 11.0,
        headline_color: Any = None,
        caption_color: Any = None,
        space_between_pt: float = 3.0,
        register: bool = True,
    ) -> Any:
        """
        Render a bold headline (`Paragraph 0`) and muted caption (`Paragraph 1`) inside a SINGLE
        unified PowerPoint `TextFrame` so they can NEVER overlap even if the headline wraps.
        """
        box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        tf = box.text_frame
        tf.word_wrap = True
        tf.margin_left = Inches(0.02)
        tf.margin_right = Inches(0.02)
        tf.margin_top = Inches(0.01)
        tf.margin_bottom = Inches(0.01)

        p0 = tf.paragraphs[0]
        p0.text = str(headline_text).strip()
        p0.font.name = "Arial"
        p0.font.size = Pt(headline_pt)
        p0.font.bold = True
        p0.font.color.rgb = self._resolve_color(headline_color) if headline_color is not None else C_INK
        p0.space_after = Pt(space_between_pt)

        if str(caption_text).strip():
            p1 = tf.add_paragraph()
            p1.text = str(caption_text).strip()
            p1.font.name = "Arial"
            p1.font.size = Pt(caption_pt)
            p1.font.bold = False
            p1.font.color.rgb = self._resolve_color(caption_color) if caption_color is not None else C_MUTED
            p1.space_before = Pt(1.0)

        if register:
            real_h = max(
                h,
                _estimate_rendered_height_in(
                    [(str(headline_text), headline_pt), (str(caption_text), caption_pt)],
                    w,
                    extra_pad_in=0.06,
                ),
            )
            self._register_box(slide_idx, label, x, y, w, real_h)
        return box

    def _add_rule(self, slide: Any, x: float, y: float, w: float, color: Any = None, thickness: float = 0.012) -> Any:
        line = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(thickness))
        line.fill.solid()
        line.fill.fore_color.rgb = color if color is not None else C_RULE
        line.line.fill.background()
        return line

    def _add_editorial_chrome(self, slide: Any, slide_idx: int, kicker: str, action_title: str) -> None:
        """
        Render the signature airy editorial header and footer:
        - Uppercase muted kicker at y=0.22"
        - Auto-scaled Single-Line Action Title at y=0.48" (never wraps into divider rule)
        - Hairline rule at y=1.02"
        - Subtle footer & slide counter at y=7.14"
        """
        self._set_bg(slide, C_BG_WHITE)
        clean_title = str(action_title).strip()
        self.slide_titles.append(clean_title)

        title_pt = _auto_scale_headline_pt(clean_title, box_width_in=12.1, max_pt=22.5, min_pt=15.5)

        self._add_textbox(
            slide,
            slide_idx,
            "kicker",
            0.55,
            0.22,
            11.5,
            0.22,
            f"{slide_idx:02d} / {kicker.upper()}",
            size=10.0,
            bold=True,
            color=C_LIGHT,
        )
        self._add_textbox(
            slide,
            slide_idx,
            "action_title",
            0.55,
            0.48,
            12.2,
            0.46,
            clean_title,
            size=title_pt,
            bold=True,
            color=C_INK,
        )
        self._add_rule(slide, 0.55, 1.02, 12.23, C_RULE)

        total_n = max(self.total_expected_slides, slide_idx)
        self._add_textbox(
            slide,
            slide_idx,
            "footer_left",
            0.55,
            7.14,
            8.5,
            0.22,
            f"{self.title} · ProteoSkill Calibrated Discovery",
            size=9.0,
            bold=False,
            color=C_LIGHT,
        )
        self._add_textbox(
            slide,
            slide_idx,
            "footer_right",
            11.20,
            7.14,
            1.58,
            0.22,
            f"{slide_idx:02d} / {total_n:02d}",
            size=9.0,
            bold=True,
            color=C_LIGHT,
            align=PP_ALIGN.RIGHT,
        )

    def _add_stat_callout(
        self,
        slide: Any,
        slide_idx: int,
        tag_id: str,
        x: float,
        y: float,
        num_headline: str,
        caption: str,
        color: Any = "blue",
        w: float = 3.45,
        h: float = 1.02,
    ) -> None:
        """
        Collision-proof unboxed typographic KPI callout:
        - Uses ONE unified `TextFrame` (`_add_stacked_textbox`) so `num_headline` and `caption`
          can NEVER overlap each other even if `num_headline` is long.
        - Auto-scales `headline_pt` (`13.5pt` to `21.5pt`) so long gene names or metrics
          (e.g. '21/22 Tubes Passed QC' or 'ATXN7 -0.70 / MCL1 +0.66') fit on a single line.
        """
        clean_num = str(num_headline).strip()
        clean_cap = str(caption).strip()
        self.slide_stat_labels.append((slide_idx, clean_num, clean_cap))

        headline_pt = _auto_scale_headline_pt(clean_num, box_width_in=w, max_pt=21.5, min_pt=13.5)
        caption_pt = 11.0 if len(clean_cap) <= 52 else 10.0

        self._add_stacked_textbox(
            slide=slide,
            slide_idx=slide_idx,
            label=tag_id,
            x=x,
            y=y,
            w=w,
            h=h,
            headline_text=clean_num,
            caption_text=clean_cap,
            headline_pt=headline_pt,
            caption_pt=caption_pt,
            headline_color=self._resolve_color(color),
            caption_color=C_MUTED,
            space_between_pt=2.5,
            register=True,
        )

    def _add_contain_image(
        self,
        slide: Any,
        slide_idx: int,
        label: str,
        img_path: Optional[str | Path],
        x: float,
        y: float,
        max_w: float,
        max_h: float,
        placeholder_title: str = "High-DPI Scientific Figure Canvas",
    ) -> None:
        """Place image centered within (x, y, max_w, max_h) preserving exact aspect ratio."""
        if img_path and Path(img_path).exists():
            resolved_p = str(Path(img_path).resolve())
            self.slide_images_used.append((slide_idx, label, resolved_p))
            try:
                from PIL import Image
                with Image.open(img_path) as im:
                    iw, ih = im.size
                aspect = iw / max(ih, 1)
                if max_w / aspect <= max_h:
                    pic_w = max_w
                    pic_h = max_w / aspect
                else:
                    pic_h = max_h
                    pic_w = max_h * aspect
                pic_x = x + (max_w - pic_w) / 2.0
                pic_y = y + (max_h - pic_h) / 2.0
                slide.shapes.add_picture(str(img_path), Inches(pic_x), Inches(pic_y), width=Inches(pic_w), height=Inches(pic_h))
                self._register_box(slide_idx, label, pic_x, pic_y, pic_w, pic_h)
                return
            except Exception:
                slide.shapes.add_picture(str(img_path), Inches(x), Inches(y), width=Inches(max_w))
                self._register_box(slide_idx, label, x, y, max_w, max_h)
                return

        canvas = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(max_w), Inches(max_h))
        canvas.fill.solid()
        canvas.fill.fore_color.rgb = C_BG_PAPER
        canvas.line.color.rgb = C_RULE
        self._register_box(slide_idx, label, x, y, max_w, max_h)
        ctf = canvas.text_frame
        ctf.word_wrap = True
        cp = ctf.paragraphs[0]
        cp.text = placeholder_title
        cp.font.name = "Arial"
        cp.font.size = Pt(13)
        cp.font.bold = True
        cp.font.color.rgb = C_MUTED
        cp.alignment = PP_ALIGN.CENTER

    def _attach_notes(self, slide: Any, slide_idx: int, kicker: str, title: str, speaker_notes: str) -> None:
        if speaker_notes:
            slide.notes_slide.notes_text_frame.text = speaker_notes
            self.bilingual_notes_log.append(
                {
                    "slide_idx": slide_idx,
                    "kicker": kicker,
                    "title": title,
                    "notes": speaker_notes,
                }
            )

    @staticmethod
    def _parse_cover_badge(badge: Any, idx: int) -> Tuple[str, str]:
        """
        Parse a cover badge without ever truncating words mid-string:
        - Accepts `(headline, caption)` tuple OR `"Label: Value"` string.
        - If `"A: B"` is passed, places whichever side starts with a digit or is shorter (`<=22 chars`)
          as the bold metric headline, and the other side as the descriptive caption.
        """
        if isinstance(badge, (tuple, list)) and len(badge) >= 2:
            return str(badge[0]).strip(), str(badge[1]).strip()
        raw = str(badge).strip()
        if ":" in raw:
            left, right = [p.strip() for p in raw.split(":", 1)]
            # If left starts with a digit or is the concise metric (e.g. "22 DIA-MS Tubes: 3 Workbooks"), keep left as headline
            if (left and left[0].isdigit()) or (len(left) <= len(right) and not (right and right[0].isdigit())):
                return left, right
            return right, left
        return raw, f"Key Study Pillar {idx + 1}"

    # ── Layout A: Cover / Full-Bleed Visual Canvas ──────────────────────────
    def add_cover_slide(
        self,
        metadata_badges: Sequence[Any],
        speaker_notes: str = "",
        graphical_abstract_path: Optional[str | Path] = None,
        kicker: str = "CALIBRATED PROTEOMICS DISCOVERY & MECHANISM REPORT",
    ) -> None:
        """
        Layout A (Cover): Warm paper background (#FAFAF7), unified title/subtitle flow on left,
        4 collision-proof metric callouts at bottom-left, and Graphical Abstract preview on right.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_A_Cover")
        self.slide_titles.append(self.title[:48])
        self._set_bg(slide, C_BG_PAPER)

        self._add_textbox(slide, s_idx, "cover_kicker", 0.65, 0.52, 6.10, 0.28, kicker.upper(), size=10.5, bold=True, color=C_CORAL)

        # Unified Title + Subtitle TextFrame so a 2-3 line title automatically pushes subtitle down
        cover_title_pt = 28.0 if len(self.title) <= 44 else (24.0 if len(self.title) <= 64 else 20.5)
        self._add_stacked_textbox(
            slide=slide,
            slide_idx=s_idx,
            label="cover_title_block",
            x=0.65,
            y=0.88,
            w=6.10,
            h=2.45,
            headline_text=self.title,
            caption_text=self.subtitle,
            headline_pt=cover_title_pt,
            caption_pt=13.0,
            headline_color=C_INK,
            caption_color=C_MUTED,
            space_between_pt=10.0,
            register=True,
        )
        self._add_rule(slide, 0.65, 3.52, 5.95, C_RULE)

        palette_cycle = ["coral", "blue", "emerald", "plum"]
        for i, badge in enumerate(metadata_badges[:4]):
            headline, sublabel = self._parse_cover_badge(badge, i)
            col_x = 0.65 + (i % 2) * 3.08
            row_y = 3.78 + (i // 2) * 1.52
            self._add_stat_callout(
                slide,
                s_idx,
                f"cover_stat_{i}",
                col_x,
                row_y,
                headline,
                sublabel,
                color=palette_cycle[i % len(palette_cycle)],
                w=2.92,
                h=1.22,
            )

        self._add_contain_image(
            slide,
            s_idx,
            "cover_graphical_abstract",
            graphical_abstract_path,
            7.05,
            0.75,
            5.65,
            5.95,
            placeholder_title="Graphical Abstract & 7-Step ProteoSkill Discovery Architecture",
        )
        self._attach_notes(slide, s_idx, kicker, self.title, speaker_notes)

    def add_full_canvas_slide(
        self,
        title: str,
        kicker: str,
        image_path: Optional[str | Path] = None,
        speaker_notes: str = "",
    ) -> None:
        """
        Layout A (Full Canvas): 12.23" x 5.82" full-bleed figure canvas for Top-Hit Spotlights
        or the Final BioRender Mechanistic Working Model.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_A_FullCanvas")
        self._add_editorial_chrome(slide, s_idx, kicker, title)
        self._add_contain_image(slide, s_idx, "full_canvas_fig", image_path, 0.55, 1.16, 12.23, 5.82, placeholder_title=title)
        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Layout B: Wide Figure Top + 3 Horizontal `stat()` Callouts Bottom ───
    def add_top_figure_bottom_stats_slide(
        self,
        title: str,
        kicker: str,
        stats: Sequence[Tuple[str, str, str]],
        image_path: Optional[str | Path] = None,
        speaker_notes: str = "",
    ) -> None:
        """
        Layout B: Wide multi-panel figure across top (12.23" x 4.58"), hairline divider at y=5.84",
        and 3 horizontal collision-proof `stat()` KPI callouts across the bottom (`h=1.12"`).
        Each item in `stats` is `(num_headline, caption, color_key)`.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_B_TopFigBottomStats")
        self._add_editorial_chrome(slide, s_idx, kicker, title)

        self._add_contain_image(slide, s_idx, "top_wide_fig", image_path, 0.55, 1.16, 12.23, 4.58, placeholder_title=title)
        self._add_rule(slide, 0.55, 5.84, 12.23, C_RULE)

        x_positions = [0.60, 4.75, 8.90]
        for i, (num_h, cap, col) in enumerate(stats[:3]):
            self._add_stat_callout(slide, s_idx, f"b_stat_{i}", x_positions[i], 5.92, num_h, cap, color=col, w=3.80, h=1.12)

        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Layout C: Large Figure Left + 4 Vertical `stat()` Callouts Right ────
    def add_left_figure_right_stats_slide(
        self,
        title: str,
        kicker: str,
        stats: Sequence[Tuple[str, str, str]],
        image_path: Optional[str | Path] = None,
        speaker_notes: str = "",
    ) -> None:
        """
        Layout C: Large figure on left (8.35" x 5.82") + 4 stacked collision-proof `stat()`
        callouts on right (`x=9.12", w=3.68"`). Zero bullet lists or gray text boxes.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_C_LeftFigRightStats")
        self._add_editorial_chrome(slide, s_idx, kicker, title)

        self._add_contain_image(slide, s_idx, "left_large_fig", image_path, 0.55, 1.16, 8.35, 5.82, placeholder_title=title)

        y_positions = [1.22, 2.62, 4.02, 5.42]
        for i, (num_h, cap, col) in enumerate(stats[:4]):
            self._add_stat_callout(slide, s_idx, f"c_stat_{i}", 9.12, y_positions[i], num_h, cap, color=col, w=3.68, h=1.26)

        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Layout D: Dual Side-by-Side Split (Data Left + Mechanism Right) ─────
    def add_dual_data_and_mechanism_slide(
        self,
        title: str,
        kicker: str,
        left_image_path: Optional[str | Path] = None,
        right_image_path: Optional[str | Path] = None,
        right_mechanism_cards: Optional[Sequence[Dict[str, str]]] = None,
        bottom_stats: Sequence[Tuple[str, str, str]] = (),
        speaker_notes: str = "",
    ) -> None:
        """
        Layout D: Side-by-side co-display of Quantitative Statistical Plot (left 5.80") and
        BioRender Organelle/Topology Schematic OR Native Editable Mechanism Comparison Cards (right 6.18"),
        paired with 3 bottom `stat()` callouts.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_D_DualDataMechanism")
        self._add_editorial_chrome(slide, s_idx, kicker, title)

        fig_h = 4.55 if bottom_stats else 5.75
        self._add_contain_image(
            slide,
            s_idx,
            "dual_left_data_fig",
            left_image_path,
            0.55,
            1.16,
            5.80,
            fig_h,
            placeholder_title="Quantitative Statistical Proof",
        )

        if right_mechanism_cards:
            n_cards = min(len(right_mechanism_cards), 3)
            card_h = (fig_h - 0.16 * (n_cards - 1)) / max(n_cards, 1)
            for i, cinfo in enumerate(right_mechanism_cards[:3]):
                cy = 1.16 + i * (card_h + 0.16)
                card = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(6.60), Inches(cy), Inches(6.18), Inches(card_h))
                card.fill.solid()
                card.fill.fore_color.rgb = C_BG_PAPER
                card.line.color.rgb = C_RULE
                self._register_box(s_idx, f"mech_card_{i}", 6.60, cy, 6.18, card_h)

                accent_col = self._resolve_color(cinfo.get("color", "blue"))
                badge = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(6.60), Inches(cy), Inches(0.10), Inches(card_h))
                badge.fill.solid()
                badge.fill.fore_color.rgb = accent_col
                badge.line.fill.background()

                head_str = str(cinfo.get("headline", "Mechanistic Axis")).strip()
                det_str = str(cinfo.get("detail", "")).strip()
                head_pt = _auto_scale_headline_pt(head_str, box_width_in=5.75, max_pt=14.5, min_pt=11.5)
                self._add_stacked_textbox(
                    slide=slide,
                    slide_idx=s_idx,
                    label=f"mech_card_text_{i}",
                    x=6.82,
                    y=cy + 0.10,
                    w=5.80,
                    h=card_h - 0.16,
                    headline_text=head_str,
                    caption_text=det_str,
                    headline_pt=head_pt,
                    caption_pt=11.0,
                    headline_color=accent_col,
                    caption_color=C_MUTED,
                    space_between_pt=3.0,
                    register=False,
                )
        else:
            self._add_contain_image(
                slide,
                s_idx,
                "dual_right_mech_fig",
                right_image_path,
                6.60,
                1.16,
                6.18,
                fig_h,
                placeholder_title="BioRender Organelle & Topology Mechanism Schematic",
            )

        if bottom_stats:
            self._add_rule(slide, 0.55, 5.84, 12.23, C_RULE)
            x_positions = [0.60, 4.75, 8.90]
            for i, (num_h, cap, col) in enumerate(bottom_stats[:3]):
                self._add_stat_callout(slide, s_idx, f"d_stat_{i}", x_positions[i], 5.92, num_h, cap, color=col, w=3.80, h=1.12)

        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Layout E: 100% Native Editable Roadmap / Experimental Design ────────
    def add_native_roadmap_slide(
        self,
        title: str,
        kicker: str,
        steps: Sequence[Tuple[str, str, str, str]],
        bottom_takeaway: str = "",
        speaker_notes: str = "",
    ) -> None:
        """
        Layout E: 100% Native Vector Numbered Experimental Design or Wet-Lab Validation Roadmap.
        Each item in `steps` is `(title_bold, method_subtitle, target_readout, color_key)`.
        Uses a unified `TextFrame` per step row so title and subtitle never collide.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_E_NativeRoadmap")
        self._add_editorial_chrome(slide, s_idx, kicker, title)

        y0 = 1.22
        row_step = 1.02
        for idx, (st_title, st_sub, st_readout, col_key) in enumerate(steps[:5], start=1):
            col_rgb = self._resolve_color(col_key)
            badge = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.60), Inches(y0 + 0.08), Inches(0.46), Inches(0.46))
            badge.fill.solid()
            badge.fill.fore_color.rgb = col_rgb
            badge.line.fill.background()
            self._register_box(s_idx, f"roadmap_badge_{idx}", 0.60, y0 + 0.08, 0.46, 0.46)
            btf = badge.text_frame
            bp = btf.paragraphs[0]
            bp.text = str(idx)
            bp.font.name = "Arial"
            bp.font.size = Pt(14)
            bp.font.bold = True
            bp.font.color.rgb = C_BG_WHITE
            bp.alignment = PP_ALIGN.CENTER

            st_head_pt = _auto_scale_headline_pt(st_title, box_width_in=7.65, max_pt=14.5, min_pt=11.5)
            self._add_stacked_textbox(
                slide=slide,
                slide_idx=s_idx,
                label=f"rm_row_text_{idx}",
                x=1.25,
                y=y0,
                w=7.65,
                h=0.82,
                headline_text=st_title,
                caption_text=st_sub,
                headline_pt=st_head_pt,
                caption_pt=11.0,
                headline_color=C_INK,
                caption_color=C_MUTED,
                space_between_pt=2.0,
                register=True,
            )
            readout_pt = _auto_scale_headline_pt(st_readout, box_width_in=3.65, max_pt=12.5, min_pt=10.0)
            self._add_textbox(
                slide,
                s_idx,
                f"rm_readout_{idx}",
                9.10,
                y0 + 0.14,
                3.65,
                0.42,
                st_readout,
                size=readout_pt,
                bold=True,
                color=col_rgb,
                align=PP_ALIGN.RIGHT,
            )
            if idx < min(len(steps), 5):
                self._add_rule(slide, 1.25, y0 + 0.90, 11.50, C_RULE)
            y0 += row_step

        if bottom_takeaway:
            self._add_textbox(
                slide,
                s_idx,
                "roadmap_takeaway",
                0.60,
                6.42,
                12.15,
                0.48,
                bottom_takeaway,
                size=11.5,
                bold=True,
                color=C_EMERALD,
            )

        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Layout F: 100% Native Editable Evidence Tier & Tractability Table ───
    def add_native_evidence_matrix_slide(
        self,
        title: str,
        kicker: str,
        headers: Sequence[str],
        rows: Sequence[Sequence[str]],
        col_widths: Optional[Sequence[float]] = None,
        bottom_stats: Sequence[Tuple[str, str, str]] = (),
        speaker_notes: str = "",
    ) -> None:
        """
        Layout F: 100% Native Editable PowerPoint Table for Study Health Scorecard or
        Multi-Layer Target Evidence & Tractability Matrix, paired with 3 bottom `stat()` callouts.
        """
        slide = self.prs.slides.add_slide(self.blank_layout)
        s_idx = len(self.prs.slides)
        self.slide_layouts_used.append("Layout_F_NativeTable")
        self._add_editorial_chrome(slide, s_idx, kicker, title)

        n_rows = len(rows) + 1
        n_cols = len(headers)
        max_table_h = 4.45 if bottom_stats else 5.60
        top, left, width = 1.22, 0.55, 12.23
        height = min(max_table_h, 0.46 * n_rows)

        table_shape = slide.shapes.add_table(n_rows, n_cols, Inches(left), Inches(top), Inches(width), Inches(height))
        self._register_box(s_idx, "native_matrix_table", left, top, width, height)
        table = table_shape.table

        if col_widths and len(col_widths) == n_cols:
            for c_i, w in enumerate(col_widths):
                table.columns[c_i].width = Inches(w)

        for c_i, h_text in enumerate(headers):
            cell = table.cell(0, c_i)
            cell.fill.solid()
            cell.fill.fore_color.rgb = C_CARD_BG
            p = cell.text_frame.paragraphs[0]
            p.text = str(h_text)
            p.font.name = "Arial"
            p.font.size = Pt(10.0)
            p.font.bold = True
            p.font.color.rgb = C_INK

        for r_i, row_vals in enumerate(rows, start=1):
            for c_i, val in enumerate(row_vals):
                sval = str(val)
                cell = table.cell(r_i, c_i)
                cell.fill.solid()
                cell.fill.fore_color.rgb = C_BG_PAPER if r_i % 2 == 1 else C_BG_WHITE
                p = cell.text_frame.paragraphs[0]
                p.text = sval
                p.font.name = "Arial"
                p.font.size = Pt(9.5)
                p.font.bold = c_i == 0 or ("E1" in sval or "E2" in sval or "E3" in sval)
                if "E1" in sval:
                    p.font.color.rgb = C_EMERALD
                elif "E2" in sval:
                    p.font.color.rgb = C_SLATE_BLUE
                elif "E3" in sval:
                    p.font.color.rgb = C_AMBER
                else:
                    p.font.color.rgb = C_INK if c_i == 0 else C_MUTED

        if bottom_stats:
            self._add_rule(slide, 0.55, 5.84, 12.23, C_RULE)
            x_positions = [0.60, 4.75, 8.90]
            for i, (num_h, cap, col) in enumerate(bottom_stats[:3]):
                self._add_stat_callout(slide, s_idx, f"f_stat_{i}", x_positions[i], 5.92, num_h, cap, color=col, w=3.80, h=1.12)

        self._attach_notes(slide, s_idx, kicker, title, speaker_notes)

    # ── Backward-Compatible Wrappers (Upgraded to New Editorial Layouts) ────
    def add_executive_scorecard_slide(
        self,
        scorecard_rows: Sequence[Dict[str, Any]],
        speaker_notes: str = "",
    ) -> None:
        """Upgraded Executive Study Health & Empirical Null FDR Scorecard using Layout F."""
        headers = [
            "Contrast / Arm",
            "Design (N)",
            "Target Audit",
            "Noise (σ)",
            "Obs vs Null",
            "Emp. FDR",
            "Max Tier",
            "Calibrated Takeaway",
        ]
        rows = [
            [
                str(r.get("contrast", "")),
                str(r.get("design_n", "")),
                str(r.get("target_status", "")),
                str(r.get("sigma_hat", "")),
                str(r.get("obs_vs_null", "")),
                str(r.get("emp_fdr", "")),
                str(r.get("max_grade", "E2")),
                str(r.get("takeaway", "")),
            ]
            for r in scorecard_rows
        ]
        col_widths = [1.80, 0.95, 1.40, 0.85, 1.25, 0.90, 0.95, 4.13]
        self.add_native_evidence_matrix_slide(
            title="Study Health & Empirical Null FDR Scorecard",
            kicker="STEP 1 · QUALITY CONTROL & CALIBRATION",
            headers=headers,
            rows=rows,
            col_widths=col_widths,
            bottom_stats=[
                ("12-Point QC Pass", "Contaminants scrubbed & missingness profiled", "emerald"),
                ("Permutation Null", "Empirical FDR calibrated before DE testing", "blue"),
                ("Evidence Capped", "Claims strictly bounded by sample power N", "coral"),
            ],
            speaker_notes=speaker_notes,
        )

    def add_figure_and_narrative_slide(
        self,
        title: str,
        kicker: str,
        bullets: Sequence[str],
        image_path: Optional[str | Path] = None,
        speaker_notes: str = "",
    ) -> None:
        """
        Upgraded legacy wrapper: automatically converts legacy `bullets` into collision-proof
        unboxed `stat()` callouts on Layout C (never renders a gray bullet-wall box!).
        """
        palette_cycle = ["coral", "blue", "emerald", "plum"]
        converted_stats: List[Tuple[str, str, str]] = []
        for i, b in enumerate(bullets[:4]):
            clean_b = str(b).lstrip("•-* ").strip()
            if ":" in clean_b:
                head, tail = clean_b.split(":", 1)
                converted_stats.append((head.strip(), tail.strip(), palette_cycle[i % 4]))
            else:
                words = clean_b.split()
                head = " ".join(words[:3])
                tail = " ".join(words[3:]) or clean_b
                converted_stats.append((head, tail, palette_cycle[i % 4]))

        self.add_left_figure_right_stats_slide(
            title=title,
            kicker=kicker,
            stats=converted_stats,
            image_path=image_path,
            speaker_notes=speaker_notes,
        )

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
        bottom_stats = [(bottom_callout[:24], bottom_callout[24:85] or "Key takeaway", "emerald")] if bottom_callout else []
        self.add_native_evidence_matrix_slide(
            title=title,
            kicker=kicker,
            headers=headers,
            rows=rows,
            col_widths=col_widths,
            bottom_stats=bottom_stats,
            speaker_notes=speaker_notes,
        )

    # ── Automated Geometry & Design Constraint Audits ───────────────────────
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

    def validate_deck_design_constraints(
        self,
        enforce_comprehensive: Optional[bool] = None,
        min_slides: Optional[int] = None,
    ) -> List[str]:
        """
        Enforce strict editorial & anti-laziness PPTX constraints:
        1. Action Title length <= 54 chars.
        2. `stat()` headline <= 28 chars & caption <= 72 chars.
        3. Zero bullet symbols (`•`) on slide canvas text frames.
        4. Layout variety: no single layout archetype repeated > 2 times consecutively.
        5. Speaker notes free of ceremonial greetings ('各位老师同学好' / 'Hello everyone').
        6. NO_DUPLICATE_IMAGES_ACROSS_SLIDES: Every figure/schematic image path may appear on at most ONE slide.
        7. MIN_COMPREHENSIVE_SLIDE_COUNT: When `enforce_comprehensive` is enabled, total slide count must be >= `min_slides` (default 20).
        """
        violations: List[str] = []
        check_comp = self.enforce_comprehensive if enforce_comprehensive is None else enforce_comprehensive
        min_s = self.min_comprehensive_slides if min_slides is None else min_slides

        for idx, t in enumerate(self.slide_titles, start=1):
            if len(t) > 54:
                violations.append(f"Slide {idx}: Action Title exceeds 54 chars ({len(t)} chars): '{t}'")

        for s_idx, num_h, cap in self.slide_stat_labels:
            if len(num_h) > 28:
                violations.append(f"Slide {s_idx}: stat() headline exceeds 28 chars ({len(num_h)} chars): '{num_h}'")
            if len(cap) > 72:
                violations.append(f"Slide {s_idx}: stat() caption exceeds 72 chars ({len(cap)} chars): '{cap}'")

        for s_idx, slide in enumerate(self.prs.slides, start=1):
            for shape in slide.shapes:
                if getattr(shape, "has_text_frame", False):
                    txt = shape.text_frame.text or ""
                    if "•" in txt:
                        violations.append(f"Slide {s_idx}: Forbidden bullet symbol ('•') detected on slide canvas.")

        for i in range(len(self.slide_layouts_used) - 2):
            if (
                self.slide_layouts_used[i]
                == self.slide_layouts_used[i + 1]
                == self.slide_layouts_used[i + 2]
            ):
                violations.append(
                    f"Slides {i+1}-{i+3}: Layout '{self.slide_layouts_used[i]}' repeated 3 times consecutively."
                )

        for entry in self.bilingual_notes_log:
            low_notes = entry["notes"].lower()
            for greet in FORBIDDEN_GREETINGS:
                if greet in low_notes:
                    violations.append(
                        f"Slide {entry['slide_idx']}: Ceremonial greeting '{greet}' found in speaker notes."
                    )

        # Gate 6: NO_DUPLICATE_IMAGES_ACROSS_SLIDES
        seen_images: Dict[str, int] = {}
        for s_idx, lbl, img_p in self.slide_images_used:
            if img_p in seen_images and seen_images[img_p] != s_idx:
                violations.append(
                    f"Slides {seen_images[img_p]} and {s_idx}: Duplicate image '{Path(img_p).name}' reused across slides (NO_DUPLICATE_IMAGES_ACROSS_SLIDES violation)."
                )
            else:
                seen_images[img_p] = s_idx

        # Gate 7: MIN_COMPREHENSIVE_SLIDE_COUNT
        if check_comp and len(self.prs.slides) < min_s:
            violations.append(
                f"Deck has only {len(self.prs.slides)} slides, below minimum comprehensive threshold of {min_s} slides (MIN_COMPREHENSIVE_SLIDE_COUNT violation)."
            )

        # Gate 8: BIOLOGIST_PLAIN_LANGUAGE_LINT (no cryptic standalone statistical jargon in slide titles or stat() headlines)
        for idx, t in enumerate(self.slide_titles, start=1):
            low_t = t.lower()
            for jargon in BANNED_OPAQUE_JARGON_HEADLINES:
                if jargon in low_t:
                    violations.append(
                        f"Slide {idx}: Opaque bioinformatic jargon '{jargon}' in Action Title '{t}'. Translate to plain biological language!"
                    )
        for s_idx, num_h, _cap in self.slide_stat_labels:
            low_h = num_h.lower()
            for jargon in BANNED_OPAQUE_JARGON_HEADLINES:
                if jargon in low_h:
                    violations.append(
                        f"Slide {s_idx}: Opaque bioinformatic jargon '{jargon}' in stat() headline '{num_h}'. Translate to plain biological language!"
                    )

        # Gate 9: PENULTIMATE_BIOLOGIST_SUMMARY_SLIDE (Comprehensive decks must include s26_biorender_single_kd_summary.png at penultimate slide)
        if check_comp:
            s26_slides = [s_idx for s_idx, _, img_p in self.slide_images_used if "s26_biorender_single_kd_summary" in Path(img_p).name]
            if not s26_slides:
                violations.append(
                    "Comprehensive deck is missing the mandatory penultimate Single-Perturbation Hero Genes, Pathways, Literature & Biological Meaning slide ('s26_biorender_single_kd_summary.png')."
                )
            elif len(self.prs.slides) >= 3 and (len(self.prs.slides) - 1) not in s26_slides:
                violations.append(
                    f"Slide 's26_biorender_single_kd_summary.png' is at slide {s26_slides[0]}, expected penultimate slide {len(self.prs.slides) - 1}."
                )

        return violations

    def export_bilingual_notes_markdown(self, md_path: str | Path) -> Path:
        """Export companion `*_Notes_Bilingual.md` for self-study and presentation."""
        out = Path(md_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            f"# {self.title} — 双语深度解析与汇报要点 (Bilingual Slide Walkthrough)",
            "",
            f"> **{self.subtitle}**",
            "> 本讲义与 PPT 页码一一对应，直接提炼每页的核心数据逻辑、统计学依据与生物学机制（省去客套开场白，兼顾个人精读复盘与组会直接汇报）。",
            "",
        ]
        for item in self.bilingual_notes_log:
            lines.append(f"## Slide {item['slide_idx']:02d} · {item['title']}")
            lines.append(f"*Section: {item['kicker']}*")
            lines.append("")
            lines.append(item["notes"])
            lines.append("")
            lines.append("---")
            lines.append("")
        out.write_text("\n".join(lines), encoding="utf-8")
        return out

    def save(self, output_path: str | Path, export_notes_md: bool = True) -> Tuple[Path, List[str]]:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        warnings = self.validate_deck_geometry() + self.validate_deck_design_constraints()
        self.prs.save(str(out))
        if export_notes_md and self.bilingual_notes_log:
            md_out = out.with_name(out.stem + "_Notes_Bilingual.md")
            self.export_bilingual_notes_markdown(md_out)
        return out, warnings


def build_reference_16_slide_blueprint(output_path: str | Path, title: str) -> Tuple[Path, List[str]]:
    """
    Build the canonical ProteoSkill synthetic demo presentation demonstrating all 6 layouts (A-F)
    with zero bullet walls, collision-proof unified TextFrames, and direct bilingual notes.
    """
    deck = EditorialSlideDeck(
        title=title[:48],
        subtitle="[SYNTHETIC DEMO BLUEPRINT] 7-Step Calibrated Proteomics Discovery, Systems Topology & Clinical Translation",
        total_expected_slides=6,
    )

    # Slide 1: Layout A (Cover)
    deck.add_cover_slide(
        metadata_badges=[
            ("6,420 Quantified Proteins", "Synthetic Layout Showcase"),
            ("Emp. FDR < 0.01", "Shared-Control Permutation Null"),
            ("125 Dual-Perturbed Hits", "Cross-Layer Bridge Verified"),
            ("6 Editorial Layouts", "Demo Blueprint Architecture"),
        ],
        speaker_notes=format_bilingual_notes(
            [
                "【演示模板说明】本 PPT 为展示 6 种版式（Layout A-F）与无碰撞双语备注系统的合成演示蓝图（Synthetic Layout Demo）。",
                "全流程在差异统计前引入共享对照置换零分布校准，确保所有入选的核心靶点具备严格的假阳性率控制。",
            ],
            [
                "[SYNTHETIC DEMO BLUEPRINT] This deck demonstrates the 6 editorial layouts (Layouts A-F) and collision-proof bilingual notes architecture.",
                "Shared-control permutation null calibration is enforced prior to differential testing to bound false discovery rates.",
            ],
        ),
    )

    # Slide 2: Layout E (Native Experimental Design & 7-Step Architecture)
    deck.add_native_roadmap_slide(
        title="Experimental Design & 7-Step Analysis Pipeline",
        kicker="ACT I · STUDY DESIGN & ARCHITECTURE",
        steps=[
            ("Step 1 · QC & Empirical Null Calibration", "Remove cRAP/keratins, audit missingness, calibrate permutation FDR", "6,420 Clean Genes", "emerald"),
            ("Step 2 · Differential & TOST Equivalence", "Peptide-stratified empirical Bayes (DEqMS) + TOST baseline proof", "142 Regulated Hits", "coral"),
            ("Step 3-4 · Complexes, Organelles & Regulators", "CORUM stoichiometry, ectodomain shedding, upstream kinase/TF ULM", "3 Core Modules", "blue"),
            ("Step 5-6 · Multi-Omics & CPTAC Clinical Cohorts", "RNA-protein decoupling, CPTAC tumor validation, ROC-AUC & ChEMBL", "AUC = 0.91 Verified", "plum"),
            ("Step 7 · BioRender Model & Wet-Lab Validation", "Evidence-graded vector schematic (E1/E2 vs E3) & orthogonal assays", "5 Follow-Up Assays", "amber"),
        ],
        bottom_takeaway="Design Principle: Every statistical claim is bounded by sample power N and cross-checked by Fact Gate.",
        speaker_notes=format_bilingual_notes(
            [
                "实验设计涵盖质控校准、分层贝叶斯差异分析、大分子复合物与跨膜拓扑解构、以及外部临床队列交叉验证。",
                "所有后续页面严格区分质谱直接证据（E1/E2）与计算推断假说（E3）。",
            ],
            [
                "The pipeline spans forensic QC, peptide-stratified empirical Bayes, macromolecular complex topology, and clinical cohort validation.",
                "Subsequent slides strictly distinguish direct MS evidence (E1/E2) from computationally inferred hypotheses (E3).",
            ],
        ),
    )

    # Slide 3: Layout B (Wide QC Figure Top + 3 Horizontal Stats Bottom)
    deck.add_top_figure_bottom_stats_slide(
        title="Rigorous Filtering Yields 6,420 Clean Proteins",
        kicker="STEP 1 · QUALITY CONTROL & SANITIZATION",
        stats=[
            ("21/22 Tubes Passed QC", "Single outlier tube (−19.8% IDs) quarantined", "emerald"),
            ("σ_null: 0.048 to 0.171", "Noise scales cleanly with perturbation mode", "blue"),
            ("Log2 Scale Auto-Proven", "Zero double-log transformation errors", "coral"),
        ],
        speaker_notes=format_bilingual_notes(
            [
                "第一步移除外源角蛋白与血清白蛋白污染并隔离单支异常管，保留 6,420 个唯一标准基因符号。",
                "生物学重复内部相关性优异，且自动核验了对数尺度。",
            ],
            [
                "Step 1 removes exogenous keratin and cRAP contaminants and quarantines 1 outlier tube, retaining 6,420 unique canonical gene symbols.",
                "Intra-group biological replicate reproducibility and log2 scale integrity are automatically verified.",
            ],
        ),
    )

    # Slide 4: Layout C (Large Volcano Left + 4 Vertical Stats Right)
    deck.add_left_figure_right_stats_slide(
        title="Stratified DEqMS Identifies 142 Core Effectors",
        kicker="STEP 2 · QUANTITATIVE DIFFERENTIAL REMODELING",
        stats=[
            ("86 Upregulated Hits", "Significant induction (FDR < 0.01, |log2FC| > 1)", "coral"),
            ("56 Downregulated Hits", "Coordinated depletion of metabolic enzymes", "blue"),
            ("4,120 TOST Stable", "Equivalence test confirms unperturbed baseline", "emerald"),
            ("Emp. FDR = 0.008", "Verified against shared-control null distribution", "plum"),
        ],
        speaker_notes=format_bilingual_notes(
            [
                "采用按肽段数分层的经验贝叶斯模型（DEqMS），共鉴定出 86 个显著上调蛋白和 56 个显著下调蛋白。",
                "同时通过双单侧等价性检验（TOST）严格证明 4,120 个背景蛋白未发生显著偏移。",
            ],
            [
                "Peptide-stratified empirical Bayes (DEqMS) identifies 86 upregulated and 56 downregulated core effectors.",
                "Two One-Sided Equivalence Testing (TOST) confirms 4,120 housekeeping proteins remain stably invariant.",
            ],
        ),
    )

    # Slide 5: Layout D (Dual Quantitative Plot Left + Native Mechanism Cards Right)
    deck.add_dual_data_and_mechanism_slide(
        title="Organelle Redistribution & Paralog Compensation",
        kicker="STEP 3 & 4 · SUBCELLULAR TOPOLOGY & COMPLEXES",
        right_mechanism_cards=[
            {
                "headline": "1. Secretory & Vesicle Cargo Enrichment (OR = 4.8x)",
                "detail": "Fisher exact test & Compartment Mapping confirm directional sorting toward extracellular vesicles.",
                "color": "coral",
            },
            {
                "headline": "2. Paralog Family Buffering (+1.82 log2FC)",
                "detail": "Primary target depletion triggers selective upregulation of its closest sequence paralog.",
                "color": "blue",
            },
            {
                "headline": "3. Ectodomain Shedding Verified by Peptide Topology",
                "detail": "100% of detected peptides map N-terminal to the transmembrane helix (residues 28-185).",
                "color": "emerald",
            },
        ],
        bottom_stats=[
            ("OR = 4.82 (p = 2e-7)", "Vesicle & endosomal compartment enrichment", "coral"),
            ("Cohesion Z = −3.4", "Coordinated complex stoichiometry shift", "blue"),
            ("N-Term Shedding", "Transmembrane domain boundary verified", "emerald"),
        ],
        speaker_notes=format_bilingual_notes(
            [
                "左侧定量富集数据与右侧细胞生物学机制相互印证：差异蛋白显著富集于囊泡运输通路（OR = 4.82）。",
                "同源家族扫描发现主靶点下调后其同源蛋白发生代偿性上调，且跨膜肽段拓扑证实了外域切割脱落。",
            ],
            [
                "Left quantitative enrichment pairs directly with right mechanistic modules: regulated hits enrich in vesicle sorting (OR = 4.82).",
                "Paralog scanning reveals compensatory upregulation, while peptide coordinate mapping confirms ectodomain shedding.",
            ],
        ),
    )

    # Slide 6: Layout F (Native Editable Target Evidence & Tractability Matrix)
    deck.add_native_evidence_matrix_slide(
        title="Multi-Layer Target Evidence & Clinical Matrix",
        kicker="STEP 5 & 6 · CLINICAL COHORTS & DRUG TRACTABILITY",
        headers=["Gene", "MS log2FC", "Peptides", "RNA-Prot Mode", "CPTAC Cohort", "ROC AUC", "Evidence Tier", "Tractability"],
        rows=[
            ["EGFR", "+2.14", "18 pep", "Concordant Up", "Up (p < 1e-6)", "0.93", "E1 Confirmatory", "Approved TKI / ADC"],
            ["KRAS", "+1.65", "9 pep", "Post-Trans Stab.", "Up (p = 2e-4)", "0.88", "E1 Confirmatory", "Phase III Covalent"],
            ["RAB27A", "+1.48", "11 pep", "Concordant Up", "Up (p = 8e-4)", "0.85", "E2 Moderated", "Tool Compound"],
            ["MAPK3", "+0.92", "14 pep", "Paralog Buffer", "Stable", "0.74", "E2 Moderated", "Selective Inhibitor"],
            ["CAND_X", "+1.88", "1 pep (Exp)", "Decoupled", "Pending", "0.79", "E3 Exploratory", "Novel Surface Antigen"],
        ],
        col_widths=[1.15, 1.05, 1.15, 1.65, 1.55, 1.00, 1.65, 3.03],
        bottom_stats=[
            ("CPTAC Concordance", "84% directional agreement in patient tumors", "emerald"),
            ("Panel AUC = 0.93", "5-gene signature stratified with 1000x bootstrap", "coral"),
            ("2 Approved Targets", "Direct clinical repurposing tractability", "blue"),
        ],
        speaker_notes=format_bilingual_notes(
            [
                "本表为 100% 原生可编辑 PowerPoint 表格，汇总了核心靶点的质谱变化倍数、肽段数、RNA-蛋白耦合模式、CPTAC 临床队列验证及成药性。",
                "单肽段蛋白（CAND_X）被强制标注为 E3 Exploratory，杜绝过度解读。",
            ],
            [
                "This native editable PowerPoint table integrates MS log2FC, peptide count, RNA-protein coupling, CPTAC validation, and drug tractability.",
                "Single-peptide detections (CAND_X) are strictly capped at E3 Exploratory to prevent overclaiming.",
            ],
        ),
    )

    return deck.save(output_path)


build_demo_layout_showcase = build_reference_16_slide_blueprint


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 7A: Editorial Minimalist PPTX Slide Deck Builder")
    parser.add_argument("--out", "-o", default="proteomics_discovery_deck.pptx", help="Output PPTX path")
    parser.add_argument("--title", default="Proteomics Discovery & Mechanistic Report", help="Deck title")
    args = parser.parse_args()

    out_p, warns = build_reference_16_slide_blueprint(args.out, args.title)
    print(f"[Step 7A] Saved Editorial PPTX to {out_p} (Design/Geometry Warnings: {len(warns)})")
    for w in warns:
        print(f"  [WARN] {w}")


if __name__ == "__main__":
    main()

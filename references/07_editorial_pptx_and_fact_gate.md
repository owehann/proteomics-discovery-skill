<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 07: Biologist-First Editorial PPTX Architecture & Pre-Delivery Fact Gate

> [!IMPORTANT]
> **Stage Scope (`Stage 6 & Stage 7` — Presentation & Verification Gates)**
> This playbook governs [`scripts/08_editorial_pptx_builder.py`](../scripts/08_editorial_pptx_builder.py) and [`scripts/09_fact_and_citation_gate.py`](../scripts/09_fact_and_citation_gate.py), enforcing publication-grade minimalist slide design, **Biologist-First plain-language storytelling**, the mandatory **Penultimate Core Genes & Literature Summary Slide (`s26`)**, and 100% numeric/bibliographic verification before delivery.

---

## 1. Minimalist Editorial 16:9 Slide Deck (`08_editorial_pptx_builder.py`)

### 1.1 Five Non-Negotiable Visual & Storytelling Rules

1. **Airy Pure-White Editorial Header (No Heavy Dark Top Banner)**:
   - Canvas is widescreen `16:9` (`13.333 in × 7.500 in`) on pure white (`#FFFFFF`) or warm paper (`#FAFAF7`). **Never place a heavy dark navy rectangle across the top of content slides.**
   - Muted uppercase section kicker at `y = 0.24"` (`10.5pt Arial Bold`, `#8A8984`), e.g., `04 / ACT I · SAMPLE QUALITY CONTROL & SANITIZATION`.
   - Single-line concise **Action Title** at `y = 0.50"` (`20–23pt Arial Bold`, `#1F1F1D`, **strictly `<= 54 characters`** so it never wraps to 2 lines or collides with the divider rule).
   - Ultra-thin hairline divider (`0.012"` height, `#D9D8D3`) at `y = 1.02"` and subtle bottom footer at `y = 7.14"`.
2. **Zero Bullet-Wall Anti-Pattern & Collision-Proof Single-`TextFrame` `stat()` Callouts**:
   - **Never place bulleted paragraph lists (`•`) or giant gray rounded text boxes on the slide face.**
   - Always use `_add_stat_callout()` / `_add_stacked_textbox()`, which places the bold colored metric (`Paragraph 0`, auto-scaled between `13.5pt` and `21.5pt` via `_auto_scale_headline_pt()` to fit on a single line) and the muted caption (`Paragraph 1`, `10–11pt Arial`, `#5A5955`) inside **one unified PowerPoint `TextFrame`**, preventing headline-caption collisions.
3. **Biologist-First Plain-Language Storytelling (`BIOLOGIST_PLAIN_LANGUAGE_LINT`)**:
   - **Never use opaque bioinformatic jargon** (`"Kolde Beta Order-Statistic"`, `"Dirichlet Stability"`, `"Empirical Null Integral"`, `"Archetype A4"`) in slide titles or `stat()` headlines.
   - Every slide's `stat()` callouts and bilingual speaker notes (`format_bilingual_notes()`) must clearly answer **4 biologist questions**:
     1. **Why we did this step** (the biological or technical question).
     2. **How to read the figure & what stands out** (visual reading guide: axes, colors, and highlighted genes).
     3. **Which specific core genes & pathways changed and their biological meaning**.
     4. **Why this discovery motivates the next slide**.
4. **Mandatory Penultimate Slide — Core Changed Genes, Pathways, Literature & Biological Meaning (`s26_biorender_single_kd_summary.png`)**:
   - In every comprehensive presentation (`>= 20` slides), the **second-to-last slide (`Slide N-1`)** right before the final Wet-Lab Validation Roadmap must display `s26_biorender_single_kd_summary.png`, summarizing each perturbation's **Top Upregulated & Downregulated Hero Genes (with exact fold-changes)**, **Core Dysregulated Pathways**, **Landmark Literature Citations (`Author et al., Journal Year, PMID`)**, and **Plain-Language Biological Meaning**.
5. **Direct Bilingual Notes (Zero Ceremonial Greetings & Zero Duplicate Images)**:
   - Every image embedded across the deck must be unique (`NO_DUPLICATE_IMAGES_ACROSS_SLIDES`).
   - Detailed walkthrough notes are written into both `slide.notes_slide` and `*_Notes_Bilingual.md` (`【中文深度导读：为什么要做这一步 + 图表怎么看 + 核心生物学发现】` + `[English Visual & Biological Walkthrough]`).
   - **Never start notes with ceremonial greetings** (`"各位老师同学好"`, `"大家好"`, `"Hello everyone"`).

---

### 1.2 Six Alternating Layout Archetypes (`Layouts A–F`)

**Anti-Monotony Rule**: Never repeat the same layout archetype more than **2 times consecutively** (`validate_deck_design_constraints()` automatically audits this).

| Layout ID | Method in `EditorialSlideDeck` | Visual Structure & Ideal Use Case |
| :--- | :--- | :--- |
| **Layout A** | `add_cover_slide()` / `add_full_canvas_slide()` | **Cover & Full-Bleed Visual Canvas (`12.23" × 5.82"`)**: Used for P1 Cover + Graphical Abstract, Multi-Panel Top-Hit Spotlights, P25 Working Model, and P26 Penultimate Biologist Literature Cheat Sheet. |
| **Layout B** | `add_top_figure_bottom_stats_slide()` | **Wide Figure Top (`12.23" × 4.68"`) + 3 Horizontal `stat()` Callouts Bottom**: Used for 4-panel sample QC overviews, Multi-condition rescue grids, and RNA–Protein coupling quadrants. |
| **Layout C** | `add_left_figure_right_stats_slide()` | **Large Figure Left (`8.55" × 5.82"`) + 4 Stacked `stat()` Callouts Right (`x = 9.35"`)**: Used for PCA/Null calibration, DEqMS Volcano + TOST equivalence, Jaccard pathways, and Upstream Kinase/TF inference. |
| **Layout D** | `add_dual_data_and_mechanism_slide()` | **Dual Side-by-Side Split (`5.80"` Plot Left + `6.18"` BioRender Schematic or Native Cards Right + 3 Bottom `stat()` Callouts)**: Pairs statistical enrichment directly with subcellular/topological mechanism on the same slide. |
| **Layout E** | `add_native_roadmap_slide()` | **100% Native Vector Numbered Workflow / Wet-Lab Validation Roadmap**: 5 numbered colored oval badges (`1`–`5`), bold assay title + method subtitle on left, expected readout on right. |
| **Layout F** | `add_native_evidence_matrix_slide()` | **100% Native Editable PowerPoint Evidence & Tractability Matrix Table + 3 Bottom `stat()` Callouts**: Color-coded `E1 / E2 / E3` evidence tiers, peptide counts, cohort concordance, and ChEMBL drug status. |

---

## 2. Pre-Delivery Fact & Citation Gate (`09_fact_and_citation_gate.py`)

No report (`.md`), schematic (`.svg`), or presentation (`.pptx`) may be finalized until `09_fact_and_citation_gate.py` and `validate_deck_design_constraints()` return **`PASS (0 Critical Errors)`**.

### 2.1 Single Source of Truth (`facts.json`)
`09_fact_and_citation_gate.py` compiles a canonical `facts.json` dictionary directly from the processed CSV/TSV matrices and QC outputs:
- **Gene-Level Facts**: `accession`, `unique_peptides`, `psm_count`, raw/normalized sample values (excluding quarantined outlier tubes), group means, $\Delta\log_2\mathrm{FC}$, $p$-value, BH $q$-value, `TOST` $p$-value, quadrant/sector label, and evidence grade (`E1/E2/E3`).
- **Contrast-Level Facts**: `n_tested`, `sigma_null`, `expected_null_hits`, `observed_hits`, `empirical_null_fdr`, and `max_evidence_grade`.

### 2.2 Nine Automated Gate Checks

| Gate ID | Verification Check | Failure Condition |
| :--- | :--- | :--- |
| **Check 1** | **Gene–Number Pair & Quadrant Reconciliation** | Any numeric value or quadrant label in `.md`/`.svg`/`.pptx` deviates by `> 0.02` from `facts.json` or quotes a quarantined outlier sample. |
| **Check 2** | **Single-Peptide Promotion Audit** | A `unique_peptides == 1` protein appears in a Core Tier table or headline without an explicit `[1-pep]` / `exploratory` badge. |
| **Check 3** | **Overclaim & Evidence-Grade Linting** | Causal verbs (`"proves"`, `"confirms"`, `"demonstrates definitively"`) applied to contrasts with `empirical_null_fdr > 0.25` or `E3` grade. |
| **Check 4** | **PMID / DOI Bibliographic Verification** | Any cited `PMID` or `DOI` fails live metadata resolution against NCBI E-utilities / EuropePMC. |
| **Check 5** | **Geometry & Typography Linting** | Any bounding-box overlap, Action Title `> 54 chars`, `stat()` caption `> 90 chars`, `•` bullet character on slide face, or layout repeated `> 2` times consecutively. |
| **Check 6** | **No Ceremonial Greetings Lint** | Speaker notes contain ceremonial opening phrases (`"各位老师同学好"`, `"大家好"`, `"Hello everyone"`). |
| **Check 7** | **No Duplicate Images Across Slides** | The same `.png` figure file is embedded on more than one slide (`NO_DUPLICATE_IMAGES_ACROSS_SLIDES`). |
| **Check 8** | **Biologist Plain-Language Lint** | Slide titles or `stat()` headlines contain unexplained statistical jargon (`BIOLOGIST_PLAIN_LANGUAGE_LINT`). |
| **Check 9** | **Penultimate Biologist Summary Slide Gate** | Comprehensive deck (`>= 20` slides) is missing `s26_biorender_single_kd_summary` on the second-to-last slide (`PENULTIMATE_BIOLOGIST_SUMMARY_SLIDE`). |

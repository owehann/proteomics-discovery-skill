# Reference Playbook 07: Editorial PPTX Deck Architecture & Pre-Delivery Fact/Citation Gate

This playbook governs **Stage 7 (`scripts/08_editorial_pptx_builder.py` and `scripts/09_fact_and_citation_gate.py`)**, ensuring zero-overlap editorial presentations and 100% numeric and bibliographic integrity before delivery.

---

## 1. Minimalist Editorial 16:9 Slide Deck (`08_editorial_pptx_builder.py`)

### 1.1 Slide Layout & Design Tokens
- **Dimensions**: Widescreen 16:9 (`13.333 in × 7.500 in`).
- **Editorial Palette**:
  - Background: Warm Alabaster `#F8FAFC` / Pure White `#FFFFFF` cards with `1 pt` subtle border `#E2E8F0`.
  - Header Band: Deep Navy Ink `#0F172A` with Slate subtitle `#475569` and Crimson accent bar `#E64B35`.
  - Evidence Tier Badges:
    - `E1 Confirmatory`: Emerald `#059669` pill
    - `E2 Moderated Hit`: Royal Blue `#2563EB` pill
    - `E3 Exploratory`: Amber `#D97706` pill
- **Zero-Overlap Bounding Box Validator**: Every text frame, KPI card, table, and image container is checked for pairwise bounding-box intersection (`check_slide_overlaps()`). Any overlap $> 0.01\text{ in}^2$ raises a layout warning and auto-adjusts vertical spacing.

### 1.2 Mandatory Slide Sequence (10–14 Slides)
1. **Slide 1 — Title & Study Architecture**: Project title, modality (`M1–M6`), archetype (`A1–A10`), sample cohort structure, and date.
2. **Slide 2 — Mandatory Executive Scorecard (One-Slide Study Health Summary)**:
   - Displays a structured table of every experimental contrast with:
     - **Sample Size ($N$)** & **Target Verification Status** (Detected / `NaN` / Peptide count)
     - **Background Noise ($\hat{\sigma}$)** & **Expected Null Hits ($E_{\text{null}}$)** vs **Observed Hits ($N_{\text{obs}}$)**
     - **Empirical Null FDR ($\widehat{\text{FDR}}_{\text{emp}}$)**
     - **Highest Evidence Grade Allowed (`E1` / `E2` / `E3`)**
     - **Headline Biological Conclusion**
3. **Slide 3 — Forensic QC & Normalization Diagnostics**: 6-panel QC figure + outlier tube audit + MISEV2023 / cRAP summary.
4. **Slides 4–7 — Core Quantitative & Archetype Discoveries**: Paired left-side narrative cards (exact numbers from `facts.json`) + right-side high-resolution NPG figures.
5. **Slide 8 — Systems Biology, CORUM Complexes, & Network Topology**: Pathway GSEA/ORA (detected background), complex coherence, and upstream kinase/TF/biophysical enrichment.
6. **Slide 9 — BioRender-Style Mechanistic Working Model**: Embedded `300 DPI` rendering of `fig_mechanism_biorender.png` with solid (`E1/E2`) vs dashed (`E3`) arrows explained.
7. **Slide 10 — Multi-Cohort Clinical Validation & Literature Synthesis**: CPTAC/TCGA concordance, Orthogonalized `RRA` top drivers, and balanced clinical trial context (including negative trials).
8. **Slide 11 — Actionable Wet-Lab Validation Roadmap**: Concrete, prioritized follow-up assays (`WB`, `PRM/MRM`, `Co-IP`, `CETSA`, `Site-Directed Mutagenesis`) required to upgrade `E2/E3` candidates to `E1`.
9. **Bilingual Speaker Notes**: Every slide includes detailed presenter notes in English, Chinese, or Bilingual mode as configured in `analysis_spec.yaml`.

---

## 2. Pre-Delivery Fact & Citation Gate (`09_fact_and_citation_gate.py`)

No report (`.md`), schematic (`.svg`), or presentation (`.pptx`) may be finalized until `09_fact_and_citation_gate.py` returns **`PASS (0 Critical Errors)`**.

### 2.1 Single Source of Truth (`facts.json`)
`09_fact_and_citation_gate.py` compiles a canonical `facts.json` dictionary directly from the processed CSV/TSV matrices and QC outputs:
- For every gene $g$: `accession`, `unique_peptides`, `psm_count`, raw/normalized sample values (excluding quarantined outlier tubes!), group means, $\Delta\log_2\text{FC}$ (both uncentered and patient-centered where applicable), $p$-value, BH $q$-value, `TOST` $p$-value, quadrant/sector label, and evidence grade (`E1/E2/E3`).
- For every contrast $c$: `n_tested`, `sigma_null`, `expected_null_hits`, `observed_hits`, `empirical_null_fdr`, `max_evidence_grade`.

### 2.2 Four Automated Gate Checks
1. **Check 1 — Gene-Number Pair & Quadrant Reconciliation**:
   - Scans `.md`, `.svg`, and `.pptx` text for gene-adjacent numbers (e.g., `CARGO_X (18.42)`, `EGFR (+3.12)`, `MARKER_A (+2.40)`) and quadrant/category labels (e.g., `Q2_hyper_secreted` vs `Q1_concordant_up`).
   - Flags any discrepancy $> 0.02$ against `facts.json` or any value accidentally pulled from a quarantined outlier tube (such as `CTRL_Rep2`).
2. **Check 2 — Single-Peptide Promotion Audit**:
   - Scans all headline tables, Executive Summaries, and SVG nodes for proteins with `unique_peptides == 1` in `facts.json`.
   - Fails if a `unique_peptides == 1` protein appears in a Core Tier table without an explicit `1-peptide` / `exploratory` annotation within the same line or paragraph.
3. **Check 3 — Overclaim & Evidence-Grade Linting**:
   - Scans for causal/definitive verbs (`"proves"`, `"confirms"`, `"demonstrates definitively"`, `"drives"`, `"indispensable"`) applied to contrasts where `empirical_null_fdr > 0.25` or `max_evidence_grade == "E3"`.
   - Verifies that any claim of `"preserved"`, `"unperturbed"`, or `"unaffected"` is backed by `tost_p < 0.05`.
4. **Check 4 — PMID / DOI Bibliographic Verification**:
   - Extracts all `PMID: \d+` and `10.\d{4,9}/[-._;()/:A-Za-z0-9]+` strings and verifies their existence and title/author match against NCBI E-utilities / EuropePMC.

# Reference Playbook 06: BioRender-Style Vector SVG Canvas & NPG Statistical Figures

This playbook specifies the visual standards and programmatic primitives implemented in **Stage 6 (`scripts/06_npg_statistical_figures.py` and `scripts/07_biorender_canvas_engine.py`)**.

---

## 1. NPG / Nature Cell Biology Color Palette & Typography Standards

### 1.1 Standard Chromatic Encoding
- **Up-Regulated / Activated / Tumor / OE Gain**: Crimson Vermilion `#E64B35` (Fill tint `#FDE8E4`, stroke `#B8321F`)
- **Down-Regulated / Suppressed / Non-Tumor / KD Loss**: Cerulean Blue `#4DBBD5` (Fill tint `#E3F5FA`, stroke `#1F829C`)
- **Reciprocal Core / Dual-Compartment / Epistatic Shared**: Emerald Teal `#00A087` (Fill tint `#E0F5F1`, stroke `#006B5A`)
- **Specific / Divergent / PTM Switch**: Royal Indigo `#3C5488` (Fill tint `#E6EAF2`, stroke `#233459`)
- **Warning / Hairpin Off-Target / Single-Peptide Exploratory**: Amber Gold `#F39B7F` / `#D97706` (Fill tint `#FEF3C7`)
- **Neutral / Statistically Equivalent (`TOST_Equiv`) / Background**: Warm Slate `#8491B4` / `#94A3B8` (Fill tint `#F1F5F9`)

### 1.2 Publication Plot Rules (`06_npg_statistical_figures.py`)
- **Resolution & Formats**: Export both vector `.pdf` / `.svg` and `300 DPI` `.png` (`bbox_inches='tight'`, `facecolor='white'`).
- **Typography**: `Arial` / `Helvetica` / `DejaVu Sans`, Title `12–13 pt bold`, Axis labels `10.5 pt`, Tick labels `9 pt`, Gene callouts `8.5 pt bold italic` with repulsive non-overlapping leader lines (`adjustText`-style force-directed relaxation).
- **Mandatory Visual Honesty**:
  - Every volcano / 4-quadrant / 9-sector plot must display the exact sample size $N$, threshold lines ($\pm\log_2\text{FC}$, FDR or Stage 1.5 Empirical Null $\widehat{\text{FDR}}_{\text{emp}}$), and distinguish `unique_peptides >= 2` (solid circles) from `unique_peptides == 1` (open triangles).

---

## 2. BioRender-Style Mechanistic Vector SVG Engine (`07_biorender_canvas_engine.py`)

When generating pathway & organelle mechanism schematics (`fig_mechanism_biorender.svg` and `.png`), use the programmatic `BioRenderCanvas` engine rather than generic flowcharts.

### 2.1 Organelle & Cellular Primitives Supported by `BioRenderCanvas`
1. **Lipid Bilayer Plasma Membrane (`draw_bilayer_membrane`)**: Dual-leaflet phospholipid headgroup array with hydrophobic acyl tails and extracellular/cytosolic shading gradients.
2. **Multivesicular Body (`MVB`) & Intraluminal Vesicles (`ILVs`) (`draw_mvb_and_ev`)**: Double-membrane late endosome encapsulating budding/secreted EVs (`30–150 nm` exosomes) studded with transmembrane tetraspanins (`CD9/CD63/CD81`), integrins (`ITGA3/ITGB1`), and Rab GTPase tethers (`RAB27A`) + SNARE complexes (`VAMP7`).
3. **Nucleus, ER-Golgi Secretory Network, & Mitochondria (`draw_organelle`)**: Soft radial-gradient SVG `<defs>` with drop shadows (`<filter id="drop-shadow">`) and inner cristae/cisternae paths.
4. **Antibody / Biologic Y-Shape (`draw_antibody_igg`)**: Fab variable domains (parental vs humanized CDR loops highlighted) and Fc constant stalk.
5. **Protein / Complex Pill Badges (`draw_protein_node`)**:
   - Displays exact `Gene Symbol` + quantitative badge (`Δlog2FC = +1.42` or `pEC50 = 7.4`) + peptide tier indicator (`*` if `1-pep`).
6. **Evidence-Graded Causal Arrows (`draw_causal_arrow`)**:
   - **`E1 / E2` Validated / Statistically Calibrated Interaction**: **Solid Bezier curve (`stroke-width="2.6"`)** with arrowhead (`->`) or inhibitory T-bar (`-|`).
   - **`E3` Exploratory / Single-Tube / Unverified Hypothesis**: **Dashed Bezier curve (`stroke-dasharray="6,4"`, `stroke-width="2.0"`)** with mandatory `[E3 Hypothesis]` tag.
   - **Never draw a solid causal arrow for an `E3` exploratory finding or a contrast exceeding the Stage 1.5 empirical null FDR threshold!**

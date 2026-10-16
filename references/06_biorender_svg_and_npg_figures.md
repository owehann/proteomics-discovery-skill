<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 06: 14-Figure NPG Statistical Suite & 7-Schematic BioRender Vector Engine

> [!IMPORTANT]
> **Stage Scope (`Stage 6` — Visual Synthesis)**
> This playbook specifies the publication visual standards and programmatic primitives implemented in [`scripts/06_npg_statistical_figures.py`](../scripts/06_npg_statistical_figures.py) and [`scripts/07_biorender_canvas_engine.py`](../scripts/07_biorender_canvas_engine.py).

---

## 1. NPG / Nature Cell Biology Color Palette & Typography Standards

### 1.1 Standard Chromatic Encoding

| Biological Direction / Tier | Primary Hex | Fill Tint | Stroke Hex | Semantic Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Up-Regulated / Activated / Tumor / OE Gain** | `#E64B35` | `#FDE8E4` | `#B8321F` | Positive fold-change, oncogene activation, hyper-secretion |
| **Down-Regulated / Suppressed / Non-Tumor / KD Loss** | `#4DBBD5` | `#E3F5FA` | `#1F829C` | Negative fold-change, knockdown depletion, tumor suppressor |
| **Reciprocal Core / Shared Co-Regulated / Dual-Up** | `#00A087` | `#E0F5F1` | `#006B5A` | True `KD × OE` rescue (`R-Score`), convergent `G1/G2` core |
| **Specific / Divergent / PTM Stoichiometric Switch** | `#3C5488` | `#E6EAF2` | `#233459` | Perturbation-specific sector (`G3–G6`), net PTM switch |
| **Warning / Off-Target / Single-Peptide (`1-pep`)** | `#D97706` | `#FEF3C7` | `#92400E` | Hairpin divergence, `unique_peptides == 1`, `E3` exploratory |
| **Neutral / `TOST_Equiv` / Homeostatic Baseline** | `#8491B4` | `#F1F5F9` | `#64748B` | Statistically equivalent (`TOST p < 0.05`), background |

### 1.2 Publication Plot Rules (`06_npg_statistical_figures.py`)
- **Resolution & Formats**: Export both vector `.pdf` / `.svg` and `600 DPI` `.png` (`bbox_inches="tight"`, `facecolor="white"`).
- **Typography & Non-Overlapping Labels**: `Arial` / `Helvetica` / `DejaVu Sans`, Title `12–13 pt bold`, Axis labels `10.5 pt`, Tick labels `9 pt`, Gene callouts `8.5 pt bold` with force-directed repulsive leader lines (`place_non_overlapping_labels()`).
- **Mandatory Complete 14-Figure Suite (`plot_01` through `plot_14`)**:
  - Every comprehensive analysis executes `generate_all_14_publication_figures()` to generate all `20+` distinct statistical figures (`s04` through `s24`), ensuring **zero duplicate images across slides** in a 20–28 slide presentation.

---

## 2. The 7-Schematic BioRender Vector SVG Engine (`07_biorender_canvas_engine.py`)

When generating cellular, organelle, and pathway schematics, `ProteoSkill` uses the programmatic `BioRenderCanvas` vector engine (`generate_all_biorender_schematics()`) rather than generic flowcharts.

### 2.1 Organelle & Cellular Primitives in `BioRenderCanvas`
1. **Lipid Bilayer Plasma Membrane (`draw_bilayer_membrane`)**: Dual-leaflet phospholipid headgroup array with hydrophobic acyl tails and extracellular/cytosolic shading gradients.
2. **Multivesicular Body (`MVB`) & Intraluminal Vesicles (`ILVs`) (`draw_mvb_and_ev`)**: Double-membrane late endosome encapsulating budding/secreted EVs (`30–150 nm` exosomes) studded with 4-TM tetraspanins/MARVEL scaffolds (`draw_4tm_protein`), $\alpha/\beta$ integrin heterodimers (`draw_integrin_dimer`), GPI-anchored proteoglycans (`draw_gpi_protein`), ECM fibrils (`draw_ecm_fibril`), and Rab/SNARE tethers.
3. **Nucleus, ER–Golgi Secretory Network, & Mitochondria (`draw_organelle`, `draw_golgi_stack`)**: Soft radial-gradient SVG `<defs>` with drop shadows (`<filter id="soft-shadow">`) and inner cisternae paths.
4. **Antibody / Biologic Y-Shape (`draw_antibody_igg`)**: Fab variable domains and Fc constant stalk.
5. **Quantitative Protein Pill Badges (`draw_protein_node`)**: Displays `Gene Symbol` + quantitative badge (`Δlog2FC = +1.42` or `R-Score = +2.1`) + automatic `[1-pep*]` dashed border if `unique_peptides == 1`.
6. **Evidence-Graded Causal Arrows (`draw_causal_arrow`)**:
   - **`E1 / E2` Validated / Calibrated Interaction**: **Solid Bezier curve (`stroke-width="2.5"`)** with arrowhead (`->`) or inhibitory T-bar (`-|`).
   - **`E3` Exploratory / Single-Tube / Unverified Hypothesis**: **Dashed Bezier curve (`stroke-dasharray="7,5"`)** with mandatory `[E3]` tag.

---

### 2.2 The 7 Mandatory BioRender Schematics Generated per Run

| Schematic ID | Builder Function | Slide Placement & Biological Role |
| :--- | :--- | :--- |
| **`s01_biorender_hero`** | `build_s01_biorender_hero()` | **Cover / Graphical Abstract (`Slide 1`)**: 3-stage end-to-end visual summary of perturbation, mechanism, and clinical cohort translation. |
| **`s03_biorender_exp_design`** | `build_s03_biorender_exp_design()` | **Experimental Architecture (`Slide 3`)**: Perturbation arms, Dual-Track `DEqMS + <LOD` quantification, and 3-Lock causal filtering. |
| **`s14_biorender_topology`** | `build_s14_biorender_topology()` | **Subcellular & Membrane Topology (`Slide 15`)**: Classical ER–Golgi `1-TM / Signal Peptide / ECM` route vs Non-Classical `4-TM MARVEL / Tetraspanin / Lipid-Raft` route. |
| **`s17_biorender_lock_and_key`** | `build_s17_biorender_lock_and_key()` | **Physical Lock-and-Key Subnetwork (`Slide 19`)**: Direct `STRING` physical/functional bridges linking perturbed baits to `KDA` hub drivers and effector modules. |
| **`s24_biorender_spotlight`** | `build_s24_biorender_spotlight()` | **Top-6 Flagship Subcellular Cards (`Slide 24`)**: 6-card grid showing each flagship target's subcellular niche, quantitative metrics, and translational hook. |
| **`s25_biorender_working_model`** | `build_s25_biorender_working_model()` | **Integrated Working Model (`Slide 25`)**: Complete Inside-Out cellular mechanism integrating `E1/E2` solid vs `E3` dashed pathways. |
| **`s26_biorender_single_kd_summary`** | `build_s26_biorender_single_kd_summary()` | **Mandatory Penultimate Slide (`Slide N-1`)**: Biologist-First Cheat Sheet of **Core Changed Genes (Up/Down), Dysregulated Pathways, Landmark Literature (`Author et al., Year, PMID`), and Biological Meaning**. |

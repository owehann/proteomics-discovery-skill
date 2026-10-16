<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 03: PTM, Interactomics, Thermal Stability (TPP/CETSA), & Time/Dose Models (M2–M6)

> [!IMPORTANT]
> **Stage Scope (`Modality M2–M6` / `Archetypes A7–A10`)**
> This playbook governs specialized proteomics modalities executed in [`scripts/02_statistical_engine.R`](../scripts/02_statistical_engine.R), [`scripts/02_statistical_engine_py.py`](../scripts/02_statistical_engine_py.py), and [`scripts/03_archetype_and_modality_solver.py`](../scripts/03_archetype_and_modality_solver.py).

---

## 1. Modality M2 / Archetype A7: Post-Translational Modifications (`PTM`)

### 1.1 Localization Hygiene & Multiplicity Disambiguation
- Retain **Class I localized sites** (`Localization Probability >= 0.75` in MaxQuant/Spectronaut or `PTM.Site.Confidence >= 0.75` in DIA-NN).
- Track site multiplicity (`_M1` single-phosphorylated vs `_M2/_M3` multiphosphorylated peptides) separately, as multiphosphorylated sites often reflect processive kinase cascades (`GSK3`, `CK1`, `CDK1`).
- Key every row by `UniProt_Accession | Gene_Symbol | Site_Code | Multiplicity` (e.g., `P04637|TP53|S15|M1`).

### 1.2 Mandatory Uncoupling of PTM Stoichiometry from Host Protein Abundance

A 4-fold rise in a phosphopeptide is not a signaling switch if the parent protein abundance also rose 4-fold. Whenever paired Global Proteome (`WCL`) data exist:

1. **Net Stoichiometric Shift (`MSstatsPTM` Contrast Difference)**:

$$
\Delta_{\mathrm{net,PTM}, s} = \Delta\log_2(\mathrm{PTM}_s) - \Delta\log_2(\mathrm{Protein}_{g(s)})
$$

   with propagated standard error:

$$
\mathrm{SE}_{\mathrm{net,PTM}, s} = \sqrt{\mathrm{SE}_{\mathrm{PTM}, s}^2 + \mathrm{SE}_{\mathrm{Prot}, g(s)}^2}
$$

2. **Sample-Wise Linear Residualization**: When samples are paired tube-for-tube, regress site $\log_2(\text{PTM})$ on host $\log_2(\text{Protein})$ across tubes $j$ and test differential residuals.
3. **Biological Classification**:
   - **`True Stoichiometric Switch`**: $|\Delta_{\mathrm{net,PTM}}| \ge c$ and $q_{\mathrm{net,PTM}} < 0.05$.
   - **`Abundance-Driven PTM Shift`**: $|\Delta_{\mathrm{raw,PTM}}| \ge c$ but $|\Delta_{\mathrm{net,PTM}}| < 0.322$.

### 1.3 Kinase & E3 Ligase Upstream Inference
- **Phosphoproteomics (`STY`)**:
  - Extract flanking $\pm 7\text{ aa}$ sequence windows around each modified residue (`Position - 7` to `Position + 7`).
  - Score against **The Kinase Library** position-specific scoring matrices (Johnson et al., *Nature* 2023; Yaron-Barir et al., *Nature* 2024), **PTM-SEA / PTMsigDB**, and **KSEA / RoKAI**:

$$
Z_K = \frac{(\bar{\Delta}_K - \bar{\Delta}_{\mathrm{all}})\sqrt{m_K}}{\sigma_{\mathrm{all}}}
$$

- **Ubiquitinomics (`K-ε-GG`)**:
  - Score E3 ligase–substrate pairs (`UbiBrowser 2.0` / `OmniPath` E3 interactions) and distinguish **degradative ubiquitylation** ($\Delta\mathrm{Ub} > +c$ and $\Delta\mathrm{WCL} < -c$) from **non-degradative signaling ubiquitylation** ($\Delta\mathrm{Ub} > +c$ and `TOST_Equiv`($\Delta\mathrm{WCL}$)).

---

## 2. Modality M3 / Archetype A8: Interactomics & Proximity Labeling (`AP-MS / BioID / TurboID / APEX`)

### 2.1 Bait-Normalized Stoichiometry & Control Subtraction
1. **Forbid Global Median Normalization**: Never apply global median normalization when comparing `Bait IP` against `IgG / Empty Vector / BirA*-GFP` negative controls (which naturally have much lower total protein yield).
2. **Bait-Normalized Stoichiometry**: For each prey $i$ relative to bait $B$ in sample $j$:

$$
\mathrm{Stoich}_{i, j} = \log_2(I_{i, j}) - \log_2(I_{B, j}) - \log_2\!\left(\frac{L_i}{L_B}\right)
$$

   where $L_i$ is the number of theoretical tryptic peptides (`iBAQ`-style observable length).

### 2.2 Specific Interactor Scoring (`SAINT`-Style Specificity / `CompPASS` / `MiST`)
- **Fold-Change over Control (`FC_A` & `FC_B`)**: Ratio against mean control intensity (`FC_A`) and 90th-percentile control intensity (`FC_B`) with background floor $\tau_{\mathrm{floor}}$.
- **Empirical Specificity Score (`heuristic_specificity_score` / `saint_lite_prob`)**: Posterior-like specificity score $\in [0, 1]$ modeling prey enrichment against `IgG/Empty` controls and `CRAPome` contaminant frequency priors.
- **MiST Composite Score** (for multi-bait networks): Weighted sum of Reproducibility ($0.3085$), Abundance ($0.0059$), and Specificity across baits ($0.6855$).
- **High-Confidence Interactor Gate**: `heuristic_specificity_score >= 0.80` (or `BH FDR < 0.05` and $\Delta\log_2(\mathrm{Bait/IgG}) \ge 2.0$), `unique_peptides >= 2`, and `CRAPome_freq < 0.20`.

---

## 3. Modality M4 / Archetype A9: Thermal & Structural Proteome Profiling (`TPP / CETSA / PISA / LiP-MS`)

### 3.1 TPP-TR Sigmoidal Melting & `NPARC` Functional $F$-Test

Given soluble fraction measurements $y_i(T) \in [0, 1.2]$ across temperatures $T_1$, $T_2$, $\dots$, $T_m$ (normalized to $T_1 = 37^\circ\mathrm{C}$):

1. **Sigmoidal Denaturation Curve & Midpoint Melting Temperature ($T_{m, i}$ where $f(T_{m, i}) = 0.5$)**:

$$
f(T;\ p_i, a_i, b_i) = p_i + \frac{1 - p_i}{1 + \exp\!\left(b_i - \frac{a_i}{T}\right)}, \qquad T_{m, i} = \frac{a_i}{b_i - \ln\!\left(\frac{0.5 - p_i}{0.5}\right)}
$$

2. **NPARC (Non-Parametric Analysis of Response Curves, Childs et al., *MCP* 2019)**:
   Fit a single null curve $f_0(T)$ to combined `Vehicle + Drug` measurements ($\mathrm{RSS}_{0, i}$) and separate condition-specific curves $f_{\mathrm{Veh}}(T)$ and $f_{\mathrm{Drug}}(T)$ with combined residual sum of squares $\mathrm{RSS}_{1, i}$. Compute the empirical-Bayes moderated $F$-statistic:

$$
F_i = \frac{(\mathrm{RSS}_{0, i} - \mathrm{RSS}_{1, i}) / (d_1 - d_0)}{\tilde{s}_i^2}
$$

   This captures both classic melting point shifts ($\Delta T_m$) and non-canonical plateau or cooperativity shifts.

### 3.2 PISA & LiP-MS
- **PISA (Proteome Integral Solubility Alteration)**: Test the pooled temperature integral shift $\Delta S_m$ between Drug and Vehicle via `DEqMS` moderated $t$-tests after subtracting $37^\circ\mathrm{C}$ abundance shifts.
- **LiP-MS (Limited Proteolysis-Coupled MS)**: Test half-tryptic conformational cleavage peptides vs fully tryptic peptides (`MSstatsLiP`) adjusting for host protein abundance to map structural binding pockets.

---

## 4. Modality M5 / Archetype A10: Dose-Response (`EC50 / IC50`) & Time-Course

### 4.1 Dose-Response (`CurveCurator` 4-Parameter Log-Logistic Calibration)

For log-transformed concentrations $x = \log_{10}(\mathrm{Dose})$:

$$
y(x) = \mathrm{Front} + \frac{\mathrm{Back} - \mathrm{Front}}{1 + 10^{\mathrm{Slope} \cdot (x - \mathrm{pEC}_{50})}}
$$

- Compute the `CurveCurator` calibrated $F$-statistic and relevance score (Bayer et al., *Nature Communications* 2023) comparing the 4-parameter Hill curve against a flat horizontal null line ($y(x) = \bar{y}$), enforcing biologically plausible slopes ($0.3 \le |\mathrm{Slope}| \le 10$) and $R^2 \ge 0.70$.
- Report $\mathrm{pEC}_{50}$ ($-\log_{10}$ of molar $\mathrm{EC}_{50}$), maximal fold-change, and potency selectivity index across off-targets.

### 4.2 Longitudinal Time-Course (`limma` Natural Cubic Splines)
- Model time $t$ via natural cubic splines (`ns(t, df = 3)`) with `Time * Condition` interaction $F$-test.
- Classify significant trajectories into **`Immediate-Early Transient Impulse`**, **`Sustained Monotonic Shift`**, or **`Delayed Secondary Wave`**.

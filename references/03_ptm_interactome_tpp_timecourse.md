# Reference Playbook 03: PTM, Interactomics, Thermal Stability (TPP/CETSA), & Time/Dose Models (M2–M6 / A7–A10)

This playbook governs specialized proteomics modalities **M2–M6** and archetypes **A7–A10** executed in `scripts/02_statistical_engine.R`, `scripts/02_statistical_engine_py.py`, and `scripts/03_archetype_and_modality_solver.py`.

---

## 1. Modality M2 / Archetype A7: Post-Translational Modifications (`PTM`)

### 1.1 Localization Hygiene & Multiplicity Disambiguation
- Retain **Class I localized sites** (`Localization Probability >= 0.75` in MaxQuant/Spectronaut or `PTM.Site.Confidence >= 0.75` in DIA-NN).
- Track site multiplicity (`_M1` single-phosphorylated vs `_M2/_M3` multiphosphorylated peptides) separately, as multiphosphorylated sites often reflect processive kinase cascades (`GSK3`, `CK1`, `CDK1`).
- Key every row by `UniProt_Accession | Gene_Symbol | Site_Code | Multiplicity` (e.g., `P04637|TP53|S15|M1`).

### 1.2 Mandatory Uncoupling of PTM Stoichiometry from Total Protein Abundance
A 4-fold rise in a phosphopeptide is biologically meaningless as a signaling event if the parent protein abundance also rose 4-fold!Whenever paired Global Proteome (`WCL`) data exist:
1. **MSstatsPTM / Contrast Difference**:
   $$\Delta_{\text{net PTM}, s} = \Delta\log_2(\text{PTM}_{s}) - \Delta\log_2(\text{Parent Protein}_{g(s)})$$
   with propagated standard error:
   $$\text{SE}_{\text{net PTM}, s} = \sqrt{\text{SE}_{\text{PTM}, s}^2 + \text{SE}_{\text{Prot}, g(s)}^2}$$
2. **Sample-Wise Linear Residualization** (when samples are paired tube-for-tube): Regress $\log_2(\text{PTM}_{s, j})$ on $\log_2(\text{Prot}_{g(s), j})$ across tubes $j$ and test differential residuals.
3. **Classification**:
   - **`True Stoichiometric Switch`**: $|\Delta_{\text{net PTM}}| \ge c$ and $q_{\text{net PTM}} < 0.05$.
   - **`Abundance-Driven PTM Shift`**: $|\Delta_{\text{raw PTM}}| \ge c$ but $|\Delta_{\text{net PTM}}| < 0.322$.

### 1.3 Kinase & E3 Ligase Upstream Inference
- **Phosphoproteomics (`STY`)**:
  - Extract flanking $\pm 7\text{ aa}$ sequence windows around each modified residue (`Position - 7` to `Position + 7`).
  - Score against **The Kinase Library** position-specific scoring matrices (Johnson et al., *Nature* 2023; Yaron-Barir et al., *Nature* 2024), **PTM-SEA / PTMsigDB**, and **KSEA / RoKAI** ($Z$-score of substrate mean shift: $Z_K = \frac{(\bar{\Delta}_K - \bar{\Delta}_{\text{all}})\sqrt{m_K}}{\sigma_{\text{all}}}$).
- **Ubiquitinomics (`K-ε-GG`)**:
  - Score E3 ligase–substrate pairs (`UbiBrowser 2.0` / `OmniPath` E3 interactions) and distinguish **degradative ubiquitylation** ($\Delta\text{Ub} > +c$ AND $\Delta\text{WCL} < -c$) from **non-degradative signaling ubiquitylation** ($\Delta\text{Ub} > +c$ AND $\text{TOST\_Equiv}(\Delta\text{WCL})$).

---

## 2. Modality M3 / Archetype A8: Interactomics & Proximity Labeling (`AP-MS / BioID / TurboID / APEX`)

### 2.1 Bait-Normalized Stoichiometry & Control Subtraction
1. **Never use global median normalization** when comparing `Bait IP` against `IgG / Empty Vector / BirA*-GFP` negative controls (which naturally have much lower total protein yield).
2. Compute **Bait-Normalized Stoichiometry** for each prey $i$ relative to bait $B$:
   $$\text{Stoich}_{i, j} = \log_2(I_{i, j}) - \log_2(I_{B, j}) - \log_2(L_i / L_B)$$
   where $L_i$ is the number of theoretical tryptic peptides (or observable sequence length, `iBAQ`-style).

### 2.2 Specific Interactor Scoring (`SAINT`-lite / `CompPASS` / `MiST`)
For each prey $i$ in bait $B$ vs negative controls $C$:
1. **Fold-Change over Control (`FC_A` & `FC_B`)**: Ratio against mean control intensity (`FC_A`) and maximum/90th-percentile control intensity (`FC_B`) with background floor $\tau_{\text{floor}}$.
2. **Bayesian Empirical Probability (`SAINT`-lite Score)**: Posterior probability $P(\text{True Interactor}_i \mid \mathbf{y}_{i, B}, \mathbf{y}_{i, C}) \in [0, 1]$ modeling prey intensity via a two-component Gaussian/Poisson mixture (background contaminant distribution parameterized by `IgG/Empty` controls + `cRAP/CRAPome` frequency prior vs specific bait capture).
3. **MiST Composite Score** (for multi-bait networks):Weighted sum of Reproducibility ($0.3085$), Abundance ($0.0059$), and Specificity across baits ($0.6855$).
4. **High-Confidence Interactor Gate**: `SAINT_prob >= 0.80` (or `BH FDR < 0.05` & $\Delta\log_2(\text{Bait/IgG}) \ge 2.0$), `unique_peptides >= 2`, and `CRAPome_freq < 0.20`.

---

## 3. Modality M4 / Archetype A9: Thermal & Structural Proteome Profiling (`TPP / CETSA / PISA / LiP-MS`)

### 3.1 TPP-TR (Temperature Range) Sigmoidal Melting & `NPARC` Functional $F$-Test
Given soluble fraction measurements $y_i(T) \in [0, 1.2]$ across temperatures $T_1 < T_2 < \dots < T_m$ (normalized to lowest temperature $T_1 = 37^\circ\text{C}$):
1. **4-Parameter / 3-Parameter Sigmoidal Denaturation Curve**:
   $$f(T; p_i, a_i, b_i) = p_i + \frac{1 - p_i}{1 + \exp\left(b_i - \frac{a_i}{T}\right)}, \qquad T_{m, i} = \frac{a_i}{b_i}$$
2. **NPARC (Non-Parametric Analysis of Response Curves, Childs et al. *MCP* 2019)**:
   Fit a single null curve $f_0(T)$ to combined `Vehicle + Drug` measurements ($\text{RSS}_0$) and separate condition-specific curves $f_{\text{Veh}}(T), f_{\text{Drug}}(T)$ ($\text{RSS}_1 = \text{RSS}_{\text{Veh}} + \text{RSS}_{\text{Drug}}$). Compute the empirical-Bayes moderated $F$-statistic:
   $$F_i = \frac{(\text{RSS}_{0, i} - \text{RSS}_{1, i}) / (d_1 - d_0)}{\tilde{s}_i^2}$$
   This captures both classic $\Delta T_m = T_{m, \text{Drug}} - T_{m, \text{Veh}}$ shifts and non-canonical plateau/slope shifts.

### 3.2 PISA (Proteome Integral Solubility Alteration)
When soluble fractions across temperatures are pooled prior to MS injection (`PISA` assay), test the integral solubility shift $\Delta S_m = \mu_{\text{PISA, Drug}} - \mu_{\text{PISA, Veh}}$ using `limma`/`DEqMS` moderated $t$-tests, subtracting any total abundance change measured at $37^\circ\text{C}$.

### 3.3 LiP-MS (Limited Proteolysis-Coupled MS)
Test half-tryptic (conformational cleavage) vs fully tryptic peptides at the peptide level (`MSstatsLiP` / linear mixed model adjusting for parent protein abundance) to pinpoint the exact structural binding pocket (residue window) on the target protein.

---

## 4. Modality M5 / Archetype A10: Dose-Response ($EC_{50}/IC_{50}$) & Longitudinal Time-Course

### 4.1 Dose-Response (`CurveCurator` 4-Parameter Log-Logistic Calibration)
For concentrations $x = \log_{10}(\text{Dose})$:
$$y(x) = \text{Front} + \frac{\text{Back} - \text{Front}}{1 + 10^{\text{Slope} \cdot (x - \text{pEC}_{50})}}$$
- Compute the `CurveCurator` calibrated $F$-statistic and relevance score (Bayer et al., *Nature Communications* 2023) comparing the 4-parameter hill curve against a flat horizontal null line ($y(x) = \bar{y}$), enforcing biologically plausible hill slopes ($0.3 \le |\text{Slope}| \le 10$) and $R^2 \ge 0.70$.
- Report $\text{pEC}_{50} = -\log_{10}(\text{EC}_{50}\text{ [M]})$, maximal fold change, and potency selectivity index across off-targets.

### 4.2 Longitudinal Time-Course (`limma` Natural Cubic Splines & Impulse Classification)
- Model time $t$ via natural cubic splines (`ns(t, df = 3)`) with `Time * Condition` interaction $F$-test.
- Classify significant trajectories into **`Immediate-Early Transient Impulse`** (peaks at early $t$, returns to baseline), **`Sustained Monotonic Shift`**, or **`Delayed Secondary Wave`**.

<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 02: Abundance & Experimental Archetype Models (A1–A6)

> [!IMPORTANT]
> **Stage Scope (`Stage 2` and `Stage 3`)**
> This playbook documents the statistical models and algebraic scoring formulas implemented in [`scripts/02_statistical_engine.R`](../scripts/02_statistical_engine.R), [`scripts/02_statistical_engine_py.py`](../scripts/02_statistical_engine_py.py), and [`scripts/03_archetype_and_modality_solver.py`](../scripts/03_archetype_and_modality_solver.py).

---

## 1. Sample-Size & Modality Routing Matrix

| Sample Regime | Recommended Primary Engine | Missing Value Strategy | Hypothesis Test & Multiple Testing |
| :--- | :--- | :--- | :--- |
| **<i>N</i> ≥ 5 per group**<br/>*(Clinical / Large Cohort)* | `limpa` (`dpc` + `dpcDE`) or `limma` + `DEqMS` + `msqrob2` Hurdle | Probabilistic DPC recovery (`limpa`) or two-part Hurdle (`msqrob2`) | Moderated <i>t</i> / peptide-count eBayes + Benjamini–Hochberg FDR + Wilcoxon sensitivity |
| **<i>N</i> = 3–4 per group**<br/>*(Standard Wet-Lab)* | `limpa` (`dpcDE`) or `DEqMS` (PSM-weighted empirical Bayes) | DPC posterior variance weighting (`limpa`) + Track B `<LOD` Binary Switches; zero random imputation | Peptide-count moderated <i>t</i>-test (`DEqMS` / `squeezeVar`) + BH FDR |
| **<i>N</i> = 2 per group**<br/>*(Small Replicate / Multi-Clone)* | Empirical Bayes `squeezeVar` + Shared-Control Covariance Adjustment | Separate continuous moderated test from Hurdle / `<LOD` presence-absence switches | Moderated <i>p</i>-value + Stage 1.5 Empirical Null FDR calibration |
| **<i>N</i> = 1 per group**<br/>*(Multi-Clone Screen / Pilot)* | Stage 1.5 Reference-Noise <i>Z</i>-score against `sigma_null` | Exact binomial / DPC probability against empirical null switch rate | Empirical Null FDR (`FDR_emp`); all single-tube claims capped at **`E3`** |

---

## 2. Analytical Empirical Bayes (`squeezeVar`) & DPC Missingness Recovery

In both the R (`limma::squeezeVar` / `DEqMS::spectraCounteBayes`) and Python (`02_statistical_engine_py.py`) engines:

1. **Peptide-Stratified Variance Shrinkage**: For each protein $i$ with residual variance $s_i^2$ and residual degrees of freedom $d_i$, we fit the prior variance $s_0^2$ and prior degrees of freedom $d_0$ (trended against peptide count or mean intensity $\bar{y}_i$) using trigamma/digamma moment matching (Smyth 2004):

$$
\tilde{s}_i^2 = \frac{d_0 s_0^2 + d_i s_i^2}{d_0 + d_i}, \qquad \tilde{t}_i = \frac{\hat{\beta}_i}{\tilde{s}_i \sqrt{v_i}}
$$

   where $\tilde{t}_i$ is evaluated against a Student's $t$-distribution with $d_0$ + $d_i$ degrees of freedom.

2. **Detection Probability Curve (`DPC`) & Track B `<LOD` Switches**: Rather than fabricating numbers via random Gaussian tail (`MinProb`) imputation, we fit the logistic detection curve across mean intensities $\bar{y}_i$:

$$
\mathrm{logit}\,P(\mathrm{observed}_{ij}) = \beta_0 + \beta_1 \bar{y}_i
$$

   Proteins completely absent (`<LOD` in $0/3$ replicates) in one group and consistently detected ($3/3$ or $\ge 2/3$ replicates) in the compared group are evaluated on **Track B (All-or-None Binary Switches)** with an adaptive intensity floor for unified volcano visualization.

---

## 3. Two One-Sided Tests (`TOST`) for Formal Equivalence Proof

> [!TIP]
> **Eliminating the "Absence-of-Evidence" Fallacy**
> Never claim a protein, complex, or pathway is **"preserved"**, **"unperturbed"**, or **"Target-A Specific"** merely because $p > 0.05$ in Condition B. Require formal **Two One-Sided Tests (`TOST`)** within an equivalence margin [−$\Delta_{\mathrm{eq}}$, +$\Delta_{\mathrm{eq}}$] (default $\Delta_{\mathrm{eq}} \approx 0.322$ for $1.25\times$, or $\pm 0.50$ in small-$N$ screens).

For each protein $i$ with estimated fold-change $\hat{\beta}_i$ and moderated standard error $\mathrm{SE}_i$:

$$
t_{\mathrm{lower}, i} = \frac{\hat{\beta}_i + \Delta_{\mathrm{eq}}}{\mathrm{SE}_i}, \qquad t_{\mathrm{upper}, i} = \frac{\Delta_{\mathrm{eq}} - \hat{\beta}_i}{\mathrm{SE}_i}
$$

$$
p_{\mathrm{TOST}, i} = \max\left(1 - F_{t(\nu)}(t_{\mathrm{lower}, i}),\ 1 - F_{t(\nu)}(t_{\mathrm{upper}, i})\right)
$$

A protein or macromolecular complex is certified as **Statistically Equivalent (`TOST_Equiv`)** only when $p_{\mathrm{TOST}, i} < 0.05$ (meaning its $90\%$ confidence interval lies strictly inside [−$\Delta_{\mathrm{eq}}$, +$\Delta_{\mathrm{eq}}$]).

---

## 4. Formal Definitions of Experimental Archetypes (A1–A5)

### 4.1 Archetype A1: Single-Gene Perturbation (`KD / KO / OE / CRISPRa`)
- **Primary Test**: Moderated contrast $\Delta_{\mathrm{KD}} = \mu_{\mathrm{KD}} - \mu_{\mathrm{Ctrl}}$ with peptide-count variance shrinkage.
- **Multi-Clone Consistency Filter**: When $k \ge 2$ independent hairpins or sgRNAs (`Clone_1, Clone_2, ...`) share a control (`NTC`), compute both:
  1. Pooled linear contrast with clone block effect.
  2. Minimum effect across clones ($\min |\Delta_j|$) with exact directional agreement, calibrated against the Stage 1.5 shared-control null integral.
- **Clone Divergence Warning**: Any protein with clone difference exceeding $1.0$ ($\log_2$ scale) or opposite signs between hairpins is flagged as **`Clone-Specific Off-Target / Seed Effect`**.

---

### 4.2 Archetype A2: Reciprocal Perturbation (`KD` × `OE` or `KO` + `Rescue`)
- **Reciprocity Score (`R-Score`)**: Quantifies strict bidirectional mirror rescues between knockdown ($\bar{\Delta}_{\mathrm{KD}, i}$) and overexpression ($\Delta_{\mathrm{OE}, i}$):

$$
R_{\mathrm{score}, i} = -\mathrm{sign}(\Delta_{\mathrm{OE}, i}) \cdot \mathbb{I}\!\left(\mathrm{sign}(\Delta_{\mathrm{OE}, i}) = -\mathrm{sign}(\bar{\Delta}_{\mathrm{KD}, i})\right) \cdot \sqrt{|\Delta_{\mathrm{OE}, i}| \cdot |\bar{\Delta}_{\mathrm{KD}, i}|}
$$

- **Directional Stouffer Meta-$Z$**: Combines $Z_{\mathrm{OE}, i}$ and $-Z_{\mathrm{KD}, i}$ across orthogonal perturbations:

$$
Z_{\mathrm{recip}, i} = \frac{w_{\mathrm{OE}} Z_{\mathrm{OE}, i} - w_{\mathrm{KD}} Z_{\mathrm{KD}, i}}{\sqrt{w_{\mathrm{OE}}^2 + w_{\mathrm{KD}}^2}}
$$

- **4-Quadrant Classification**:
  - **True Reciprocal Axis (`Q2 / Q4`)**: `OE Up & KD Down` (Positive Regulated Cargo/Target) or `OE Down & KD Up` (Negative Regulated Target).
  - **Concordant Stress Artifact (`Q1 / Q3`)**: `OE Up & KD Up` or `OE Down & KD Down` — reflects general viral transduction or selection stress rather than specific target activity.

---

### 4.3 Archetype A3: Multi-Perturbation Epistasis & 9-Sector Divergence (`Perturbation A` vs `Perturbation B`)

Given two perturbations $A$ and $B$, every protein's composite effective fold-change ($\Delta_{A, i}$ and $\Delta_{B, i}$) is mapped into a **9-Sector Cartesian Grid** combined with `TOST` equivalence proof and divergence index $\mathrm{DI}_i$ = $\Delta_{A, i}$ − $\Delta_{B, i}$:

| Sector Class | Mathematical Definition | Biological Interpretation |
| :--- | :--- | :--- |
| **`S1 / S2` (`G1 / G2`) — Shared Co-Regulated Core** | `Δ_A >= c` and `Δ_B >= c` (or both `<= -c`) | Convergent epistatic axis co-regulated by both perturbations |
| **`S3 / S4` (`G3 / G4`) — Perturbation-A Specific** | `\|Δ_A\| >= c` AND `TOST_Equiv(Δ_B)` | Strictly specific to Perturbation A (proven invariant in B) |
| **`S5 / S6` (`G5 / G6`) — Perturbation-B Specific** | `\|Δ_B\| >= c` AND `TOST_Equiv(Δ_A)` | Strictly specific to Perturbation B (proven invariant in A) |
| **`S7 / S8` (`G7`) — Seesaw Antagonistic Flux** | `Δ_A >= c` with `Δ_B <= -c` (or vice versa) | Opposing regulatory tug-of-war / route competition between A and B |
| **`S0` (`G0`) — Dual Equivalent Baseline** | `TOST_Equiv(Δ_A)` AND `TOST_Equiv(Δ_B)` | Homeostatic background proteome |

---

### 4.4 Archetype A4: Biologics & Antibody Engineering (`Parental / Chimeric / Humanized / Isotype`)
- **3-Axis Quantitative Decomposition**:

$$
\Delta_{\mathrm{Efficacy}} = \frac{1}{2}\left(\mu_{\mathrm{Parental}} + \mu_{\mathrm{Humanized}}\right) - \mu_{\mathrm{Isotype}}
$$

$$
\Delta_{\mathrm{Fc}} = \mu_{\mathrm{Isotype}} - \mu_{\mathrm{Mock}}, \qquad \Delta_{\mathrm{Hum}} = \mu_{\mathrm{Humanized}} - \mu_{\mathrm{Parental}}
$$

  where $\Delta_{\mathrm{Hum}}$ is evaluated via **`TOST`** equivalence ($|\Delta_{\mathrm{Hum}}| < 0.322$) vs humanization-induced drift ($|\Delta_{\mathrm{Hum}}| \ge 0.585$).

---

### 4.5 Archetype A5: Clinical Paired Multi-Compartment & Secretion Gradient Index (`SGI`)

For paired biofluid (interstitial fluid, EV, or plasma) and solid tumor/tissue proteomes across patients:

$$
\mathrm{SGI}_i = \Delta\log_2\mathrm{FC}_{\mathrm{Fluid}, i} - \Delta\log_2\mathrm{FC}_{\mathrm{Tissue}, i}
$$

*(Defined strictly as the difference of within-compartment Tumor-vs-Normal $\log_2$ fold-changes, never raw cross-compartment intensity differences).*

| Clinical Quadrant | Rule (`c = log2(1.5) ≈ 0.585`) | Translational Meaning |
| :--- | :--- | :--- |
| **`Q1_concordant_up`** | `Δ_Fluid >= c` AND `Δ_Tissue >= c` | Dual Tissue + Non-Invasive Liquid Biopsy Biomarker |
| **`Q2_hyper_secreted`** | `Δ_Fluid >= c` AND `Δ_Tissue < c` | Active Vesicular Export / Ectodomain Shedding Candidate |
| **`Q3_concordant_down`** | `Δ_Fluid <= -c` AND `Δ_Tissue <= -c` | Concordantly Depleted Tumor Suppressor / Metabolic Axis |
| **`Q4_intracellular_retained`** | `Δ_Tissue >= c` AND `Δ_Fluid < c` | Intracellularly Trapped Oncoprotein |

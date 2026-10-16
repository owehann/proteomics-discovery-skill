# Reference Playbook 02: Abundance & Experimental Archetype Models (A1–A10)

This playbook documents the statistical models and algebraic scoring formulas used in **Stage 2 (`scripts/02_statistical_engine.R` & `scripts/02_statistical_engine_py.py`)** and **Stage 3 (`scripts/03_archetype_and_modality_solver.py`)**.

---

## 1. Sample-Size & Modality Routing Matrix

| Sample Regime | Recommended Primary Engine | Missing Value Strategy | Hypothesis Test & Multiple Testing |
| :--- | :--- | :--- | :--- |
| **$N \ge 5$ per group** (Clinical / Large Cohort) | `limpa` (`dpc` + `dpcDE`) or `limma` + `DEqMS` + `msqrob2` Hurdle | Probabilistic DPC recovery (`limpa`) or two-part Hurdle (`msqrob2`) | Moderated $t$ / peptide-count eBayes + Benjamini-Hochberg FDR + Wilcoxon sensitivity |
| **$N = 3\text{–}4$ per group** (Standard Wet-Lab) | `limpa` (`dpcDE`) or `DEqMS` (PSM-weighted empirical Bayes) | DPC posterior variance weighting (`limpa`); zero random imputation | Peptide-count moderated $t$-test (`DEqMS` / `squeezeVar`) + BH FDR |
| **$N = 2$ per group** (Small Replicate / Multi-Clone) | Empirical Bayes `squeezeVar` + Shared-Control Covariance Adjustment | Separate continuous moderated test from `MSqRob` Hurdle presence/absence | Moderated $p$-value + Stage 1.5 Empirical Null FDR calibration |
| **$N = 1$ per group** (Multi-Clone Screen / Pilot) | Stage 1.5 Reference-Noise Z-score: $Z_i = \Delta_i / (\sqrt{2}\hat{\sigma}_{\text{null}})$ | Exact binomial / DPC probability against null switch rate $p_{\text{switch,null}}$ | Empirical Null FDR $\widehat{\text{FDR}}_{\text{emp}}$; all single-tube claims capped at **`E3`** |

---

## 2. Analytical Empirical Bayes (`squeezeVar`) & DPC Equivalence

In both the R (`limma::squeezeVar` / `DEqMS::spectraCounteBayes`) and Python (`02_statistical_engine_py.py`) engines:
1. For each protein $i$ with residual variance $s_i^2$ and residual degrees of freedom $d_i$, fit the log-variance distribution $\log s_i^2$ (optionally trended against $\log_2(\text{unique\_peptides}_i)$ or mean intensity $\bar{y}_i$) using trigamma/digamma moment matching (Smyth 2004):
   $$\tilde{s}_i^2 = \frac{d_0 s_0^2 + d_i s_i^2}{d_0 + d_i}, \qquad \tilde{t}_i = \frac{\hat{\beta}_i}{\tilde{s}_i \sqrt{v_i}}$$
   evaluated against a Student's $t$-distribution with $d_0 + d_i$ degrees of freedom.
2. **Detection Probability Curve (DPC)**: Fit $\text{logit}\,P(\text{observed}_{ij}) = \beta_0 + \beta_1 \bar{y}_i$. For proteins completely missing (`NaN`) in one condition and observed in $m$ replicates of the other condition with mean $\bar{y}_{\text{obs}}$, compute the exact compound probability that $n_{\text{miss}}$ replicates fell below detection given $\bar{y}_{\text{obs}}$ vs the empirical background switch rate.

---

## 3. Two One-Sided Tests (`TOST`) for Equivalence / "Non-Perturbation"

Never claim a protein, complex, or pathway is **"preserved"**, **"unperturbed"**, or **"safe"** based on a non-significant $p > 0.05$ ("absence of evidence is not evidence of absence").
Instead, require **Two One-Sided Tests (`TOST`)** within an equivalence margin $[-\Delta_{\text{eq}}, +\Delta_{\text{eq}}]$ (default $\Delta_{\text{eq}} = \log_2(1.25) \approx 0.322$):
$$t_{\text{lower}, i} = \frac{\hat{\beta}_i - (-\Delta_{\text{eq}})}{\text{SE}_i}, \qquad t_{\text{upper}, i} = \frac{\Delta_{\text{eq}} - \hat{\beta}_i}{\text{SE}_i}$$
$$p_{\text{TOST}, i} = \max\left(1 - F_{t(\nu)}(t_{\text{lower}, i}), \; 1 - F_{t(\nu)}(t_{\text{upper}, i})\right)$$
A protein or complex is certified as **Statistically Equivalent (`TOST_Equiv`)** only when $p_{\text{TOST}, i} < 0.05$ (or its $90\%$ confidence interval lies strictly inside $[-\Delta_{\text{eq}}, +\Delta_{\text{eq}}]$).

---

## 4. Formal Definitions of Archetypes A1–A6

### Archetype A1: Single-Gene Perturbation (`KD / KO / OE / Crispra`)
- **Primary Test**: Moderated contrast $\Delta_{\text{KD}} = \mu_{\text{KD}} - \mu_{\text{Ctrl}}$ with peptide-count variance shrinkage.
- **Clone Consistency Filter**: When $k \ge 2$ independent hairpins/sgRNAs (`sh1, sh2, ...`) share a control (`NTC`), compute both:
  1. Pooled linear contrast with hairpin random/block effect.
  2. Minimum effect across clones $\text{min\_abs\_FC} = \min_j |\Delta_j|$ with exact directional agreement, calibrated against the Stage 1.5 shared-control null integral.
- **Clone Divergence Warning**: Any protein with $|\Delta_{\text{sh1}} - \Delta_{\text{sh2}}| > 1.0\,\log_2$ or opposite signs is flagged as **`Hairpin-Specific Off-Target / Seed Effect`** (e.g., divergent shift between `shRNA_1` and `shRNA_2`).

### Archetype A2: Reciprocal Perturbation (`KD` $\times$ `OE` or `KO` + `Rescue`)
- **Reciprocity Score (`R-Score`)**:
  $$\text{R-Score}_i = -\text{sign}(\Delta_{\text{OE}, i}) \cdot \mathbb{I}(\text{sign}(\Delta_{\text{OE}, i}) = -\text{sign}(\bar{\Delta}_{\text{KD}, i})) \cdot \sqrt{|\Delta_{\text{OE}, i}| \cdot |\bar{\Delta}_{\text{KD}, i}|}$$
- **4-Quadrant Classification**:
  - **True Reciprocal Axis (Q2/Q4)**: `OE Up & KD Down` (Positive Driver Cargo/Target) or `OE Down & KD Up` (Negative Regulated Target).
  - **Concordant Stress Artifact (Q1/Q3)**: `OE Up & KD Up` or `OE Down & KD Down` — reflects general lentiviral/transfection/selection stress rather than target activity.
- **Combined Meta-P-Value**: Directional Stouffer's $Z$ combining $Z_{\text{OE}}$ and $-Z_{\text{KD}}$:
  $$Z_{\text{recip}, i} = \frac{w_{\text{OE}} Z_{\text{OE}, i} - w_{\text{KD}} Z_{\text{KD}, i}}{\sqrt{w_{\text{OE}}^2 + w_{\text{KD}}^2}}$$

### Archetype A3: Multi-Perturbation Epistasis & 9-Sector Divergence (`Gene A` vs `Gene B`)
Given two perturbations $A$ and $B$ (e.g., `Target_A` vs `Target_B`, or Drug A vs Drug B), standardize or threshold $\Delta_{A, i}$ and $\Delta_{B, i}$ into the **9-Sector Cartesian Grid**:
1. **Shared Up (`A_Up & B_Up`)** / **Shared Down (`A_Down & B_Down`)**: Convergent epistatic axis.
2. **A-Specific Up / Down (`|ΔA| >= c` & `TOST_Equiv(ΔB)`)**: Strictly specific to pathway A.
3. **B-Specific Up / Down (`|ΔB| >= c` & `TOST_Equiv(ΔA)`)**: Strictly specific to pathway B.
4. **Antagonistic (`A_Up & B_Down` or `A_Down & B_Up`)**: Opposing regulatory control.
5. **Dual-Equiv (`TOST_Equiv(ΔA)` & `TOST_Equiv(ΔB)`)**: Unperturbed core proteome.

### Archetype A4: Biologics & Antibody Engineering (`Parental / Chimeric / Humanized / Isotype`)
- **3-Axis Decomposition**:
  1. **Target Efficacy Axis**: $\Delta_{\text{Efficacy}} = \frac{1}{2}(\mu_{\text{Parental}} + \mu_{\text{Humanized}}) - \mu_{\text{Isotype\_IgG}}$.
  2. **Fc / IgG Non-Specific Background Axis**: $\Delta_{\text{Fc}} = \mu_{\text{Isotype\_IgG}} - \mu_{\text{Mock/Untreated}}$.
  3. **Humanization Equivalence & Liability Axis**: $\Delta_{\text{Hum}} = \mu_{\text{Humanized}} - \mu_{\text{Parental}}$, evaluated via **`TOST`** equivalence ($|\Delta_{\text{Hum}}| < 0.322$) vs humanization-induced drift ($|\Delta_{\text{Hum}}| \ge 0.585$).

### Archetype A5: Clinical Paired Multi-Compartment & Secretion Gradient Index (`SGI`)
- **Paired Patient Linear Mixed Model**: $y_{p, c, t} = \mu + \beta_{\text{Tumor}} + \gamma_{\text{Compartment}} + \delta_{\text{Tumor}\times\text{Compartment}} + u_{\text{Patient } p} + \epsilon$.
- **Secretion Gradient Index (`SGI`)**:
  $$\text{SGI}_i = \Delta\log_2\text{FC}_{\text{Fluid (TIF/EV/Plasma)}, i} - \Delta\log_2\text{FC}_{\text{Tissue}, i}$$
  *(Explicitly defined as the difference of within-compartment Tumor-vs-NonTumor $\log_2$ fold-changes, NOT raw concentration difference between Fluid and Tissue).*
- **4-Quadrant Clinical Partition** (using strict boundary $c = \log_2(1.5) = 0.58496$ and `unique_peptides >= 2` for Core Clinical Candidates):
  - **`Q1_concordant_up`**: $\Delta_{\text{Fluid}} \ge c$ AND $\Delta_{\text{Tissue}} \ge c$ (Dual Tissue + Liquid Biopsy Biomarker).
  - **`Q2_hyper_secreted`**: $\Delta_{\text{Fluid}} \ge c$ AND $\Delta_{\text{Tissue}} < c$ (Active Secretion / Ectodomain Shedding Candidate).
  - **`Q3_concordant_down`**: $\Delta_{\text{Fluid}} \le -c$ AND $\Delta_{\text{Tissue}} \le -c$.
  - **`Q4_intracellular_retained`**: $\Delta_{\text{Tissue}} \ge c$ AND $\Delta_{\text{Fluid}} < c$ (Intracellularly Trapped Oncoprotein).
- **Scale Transparency Rule**: Whenever patient-level single-sample profiles are reported, explicitly label whether values represent **Uncentered Compartment $\log_2(\text{T/N})$** or **Patient-Centered Deviation ($\Delta\log_2(\text{T/N}) - \text{median}_p$)**.

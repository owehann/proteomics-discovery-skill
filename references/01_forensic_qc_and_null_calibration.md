<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 01: 12-Point Sample & Data QC & Empirical Null FDR Calibration

> [!IMPORTANT]
> **Stage Scope (`Stage 0`, `Stage 1`, and `Stage 1.5`)**
> This playbook defines the mandatory **Schema & Power Audit**, **12-Point Sample & Data QC**, and **Empirical Null FDR & Shared-Control Variance Inflation** rules executed by [`scripts/00_forensic_schema_indexer.py`](../scripts/00_forensic_schema_indexer.py) and [`scripts/01_forensic_qc_and_null_audit.py`](../scripts/01_forensic_qc_and_null_audit.py).

---

## 1. Stage 0: Multi-Format Schema & Power Inference

### 1.1 Supported Search Engine & Matrix Formats

| Platform | Supported Export Files | Mandatory Pre-Filtering Rule |
| :--- | :--- | :--- |
| **DIA-NN** | `report.pg_matrix.tsv`, `report.pr_matrix.tsv`, `report.parquet` | Retain `Q.Value <= 0.01`, `PG.Q.Value <= 0.01`, `Lib.Q.Value <= 0.01` |
| **MaxQuant** | `proteinGroups.txt`, `Phospho (STY)Sites.txt`, `GlyGly (K)Sites.txt`, `evidence.txt` | Exclude `Potential contaminant == '+'`, `Reverse == '+'`, `Only identified by site == '+'`; require `Localization prob >= 0.75` for PTM sites |
| **FragPipe / IonQuant** | `combined_protein.tsv`, `combined_modified_peptide.tsv`, `abundance_protein_MD.tsv` | Retain `Protein Probability >= 0.99` and strip decoy/contaminant prefixes |
| **Proteome Discoverer** | `.xlsx` / `.tsv` exports (`Abundances (Normalized)`, `Abundances (Scaled)`) | Require `Master == 'IsMasterProtein'` and `Protein FDR Confidence == 'High'` |
| **Spectronaut** | `PG.Quantity`, `PG.Genes`, `PG.ProteinAccessions`, `EG.PTMAssayProbability` | Retain `PG.Qvalue <= 0.01` and `EG.PTMAssayProbability >= 0.75` for PTMs |
| **Olink / SomaScan** | `NPX` (log₂ Normalized Protein Expression) or `RFU` matrices | Audit `LOD` and `Assay_Warning` flags per panel |

### 1.2 Mandatory Identifier & Duplicate Hygiene

1. **Excel Date-Corruption Recovery**: Before any gene-level operation, sanitize strings matching Excel auto-converted dates (`1-Mar` → `MARCHF1`, `1-Sep` → `SEPTIN1`, `1-Dec` → `DELEC1`, etc., using [`resources/curated_proteomics_signatures.json`](../resources/curated_proteomics_signatures.json)).
2. **Primary Key Policy**: Never aggregate or merge tables on `Gene Symbol` alone. Use the composite key `UniProt_Accession | Gene_Symbol` (or `UniProt_Accession | Gene_Symbol | Mod_Residue_Pos` for PTMs).
3. **Isoform Disambiguation**: When multiple protein groups map to the same gene symbol, retain both rows keyed by `Accession|Gene` and mark the canonical or highest-peptide-coverage entry as `is_canonical_rep = True`.

---

## 2. Stage 1: The 12-Point Sample & Data QC Checklist

Every proteomics dataset is audited against all 12 quality control checks before statistical modeling:

| Check ID | Quality Control Check | Quantitative Threshold / Rule | Automatic Action / Gate 2 Trigger |
| :--- | :--- | :--- | :--- |
| **QC-01** | **Contaminant Audit** | Sum of `cRAP` (`KRT1/2/9/10/14`, `TRY1`, `ALB`, `CSN2`, `IGHG*`, BSA `P02769`, trypsin `P00761`, streptavidin `P22629`) > 5% of total intensity or differential (\|Δlog₂\| > 1.0) between conditions. | Flag culture serum (`ALB/FETUA`), exogenous antibody (`IGHG*`), or streptavidin contamination; remove non-sample `cRAP` rows before normalization. |
| **QC-02** | **Peptide Support Tiering** | Partition rows into **Core Quantitative Tier** (`unique_peptides >= 2` and `PSM/precursors >= 2`) vs **Single-Peptide Exploratory Tier** (`unique_peptides == 1`). | Never allow `unique_peptides == 1` rows into headline claims, Executive Scorecards, or mechanistic core models without an explicit `[1-Peptide Exploratory]` badge. |
| **QC-03** | **Perturbation Target Audit** | Verify the perturbed target (`Gene_X` in KD/OE/KO/Bait/Drug) directly in the raw matrix. | If target is undetected (`NaN`), discordant with perturbation direction, or `unique_peptides == 1`, trigger **Gate 2** and cap target-direct claims at `E3` unless orthogonal WB/qPCR is provided. |
| **QC-04** | **Paralog & Isoform Compensation** | Inspect paralogs of the perturbed target (e.g., `MAPK3` for `MAPK1`; `RAB27B` for `RAB27A`). | Report compensatory upregulation or shared-peptide cross-interference (off-target shRNA/siRNA vs paralog buffering). |
| **QC-05** | **Tube-Level Anomaly Detection** | Compute sample-wise ID count, median log₂ intensity, IQR, pairwise Pearson <i>r</i>, PCA Mahalanobis distance, and Hierarchical Clustering. | Trigger **Gate 2** if any tube has < 65% of group median IDs, PCA > 3 SD, or <i>r</i> < 0.80 with biological replicates. |
| **QC-06** | **Sample Sex / Line Identity Check** | Cross-check sex-linked (`RPS4Y1`, `DDX3Y`, `XIST`) and lineage markers (`ALB`, `AFP`, `EPCAM`, `VIM`) across samples. | Flag tube swaps, cross-line contamination, or clinical metadata mismatches. |
| **QC-07** | **TMT / Isobaric Ratio Compression** | Check `Isolation Interference (%) > 30%`, co-isolation, and bridge-channel CV across plexes. | Recommend `IRS + TMM` + `MSstatsTMT` and warn that fold-changes are compressed toward 0. |
| **QC-08** | **DIA-NN / Spectronaut PTM & Lib Audit** | Check `PTM.Site.Confidence >= 0.75`, RT calibration drift, and `Lib.Q.Value <= 0.01`. | Separate localized (`Class I >= 0.75`) from unlocalized PTM sites. |
| **QC-09** | **Normalization Suitability Audit** | Compare global median shift across biological groups (\|Δmedian\| > 0.15 log₂ or > 15% of proteins shifting unidirectionally). | **Forbid naive median/quantile normalization** when global secretion or AP-MS pulldown asymmetry exists. Use housekeepers, spike-in standards, bait-normalization, or `TMM` / `VSN`. |
| **QC-10** | **MISEV2023 EV Purity Audit** | (EV / Secretome only) Quantify Category 1 (`CD9/63/81`, `TSG101`, `SDCBP`, `PDCD6IP`, `FLOT1`) vs Category 3 (`CANX`, `HSP90B1`, `GOLGA2`, `HSPA5`, `CYC1`, `ALB`, `APOA1`). | Report log₂(Cat 1+2 / Cat 3) purity index and flag co-isolated ER/Golgi/nuclear debris. |
| **QC-11** | **Missingness Mechanism (MNAR vs MAR)** | Fit the **Detection Probability Curve (DPC)** across mean protein intensities (`mean_log2`). | Estimate logistic slope `beta_1`. Use DPC (`limpa`) or Hurdle (`msqrob2`) + Track B `<LOD` switches rather than random left-censored (`MinProb`) imputation when replicate count is < 3. |
| **QC-12** | **Confounder Signature Scan** | Score cell-cycle (`MKI67`, `TOP2A`, `PCNA`, `CDK1`), stress/UPR (`HSPA5`, `DDIT3`), and IFN (`ISG15`, `MX1`, `STAT1`) indices per sample. | Alert when downstream enrichment is dominated by proliferation or stress arrest rather than pathway-specific wiring. |

---

## 3. Stage 1.5: Empirical Null FDR & Shared-Control Variance Inflation Audit

### 3.1 Why Empirical Null Calibration is Mandatory

In small-$N$ functional proteomics ($n = 1\text{--}2$ per clone/condition), researchers frequently use ad-hoc multi-clone consistency rules (e.g., *"all $k$ shRNAs shift $\ge 1.5\times$ in the same direction relative to a shared control (`NTC`)"* or *"detected in treatment, `NaN` in control"*).

Because all $k$ knockdowns are subtracted against the **same shared control tube** $C$, the fold-changes $D_j$ = $X_j - C$ ($j = 1, \dots, k$) are **not independent**. Under the null hypothesis where $X_j$ and $C$ are independent and identically distributed as $\mathcal{N}(0, \sigma^2)$:

$$
\mathrm{Var}(D_j) = 2\sigma^2, \qquad \mathrm{Cov}(D_j, D_l) = \mathrm{Var}(C) = \sigma^2 \implies \mathrm{Corr}(D_j, D_l) = 0.50
$$

Even for pure Gaussian noise, any random downward fluctuation in the single control tube $C$ simultaneously inflates $D_1$, $D_2$, $\dots$, $D_k$ upward across all $k$ clones.

---

### 3.2 Exact Shared-Control Null Probability Formula

For a fold-change cutoff $c$ = $\log_2(\text{FC cutoff})$ (e.g., $c = 0.585$ for $1.5\times$) and robust background noise standard deviation $\hat{\sigma}$ estimated from negative-control pairs (`VC vs NTC`, `IgG vs Mock`, or replicate splits):

$$
\hat{\sigma} = \frac{1.4826 \cdot \mathrm{MAD}(\Delta_{\mathrm{null}})}{\sqrt{2}}
$$

The exact probability that all $k$ clones exceed cutoff $c$ in the same direction by chance alone is given by integrating over the shared control $C = \hat{\sigma} z$:

$$
P_{\mathrm{null}}\left(\min_{1 \le j \le k} D_j > c\right) = \int_{-\infty}^{\infty} \phi(z) \left[ 1 - \Phi\left(\frac{c + \hat{\sigma} z}{\hat{\sigma}}\right) \right]^k dz
$$

where $\phi(z)$ and $\Phi(z)$ denote the standard normal probability density function (PDF) and cumulative distribution function (CDF).

> [!WARNING]
> **20x–29x False-Positive Inflation Under Naive Independence**
> Assuming the $k$ clone comparisons are statistically independent underestimates the true false-positive rate by **5-fold to 29-fold**. `ProteoSkill` always evaluates the exact shared-control integral above.

---

### 3.3 Empirical Null FDR Metrics (`01_forensic_qc_and_null_audit.py`)

1. **Robust Null Noise ($\hat{\sigma}$)**: Estimated via $1.4826 \cdot \mathrm{MAD}(\Delta_{\mathrm{null}}) / \sqrt{2}$ across designated negative-control or replicate pairs.
2. **Expected Null Hits ($E_{\mathrm{null}}$)**: Computed both via exact Gaussian quadrature of the shared-control integral and via empirical resampling of observed null differences.
3. **Presence/Absence Switch Null Rate**: Number of proteins toggling between `Detected` and `NaN` across negative-control tubes (`VC <-> NTC`).
4. **Empirical Null False Discovery Rate ($\hat{\mathrm{FDR}}_{\mathrm{emp}}$)**:

$$
\hat{\mathrm{FDR}}_{\mathrm{emp}} = \min\left(1.0,\ \frac{E_{\mathrm{null}}(c, k)}{N_{\mathrm{observed}}(c, k)}\right)
$$

---

### 3.4 Automatic Evidence Tiering from Empirical Null FDR

| Empirical Null FDR Range | Assigned Evidence Grade | Reporting & Visual Rule |
| :--- | :--- | :--- |
| **`FDR_emp <= 0.05`** | **`E1 / E2` (Confirmatory / Calibrated Hit)** | Eligible for solid causal arrows in BioRender schematics and headline claims (requires `unique_peptides >= 2`). |
| **`0.05 < FDR_emp <= 0.25`** | **`E2` (Candidate Hit — Moderate Noise)** | Eligible for candidate prioritization when supported by multi-condition concordance or external cohorts. |
| **`FDR_emp > 0.25`** | **`E3` (Exploratory Hypothesis)** | Automatically downgraded to exploratory tier (`dashed arrows` in schematics); if `FDR_emp >= 0.80`, flagged as **Below Shared-Control Noise Floor**. |

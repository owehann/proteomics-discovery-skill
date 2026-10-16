# Reference Playbook 01: Forensic QC Checklist & Empirical Null FDR Calibration

This playbook specifies the mandatory **Stage 0 (Schema & Power Audit)**, **Stage 1 (12-Point Forensic QC)**, and **Stage 1.5 (Empirical Null FDR & Shared-Control Variance Inflation)** rules executed by `scripts/00_forensic_schema_indexer.py` and `scripts/01_forensic_qc_and_null_audit.py`.

---

## 1. Stage 0: Multi-Format Schema & Power Inference

### 1.1 Supported Search Engine & Matrix Formats
- **DIA-NN**: `report.pg_matrix.tsv`, `report.pr_matrix.tsv`, `report.parquet` (Filter `Q.Value <= 0.01`, `PG.Q.Value <= 0.01`, `Lib.Q.Value <= 0.01`).
- **MaxQuant**: `proteinGroups.txt`, `Phospho (STY)Sites.txt`, `GlyGly (K)Sites.txt`, `evidence.txt` (Exclude `Potential contaminant == '+'`, `Reverse == '+'`, `Only identified by site == '+'`; retain `Localization prob >= 0.75` for PTM sites).
- **FragPipe / Philosopher / IonQuant**: `combined_protein.tsv`, `combined_modified_peptide.tsv`, `abundance_protein_MD.tsv`.
- **Proteome Discoverer (PD)**: Excel/TSV exports (`Abundances (Normalized)`, `Abundances (Scaled)`, `# Unique Peptides`, `Master == 'IsMasterProtein'`, `Protein FDR Confidence == 'High'`).
- **Spectronaut**: `PG.Quantity`, `PG.Genes`, `PG.ProteinAccessions`, `EG.PTMAssayProbability`.
- **Olink / SomaScan**: `NPX` (Log2 Normalized Protein Expression) or `RFU` matrices with `LOD` and `Assay_Warning` columns.

### 1.2 Mandatory Identifier & Duplicate Hygiene
1. **Excel Date-Corruption Recovery**: Before any gene-level operation, sanitize strings matching Excel auto-converted dates (`1-Mar` $\rightarrow$ `MARCHF1`, `1-Sep` $\rightarrow$ `SEPTIN1`, `1-Dec` $\rightarrow$ `DELEC1`, etc. using `resources/curated_proteomics_signatures.json`).
2. **Primary Key Policy**: Never aggregate or merge tables on `Gene Symbol` alone. Use composite key `UniProt_Accession | Gene_Symbol` (or `UniProt_Accession | Gene_Symbol | Mod_Residue_Pos` for PTMs).
3. **Isoform Disambiguation**: When multiple protein groups map to the same gene symbol, retain both rows keyed by `Accession|Gene` and mark the canonical / highest-peptide-coverage entry as `is_canonical_rep = True`.

---

## 2. Stage 1: The 12-Point Forensic QC Checklist

Every proteomics dataset must be evaluated against all 12 diagnostic checks before statistical modeling:

| Check ID | Diagnostic Check | Quantitative Threshold / Rule | Automatic Action / Gate 2 Trigger |
| :--- | :--- | :--- | :--- |
| **QC-01** | **Contaminant Audit** | Sum of `cRAP` (`KRT1/2/9/10/14`, `TRY1`, `ALB`, `CSN2`, `IGHG*`, `P02769` BSA, `P00761` pig trypsin, `P22629` streptavidin) $>5\%$ of total intensity or differential ($|\Delta\log_2| > 1.0$) between conditions. | Flag culture serum (`ALB/FETUA`), antibody (`IGHG*`), or streptavidin contamination; remove non-human cRAP rows before normalization. |
| **QC-02** | **Peptide Support Tiering** | Partition rows into **Core Quantitative Tier** (`unique_peptides >= 2` AND `PSM/precursors >= 2`) vs **Single-Peptide Exploratory Tier** (`unique_peptides == 1`). | Never allow `unique_peptides == 1` rows into headline claims, Executive Scorecards, or mechanistic core models without an explicit `[1-Peptide Exploratory]` badge. |
| **QC-03** | **Perturbation Target Audit** | Verify the perturbed target (`Gene_X` in KD/OE/KO/Bait/Drug) in the raw matrix. | If target is undetected (`NaN`), discordant with perturbation direction, or `unique_peptides == 1`, trigger **Gate 2** and cap target-direct claims at `E3` unless orthogonal WB/qPCR is provided. |
| **QC-04** | **Paralog & Isoform Compensation** | Inspect paralogs of perturbed target (e.g., `MAPK3` for `MAPK1`; `RAB27B` for `RAB27A`). | Report compensatory upregulation or shared-peptide cross-interference (off-target shRNA/siRNA vs paralog compensation). |
| **QC-05** | **Tube-Level Anomaly Detection** | Compute sample-wise ID count, median $\log_2$ intensity, IQR, pairwise Pearson $r$, PCA Mahalanobis distance, and Hierarchical Clustering. | Trigger **Gate 2** if any tube has $<65\%$ of group median IDs, PCA $>3\,\text{SD}$, or $r < 0.80$ with biological replicates (e.g., `CTRL_Rep2` with $602$ IDs vs $2,095$). |
| **QC-06** | **Sample Sex / Line Identity Check** | Cross-check sex-linked (`RPS4Y1`, `DDX3Y`, `XIST`) and lineage markers (`ALB`, `AFP`, `EPCAM`, `VIM`) across samples. | Flag tube swaps, cross-line contamination, or clinical metadata mismatches. |
| **QC-07** | **TMT / Isobaric Ratio Compression** | Check `Isolation Interference (%) > 30%`, co-isolation, and bridge-channel CV across plexes. | Recommend `IRS + TMM` + `MSstatsTMT` and warn that fold-changes are compressed toward $0$. |
| **QC-08** | **DIA-NN / Spectronaut PTM & Lib Audit** | Check `PTM.Site.Confidence >= 0.75`, RT calibration drift, and `Lib.Q.Value <= 0.01`. | Separate localized (`Class I >= 0.75`) from unlocalized PTM sites. |
| **QC-09** | **Normalization Suitability Audit** | Compare global median shift across biological groups. If $|\Delta\text{median}| > 0.15\,\log_2$ or $>15\%$ of proteins shift unidirectionally (e.g., EV cargo collapse, AP-MS vs IgG, global translation arrest): | **FORBID naive median/quantile normalization**. Use housekeepers, spike-in standards, bait-normalization, or `TMM` / `VSN`. |
| **QC-10** | **MISEV2023 EV Purity Audit** | (EV / Secretome only) Quantify Cat 1 (`CD9/63/81`, `TSG101`, `SDCBP`, `PDCD6IP`, `FLOT1`) vs Cat 3 (`CANX`, `HSP90B1`, `GM130/GOLGA2`, `HSPA5`, `CYC1`, `HIST1H*`, `ALB`, `APOA1`). | Report $\log_2(\text{Cat 1+2}/\text{Cat 3})$ purity index and flag co-isolated ER/Golgi/nuclear debris. |
| **QC-11** | **Missingness Mechanism (MNAR vs MAR)** | Fit the **Detection Probability Curve (DPC)**: $\text{logit}\,P(\text{observed}_{ij}) = \beta_0 + \beta_1 \bar{y}_i$. | Estimate $\hat{\beta}_1$. Use DPC (`limpa`) or Hurdle (`msqrob2`) rather than random left-censored (`MinProb`) imputation when $N_{\text{rep}} < 3$. |
| **QC-12** | **Confounder Signature Scan** | Score cell-cycle (`MKI67`, `TOP2A`, `PCNA`, `CDK1`), stress/UPR (`HSPA5`, `DDIT3`), and IFN (`ISG15`, `MX1`, `STAT1`) indices per sample. | Alert when downstream enrichment is dominated by proliferation or stress arrest rather than specific pathway wiring. |

---

## 3. Stage 1.5: Empirical Null FDR & Shared-Control Variance Inflation Audit

### 3.1 Why Empirical Null Calibration is Mandatory
In small-$N$ functional proteomics ($n=1\text{–}2$ per clone/condition), researchers frequently use ad-hoc multi-clone consistency rules (e.g., *"all $k$ shRNAs shift $\ge 1.5\times$ in the same direction relative to a shared control (`NTC`)"* or *"detected in treatment, `NaN` in control"*).

Because all $k$ knockdowns are subtracted against the **same shared control tube** $C$, the fold-changes $D_j = X_j - C$ ($j=1,\dots,k$) are **not independent**! Under the null hypothesis $X_j, C \stackrel{\text{i.i.d.}}{\sim} \mathcal{N}(0, \sigma^2)$:
$$\text{Var}(D_j) = 2\sigma^2, \qquad \text{Cov}(D_j, D_l) = \text{Var}(C) = \sigma^2 \implies \text{Corr}(D_j, D_l) = 0.50$$
Even for pure Gaussian noise, any random downward fluctuation in the single control tube $C$ simultaneously inflates $D_1, D_2, \dots, D_k$ upward across all $k$ clones!

### 3.2 Exact Shared-Control Null Probability Formula
For a fold-change cutoff $c = \log_2(\text{FC}_{\text{cut}})$ (e.g., $c = 0.585$ for $1.5\times$) and robust background noise SD $\hat{\sigma} = \text{MAD}(\Delta_{\text{null}}) / (1.4826 \sqrt{2})$ estimated from negative-control pairs (`VC vs NTC`, `IgG vs Mock`, or replicate splits):
$$P_{\text{null}}\left(\min_{1\le j \le k} D_j > c\right) = \int_{-\infty}^{\infty} \phi(z) \left[ 1 - \Phi\left(\frac{c + \hat{\sigma} z}{\hat{\sigma}}\right) \right]^k dz$$
where $\phi$ and $\Phi$ are the standard normal PDF and CDF (integrated over the shared control $C = \hat{\sigma} z$).

By contrast, if someone mistakenly assumed the $k$ comparisons were independent, they would underestimate the false-positive rate by $5\times\text{ to }25\times$!

### 3.3 Empirical Null FDR Metrics Computed by `01_forensic_qc_and_null_audit.py`
1. **Robust Null Noise $\hat{\sigma}$**: Estimated via $\text{MAD} / (1.4826 \sqrt{2})$ across all designated negative-control or replicate pairs.
2. **Expected Null Hits $E_{\text{null}}(c, k)$**: Computed both via exact Gaussian quadrature of the shared-control integral and via empirical resampling of observed null differences.
3. **Presence/Absence Switch Null Rate $E_{\text{switch, null}}$**: Number of proteins toggling between `Detected` and `NaN` across negative-control tubes (`VC <-> NTC`).
4. **Empirical Null FDR**:
   $$\widehat{\text{FDR}}_{\text{emp}} = \min\left(1.0, \; \frac{E_{\text{null}}(c, k)}{N_{\text{observed}}(c, k)}\right)$$

### 3.4 Automatic Evidence Tiering from Empirical Null FDR
- $\widehat{\text{FDR}}_{\text{emp}} \le 0.05$: Eligible for **`E1` / `E2` (Confirmatory / Moderated Hit)** if `unique_peptides >= 2`.
- $0.05 < \widehat{\text{FDR}}_{\text{emp}} \le 0.25$: **`E2` (Candidate Hit — moderate noise background)**.
- $\widehat{\text{FDR}}_{\text{emp}} > 0.25$: Automatically downgraded to **`E3 (Exploratory Hypothesis)`**. If $\widehat{\text{FDR}}_{\text{emp}} \ge 0.80$ (i.e., $N_{\text{observed}} \le E_{\text{null}}$), the contrast is flagged as **`Below Shared-Control Noise Floor`** and cannot be presented as a biological discovery without orthogonal validation.

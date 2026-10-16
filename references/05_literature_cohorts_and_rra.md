# Reference Playbook 05: Live Literature Verification, Clinical Cohorts, & Orthogonalized RRA

This playbook governs **Stage 5 (`scripts/05_literature_cohorts_and_rra.py`)**, enforcing live citation verification, negative trial awareness, multi-cohort validation, and Orthogonalized Robust Rank Aggregation (`RRA`).

---

## 1. Live Literature Verification & Anti-Hallucination Policy

1. **Zero Fabricated Citations**: Every `PMID` or `DOI` cited in a report or slide deck must be verified via the NCBI Entrez E-Utilities API (`esearch.fcgi` / `esummary.fcgi`) or EuropePMC REST API (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`) and cached in `facts.json`.
2. **Bidirectional Clinical Trial Query**: Whenever a therapeutic target or metabolic vulnerability is nominated (e.g., `IDO1` tryptophan catabolism $\rightarrow$ epacadostat; `ERBB2` $\rightarrow$ ADCs/bispecifics), the literature engine **MUST** query both positive preclinical studies AND **negative/failed Phase II/III clinical trials** (e.g., ECHO-301 Phase III trial of epacadostat in melanoma, Long et al. *Lancet Oncol* 2019, `PMID: 31026679`).
3. **Mechanistic Novelty Classification**:
   - **`Known Canonical Target`**: Direct literature evidence in the same tissue/disease and modality.
   - **`Cross-Lineage / Cross-Compartment Extension`**: Known in another cancer/organelle, newly demonstrated in the current compartment.
   - **`Novel Unreported Axis`**: Zero prior co-mentions of `(Perturbation Target, Candidate Cargo/Substrate)` in PubMed/EuropePMC.

---

## 2. External Multi-Omics & Clinical Cohort Alignment (Gate 3)

Depending on the organ/disease preset in `analysis_spec.yaml`, cross-validate top `E1/E2` candidates against:
- **CPTAC / ICPC Pan-Cancer Proteomics**: Paired Tumor-vs-NAT protein abundance and phosphoproteome concordance ($\%$ directional agreement + Spearman $\rho$).
- **TCGA / ICGC / GEO Transcriptomics**: RNA-vs-Protein coupling ($r_{\text{RNA, Prot}}$) to distinguish transcriptional induction from post-transcriptional protein stabilization or secretion.
- **DepMap CRISPR/RNAi Essentiality**: Lineage-specific Chronos dependency scores.
- **Human Protein Atlas (`HPA`), UniProt Subcellular Location, Vesiclepedia / ExoCarta, & ChEMBL / OpenTargets**: Secretome classification, druggability tractability (Clinical / Tool Compound / Biologic Accessible / Undruggable).

*If external cohort files are absent locally and cannot be fetched automatically via open API, trigger **Gate 3** with exact download links and target file paths.*

---

## 3. Orthogonalized Robust Rank Aggregation (`RRA`)

### 3.1 Avoiding Double-Counting Collinear Metrics
A common mistake in multi-stage proteomics integration is feeding 5–8 collinear metrics (e.g., $\log_2\text{FC}$, moderated $p$-value, ROC AUC, ElasticNet frequency, Random Forest importance—all derived from the exact same differential contrast!) into an unweighted rank average. That artificially inflates fold-change outliers while drowning out orthogonal biological evidence.

### 3.2 Orthogonal Evidence Pillars
Group all candidate evidence into **strictly orthogonal pillars** before aggregation:
1. **Pillar 1 — Internal Quantitative Effect & Statistical Precision**: Moderated effect size & empirical-null calibrated FDR (plus `R-Score` in reciprocal `KD x OE` or `SGI` in multi-compartment designs).
2. **Pillar 2 — Multi-Clone / Multi-Cohort / Prognostic Reproducibility**: External cohort concordance (e.g., CPTAC/TCGA HR) or independent biological clone agreement.
3. **Pillar 3 — Systems & Topological Centrality**: CORUM complex membership, degree-corrected `RWR` / `KDA` bridge score, or PTM kinase/E3 substrate specificity.
4. **Pillar 4 — Biophysical / Translational Tractability**: Peptide support (`unique_peptides >= 2`), surface/extracellular accessibility (`HPA`/`Vesiclepedia`), or small-molecule/biologic tractability.

### 3.3 Beta-Order Statistic RRA + Bootstrap Confidence Intervals
1. Within each pillar $m \in \{1,\dots,M\}$, convert raw metrics into normalized rank percentiles $u_{i, m} \in (0, 1]$.
2. Compute the Kolde et al. (2012) **Robust Rank Aggregation (`RRA`)** $\rho$-score from the minimum order statistics under a uniform $[0,1]$ null:
   $$\rho_i = \min_{1 \le k \le M} P\left(U_{(k)} \le u_{i, (k)}\right) = \min_{1 \le k \le M} \sum_{j=k}^M \binom{M}{j} u_{i, (k)}^j (1 - u_{i, (k)})^{M - j}$$
   adjusted via Bonferroni/Benjamini-Hochberg across $M$ pillars.
3. **Dirichlet Weight Sensitivity Bootstrap ($B = 1,000$)**: Draw pillar weights $\mathbf{w}^{(b)} \sim \text{Dirichlet}(\alpha_1, \dots, \alpha_M)$ across $1,000$ iterations and compute each protein's **Top-$K$ Rank Stability Frequency** ($\text{Stability}_i \in [0\%, 100\%]$) and $95\%$ rank confidence interval. Only proteins with $\text{Stability}_i \ge 80\%$ and `unique_peptides >= 2` are crowned **`Tier-1 Core Drivers`**.

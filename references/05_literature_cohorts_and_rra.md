<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 05: Live Literature Verification, Clinical Cohorts, & 5-Layer Multi-Lens RRA

> [!IMPORTANT]
> **Stage Scope (`Stage 5` — Engines `[09]` through `[12]`)**
> This playbook governs [`scripts/05_literature_cohorts_and_rra.py`](../scripts/05_literature_cohorts_and_rra.py), enforcing live citation verification, negative clinical trial awareness, multi-cohort validation (`CPTAC`, `TCGA`, biofluid `SGI`, `DepMap`, `HPA`), and 5-Layer Multi-Lens Robust Rank Aggregation (`RRA`).

---

## 1. Live Literature Verification & Anti-Hallucination Policy

1. **Zero Fabricated Citations**: Every `PMID` or `DOI` cited in a report or slide deck must be verified via the NCBI Entrez E-Utilities API (`esearch.fcgi` / `esummary.fcgi`) or EuropePMC REST API (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`) and cached in `verified_citations.json`.
2. **Bidirectional Clinical Trial Query (`Engine [12]`)**: Whenever a therapeutic target or metabolic vulnerability is nominated (e.g., `IDO1` tryptophan catabolism → epacadostat; `ERBB2` → ADCs/bispecifics), the literature engine queries both positive preclinical studies **and** negative or terminated Phase II/III clinical trials (e.g., ECHO-301 Phase III trial of epacadostat in melanoma, Long et al., *Lancet Oncol* 2019, `PMID: 31026679`).
3. **Mechanistic Novelty Classification**:
   - **`Known Canonical Target`**: Direct literature evidence in the same tissue/disease and modality.
   - **`Cross-Lineage / Cross-Compartment Extension`**: Known in another cancer or organelle, newly demonstrated in the current compartment.
   - **`Novel Unreported Axis`**: Zero prior co-mentions of `(Perturbation Target, Candidate Cargo/Substrate)` in PubMed/EuropePMC.

---

## 2. External Multi-Omics & Clinical Cohort Alignment (`Engines [09] & [10]`)

Depending on the organ/disease preset in `analysis_spec.yaml`, top `E1/E2` candidates are cross-validated against:

| External Resource | Quantitative Metrics Extracted | Biological & Translational Purpose |
| :--- | :--- | :--- |
| **CPTAC / ICPC Pan-Cancer Proteomics** | Paired Tumor-vs-NAT protein $\log_2\mathrm{FC}$, paired $t$-statistic, and $1,000\times$ bootstrap diagnostic `ROC-AUC` | Confirms in-vivo human tumor dysregulation and clinical biomarker performance |
| **TCGA / Paired Transcriptomics (`Engine [09]`)** | Paired mRNA $\log_2\mathrm{FC}$ vs protein $\log_2\mathrm{FC}$ (`RNA–Protein Coupling Quadrant`) | Separates transcriptional induction from post-translational stabilization, degradation, or vesicular export |
| **Paired Clinical Biofluid / Secretome (`SGI`)** | Secretion Gradient Index ($\mathrm{SGI} = \Delta\mathrm{Fluid} - \Delta\mathrm{Tissue}$) | Identifies active tumor-to-fluid hyper-secretion vs passive tissue overflow |
| **DepMap CRISPR Chronos (`Engine [10]`)** | Lineage-specific cell-autonomous dependency count ($n$ dependent lines) | Distinguishes **non-autonomous paracrine/TME remodeling cargos** (`0` dependent lines in 2D monoculture) from **cell-autonomous survival hubs** |
| **HPA, UniProt, Vesiclepedia & ChEMBL (`Engine [12]`)** | Plasma/secretome detectability, subcellular accessibility, and approved/clinical drug status | Guides liquid-biopsy assay design and pharmacological tractability |

---

## 3. Orthogonalized 5-Layer Multi-Lens Robust Rank Aggregation (`Engine [11]`)

### 3.1 Avoiding Double-Counting Collinear Metrics

A common flaw in multi-stage proteomics integration is feeding 5–8 collinear metrics (e.g., $\log_2\mathrm{FC}$, moderated $p$-value, single-dataset ROC AUC, and Lasso frequency—all derived from the exact same differential contrast) into an unweighted rank average. That artificially inflates fold-change outliers while drowning out orthogonal biological evidence.

### 3.2 Five Strictly Orthogonal Evidence Layers (`L1`–`L5`)

`ProteoSkill` groups all candidate evidence into **5 orthogonal evidence layers** normalized to $[0, 1]$ before aggregation:

| Layer | Evidence Dimension | Underlying Quantitative Input |
| :--- | :--- | :--- |
| **`L1_MS_Perturbation`** | Internal MS Perturbation & Rescue Strength | Moderated effect size, empirical-null calibrated FDR, Track B `<LOD` switch tier, and `KD × OE` Reciprocity Score (`R-Score`) |
| **`L2_Pathway_Topology`** | Biophysical & Pathway Module Specificity | De-duplicated Jaccard pathway membership, CORUM complex co-shift, and UniProt membrane/secretory topology enrichment |
| **`L3_Network_Hub`** | PPI Network Centrality & Lock-and-Key Bridge | Degree-binned `STRING` Random Walk with Restart (`RWR`) score ($Z_{\mathrm{deg\text{-}null}}$) and Key Driver (`KDA`) bridge status |
| **`L4_Clinical_Cohort`** | Human Clinical Cohort & Biofluid Concordance | `CPTAC` paired Tumor-vs-Normal protein concordance and clinical biofluid Secretion Gradient Index (`SGI`) |
| **`L5_Translational`** | Druggability, Liquid-Biopsy & CRISPR Tractability | Peptide support (`unique_peptides >= 2`), `HPA` plasma detectability, `ChEMBL` tractability, and `DepMap` functional profile |

---

### 3.3 Kolde Beta-Order Statistic RRA + $1,000\times$ Dirichlet Stability Bootstrap

1. **Normalized Rank Percentiles**: Within each evidence layer $m \in \{1, \dots, M\}$ ($M = 5$), raw scores are converted into rank percentiles $u_{i, m} \in (0, 1]$.
2. **Order-Statistic Beta Probability (Kolde et al., *Bioinformatics* 2012)**: From the sorted rank percentiles $u_{i, (1)}$, $u_{i, (2)}$, $\dots$, $u_{i, (M)}$, compute the exact binomial tail probability under a uniform $[0, 1]$ null:

$$
\rho_i = \min_{1 \le k \le M} P\!\left(U_{(k)} \le u_{i, (k)}\right) = \min_{1 \le k \le M} \sum_{j=k}^{M} \binom{M}{j} u_{i, (k)}^j \left(1 - u_{i, (k)}\right)^{M - j}
$$

3. **Three-Lens Multi-Perspective Prioritization (`3 × Top 12 = 36 Core Targets` → `Top 6 Flagship Heroes`)**:
   - **Perspective I (`Shared Co-Regulated Core`)**: Top 12 targets driving the convergent epistatic core (`G1 / G2`).
   - **Perspective II (`Perturbation-A Specific & Reciprocal`)**: Top 12 targets specific to Perturbation A (`G3 / G4`) or reciprocal rescue.
   - **Perspective III (`Perturbation-B Specific & Seesaw Flux`)**: Top 12 targets specific to Perturbation B (`G5 / G6`) or antagonistic Seesaw (`G7`) routing.
4. **Dirichlet Weight Sensitivity Bootstrap ($B = 1,000$)**: Draw layer weight vectors $\mathbf{w}^{(b)}$ from a 5-dimensional symmetric Dirichlet distribution across $1,000$ iterations:

$$
S_i^{(b)} = \sum_{m=1}^{5} w_m^{(b)} L_{i, m}
$$

   Only proteins achieving high Top-$K$ Bootstrap Stability ($\ge 80\%$) and multi-peptide support (`unique_peptides >= 2`) are promoted to **`Top-6 Flagship Validation Targets`**.

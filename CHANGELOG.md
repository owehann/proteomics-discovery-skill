# Changelog

All notable changes to `proteomics-discovery-skill` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.1.0] - 2026-10-16

### Added & Improved
- **Calibrated Statistical Fixes**:
  - Corrected Stage 1.5 empirical-null noise estimator to $\hat{\sigma} = 1.4826 \cdot \mathrm{MAD}(\Delta_{\mathrm{null}}) / \sqrt{2}$ (`scripts/01_forensic_qc_and_null_audit.py` and `references/01_forensic_qc_and_null_calibration.md`).
  - Corrected prior variance sign in Python `squeeze_var_ebayes` (`scripts/02_statistical_engine_py.py`), added `s2_raw`/`s2_post` export columns, and required `sigma_null_override` for unreplicated `1v1` designs.
  - Added `limma::eBayes` call prior to `DEqMS::spectraCounteBayes` in `scripts/02_statistical_engine.R`.
  - Updated PTM stoichiometry decoupling (`scripts/03_archetype_and_modality_solver.py`) to classify missing parent proteins as `Parent_Unquantified` rather than imputing `0.0`.
  - Eliminated synthetic fallback numbers in clinical cohort integration (`scripts/05_literature_cohorts_and_rra.py`), recording `NaN` with explicit `clinical_data_source = "missing"` when external cohort tables are absent.
- **Pre-Delivery Fact & Citation Gate (`scripts/09_fact_and_citation_gate.py`)**:
  - Added signed gene–metric reconciliation ($\pm 0.02$ tolerance), live EuropePMC REST API PMID/DOI verification, E3 causal-verb overclaim linting, and native `.pptx`/`.svg` text inspection.
- **Biologist-First Visual & Presentation Suite (`scripts/06`–`08`)**:
  - Added `generate_all_publication_figures()` orchestrator (`scripts/06_npg_statistical_figures.py`), 7-module BioRender vector suite (`s01`, `s03`, `s14`, `s17`, `s24`, `s25`, `s26`), and penultimate literature/pathway summary slide verification (`scripts/08_editorial_pptx_builder.py`).

## [4.0.0] - 2026-10-16

### Added
- **7-Stage End-to-End Proteomics Discovery Architecture (`SKILL.md`)** supporting **6 Proteomics Modalities (`M1`–`M6`)** and **10 Experimental Archetypes (`A1`–`A10`)**.
- **Stage 0 Zero-Assumption Schema Indexer (`scripts/00_forensic_schema_indexer.py`)** with automatic detection of DIA-NN, Spectronaut, MaxQuant, FragPipe, Proteome Discoverer, and Olink NPX exports, plus interactive Gate 1 Markdown confirmation generation.
- **Stage 1 & 1.5 Sample & Data QC & Empirical Null Calibration (`scripts/01_forensic_qc_and_null_audit.py`)**:
  - Automatic recovery of 35+ Excel date-corrupted HGNC symbols (`MARCHF1`, `SEPTIN1`, `DELEC1`, etc.).
  - Composite key construction (`Gene_Symbol|UniProt_Accession|PTM_Site`).
  - Automatic flagging and quarantine of exogenous spike-ins (e.g., therapeutic IgG chains `IGHG1`/`IGKC`), bovine FBS carry-over, and keratins.
  - MISEV2023 purity auditing for EV/secretome/biofluid datasets.
  - 4-metric composite outlier tube detector (ID drop, median shift, pairwise correlation, and PCA Mahalanobis distance).
  - Stage 1.5 analytical & permutation-based **Shared-Control Empirical Null FDR Calibration** and **`E1`–`E3` Evidence Grading**.
- **Stage 2 Dual Statistical Engine (`scripts/02_statistical_engine.R` & `scripts/02_statistical_engine_py.py`)**:
  - `limpa` Detection Probability Curve (DPC) and `DEqMS` peptide-stratified empirical Bayes variance moderation (`unique_peptides >= 2` Tier-1 Core vs `1-peptide` Tier-2 Exploratory).
  - Two One-Sided Tests (`TOST`) equivalence testing for proving functional/humanization conservation.
  - Specialized solvers for `AP-MS` interactome specificity scoring, `NPARC` (TPP/CETSA thermal stability), and `CurveCurator` (dose/time-course $\mathrm{pEC}_{50}$).
- **Stage 3–5 Biological Mechanism & Literature Triangulation (`scripts/03`–`05`)**:
  - `KD × OE` 9-sector reciprocal rescue vs hairpin divergence solver and paralog compensation scanner.
  - Paired fluid-vs-tissue Secretion Gradient Index (`SGI`) 4-quadrant classifier.
  - `CORUM` complex stoichiometry & orphan-subunit degradation, `ULM` signed regulator Z-scores, and `UniProt` ectodomain shedding topology mapper.
  - Live `EuropePMC`/`PubMed`/`UniProt`/`OpenTargets` REST API verification, `Robust Rank Aggregation (RRA)`, and negative clinical trial guardrails.
- **Stage 6–7B Visual & Presentation Engines (`scripts/06`–`09`)**:
  - NPG statistical figure suite (300-DPI PNG + editable PDF/SVG).
  - Programmatic BioRender-style vector SVG schematic engine (`scripts/07_biorender_canvas_engine.py`) with evidence-graded stroke styling (`E1`/`E2`/`E3`).
  - 16:9 Editorial PPTX generator (`scripts/08_editorial_pptx_builder.py`).
  - Mandatory Stage 7B Pre-Delivery Fact, Peptide-Tier, Overclaim, & Citation Gate (`scripts/09_fact_and_citation_gate.py`).
- **Automated 8-Suite Validation Harness (`tests/run_validation_suite.py`)**.

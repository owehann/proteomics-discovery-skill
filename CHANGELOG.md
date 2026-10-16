# Changelog

All notable changes to `proteomics-discovery-skill` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-16

### Added
- **7-Stage End-to-End Proteomics Discovery Architecture (`SKILL.md`)** supporting **6 Proteomics Modalities (`M1`–`M6`)** and **10 Experimental Archetypes (`A1`–`A10`)**.
- **Stage 0 Zero-Assumption Schema Indexer (`scripts/00_forensic_schema_indexer.py`)** with automatic detection of DIA-NN, Spectronaut, MaxQuant, FragPipe, Proteome Discoverer, and Olink NPX exports, plus interactive Gate 1 Markdown confirmation generation.
- **Stage 1 & 1.5 Forensic QC & Empirical Null Calibration (`scripts/01_forensic_qc_and_null_audit.py`)**:
  - Automatic recovery of 35+ Excel date-corrupted HGNC symbols (`MARCHF1`, `SEPTIN1`, `DELEC1`, etc.).
  - Composite key construction (`Gene_Symbol|UniProt_Accession|PTM_Site`).
  - Automatic flagging and quarantine of exogenous spike-ins (e.g., therapeutic IgG chains `IGHG1`/`IGKC`), bovine FBS carry-over, and keratins.
  - MISEV2023 purity auditing for EV/secretome/TIF datasets.
  - 4-metric composite outlier tube detector (ID drop, median shift, pairwise correlation, and PCA Mahalanobis distance).
  - Stage 1.5 analytical & permutation-based **Shared-Control Empirical Null FDR Calibration** and **`E1`–`E4` Evidence Grading**.
- **Stage 2 Dual Statistical Engine (`scripts/02_statistical_engine.R` & `scripts/02_statistical_engine_py.py`)**:
  - `limpa` Detection Probability Curve (DPC) and `DEqMS` peptide-stratified empirical Bayes variance moderation (`unique_peptides >= 2` Tier-1 Core vs `1-peptide` Tier-2 Exploratory).
  - Two One-Sided Tests (`TOST`) equivalence testing for proving functional/humanization conservation.
  - Specialized solvers for `SAINT`-lite (AP-MS), `NPARC` (TPP/CETSA thermal stability), and `CurveCurator` (dose/time-course $\text{pEC}_{50}$).
- **Stage 3–5 Biological Mechanism & Literature Triangulation (`scripts/03`–`05`)**:
  - `KD × OE` 9-sector reciprocal rescue vs hairpin divergence solver and paralog compensation scanner.
  - Paired fluid-vs-tissue Secretion Gradient Index (`SGI`) 4-quadrant classifier.
  - `CORUM` complex stoichiometry & orphan-subunit degradation, `ULM` signed regulator Z-scores, and `UniProt` ectodomain shedding topology mapper.
  - Live `EuropePMC`/`PubMed`/`UniProt`/`OpenTargets` REST API verification, `Robust Rank Aggregation (RRA)`, and negative clinical trial guardrails.
- **Stage 6–7B Visual & Presentation Engines (`scripts/06`–`09`)**:
  - NPG statistical figure suite (300-DPI PNG + editable SVG).
  - Programmatic BioRender-style vector SVG schematic engine (`scripts/07_biorender_canvas_engine.py`) with evidence-graded stroke styling (`E1`/`E2`/`E3`).
  - 16:9 Editorial PPTX generator (`scripts/08_editorial_pptx_builder.py`).
  - Mandatory Stage 7B Pre-Delivery Fact, Peptide-Tier, Overclaim, & Citation Gate (`scripts/09_fact_and_citation_gate.py`).
- **Automated 5-Suite Validation Harness (`tests/run_validation_suite.py`)**.

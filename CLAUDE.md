# ProteoSkill — Agent Directives (`CLAUDE.md`)

When operating inside this repository or when `proteomics-discovery-skill` is invoked on a proteomics dataset (`.xlsx`, `.tsv`, `.csv`, `.txt`, `.parquet`), always load and follow [`SKILL.md`](./SKILL.md).

## Quick Execution Rules
1. **Gate 1 First (Never Guess Sample Groups Silently)**:
   - Run `python3 scripts/00_forensic_schema_indexer.py --input <path>` first.
   - Present the inferred experimental design, sample sizes, and Evidence Grade ceiling (`E1–E3`) to the user and wait for confirmation before running downstream steps.
2. **12 Deep Mechanistic & Clinical Engines**:
   - Execute the 7-step pipeline sequentially (`scripts/00` through `scripts/09`), covering quality control, empirical null FDR calibration, stratified `DEqMS`/`limma` + `TOST` equivalence, paralog compensation, Jaccard pathway clustering, `CORUM` complex stoichiometry, degree-corrected `STRING` PPI hubs, subcellular & ectodomain shedding topology, upstream kinase/TF/E3 `ULM` inference, RNA-protein coupling, `CPTAC`/`HPA` clinical cohort validation, patient subtyping + bootstrap `ROC-AUC`, and `ChEMBL`/`ClinicalTrials.gov` tractability.
3. **Editorial PPTX & Fact Gate**:
   - Build slides exclusively via `scripts/08_editorial_pptx_builder.py` (airy white header, `<= 48 char` action titles, unboxed `stat()` KPI callouts, 6 alternating layouts `A–F`, native editable tables/roadmaps, and greeting-free bilingual notes).
   - Run `python3 scripts/09_fact_and_citation_gate.py` before delivering final outputs.
4. **Windows Fast-Track Channel (`win64` / `x64` / `ARM64` / `win32`)**:
   - If the host OS is Windows (`platform.system() == "Windows"` or PowerShell/CMD on 64-bit `win64`/`x64`, `ARM64`, or `win32`), immediately read and follow [`references/windows_setup.md`](./references/windows_setup.md) for the 30-second binary-only installation and execution path (skipping `conda`/`environment.yml` and `libcairo` dependencies).

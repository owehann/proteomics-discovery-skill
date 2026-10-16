# Contributing to `proteomics-discovery-skill`

Thank you for your interest in contributing to **`proteomics-discovery-skill`**!

## How to Contribute

1. **Fork & Branch**: Fork the repository and create a feature branch (`git checkout -b feat/your-feature`).
2. **Maintain Zero External Skill Dependencies**: All statistical, visualization, and verification capabilities must remain self-contained inside `scripts/`, `references/`, and `resources/`.
3. **Preserve Strict Statistical & Fact-Checking Guardrails**:
   - Any new differential abundance solver must propagate `unique_peptides` and support `E1`–`E4` evidence grading.
   - Any new visualization or report template must pass `scripts/09_fact_and_citation_gate.py` with `0` critical errors.
4. **Run the Validation Suite**:
   Before opening a Pull Request, verify that all 5 end-to-end validation test suites pass:
   ```bash
   python3 tests/run_validation_suite.py
   ```
5. **Adding Organ/Disease Signatures**:
   To add new tissue/organ marker panels or contaminant signatures, edit `resources/curated_proteomics_signatures.json` and include the primary literature `PMID` or consensus guideline reference.

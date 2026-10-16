#!/usr/bin/env python3
"""
Validation & Regression Test Suite for `ProteoSkill` (`proteomics-discovery-skill`)
===================================================================================
Executes 5 rigorous end-to-end synthetic tests:
  1. Stage 0 Excel Date-Corruption Sanitizer & Schema Inference Test
     (`1-Mar` -> `MARCHF1`, `1-Sep` -> `SEPTIN1`, composite ID construction).
  2. Stage 1.5 Shared-Control Variance Inflation & Empirical Null Integral Test
     (Verifies that shared-control k-clone covariance `Cov(X1-C, X2-C) = sigma^2`
     is accurately captured vs naive independence).
  3. Stage 2 Analytical Empirical Bayes (`squeezeVar`), DPC Missingness, & TOST Equivalence
     Synthetic Ground-Truth Test.
  4. Stage 3 Archetype Solver Test (verifying synthetic `MARKER_A` with `tissue_log2FC = 0.62 >= 0.585`
     is strictly classified into `Q1_concordant_up`, never `Q2_hyper_secreted`).
  5. Stage 7B Fact Gate Adversarial Injection Test (verifies that `09_fact_and_citation_gate.py`
     catches injected `CTRL_Rep2` outlier leaks, `Q2_hyper_secreted` misclassifications,
     and untagged `unique_peptides == 1` promotions).
"""

from __future__ import annotations
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


SKILL_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = SKILL_ROOT / "scripts"


def _load_script_module(name: str, filename: str):
    path = SCRIPTS_DIR / filename
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_1_stage0_sanitizer_and_schema() -> None:
    s0 = _load_script_module("s0", "00_forensic_schema_indexer.py")
    assert s0.sanitize_gene_symbol("1-Mar") == "MARCHF1"
    assert s0.sanitize_gene_symbol("02-Sep") == "SEPTIN2"
    assert s0.sanitize_gene_symbol("1-Dec") == "DELEC1"
    assert s0.sanitize_gene_symbol("TP53;TP53TG1") == "TP53"
    print("  [PASS] Test 1: Stage 0 HGNC Excel Date-Corruption Recovery & Composite Key Builder")


def test_2_stage15_shared_control_null_integral() -> None:
    s1 = _load_script_module("s1", "01_forensic_qc_and_null_audit.py")
    res = s1.shared_control_null_probability(sigma=0.269, k_clones=3, fc_cutoff_log2=0.5849625)
    assert res["p_one_sided_shared_control"] > 0.007
    assert res["variance_inflation_ratio"] > 20.0
    print(
        f"  [PASS] Test 2: Stage 1.5 Shared-Control Null Integral "
        f"(p_shared={res['p_one_sided_shared_control']:.5f} vs p_naive={res['p_one_sided_naive_independent']:.5f}, "
        f"Inflation={res['variance_inflation_ratio']}x)"
    )


def test_3_stage2_statistical_engines() -> None:
    s2 = _load_script_module("s2", "02_statistical_engine_py.py")
    rng = np.random.default_rng(123)
    n_prot = 500
    ctrl = rng.normal(20.0, 0.20, size=(n_prot, 4))
    treat = rng.normal(20.0, 0.20, size=(n_prot, 4))
    treat[:30, :] += 1.8

    df = pd.DataFrame(
        {
            "composite_id": [f"P{i:04d}|GENE{i}" for i in range(n_prot)],
            "std_accession": [f"P{i:04d}" for i in range(n_prot)],
            "std_gene": [f"GENE{i}" for i in range(n_prot)],
            "std_unique_peptides": rng.integers(1, 12, size=n_prot),
        }
    )
    for j in range(4):
        df[f"T_{j+1}"] = treat[:, j]
        df[f"C_{j+1}"] = ctrl[:, j]

    res = s2.run_moderated_differential_analysis(
        df, [f"T_{j+1}" for j in range(4)], [f"C_{j+1}" for j in range(4)], delta_eq_log2=0.50
    )
    sig_hits = set(res.loc[res["fdr_q_value"] < 0.05, "std_gene"])
    true_positives = len(sig_hits & {f"GENE{i}" for i in range(30)})
    false_positives = len(sig_hits - {f"GENE{i}" for i in range(30)})
    assert true_positives == 30
    assert false_positives <= 3
    assert int(res["is_tost_equivalent"].iloc[30:].sum()) > 300
    print(
        f"  [PASS] Test 3: Stage 2 Moderated eBayes DE + TOST Equivalence "
        f"(TP={true_positives}/30, FP={false_positives}, TOST_Equiv={int(res['is_tost_equivalent'].sum())})"
    )


def test_4_stage3_archetype_solvers() -> None:
    s3 = _load_script_module("s3", "03_archetype_and_modality_solver.py")
    df_sgi = pd.DataFrame(
        {
            "std_gene": ["MARKER_A", "MARKER_B", "MARKER_C"],
            "fluid_log2FC": [2.40, 3.10, 2.15],
            "tissue_log2FC": [0.62, 0.20, 2.50],
            "std_unique_peptides": [1, 1, 12],
        }
    )
    out = s3.solve_archetype_a5_clinical_sgi(df_sgi, "fluid_log2FC", "tissue_log2FC")
    row_a = out.loc[out["std_gene"] == "MARKER_A"].iloc[0]
    row_b = out.loc[out["std_gene"] == "MARKER_B"].iloc[0]
    row_c = out.loc[out["std_gene"] == "MARKER_C"].iloc[0]

    assert row_a["sgi_quadrant"] == "Q1_concordant_up"
    assert "1-Peptide_Exploratory" in row_a["sgi_clinical_tier"]
    assert row_b["sgi_quadrant"] == "Q2_hyper_secreted"
    assert "Core_ge2_peptides" in row_c["sgi_clinical_tier"]
    print("  [PASS] Test 4: Stage 3 Multi-Compartment SGI Quadrant Boundary & Single-Peptide Tiering")


def test_5_stage7b_adversarial_fact_gate() -> None:
    s9 = _load_script_module("s9", "09_fact_and_citation_gate.py")
    facts = {
        "genes": {
            "CARGO_X": {
                "unique_peptides": 5,
                "metrics": {"CTRL_Rep1": 18.42, "TRT_Rep1": 18.45},
                "categories": {},
                "quarantined_outlier_values": [19.88],
            },
            "MARKER_A": {
                "unique_peptides": 1,
                "metrics": {"fluid_log2FC": 2.40, "tissue_log2FC": 0.62},
                "categories": {"quadrant": "Q1_concordant_up"},
                "quarantined_outlier_values": [],
            },
            "CARGO_Y": {
                "unique_peptides": 1,
                "metrics": {"r_score": 1.15},
                "categories": {},
                "quarantined_outlier_values": [],
            },
        }
    }

    bad_md = """# Flawed Draft Report
1. CARGO_X was strongly enriched in TRT_Rep1 (19.88) compared to CTRL_Rep1 (18.42).
2. MARKER_A is a top Q2_hyper_secreted target with massive fluid enrichment.
3. CARGO_Y is our #1 core reciprocal cargo driving vesicle biogenesis.
"""
    good_md = """# Calibrated Report
1. CARGO_X was detected in CTRL_Rep1 (18.42) and TRT_Rep1 (18.45), whereas 19.88 occurred in the excluded outlier tube CTRL_Rep2.
2. MARKER_A falls in Q1_concordant_up (1-peptide exploratory detection, unique_peptides == 1).
3. CARGO_Y shows a reciprocal shift but is a single-peptide exploratory candidate (1-pep).
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        bad_path = Path(tmpdir) / "bad.md"
        good_path = Path(tmpdir) / "good.md"
        bad_path.write_text(bad_md, encoding="utf-8")
        good_path.write_text(good_md, encoding="utf-8")

        res_bad = s9.audit_markdown_against_facts(bad_path, facts)
        res_good = s9.audit_markdown_against_facts(good_path, facts)

        assert res_bad["status"] == "FAIL"
        assert res_bad["n_critical_errors"] == 4
        assert res_good["status"] == "PASS"
        assert res_good["n_critical_errors"] == 0
    print(
        f"  [PASS] Test 5: Stage 7B Adversarial Fact Gate "
        f"(Caught {res_bad['n_critical_errors']}/4 injected errors in bad.md; 0 errors in calibrated good.md)"
    )


def main() -> None:
    print("=" * 76)
    print("Running `ProteoSkill` End-to-End Validation Suite")
    print("=" * 76)
    test_1_stage0_sanitizer_and_schema()
    test_2_stage15_shared_control_null_integral()
    test_3_stage2_statistical_engines()
    test_4_stage3_archetype_solvers()
    test_5_stage7b_adversarial_fact_gate()
    print("=" * 76)
    print("ALL 5 VALIDATION SUITES PASSED (0 ERRORS)")
    print("=" * 76)


if __name__ == "__main__":
    main()

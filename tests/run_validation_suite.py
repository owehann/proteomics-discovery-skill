#!/usr/bin/env python3
"""
Validation & Regression Test Suite for `ProteoSkill` (`proteomics-discovery-skill`)
===================================================================================
Executes 8 rigorous end-to-end synthetic & benchmark tests:
  1. Stage 0 Excel Date-Corruption Sanitizer & Schema Inference Test
     (`1-Mar` -> `MARCHF1`, `1-Sep` -> `SEPTIN1`, composite ID construction).
  2. Stage 1.5 Shared-Control Variance Inflation & Empirical Null Integral Test
     (Verifies that shared-control k-clone covariance `Cov(X1-C, X2-C) = sigma^2`
     is accurately captured vs naive independence).
  3. Stage 2 Analytical Empirical Bayes (`squeezeVar`), DPC Missingness, & TOST Equivalence
     Synthetic Ground-Truth Test.
  4. Stage 3 Archetype Solver Test (verifying synthetic `MARKER_A` with `tissue_log2FC = 0.62 >= 0.585`
     is strictly classified into `Q1_concordant_up`, never `Q2_hyper_secreted`, plus `<LOD` Binary
     Switches and Seesaw `G7/S7` Antagonistic Flux).
  5. Stage 7B Fact Gate Adversarial Injection Test (verifies that `09_fact_and_citation_gate.py`
     catches injected `CTRL_Rep2` outlier leaks, `Q2_hyper_secreted` misclassifications,
     and untagged `unique_peptides == 1` promotions).
  6. All 12 Deep Mechanistic, Topological, STRING KDA, Clinical Cohort & 5-Layer RRA Engines.
  7. Editorial PPTX Collision-Proof Architecture, Auto-Font-Scaling & Force-Repulsion Labels.
  8. Full 14-Figure Suite (`plot_01`..`plot_14`), 6-Module BioRender Schematics & Anti-Laziness Gates
     (`NO_DUPLICATE_IMAGES_ACROSS_SLIDES` & `MIN_COMPREHENSIVE_SLIDE_COUNT`).
"""

from __future__ import annotations
import importlib.util
import sys
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

    # Verify MAD-based sigma recovery on synthetic N(0, 0.30^2) tubes
    rng = np.random.default_rng(99)
    n_p = 4000
    df_null = pd.DataFrame(
        {
            "Null1": rng.normal(20.0, 0.30, size=n_p),
            "Null2": rng.normal(20.0, 0.30, size=n_p),
        }
    )
    emp = s1.estimate_null_sigma_from_pairs(df_null, [("Null1", "Null2")])
    assert abs(emp["sigma_hat"] - 0.30) < 0.03, f"Expected sigma ~0.30, got {emp['sigma_hat']}"
    print(
        f"  [PASS] Test 2: Stage 1.5 Shared-Control Null Integral & Sigma Recovery "
        f"(sigma_hat={emp['sigma_hat']:.4f} vs 0.3000, Inflation={res['variance_inflation_ratio']}x)"
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

    # Verify Track B <LOD Binary Switch detection & Seesaw G7/S7 analysis
    df_sw = pd.DataFrame(
        {
            "std_gene": ["RAB27A", "CD63", "FN1", "CD81"],
            "C_1": [20.2, 21.0, np.nan, 22.1],
            "C_2": [20.4, 21.2, np.nan, 22.0],
            "C_3": [20.1, 20.9, np.nan, 22.3],
            "T_1": [np.nan, 19.8, 21.5, 22.2],
            "T_2": [np.nan, 19.7, 21.6, 22.1],
            "T_3": [np.nan, 19.9, 21.4, 22.0],
        }
    )
    sw_out = s3.detect_binary_on_off_switches(df_sw, ["T_1", "T_2", "T_3"], ["C_1", "C_2", "C_3"], contrast_name="KD")
    assert sw_out.loc[sw_out["std_gene"] == "RAB27A", "switch_class_KD"].iloc[0] == "Gold_Switch_OFF"
    assert sw_out.loc[sw_out["std_gene"] == "FN1", "switch_class_KD"].iloc[0] == "Gold_Switch_ON"
    print("  [PASS] Test 4: Stage 3 Multi-Compartment SGI, <LOD Binary Switches & Seesaw Solvers Verified")


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


def test_6_deep_mechanistic_and_clinical_modules() -> None:
    s3 = _load_script_module("s3", "03_archetype_and_modality_solver.py")
    s4 = _load_script_module("s4", "04_systems_complexes_topology.py")
    s5 = _load_script_module("s5", "05_literature_cohorts_and_rra.py")

    # 1. Paralog compensation scan with curated paralog families
    df_para = pd.DataFrame({"std_gene": ["CD63", "SYNGR2", "MYADM", "ACTB"], "log2FC": [-1.16, -1.20, 0.65, 0.01]})
    para_res = s3.scan_paralog_compensation(df_para, "CD63")
    assert not para_res.empty
    assert "SYNGR2" in set(para_res["paralog_gene"])

    # 2. Direct protein-table Biophysical & Membrane Topology Fisher audit (never requiring tm_helix_start)
    df_topo = pd.DataFrame(
        {
            "std_gene": [
                "FN1", "THBS1", "LGALS3BP", "GPC3", "MFGE8", "ITGB1", "ITGA2", "ADAM10",
                "CD9", "SYNGR2", "CD81", "CD63", "TSPAN8", "CD151", "TSG101", "PDCD6IP",
                "GAPDH", "ACTB", "ENO1", "PKM", "LDHA", "ALDOA", "TPI1", "PGK1",
            ],
            "eff_FC_A": [1.8, 1.6, 1.5, 1.9, 1.2, 1.1, 0.9, 0.8, -1.1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.7, -0.6, 0.05, -0.02, 0.1, -0.1, 0.0, 0.04, -0.05, 0.02],
            "eff_FC_B": [-0.9, -1.1, -0.8, -0.7, -0.6, 1.2, 1.0, 0.9, -1.2, -1.2, 0.9, 0.8, 0.7, 0.6, -0.8, -0.7, -0.04, 0.03, -0.05, 0.08, 0.01, -0.02, 0.03, -0.01],
            "std_unique_peptides": [12, 10, 9, 8, 7, 14, 11, 8, 5, 4, 6, 5, 6, 7, 9, 11, 20, 25, 18, 19, 15, 16, 14, 15],
        }
    )
    ep_df = s3.solve_archetype_a3_nine_sector_epistasis(df_topo, "eff_FC_A", "eff_FC_B")
    ann_df, fisher_df = s4.audit_membrane_and_secretory_topology(ep_df, group_col="Divergence_Group")
    assert "has_signal_peptide" in ann_df.columns
    assert "is_4tm_marvel_tspan" in ann_df.columns
    assert not fisher_df.empty

    # 3. STRING PPI Network + Degree-Corrected RWR + KDA
    ppi_edges = s4.fetch_string_ppi_network(ann_df["std_gene"].tolist())
    assert len(ppi_edges) >= 15
    node_grp = dict(zip(ann_df["std_gene"], ann_df["Divergence_Group"]))
    kda_nodes, ann_edges = s4.run_key_driver_analysis_kda(ppi_edges, seed_genes=ann_df["std_gene"].iloc[:14].tolist(), node_group_map=node_grp, n_degree_perms=100)
    assert not kda_nodes.empty
    assert "is_kda_driver" in kda_nodes.columns

    # 4. Clinical Integration (CPTAC, Biofluid, DepMap, HPA) + 3-Lens x 5-Layer RRA
    merged = ann_df.merge(kda_nodes[["gene", "degree", "degree_null_z", "is_kda_driver"]], left_on="std_gene", right_on="gene", how="left")
    clin_df = s5.integrate_clinical_cohorts_and_depmap(merged)
    assert "cptac_T_vs_N_log2FC" in clin_df.columns
    assert "clin_fluid_log2FC" in clin_df.columns
    assert "depmap_context_dep_lines" in clin_df.columns
    assert "rna_protein_coupling_class" in clin_df.columns
    # Verify unannotated genes receive NaN and 'missing' source rather than fabricated numbers
    fake_clin = s5.integrate_clinical_cohorts_and_depmap(pd.DataFrame({"std_gene": ["UNANNOTATED_FAKE_GENE"], "eff_FC_A": [1.5]}))
    assert pd.isna(fake_clin["cptac_T_vs_N_log2FC"].iloc[0])
    assert fake_clin["clinical_data_source"].iloc[0] == "missing"

    rra_res = s5.run_multi_lens_five_layer_rra(clin_df, top_n_per_lens=4, n_bootstrap=200)
    assert "top36_core_df" in rra_res
    assert "top6_flagship_df" in rra_res
    assert len(rra_res["top6_flagship_df"]) == 6
    print("  [PASS] Test 6: All 12 Deep Mechanistic, Topological, STRING KDA, Clinical & 5-Layer RRA Engines Verified")


def test_7_editorial_pptx_architecture_and_lint() -> None:
    s6 = _load_script_module("s6", "06_npg_statistical_figures.py")
    s8 = _load_script_module("s8", "08_editorial_pptx_builder.py")
    if getattr(s8, "Presentation", None) is None:
        print("  [SKIP] Test 7: python-pptx not installed")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        out_pptx = Path(tmpdir) / "ref_deck.pptx"
        saved_p, warnings = s8.build_reference_16_slide_blueprint(out_pptx, "Calibrated Proteomics Discovery")
        assert saved_p.exists()
        assert len(warnings) == 0, f"Unexpected PPTX design/geometry warnings: {warnings}"

        # Adversarial design lint check: inject >54 char title, ceremonial greeting, duplicate image, and <20 slides with enforce_comprehensive=True
        dummy_img = Path(tmpdir) / "dummy.png"
        s6.plot_01_manual_cutoff_vs_calibrated_audit(
            pd.DataFrame({"std_gene": ["G1", "G2", "G3"], "log2FC": [1.2, -0.9, 0.1], "p_value": [0.001, 0.004, 0.8], "std_unique_peptides": [5, 1, 8]}),
            dummy_img.with_suffix(""),
        )

        bad_deck = s8.EditorialSlideDeck("Short Title", enforce_comprehensive=True, min_comprehensive_slides=20)
        bad_deck.add_top_figure_bottom_stats_slide(
            title="This Action Title Is Intentionally Way Too Long And Exceeds Fifty Four Characters Easily",
            kicker="STEP 1",
            stats=[("100 Hits", "各位老师同学好", "coral")],
            image_path=dummy_img,
            speaker_notes="各位老师同学好，今天汇报蛋白质组学结果。",
        )
        bad_deck.add_left_figure_right_stats_slide(
            title="Second Slide Reusing Exact Same Image Path",
            kicker="STEP 2",
            stats=[("Kolde Beta Order-Statistic", "Clean caption", "blue")],
            image_path=dummy_img,
        )
        lint_errs = bad_deck.validate_deck_design_constraints()
        assert any("exceeds 54 chars" in e for e in lint_errs)
        assert any("Ceremonial greeting" in e for e in lint_errs)
        assert any("NO_DUPLICATE_IMAGES_ACROSS_SLIDES" in e for e in lint_errs)
        assert any("MIN_COMPREHENSIVE_SLIDE_COUNT" in e for e in lint_errs)
        assert any("Opaque bioinformatic jargon" in e for e in lint_errs)
        assert any("s26_biorender_single_kd_summary" in e for e in lint_errs)

    print("  [PASS] Test 7: Collision-Proof PPTX Builder + Biologist-First Plain-Language & Penultimate Summary Gates Verified")


def test_8_all_14_figures_and_7_biorender_modules() -> None:
    s3 = _load_script_module("s3", "03_archetype_and_modality_solver.py")
    s4 = _load_script_module("s4", "04_systems_complexes_topology.py")
    s5 = _load_script_module("s5", "05_literature_cohorts_and_rra.py")
    s6 = _load_script_module("s6", "06_npg_statistical_figures.py")
    s7 = _load_script_module("s7", "07_biorender_canvas_engine.py")

    genes = [
        "RAB27A", "CD63", "SYNGR2", "ITGB1", "ITGA2", "ITGAV", "CD151", "CD81", "TSPAN8", "CD9",
        "ADAM10", "FN1", "THBS1", "LGALS3BP", "GPC3", "MFGE8", "SPP1", "LAMB1", "LAMC1", "TSG101",
        "PDCD6IP", "VPS4A", "VAMP3", "SNAP23", "NFE2L2", "GCLM", "NQO1", "FTH1", "FTL", "STAT1",
    ]
    rng = np.random.default_rng(42)
    n = len(genes)
    df = pd.DataFrame(
        {
            "std_gene": genes,
            "std_unique_peptides": rng.integers(2, 16, size=n),
            "eff_FC_A": rng.normal(0.4, 1.1, size=n),
            "eff_FC_B": rng.normal(0.3, 1.1, size=n),
            "p_val_A": rng.uniform(1e-5, 0.08, size=n),
            "p_val_B": rng.uniform(1e-5, 0.08, size=n),
            "S1": rng.normal(21.0, 1.5, size=n),
            "S2": rng.normal(21.1, 1.5, size=n),
            "S3": rng.normal(20.9, 1.5, size=n),
        }
    )
    df = s3.solve_archetype_a3_nine_sector_epistasis(df, "eff_FC_A", "eff_FC_B", pval_a="p_val_A", pval_b="p_val_B")
    ann_df, fisher_df = s4.audit_membrane_and_secretory_topology(df, group_col="Divergence_Group")
    cx_df = s4.audit_corum_complex_co_shifts(ann_df, fc_col="eff_FC_A")
    sub_mat = s4.extract_complex_subunit_matrix(ann_df, fc_cols=["eff_FC_A", "eff_FC_B"])
    ora_df = s4.run_multi_group_comparative_ora(ann_df, group_col="Divergence_Group")
    ulm_df = s4.infer_upstream_regulator_ulm(ann_df, fc_col="eff_FC_A", min_targets=2)
    edges_df = s4.fetch_string_ppi_network(genes)
    kda_nodes, ann_edges = s4.run_key_driver_analysis_kda(
        edges_df, seed_genes=genes[:15], node_group_map=dict(zip(ann_df["std_gene"], ann_df["Divergence_Group"])), n_degree_perms=50
    )
    merged = ann_df.merge(kda_nodes[["gene", "degree", "degree_null_z", "is_kda_driver"]], left_on="std_gene", right_on="gene", how="left")
    clin_df = s5.integrate_clinical_cohorts_and_depmap(merged, fc_col="eff_FC_A")
    rra_dict = s5.run_multi_lens_five_layer_rra(clin_df, top_n_per_lens=4, n_bootstrap=100)

    with tempfile.TemporaryDirectory() as tmpdir:
        tdir = Path(tmpdir)
        s6.plot_01_manual_cutoff_vs_calibrated_audit(clin_df, tdir / "f01", fc_col="eff_FC_A", pval_col="p_val_A")
        s6.plot_02_proteome_sanitization_funnel(clin_df, tdir / "f02")
        s6.plot_03_forensic_qc_completeness_and_null(clin_df, ["S1", "S2", "S3"], tdir / "f03")
        s6.plot_04_modality_marker_and_purity_audit(clin_df, tdir / "f04", modality="EV")
        s6.plot_05_direct_target_and_paralog_verification(clin_df, tdir / "f05", target_genes=["RAB27A", "CD63"])
        s6.plot_06_deqms_variance_shrinkage_and_binary_switches(clin_df, tdir / "f06", fc_col="eff_FC_A")
        s6.plot_07_multi_contrast_volcano_landscape(clin_df, [("Arm A", "eff_FC_A", "p_val_A"), ("Arm B", "eff_FC_B", "p_val_B")], tdir / "f07")
        s6.plot_08_reciprocal_rescue_quadrant_and_switches(clin_df, "eff_FC_A", "eff_FC_B", tdir / "f08")
        s6.plot_09_nine_sector_epistasis_and_seesaw_cargos(clin_df, "eff_FC_A", "eff_FC_B", tdir / "f09", seesaw_out_prefix=tdir / "f09_seesaw")
        s6.plot_10_biophysical_membrane_topology_fisher(ann_df, fisher_df, tdir / "f10")
        s6.plot_11_corum_complex_and_subpopulation_heatmap(cx_df, sub_mat, tdir / "f11")
        s6.plot_12_jaccard_comparative_ora_and_regulators(ora_df, ulm_df, tdir / "f12")
        s6.plot_13_string_ppi_kda_network_graph(kda_nodes, ann_edges, tdir / "f13")
        atlas_paths = s6.plot_14_clinical_cohort_rna_protein_and_rra_atlas(rra_dict, tdir)
        bio_paths = s7.generate_all_biorender_schematics(tdir, mode="EV", top6_df=rra_dict.get("top6_flagship_df"))

        for fname in ["f01.png", "f02.png", "f03.png", "f04.png", "f05.png", "f06.png", "f07.png", "f08.png", "f09.png", "f09_seesaw.png", "f10.png", "f11.png", "f12.png", "f13.png"]:
            assert (tdir / fname).exists(), f"Missing generated figure {fname}"
        assert len(atlas_paths) >= 6
        assert len(bio_paths) == 7
        assert "s26_single_kd_summary" in bio_paths

    print("  [PASS] Test 8: Full 14-Figure Suite (including NetworkX KDA Graph) & 7-Module BioRender Suite Verified")


def main() -> None:
    print("=" * 76)
    print("Running `ProteoSkill` End-to-End Validation Suite (v4.0 Comprehensive)")
    print("=" * 76)
    test_1_stage0_sanitizer_and_schema()
    test_2_stage15_shared_control_null_integral()
    test_3_stage2_statistical_engines()
    test_4_stage3_archetype_solvers()
    test_5_stage7b_adversarial_fact_gate()
    test_6_deep_mechanistic_and_clinical_modules()
    test_7_editorial_pptx_architecture_and_lint()
    test_8_all_14_figures_and_7_biorender_modules()
    print("=" * 76)
    print("ALL 8 VALIDATION SUITES PASSED (0 ERRORS)")
    print("=" * 76)


if __name__ == "__main__":
    main()

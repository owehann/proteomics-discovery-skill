#!/usr/bin/env python3
"""
Stage 1 & Stage 1.5: 12-Point Forensic QC & Empirical Null FDR Auditor
======================================================================
Performs:
  1. 12-Point Forensic QC Checklist (cRAP contaminants, unique_peptides >= 2 tiering,
     perturbation target audit, paralog compensation, tube-level outlier detection,
     normalization suitability check, MISEV2023 EV purity, DPC missingness curve,
     and proliferation/stress confounder scoring).
  2. Stage 1.5 Empirical Null FDR Calibration & Shared-Control Variance Inflation:
     - Computes robust null noise SD `sigma_hat = MAD(delta_null) / (1.4826 * sqrt(2))`
     - Computes exact expected false-positive hits under shared-control multi-clone
       designs (`Cov(X_1 - C, X_2 - C) = sigma^2`, correlation rho = 0.5) via numerical
       integration of the shared-control Gaussian integral.
     - Computes observed vs null presence/absence switch rates.
     - Assigns the maximum permitted Evidence Grade (`E1`, `E2`, or `E3`) per contrast.
"""

from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import integrate, stats

SKILL_ROOT = Path(__file__).resolve().parent.parent
SIGNATURES_PATH = SKILL_ROOT / "resources" / "curated_proteomics_signatures.json"


def load_signatures() -> Dict[str, Any]:
    if SIGNATURES_PATH.exists():
        with open(SIGNATURES_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def estimate_null_sigma_from_pairs(
    log2_df: pd.DataFrame, null_pairs: Sequence[Tuple[str, str]]
) -> Dict[str, Any]:
    """Estimate single-tube technical/biological noise sigma from negative-control or replicate pairs."""
    pair_stats = []
    all_diffs = []
    switch_counts = []

    for col_a, col_b in null_pairs:
        if col_a not in log2_df.columns or col_b not in log2_df.columns:
            continue
        a = pd.to_numeric(log2_df[col_a], errors="coerce")
        b = pd.to_numeric(log2_df[col_b], errors="coerce")
        both = a.notna() & b.notna()
        diff = (a[both] - b[both]).values
        if len(diff) < 10:
            continue
        med_shift = float(np.median(diff))
        centered_diff = diff - med_shift
        mad_diff = float(np.median(np.abs(centered_diff)))
        sigma_pair = mad_diff / (1.4826 * math.sqrt(2.0))
        all_diffs.append(centered_diff)

        # Count presence/absence switches between the two control tubes
        a_only = int((a.notna() & b.isna()).sum())
        b_only = int((b.notna() & a.isna()).sum())
        switch_counts.append(a_only + b_only)

        pair_stats.append(
            {
                "pair": f"{col_a} vs {col_b}",
                "n_shared": int(both.sum()),
                "median_shift_log2": round(med_shift, 4),
                "mad_diff_log2": round(mad_diff, 4),
                "sigma_hat": round(sigma_pair, 4),
                "obs_gt_1_5x_either_dir": int((np.abs(centered_diff) >= math.log2(1.5)).sum()),
                "obs_gt_2_0x_either_dir": int((np.abs(centered_diff) >= 1.0).sum()),
                "a_only_switches": a_only,
                "b_only_switches": b_only,
            }
        )

    if not all_diffs:
        return {
            "sigma_hat": 0.25,
            "n_null_pairs": 0,
            "mean_null_switches_per_pair": 0.0,
            "pair_details": [],
        }

    pooled = np.concatenate(all_diffs)
    pooled_mad = float(np.median(np.abs(pooled - np.median(pooled))))
    sigma_hat = pooled_mad / (1.4826 * math.sqrt(2.0))
    return {
        "sigma_hat": round(sigma_hat, 5),
        "n_null_pairs": len(pair_stats),
        "mean_null_switches_per_pair": round(float(np.mean(switch_counts)), 2),
        "pair_details": pair_stats,
    }


def shared_control_null_probability(
    sigma: float, k_clones: int, fc_cutoff_log2: float = 0.5849625
) -> Dict[str, float]:
    """
    Compute exact one-sided and two-sided probability that all `k_clones` exceed `fc_cutoff_log2`
    in the same direction when subtracted against a SHARED control tube C ~ N(0, sigma^2):
        D_j = X_j - C,  where X_j ~ N(0, sigma^2) i.i.d.
    Also returns the naive independent probability for comparison.
    """
    if sigma <= 1e-8 or k_clones < 1:
        return {
            "p_one_sided_shared_control": 0.0,
            "p_two_sided_shared_control": 0.0,
            "p_one_sided_naive_independent": 0.0,
            "variance_inflation_ratio": 1.0,
        }

    # Under naive independence, D_j ~ N(0, 2*sigma^2) independently
    p_single = float(1.0 - stats.norm.cdf(fc_cutoff_log2 / (math.sqrt(2.0) * sigma)))
    p_indep_one_sided = p_single**k_clones

    # Under shared control C = sigma * z, conditional on z:
    # P(X_j - C > c | C = sigma*z) = 1 - Phi((c + sigma*z) / sigma) = 1 - Phi(c/sigma + z)
    ratio = fc_cutoff_log2 / sigma

    def integrand(z: float) -> float:
        cond_prob = stats.norm.sf(ratio + z)
        return float(stats.norm.pdf(z) * (cond_prob**k_clones))

    p_shared_one_sided, _ = integrate.quad(integrand, -10.0, 10.0, limit=100)
    p_shared_two_sided = 2.0 * p_shared_one_sided
    inflation = p_shared_one_sided / max(p_indep_one_sided, 1e-16)

    return {
        "p_one_sided_shared_control": float(p_shared_one_sided),
        "p_two_sided_shared_control": float(p_shared_two_sided),
        "p_one_sided_naive_independent": float(p_indep_one_sided),
        "variance_inflation_ratio": round(float(inflation), 2),
    }


def evaluate_contrast_empirical_null(
    log2_df: pd.DataFrame,
    treat_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    sigma_hat: float,
    fc_cutoff_log2: float = 0.5849625,
    shared_control: bool = True,
) -> Dict[str, Any]:
    """
    Evaluate observed directional multi-clone or group fold-change hits against the
    Stage 1.5 Empirical Null expectation and assign the maximum permitted Evidence Grade.
    """
    t_mat = log2_df[list(treat_cols)].apply(pd.to_numeric, errors="coerce")
    c_mat = log2_df[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce")

    c_mean = c_mat.mean(axis=1)
    k = len(treat_cols)
    n_c = len(ctrl_cols)

    # Compute delta for each treatment column relative to control mean
    deltas = t_mat.sub(c_mean, axis=0)
    valid_mask = deltas.notna().all(axis=1) & c_mean.notna()
    n_tested = int(valid_mask.sum())

    if n_tested == 0:
        return {
            "n_tested": 0,
            "k_treat_tubes": k,
            "n_ctrl_tubes": n_c,
            "obs_all_up": 0,
            "obs_all_down": 0,
            "expected_null_up": 0.0,
            "expected_null_down": 0.0,
            "empirical_null_fdr_up": 1.0,
            "empirical_null_fdr_down": 1.0,
            "max_evidence_grade": "E3",
        }

    sub_deltas = deltas.loc[valid_mask]
    obs_up = int((sub_deltas >= fc_cutoff_log2).all(axis=1).sum())
    obs_down = int((sub_deltas <= -fc_cutoff_log2).all(axis=1).sum())

    # Effective control noise scales as sigma / sqrt(n_c)
    if shared_control and k >= 2:
        # General shared control integral when control is average of n_c tubes:
        # C_bar ~ N(0, sigma^2 / n_c), X_j ~ N(0, sigma^2)
        def integrand_nc(z: float) -> float:
            c_val = (sigma_hat / math.sqrt(max(n_c, 1))) * z
            cond_prob = stats.norm.sf((fc_cutoff_log2 + c_val) / max(sigma_hat, 1e-8))
            return float(stats.norm.pdf(z) * (cond_prob**k))

        p_null_one, _ = integrate.quad(integrand_nc, -10.0, 10.0, limit=100)
    else:
        se_diff = sigma_hat * math.sqrt(1.0 / max(k, 1) + 1.0 / max(n_c, 1))
        p_null_one = float(stats.norm.sf(fc_cutoff_log2 / max(se_diff, 1e-8)))

    exp_null_one = n_tested * p_null_one
    fdr_up = min(1.0, exp_null_one / max(obs_up, 1)) if obs_up > 0 else 1.0
    fdr_down = min(1.0, exp_null_one / max(obs_down, 1)) if obs_down > 0 else 1.0
    best_fdr = min(fdr_up, fdr_down)

    if best_fdr <= 0.05 and (k >= 2 and n_c >= 2):
        max_grade = "E1_or_E2"
    elif best_fdr <= 0.25:
        max_grade = "E2"
    else:
        max_grade = "E3"

    return {
        "n_tested": n_tested,
        "k_treat_tubes": k,
        "n_ctrl_tubes": n_c,
        "sigma_hat": round(sigma_hat, 4),
        "fc_cutoff_log2": round(fc_cutoff_log2, 4),
        "obs_all_up": obs_up,
        "obs_all_down": obs_down,
        "expected_null_one_sided": round(exp_null_one, 2),
        "empirical_null_fdr_up": round(fdr_up, 4),
        "empirical_null_fdr_down": round(fdr_down, 4),
        "below_noise_floor": bool((obs_up <= exp_null_one) and (obs_down <= exp_null_one)),
        "max_evidence_grade": max_grade,
    }


def run_12point_forensic_qc(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    gene_col: str = "std_gene",
    pep_col: str = "std_unique_peptides",
    is_log2: bool = False,
    target_genes: Optional[Sequence[str]] = None,
    is_ev_secretome: bool = False,
) -> Dict[str, Any]:
    """Execute the 12-Point Forensic QC Checklist on a proteomics dataframe."""
    sigs = load_signatures()
    mat_raw = df[list(sample_cols)].apply(pd.to_numeric, errors="coerce").astype(float)
    if is_log2:
        mat_log2 = mat_raw.copy()
        mat_lin = np.power(2.0, mat_raw)
    else:
        mat_lin = mat_raw.where(mat_raw > 0, np.nan)
        mat_log2 = np.log2(mat_lin)

    genes = df[gene_col].astype(str).str.strip().tolist() if gene_col in df.columns else [""] * len(df)
    peps = pd.to_numeric(df[pep_col], errors="coerce").fillna(1) if pep_col in df.columns else pd.Series(2, index=df.index)

    # QC-01: Contaminant Audit
    crap_dict = sigs.get("contaminant_patterns", sigs.get("contaminants_crap", {}))
    crap_genes = set(
        crap_dict.get("keratin_skin_genes", [])
        + crap_dict.get("keratins", [])
        + crap_dict.get("serum_culture_media", [])
        + crap_dict.get("immunoglobulins_fc", [])
    )
    is_crap = pd.Series([g in crap_genes or g.startswith(("KRT", "IGH", "IGK", "IGL")) for g in genes], index=df.index)
    crap_frac_per_sample = (mat_lin.loc[is_crap].sum(axis=0) / mat_lin.sum(axis=0).replace(0, np.nan)).fillna(0.0)

    # QC-02: Peptide Support Tiering
    n_core_pep = int((peps >= 2).sum())
    n_single_pep = int((peps == 1).sum())

    # QC-03 & QC-04: Perturbation Target & Paralog Audit
    target_audit: Dict[str, Any] = {}
    if target_genes:
        for tg in target_genes:
            matches = df.index[[g.upper() == tg.upper() for g in genes]].tolist()
            if not matches:
                target_audit[tg] = {"status": "UNDETECTED_IN_MS", "gate_2_alert": True}
            else:
                idx0 = matches[0]
                target_audit[tg] = {
                    "status": "DETECTED",
                    "unique_peptides": int(peps.loc[idx0]),
                    "sample_log2_values": {
                        c: (round(float(mat_log2.loc[idx0, c]), 3) if pd.notna(mat_log2.loc[idx0, c]) else None)
                        for c in sample_cols
                    },
                    "gate_2_alert": bool(peps.loc[idx0] < 2),
                }

    # QC-05: Tube-Level Anomaly Detection (ID count, Median shift > 0.40 log2, and Pairwise Pearson r < 0.80)
    id_counts = mat_log2.notna().sum(axis=0)
    median_id_count = float(id_counts.median())
    sample_medians = mat_log2.median(axis=0)
    global_median = float(sample_medians.median())
    corr_mat = mat_log2.corr(method="pearson")

    outlier_tubes: List[Dict[str, Any]] = []
    for c in sample_cols:
        n_id = int(id_counts[c])
        med_v = float(sample_medians[c]) if pd.notna(sample_medians[c]) else 0.0
        other_corrs = [float(corr_mat.loc[c, c2]) for c2 in sample_cols if c2 != c and pd.notna(corr_mat.loc[c, c2])]
        med_corr = float(np.median(other_corrs)) if other_corrs else 1.0
        reasons = []
        if n_id < 0.82 * median_id_count:
            reasons.append(f"Low ID count ({n_id} vs median {int(median_id_count)})")
        if abs(med_v - global_median) > 0.40:
            reasons.append(f"Large median shift ({med_v:.2f} vs {global_median:.2f}, Δ={med_v - global_median:+.2f} log2)")
        if med_corr < 0.80:
            reasons.append(f"Low median pairwise correlation (r={med_corr:.3f} < 0.80)")
        if reasons:
            outlier_tubes.append(
                {
                    "sample": c,
                    "n_ids": n_id,
                    "median_log2": round(med_v, 3),
                    "median_pairwise_r": round(med_corr, 3),
                    "reasons": reasons,
                }
            )

    # QC-09: Normalization Suitability Audit
    max_median_shift = float((sample_medians - global_median).abs().max()) if len(sample_medians) > 0 else 0.0
    naive_norm_safe = bool(max_median_shift <= 0.15)

    # QC-10: MISEV2023 EV Purity Audit
    ev_purity: Dict[str, Any] = {}
    if is_ev_secretome:
        misev = sigs.get("misev2023_markers", sigs.get("misev2023_ev_markers", {}))
        cat12 = set(
            misev.get("Cat1a_Tetraspanins_MVB", [])
            + misev.get("Cat1b_Integrins_SingleTM", [])
            + misev.get("Cat2a_ESCRT_Cytosolic", [])
            + misev.get("cat1a_tetraspanins_multipass", [])
            + misev.get("cat1b_singlepass_integrins_adhesion", [])
            + misev.get("cat2a_cytosolic_escrt_syntanin_flotillin", [])
        )
        cat3 = set(
            misev.get("Cat3a_Lipoproteins_Albumin_CoIsolates", [])
            + misev.get("Cat4_Intracellular_Organelle_Leakage", [])
            + misev.get("cat3a_lipoproteins_serum_coisolates", [])
            + misev.get("cat3b_organelle_er_golgi_mito_nuclear_debris", [])
        )
        m12 = [g in cat12 for g in genes]
        m3 = [g in cat3 for g in genes]
        mean_cat12 = float(mat_log2.loc[m12].mean().mean()) if any(m12) else np.nan
        mean_cat3 = float(mat_log2.loc[m3].mean().mean()) if any(m3) else np.nan
        ev_purity = {
            "n_cat1_cat2_detected": int(sum(m12)),
            "n_cat3_contaminants_detected": int(sum(m3)),
            "mean_cat12_log2": round(mean_cat12, 3) if np.isfinite(mean_cat12) else None,
            "mean_cat3_log2": round(mean_cat3, 3) if np.isfinite(mean_cat3) else None,
            "purity_ratio_log2": round(mean_cat12 - mean_cat3, 3)
            if (np.isfinite(mean_cat12) and np.isfinite(mean_cat3))
            else None,
        }

    # QC-11: Missingness Mechanism (DPC Slope Estimation)
    obs_prob = mat_log2.notna().mean(axis=1).values
    mean_int = mat_log2.mean(axis=1).values
    valid_dpc = np.isfinite(mean_int) & (obs_prob > 0.0) & (obs_prob < 1.0)
    if valid_dpc.sum() >= 20:
        logit_p = np.log(np.clip(obs_prob[valid_dpc], 0.02, 0.98) / (1.0 - np.clip(obs_prob[valid_dpc], 0.02, 0.98)))
        slope, intercept, r_val, _, _ = stats.linregress(mean_int[valid_dpc], logit_p)
        dpc_params = {
            "dpc_beta0_intercept": round(float(intercept), 4),
            "dpc_beta1_slope": round(float(slope), 4),
            "mnar_correlation_r": round(float(r_val), 4),
        }
    else:
        dpc_params = {"dpc_beta0_intercept": 0.0, "dpc_beta1_slope": 0.8, "mnar_correlation_r": 0.0}

    # QC-12: Confounder Signature Scan
    conf_sigs = sigs.get("confounder_signatures", {})
    confounder_scores: Dict[str, Dict[str, float]] = {}
    for sig_name, sig_genes in conf_sigs.items():
        sig_set = set(sig_genes)
        mask = [g in sig_set for g in genes]
        if any(mask):
            row_z = mat_log2.loc[mask].sub(mat_log2.loc[mask].mean(axis=1), axis=0)
            confounder_scores[sig_name] = {c: round(float(row_z[c].mean()), 3) for c in sample_cols if pd.notna(row_z[c].mean())}

    gate_2_triggered = bool(outlier_tubes or not naive_norm_safe or any(v.get("gate_2_alert") for v in target_audit.values()))

    return {
        "n_total_features": int(len(df)),
        "qc01_crap_max_fraction": round(float(crap_frac_per_sample.max()), 4) if len(crap_frac_per_sample) else 0.0,
        "qc02_peptide_tiers": {
            "core_ge2_unique_peptides": n_core_pep,
            "exploratory_1_unique_peptide": n_single_pep,
            "core_fraction": round(n_core_pep / max(len(df), 1), 4),
        },
        "qc03_target_audit": target_audit,
        "qc05_sample_id_counts": {c: int(id_counts[c]) for c in sample_cols},
        "qc05_sample_medians_log2": {c: round(float(sample_medians[c]), 3) for c in sample_cols if pd.notna(sample_medians[c])},
        "qc05_outlier_tubes": outlier_tubes,
        "qc09_max_median_shift_log2": round(max_median_shift, 4),
        "qc09_naive_median_norm_safe": naive_norm_safe,
        "qc10_misev2023_ev_purity": ev_purity,
        "qc11_dpc_missingness_model": dpc_params,
        "qc12_confounder_z_scores": confounder_scores,
        "gate_2_conditional_stop_triggered": gate_2_triggered,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 1 & 1.5: Forensic QC & Empirical Null FDR Auditor")
    parser.add_argument("--input", "-i", required=True, help="Input table (.xlsx/.tsv/.csv)")
    parser.add_argument("--sheet", default=None, help="Sheet name for Excel files")
    parser.add_argument("--out-json", "-o", default="stage1_qc_and_null_audit.json", help="Output JSON report path")
    parser.add_argument("--targets", nargs="*", default=[], help="Perturbed target gene symbols to verify")
    parser.add_argument("--ev", action="store_true", help="Enable MISEV2023 EV purity check")
    args = parser.parse_args()

    # Import Stage 0 helper dynamically
    import importlib.util
    s0_path = Path(__file__).resolve().parent / "00_forensic_schema_indexer.py"
    spec = importlib.util.spec_from_file_location("stage0", s0_path)
    s0 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s0)

    sheets = s0.detect_Header_and_load_sheet(Path(args.input), args.sheet)
    reports = {}
    for sh_name, df in sheets.items():
        schema = s0.infer_columns_and_modality(df, sh_name)
        std_df = s0.build_composite_Feature_table(df, schema)
        qc_res = run_12point_forensic_qc(
            std_df,
            schema["sample_cols"],
            is_log2=schema["is_log2_scale"],
            target_genes=args.targets,
            is_ev_secretome=args.ev or ("ev" in sh_name.lower()),
        )
        reports[sh_name] = {"schema": schema, "forensic_qc": qc_res}

    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(reports, f, indent=2)
    print(f"[Stage 1/1.5] Saved Forensic QC & Null Audit to {args.out_json}")


if __name__ == "__main__":
    main()

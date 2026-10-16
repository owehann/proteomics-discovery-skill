#!/usr/bin/env python3
"""
Stage 3: Universal Experimental Archetype (A1-A10), Binary `<LOD` Switch, Seesaw & Paralog Solver
=================================================================================================
Implements exact algebraic and statistical scoring for Engines [01]-[03]:
  - Engine [01]:
      * Track B Qualitative `<LOD` All-or-Nothing Binary Switch Detection (`detect_binary_on_off_switches`)
        with sample-specific QRILC 5th-percentile `- 0.5` left-censored bounds and conservative
        effective fold-change (`eff_FC_*`) calculation.
      * Archetype A1: Multi-Clone Consistency & Hairpin Seed Divergence Audit (`solve_archetype_a1_multi_clone`).
      * Archetype A2: Reciprocal Perturbation `KD x OE` (`solve_archetype_a2_reciprocal_kd_oe`) with
        `R-Score`, 4-Quadrant Driver vs Stress Artifact Partition, and Gold-Standard `ON-OFF-ON` Switch rescue.
  - Engine [02]:
      * Archetype A3: 9-Sector Cartesian Epistasis & Divergence Grid (`solve_archetype_a3_nine_sector_epistasis`)
        with algebraically control-free Head-to-Head Divergence Index (`DI = mean_B - mean_A`),
        moderated t-test `p_DI_divergence`, and **Seesaw Antagonistic Flux Competition (`Sector G7/S7`)**
        subtyping (`analyze_seesaw_antagonistic_flux`).
      * Archetype A4: Biologics & Antibody Engineering 3-Axis Decomposition (`solve_archetype_a4_biologics_humanization`).
      * Archetype A5: Clinical Paired Multi-Compartment Secretion Gradient Index (`SGI`, `solve_archetype_a5_clinical_sgi`).
      * Modality M2 / A7: PTM Stoichiometry Uncoupling (`solve_modality_m2_ptm_stoichiometry`).
  - Engine [03]:
      * Paralog & Isoform Compensation Scanner (`scan_paralog_compensation`) combining curated
        paralog families (`resources/curated_proteomics_signatures.json`) and HGNC gene-root prefixes,
        outputting both fold-change classifications and a cross-sample intensity matrix for plotting.
"""

from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats


SKILL_ROOT = Path(__file__).resolve().parent.parent
SIGNATURES_PATH = SKILL_ROOT / "resources" / "curated_proteomics_signatures.json"


def _load_paralog_families() -> Dict[str, List[str]]:
    fam_index: Dict[str, set] = {
        "CD63": {"CD63", "CD81", "CD9", "CD151", "CD82", "TSPAN8", "TSPAN4", "MYADM", "MAL2"},
        "RAB27A": {"RAB27A", "RAB27B", "RAB3A", "RAB3D", "RAB8A", "RAB10", "RAB11A", "RAB11B"},
        "NFE2L2": {"NFE2L2", "NFE2L1", "GCLM", "GCLC", "GSTP1", "NQO1", "TXNRD1", "HMOX1", "FTH1", "FTL"},
        "SMARCA4": {"SMARCA4", "SMARCA2", "SMARCA5", "CHD4", "DNMT1", "DNMT3A", "DNMT3B", "UHRF1"},
        "STAT1": {"STAT1", "STAT2", "STAT3", "IRF9", "ISG15", "PARP9", "DTX3L"},
    }
    if SIGNATURES_PATH.exists():
        try:
            with open(SIGNATURES_PATH, "r", encoding="utf-8") as f:
                data = json.load(f)
            raw_fams = data.get("paralog_families", {})
            if isinstance(raw_fams, dict):
                for k, members in raw_fams.items():
                    if not isinstance(members, list):
                        continue
                    m_up = {str(g).strip().upper() for g in members if str(g).strip()}
                    k_up = str(k).strip().upper()
                    fam_index.setdefault(k_up, set()).update(m_up)
                    for g in m_up:
                        fam_index.setdefault(g, set()).update(m_up | ({k_up} if not k_up.endswith("FAMILY") else set()))
        except Exception:
            pass
    # Ensure symmetric closure for built-in defaults
    for k, m_set in list(fam_index.items()):
        for g in list(m_set):
            fam_index.setdefault(g, set()).update(m_set)
    return {k: sorted(list(v)) for k, v in fam_index.items()}


def compute_sample_lod_bounds(
    df: pd.DataFrame,
    sample_cols: Sequence[str],
    quantile: float = 0.05,
    offset_log2: float = -0.50,
) -> Dict[str, float]:
    """
    Compute sample-specific left-censored Limit of Detection (LOD) bounds on log2 scale:
        LOD_s = Q_{0.05}(observed log2 intensities in sample s) - 0.50
    """
    lod_map: Dict[str, float] = {}
    for col in sample_cols:
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(vals) > 0:
                lod_map[col] = round(float(np.quantile(vals, quantile) + offset_log2), 4)
            else:
                lod_map[col] = 10.5
    return lod_map


def detect_binary_on_off_switches(
    df: pd.DataFrame,
    treat_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    oe_treat_cols: Optional[Sequence[str]] = None,
    oe_ctrl_cols: Optional[Sequence[str]] = None,
    presence_min_log2: float = 11.8,
    tier1_strict_min_log2: float = 13.0,
    lod_map: Optional[Dict[str, float]] = None,
    prefix: str = "Perturb",
    contrast_name: Optional[str] = None,
) -> pd.DataFrame:
    """
    Engine [01] Track B: Qualitative MNAR All-or-Nothing (`<LOD`) Binary Switch Detector
    and QRILC Lower-Bound Effective Fold-Change (`eff_FC`) Calculator.
    """
    if contrast_name is not None:
        prefix = contrast_name
    out = df.copy()
    all_cols = list(dict.fromkeys(list(treat_cols) + list(ctrl_cols) + list(oe_treat_cols or []) + list(oe_ctrl_cols or [])))
    if lod_map is None:
        lod_map = compute_sample_lod_bounds(out, all_cols)

    t_mat = out[list(treat_cols)].apply(pd.to_numeric, errors="coerce")
    c_mat = out[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce")

    t_all_na = t_mat.isna().all(axis=1)
    c_all_na = c_mat.isna().all(axis=1)
    t_all_obs = t_mat.notna().all(axis=1)
    c_any_obs = c_mat.notna().any(axis=1)

    t_mean_obs = t_mat.mean(axis=1, skipna=True)
    c_mean_obs = c_mat.mean(axis=1, skipna=True)

    sw_off = c_any_obs & (c_mean_obs >= presence_min_log2) & t_all_na
    sw_on = c_all_na & t_all_obs & (t_mean_obs >= presence_min_log2)

    out[f"switch_{prefix}_OFF"] = sw_off
    out[f"switch_{prefix}_ON"] = sw_on

    # Impute with sample-specific LOD bound only where at least one side is observed
    t_imp = t_mat.copy()
    for col in treat_cols:
        t_imp[col] = t_imp[col].fillna(lod_map.get(col, 10.5))
    c_imp = c_mat.copy()
    for col in ctrl_cols:
        c_imp[col] = c_imp[col].fillna(lod_map.get(col, 10.5))

    c_imp_mean = c_imp.mean(axis=1)
    eff_clone_fcs = []
    for col in treat_cols:
        valid_pair = out[col].notna() | c_any_obs
        efc = np.where(valid_pair, t_imp[col] - c_imp_mean, np.nan)
        out[f"eff_FC_{col}"] = np.round(efc, 4)
        eff_clone_fcs.append(pd.Series(efc, index=out.index))

    cont_fc = t_mat.mean(axis=1, skipna=False) - c_mat.mean(axis=1, skipna=False)
    out[f"mean_FC_{prefix}"] = np.round(cont_fc, 4)
    eff_mean = pd.concat(eff_clone_fcs, axis=1).mean(axis=1, skipna=False)
    out[f"eff_mean_FC_{prefix}"] = np.round(np.where(cont_fc.notna(), cont_fc, eff_mean), 4)
    out[f"eff_FC_{prefix}"] = out[f"eff_mean_FC_{prefix}"]

    # Optional OE switch check for Gold-Standard ON-OFF-ON
    sw_oe_on = pd.Series(False, index=out.index)
    sw_oe_off = pd.Series(False, index=out.index)
    if oe_treat_cols and oe_ctrl_cols:
        ot_mat = out[list(oe_treat_cols)].apply(pd.to_numeric, errors="coerce")
        oc_mat = out[list(oe_ctrl_cols)].apply(pd.to_numeric, errors="coerce")
        ot_mean = ot_mat.mean(axis=1, skipna=True)
        oc_mean = oc_mat.mean(axis=1, skipna=True)
        sw_oe_on = oc_mat.isna().all(axis=1) & ot_mat.notna().all(axis=1) & (ot_mean >= presence_min_log2)
        sw_oe_off = oc_mat.notna().any(axis=1) & (oc_mean >= presence_min_log2) & ot_mat.isna().all(axis=1)
        out["switch_OE_ON"] = sw_oe_on
        out["switch_OE_OFF"] = sw_oe_off

        ot_imp = ot_mat.copy()
        for col in oe_treat_cols:
            ot_imp[col] = ot_imp[col].fillna(lod_map.get(col, 10.5))
        oc_imp = oc_mat.copy()
        for col in oe_ctrl_cols:
            oc_imp[col] = oc_imp[col].fillna(lod_map.get(col, 10.5))
        valid_oe = ot_mat.notna().any(axis=1) | oc_mat.notna().any(axis=1)
        out["eff_FC_OE"] = np.round(np.where(valid_oe, ot_imp.mean(axis=1) - oc_imp.mean(axis=1), np.nan), 4)

    gold_on_off_on = (sw_oe_on & (sw_off | (out[f"eff_mean_FC_{prefix}"] <= -0.5849625))) | (
        sw_off & (out.get("eff_FC_OE", pd.Series(0.0, index=out.index)) >= 0.5849625)
    )
    out[f"gold_standard_ON_OFF_ON_{prefix}"] = gold_on_off_on

    sw_class = np.full(len(out), "Continuous_Quantified", dtype=object)
    sw_class[sw_off & (c_mean_obs >= tier1_strict_min_log2)] = "Gold_Switch_OFF"
    sw_class[sw_off & (c_mean_obs < tier1_strict_min_log2)] = "Moderate_Switch_OFF"
    sw_class[sw_on & (t_mean_obs >= tier1_strict_min_log2)] = "Gold_Switch_ON"
    sw_class[sw_on & (t_mean_obs < tier1_strict_min_log2)] = "Moderate_Switch_ON"
    sw_class[gold_on_off_on] = "Gold_Standard_Bidirectional_Switch (OE_ON & KD_OFF)"
    out[f"trackB_switch_class_{prefix}"] = sw_class
    out[f"switch_class_{prefix}"] = sw_class
    return out


def solve_archetype_a1_multi_clone(
    df: pd.DataFrame,
    clone_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    fc_cutoff_log2: float = 0.5849625,
    divergence_warn_log2: float = 1.0,
    include_binary_switches: bool = True,
) -> pd.DataFrame:
    """
    Archetype A1: Multi-clone knockdown/overexpression consistency & seed off-target audit.
    Integrates Track B `<LOD` Binary Switch detection so all-or-nothing clone silences are rescued.
    """
    out = df.copy()
    if include_binary_switches:
        out = detect_binary_on_off_switches(out, clone_cols, ctrl_cols, prefix="A1")

    c_mean = out[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    clone_deltas = pd.DataFrame(index=out.index)
    for col in clone_cols:
        raw_d = pd.to_numeric(out[col], errors="coerce") - c_mean
        if f"eff_FC_{col}" in out.columns:
            raw_d = raw_d.fillna(out[f"eff_FC_{col}"])
        clone_deltas[f"delta_{col}"] = raw_d

    out["mean_clone_log2FC"] = clone_deltas.mean(axis=1).round(4)
    out["min_clone_log2FC"] = clone_deltas.min(axis=1).round(4)
    out["max_clone_log2FC"] = clone_deltas.max(axis=1).round(4)
    out["max_inter_clone_divergence_log2"] = (clone_deltas.max(axis=1) - clone_deltas.min(axis=1)).round(4)

    all_up = (clone_deltas >= fc_cutoff_log2).all(axis=1)
    all_down = (clone_deltas <= -fc_cutoff_log2).all(axis=1)
    opp_sign = (clone_deltas.max(axis=1) > 0.3) & (clone_deltas.min(axis=1) < -0.3)

    status = np.full(len(out), "Unchanged_or_Subthreshold", dtype=object)
    status[all_up] = "Concordant_MultiClone_Up"
    status[all_down] = "Concordant_MultiClone_Down"
    status[opp_sign | (out["max_inter_clone_divergence_log2"] >= divergence_warn_log2)] = (
        "Hairpin_Divergent_Seed_Artifact_Warning"
    )
    if "switch_A1_OFF" in out.columns:
        status[out["switch_A1_OFF"].astype(bool)] = "Concordant_MultiClone_Down [Binary_Switch_OFF]"
    if "switch_A1_ON" in out.columns:
        status[out["switch_A1_ON"].astype(bool)] = "Concordant_MultiClone_Up [Binary_Switch_ON]"

    out["a1_clone_consistency_class"] = status
    return pd.concat([out, clone_deltas.round(4)], axis=1)


def solve_archetype_a2_reciprocal_kd_oe(
    df: pd.DataFrame,
    delta_oe_col: str,
    delta_kd_col: str,
    fc_cutoff_log2: float = 0.5849625,
) -> pd.DataFrame:
    """
    Archetype A2: Reciprocal Perturbation (`KD x OE` or `KO + Rescue`).
    Computes both geometric Reciprocity Score (`R-Score`) and additive span (`rec_score = OE - KD`),
    and partitions proteins into True Reciprocal Drivers (Q2/Q4) vs Concordant Stress Artifacts (Q1/Q3).
    """
    out = df.copy()
    d_oe = pd.to_numeric(out[delta_oe_col], errors="coerce")
    d_kd = pd.to_numeric(out[delta_kd_col], errors="coerce")

    opp_dir = (np.sign(d_oe) * np.sign(d_kd)) < 0
    r_mag = np.sqrt(np.abs(d_oe * d_kd))
    r_score = np.where(opp_dir, np.sign(d_oe) * r_mag, 0.0)
    out["reciprocal_r_score"] = np.round(r_score, 4)
    out["rec_score"] = np.round(d_oe.fillna(0.0) - d_kd.fillna(0.0), 4)

    quad = np.full(len(out), "Background", dtype=object)
    q_pos_driver = (d_oe >= fc_cutoff_log2) & (d_kd <= -fc_cutoff_log2)
    q_neg_repress = (d_oe <= -fc_cutoff_log2) & (d_kd >= fc_cutoff_log2)
    q_stress_up = (d_oe >= fc_cutoff_log2) & (d_kd >= fc_cutoff_log2)
    q_stress_down = (d_oe <= -fc_cutoff_log2) & (d_kd <= -fc_cutoff_log2)

    quad[q_pos_driver] = "Reciprocal_Positive_Target_Cargo (OE_Up_KD_Down)"
    quad[q_neg_repress] = "Reciprocal_Repressed_Target (OE_Down_KD_Up)"
    quad[q_stress_up] = "Concordant_Up_Stress_Artifact"
    quad[q_stress_down] = "Concordant_Down_Stress_Artifact"
    out["a2_reciprocal_quadrant"] = quad
    return out


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return q
    pv = np.clip(p[valid], 0.0, 1.0)
    n = len(pv)
    order = np.argsort(pv)
    adj = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out_v = np.empty_like(pv)
    out_v[order] = np.clip(adj, 0.0, 1.0)
    q[valid] = out_v
    return q


def solve_archetype_a3_nine_sector_epistasis(
    df: pd.DataFrame,
    delta_a_col: str,
    delta_b_col: str,
    fc_cutoff_log2: float = 0.5849625,
    equiv_margin_log2: float = 0.35,
    di_min_log2: float = 0.45,
    label_a: str = "PerturbA",
    label_b: str = "PerturbB",
    s2_post_col: Optional[str] = None,
    df_post_col: Optional[str] = None,
    pval_a: Optional[str] = None,
    pval_b: Optional[str] = None,
    se_a_col: Optional[str] = None,
    se_b_col: Optional[str] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Archetype A3: 9-Sector Cartesian Divergence & Epistasis Grid for two perturbations A vs B.
    Computes:
      1. Algebraically Control-Free Head-to-Head Divergence Index:
           `eff_DI_B_minus_A = delta_B - delta_A`
         (When A and B share a common control tube, the control cancels 100% algebraically!)
      2. Moderated 2-vs-2 Head-to-Head t-test `p_DI_divergence` & `fdr_DI_divergence`.
      3. Both canonical `a3_nine_sector_class` (`S1`..`S9`) AND `Divergence_Group` (`G1_Shared_Down`,
         `G2_{label_b}_Specific_Down`, `G3_{label_a}_Specific_Down`, `G7_Opposite_Seesaw`,
         `G4_Shared_Up`, `G5_{label_b}_Specific_Up`, `G6_{label_a}_Specific_Up`,
         `G1b_Partial_Shared_Down`, `G0_Unchanged_or_Minor`).
    """
    out = df.copy()
    da = pd.to_numeric(out[delta_a_col], errors="coerce")
    db = pd.to_numeric(out[delta_b_col], errors="coerce")

    eff_di = db - da
    out["eff_DI_B_minus_A"] = eff_di.round(4)
    out[f"eff_DI_{label_b}_minus_{label_a}"] = eff_di.round(4)

    # Compute Head-to-Head Divergence p-value & BH FDR
    if se_a_col and se_b_col and se_a_col in out.columns and se_b_col in out.columns:
        sea = pd.to_numeric(out[se_a_col], errors="coerce").fillna(0.15)
        seb = pd.to_numeric(out[se_b_col], errors="coerce").fillna(0.15)
        s2_p = sea**2 + seb**2
        df_p = pd.to_numeric(out[df_post_col], errors="coerce").fillna(6.0) if (df_post_col and df_post_col in out.columns) else 6.0
    elif s2_post_col and s2_post_col in out.columns:
        s2_p = pd.to_numeric(out[s2_post_col], errors="coerce").fillna(0.04)
        df_p = pd.to_numeric(out[df_post_col], errors="coerce").fillna(6.0) if (df_post_col and df_post_col in out.columns) else 6.0
    else:
        diff_sd = float(np.nanmedian(np.abs(eff_di - np.nanmedian(eff_di))) * 1.4826)
        s2_p = max(diff_sd ** 2, 0.025)
        df_p = 6.0

    se_di = np.sqrt(np.maximum(s2_p, 1e-6))
    t_di = eff_di / se_di
    p_di = 2.0 * stats.t.sf(np.abs(t_di), df=df_p)
    out["t_DI_divergence"] = np.round(t_di, 4)
    out["p_DI_divergence"] = p_di
    out["fdr_DI_divergence"] = _bh_adjust(p_di)

    # Canonical S1-S9 classification
    sector = np.full(len(out), "S0_Subthreshold_Intermediate", dtype=object)
    a_up, a_dn, a_eq = da >= fc_cutoff_log2, da <= -fc_cutoff_log2, da.abs() <= equiv_margin_log2
    b_up, b_dn, b_eq = db >= fc_cutoff_log2, db <= -fc_cutoff_log2, db.abs() <= equiv_margin_log2

    sector[a_up & b_up] = f"S1_Shared_Up ({label_a}_Up & {label_b}_Up)"
    sector[a_dn & b_dn] = f"S2_Shared_Down ({label_a}_Down & {label_b}_Down)"
    sector[a_up & b_eq] = f"S3_{label_a}_Specific_Up"
    sector[a_dn & b_eq] = f"S4_{label_a}_Specific_Down"
    sector[b_up & a_eq] = f"S5_{label_b}_Specific_Up"
    sector[b_dn & a_eq] = f"S6_{label_b}_Specific_Down"
    sector[a_up & b_dn] = f"S7_Antagonistic ({label_a}_Up & {label_b}_Down)"
    sector[a_dn & b_up] = f"S8_Antagonistic ({label_a}_Down & {label_b}_Up)"
    sector[a_eq & b_eq] = "S9_Dual_Equivalent_Core"
    out["a3_nine_sector_class"] = sector

    # Benchmark G0-G7 Divergence_Group partition with DI guard
    not_a_dn = (~a_dn) & (da.fillna(0) > -equiv_margin_log2) & (eff_di.fillna(-1.0) <= -di_min_log2)
    not_b_dn = (~b_dn) & (db.fillna(0) > -equiv_margin_log2) & (eff_di.fillna(1.0) >= di_min_log2)
    not_a_up = (~a_up) & (da.fillna(0) < equiv_margin_log2) & (eff_di.fillna(1.0) >= di_min_log2)
    not_b_up = (~b_up) & (db.fillna(0) < equiv_margin_log2) & (eff_di.fillna(-1.0) <= -di_min_log2)

    cond_G1 = a_dn & b_dn
    cond_G7 = (a_dn & b_up) | (a_up & b_dn)
    cond_G2 = b_dn & (~cond_G1) & (~cond_G7) & not_a_dn
    cond_G3 = a_dn & (~cond_G1) & (~cond_G7) & not_b_dn
    cond_G4 = a_up & b_up
    cond_G5 = b_up & (~cond_G4) & (~cond_G7) & not_a_up
    cond_G6 = a_up & (~cond_G4) & (~cond_G7) & not_b_up
    cond_G1b = (
        ((a_dn & (~b_dn) & (db <= -equiv_margin_log2)) | (b_dn & (~a_dn) & (da <= -equiv_margin_log2)))
        & (~cond_G1)
        & (~cond_G2)
        & (~cond_G3)
    )

    out["Divergence_Group"] = np.select(
        [cond_G1, cond_G2, cond_G3, cond_G7, cond_G4, cond_G5, cond_G6, cond_G1b],
        [
            "G1_Shared_Down",
            f"G2_{label_b}_Specific_Down",
            f"G3_{label_a}_Specific_Down",
            "G7_Opposite_Seesaw",
            "G4_Shared_Up",
            f"G5_{label_b}_Specific_Up",
            f"G6_{label_a}_Specific_Up",
            "G1b_Partial_Shared_Down",
        ],
        default="G0_Unchanged_or_Minor",
    )
    return out


def analyze_seesaw_antagonistic_flux(
    df: pd.DataFrame,
    delta_a_col: str,
    delta_b_col: str,
    clone_a_cols: Optional[Sequence[str]] = None,
    clone_b_cols: Optional[Sequence[str]] = None,
    gene_col: str = "std_gene",
    fc_cutoff_log2: float = 0.5849625,
    label_a: str = "PerturbA",
    label_b: str = "PerturbB",
) -> pd.DataFrame:
    """
    Engine [02] Seesaw Antagonistic Flux Competition Analyzer (`Sector G7/S7`).
    Extracts all proteins where Perturbation A and Perturbation B push cargo/abundance in
    strictly opposite directions (`A_Down & B_Up` vs `A_Up & B_Down`), quantifies the
    net Flux Shunting Span (`DI = delta_B - delta_A`), and returns a sorted table ready
    for `plot_09_nine_sector_epistasis_and_seesaw_cargos()` (`s13_seesaw_cargos.png`).
    """
    work = solve_archetype_a3_nine_sector_epistasis(
        df,
        delta_a_col=delta_a_col,
        delta_b_col=delta_b_col,
        fc_cutoff_log2=fc_cutoff_log2,
        label_a=label_a,
        label_b=label_b,
    )
    da = pd.to_numeric(work[delta_a_col], errors="coerce")
    db = pd.to_numeric(work[delta_b_col], errors="coerce")

    seesaw_mask = work["Divergence_Group"] == "G7_Opposite_Seesaw"
    relaxed_flag = False
    if int(seesaw_mask.sum()) < 3:
        seesaw_mask = ((da <= -0.40) & (db >= 0.40)) | ((da >= 0.40) & (db <= -0.40))
        relaxed_flag = True

    sub = work.loc[seesaw_mask].copy()
    sub["seesaw_relaxed_threshold"] = relaxed_flag
    sub["seesaw_subtype"] = np.where(
        pd.to_numeric(sub[delta_b_col], errors="coerce") > pd.to_numeric(sub[delta_a_col], errors="coerce"),
        f"G7-A: {label_a}_Down / {label_b}_Up (DI > 0)",
        f"G7-B: {label_a}_Up / {label_b}_Down (DI < 0)",
    )
    sub["seesaw_span_abs"] = sub["eff_DI_B_minus_A"].abs()
    if clone_a_cols:
        for c in clone_a_cols:
            if c in sub.columns:
                sub[c] = pd.to_numeric(sub[c], errors="coerce").round(4)
    if clone_b_cols:
        for c in clone_b_cols:
            if c in sub.columns:
                sub[c] = pd.to_numeric(sub[c], errors="coerce").round(4)

    return sub.sort_values("eff_DI_B_minus_A", ascending=True).reset_index(drop=True)


def solve_archetype_a4_biologics_humanization(
    df: pd.DataFrame,
    parental_col: str,
    humanized_col: str,
    isotype_col: str,
    mock_col: Optional[str] = None,
    fc_cutoff_log2: float = 0.5849625,
    equiv_margin_log2: float = 0.321928,
) -> pd.DataFrame:
    """
    Archetype A4: Biologics & Antibody Engineering 3-Axis Decomposition:
      1. Target Efficacy Axis: mean(Parental, Humanized) - Isotype_IgG
      2. Fc / IgG Background Axis: Isotype_IgG - Mock (if Mock provided)
      3. Humanization Drift vs Equivalence: Humanized - Parental (certified within [-equiv_margin, +equiv_margin])
    """
    out = df.copy()
    p_val = pd.to_numeric(out[parental_col], errors="coerce")
    h_val = pd.to_numeric(out[humanized_col], errors="coerce")
    i_val = pd.to_numeric(out[isotype_col], errors="coerce")

    d_par_vs_igg = p_val - i_val
    d_hum_vs_igg = h_val - i_val
    d_efficacy = 0.5 * (p_val + h_val) - i_val
    d_hum_drift = h_val - p_val

    out["delta_parental_vs_isotype_log2"] = d_par_vs_igg.round(4)
    out["delta_humanized_vs_isotype_log2"] = d_hum_vs_igg.round(4)
    out["axis1_target_efficacy_log2"] = d_efficacy.round(4)
    out["axis3_humanization_drift_log2"] = d_hum_drift.round(4)

    if mock_col and mock_col in out.columns:
        m_val = pd.to_numeric(out[mock_col], errors="coerce")
        out["axis2_fc_igg_background_log2"] = (i_val - m_val).round(4)

    both_up = (d_par_vs_igg >= fc_cutoff_log2) & (d_hum_vs_igg >= fc_cutoff_log2)
    both_dn = (d_par_vs_igg <= -fc_cutoff_log2) & (d_hum_vs_igg <= -fc_cutoff_log2)
    hum_equiv = d_hum_drift.abs() <= equiv_margin_log2

    cls = np.full(len(out), "Background", dtype=object)
    cls[both_up & hum_equiv] = "Conserved_Ab_Target_Up (Humanization_Equivalent)"
    cls[both_dn & hum_equiv] = "Conserved_Ab_Target_Down (Humanization_Equivalent)"
    cls[(both_up | both_dn) & (~hum_equiv)] = "Target_Hit_With_Humanization_Potency_Shift"
    cls[(~both_up) & (~both_dn) & (d_hum_drift.abs() >= fc_cutoff_log2)] = "Humanization_Specific_OffTarget_Drift"
    out["a4_biologics_class"] = cls
    return out


def solve_archetype_a5_clinical_sgi(
    df: pd.DataFrame,
    fluid_log2fc_col: str,
    tissue_log2fc_col: str,
    pep_col: str = "std_unique_peptides",
    fc_cutoff_log2: float = 0.5849625,
) -> pd.DataFrame:
    """
    Archetype A5: Clinical Paired Multi-Compartment Secretion Gradient Index (`SGI`).
    SGI = log2FC_Fluid(T/N) - log2FC_Tissue(T/N).
    """
    out = df.copy()
    f_fc = pd.to_numeric(out[fluid_log2fc_col], errors="coerce")
    t_fc = pd.to_numeric(out[tissue_log2fc_col], errors="coerce")
    peps = pd.to_numeric(out[pep_col], errors="coerce").fillna(1) if pep_col in out.columns else pd.Series(2, index=out.index)

    out["SGI_fluid_minus_tissue"] = (f_fc - t_fc).round(4)

    quad = np.full(len(out), "Unclassified", dtype=object)
    q1 = (f_fc >= fc_cutoff_log2) & (t_fc >= fc_cutoff_log2)
    q2 = (f_fc >= fc_cutoff_log2) & (t_fc < fc_cutoff_log2)
    q3 = (f_fc <= -fc_cutoff_log2) & (t_fc <= -fc_cutoff_log2)
    q4 = (t_fc >= fc_cutoff_log2) & (f_fc < fc_cutoff_log2)

    quad[q1] = "Q1_concordant_up"
    quad[q2] = "Q2_hyper_secreted"
    quad[q3] = "Q3_concordant_down"
    quad[q4] = "Q4_intracellular_retained"

    out["sgi_quadrant"] = quad
    out["sgi_clinical_tier"] = np.where(
        peps >= 2,
        quad + " [Core_ge2_peptides]",
        quad + " [1-Peptide_Exploratory]",
    )
    return out


def solve_modality_m2_ptm_stoichiometry(
    ptm_df: pd.DataFrame,
    prot_df: pd.DataFrame,
    ptm_fc_col: str = "log2FC",
    prot_fc_col: str = "log2FC",
    ptm_se_col: str = "se_moderated",
    prot_se_col: str = "se_moderated",
    gene_col: str = "std_gene",
    flanking_seq_col: Optional[str] = None,
    fc_cutoff_log2: float = 0.5849625,
) -> pd.DataFrame:
    """
    Modality M2 / Archetype A7: Uncouples PTM site abundance from parent protein abundance.
    """
    prot_lookup = (
        prot_df.dropna(subset=[gene_col])
        .drop_duplicates(subset=[gene_col])
        .set_index(gene_col)[[prot_fc_col, prot_se_col]]
        .rename(columns={prot_fc_col: "parent_prot_log2FC", prot_se_col: "parent_prot_SE"})
    )

    merged = ptm_df.merge(prot_lookup, left_on=gene_col, right_index=True, how="left")
    ptm_fc = pd.to_numeric(merged[ptm_fc_col], errors="coerce")
    prot_fc_raw = pd.to_numeric(merged["parent_prot_log2FC"], errors="coerce")
    parent_missing = prot_fc_raw.isna()
    prot_fc = prot_fc_raw
    ptm_se = pd.to_numeric(merged[ptm_se_col], errors="coerce").fillna(0.25)
    prot_se = pd.to_numeric(merged["parent_prot_SE"], errors="coerce").fillna(0.15)

    net_fc = ptm_fc - prot_fc
    net_se = np.sqrt(ptm_se**2 + prot_se**2)
    net_z = net_fc / np.maximum(net_se, 1e-6)
    net_p = 2.0 * stats.norm.sf(np.abs(net_z))

    merged["net_stoichiometry_log2FC"] = np.round(net_fc, 4)
    merged["net_stoichiometry_SE"] = np.round(net_se, 4)
    merged["net_stoichiometry_p_value"] = net_p

    stoich_class = np.full(len(merged), "Unchanged", dtype=object)
    true_up = (~parent_missing) & (net_fc >= fc_cutoff_log2) & (net_p < 0.05)
    true_dn = (~parent_missing) & (net_fc <= -fc_cutoff_log2) & (net_p < 0.05)
    abund_driven = (~parent_missing) & (ptm_fc.abs() >= fc_cutoff_log2) & (net_fc.abs() < 0.322)

    stoich_class[true_up] = "True_Stoichiometric_PTM_Gain"
    stoich_class[true_dn] = "True_Stoichiometric_PTM_Loss"
    stoich_class[abund_driven] = "Parent_Abundance_Driven_Shift"
    stoich_class[parent_missing] = "Parent_Unquantified"
    merged["ptm_stoichiometry_class"] = stoich_class

    if flanking_seq_col and flanking_seq_col in merged.columns:
        motifs = []
        for seq in merged[flanking_seq_col].fillna("").astype(str):
            s = seq.upper()
            center = len(s) // 2
            if len(s) >= 9 and center + 1 < len(s) and s[center + 1] == "P":
                motifs.append("Proline-Directed (CDK/MAPK/ERK/GSK3: [S/T]P)")
            elif len(s) >= 9 and center >= 3 and (s[center - 3] in {"R", "K"} and s[center - 2] in {"R", "K"}):
                motifs.append("Basophilic (AKT/PKA/PKC/RSK: R-R-X-[S/T])")
            elif len(s) >= 9 and center + 3 < len(s) and (s[center + 1 : center + 4].count("D") + s[center + 1 : center + 4].count("E") >= 2):
                motifs.append("Acidophilic (CK2/PLK: [S/T]-X-X-[D/E])")
            elif len(s) >= 9 and center + 1 < len(s) and s[center + 1] == "Q":
                motifs.append("DDR PIKK (ATM/ATR/DNA-PK: [S/T]Q)")
            else:
                motifs.append("Other / Non-Canonical Kinase Motif")
        merged["inferred_kinase_motif_family"] = motifs

    return merged


def scan_paralog_compensation(
    df: pd.DataFrame,
    target_gene: str | Sequence[str],
    fc_col: str = "log2FC",
    gene_col: str = "std_gene",
    sample_cols: Optional[Sequence[str]] = None,
    compensation_cutoff_log2: float = 0.40,
) -> pd.DataFrame:
    """
    Engine [03]: Paralog & Isoform Compensation Scanner.
    Combines:
      1. Curated functional paralog & pathway buffering families from `resources/curated_proteomics_signatures.json`
         (e.g., `NFE2L2 -> NFE2L1/FTH1/FTL`, `RAB27A -> RAB27B/RAB3A/PDCD6IP`,
         `CD63 -> CD9/CD81/CD151`).
      2. HGNC gene-symbol root prefix matching (`MAPK*`, `RAB*`, `SLC*`, `ITG*`, `S100*`).
    Also retains cross-sample intensity columns (`sample_cols`) with `<LOD` floor indicators so
    `plot_05_direct_target_and_paralog_verification()` (`s08_target_verification.png`) can plot
    both sample-by-sample trajectories and fold-change bars directly from this DataFrame!
    """
    targets = [target_gene] if isinstance(target_gene, str) else list(target_gene)
    targets_up = [str(t).strip().upper() for t in targets if str(t).strip()]
    curated_map = _load_paralog_families()

    candidate_set = set(targets_up)
    prefixes = []
    for t in targets_up:
        if t in curated_map:
            candidate_set.update(g.upper() for g in curated_map[t])
        m = re.match(r"^([A-Za-z]+)", t)
        pref = m.group(1) if m else t[:4]
        if len(pref) >= 3:
            prefixes.append(pref)

    sub = df.copy()
    if gene_col not in sub.columns:
        sub[gene_col] = sub.index.astype(str)
    sub[gene_col] = sub[gene_col].astype(str).str.strip().str.upper()
    if fc_col in sub.columns:
        sub[fc_col] = pd.to_numeric(sub[fc_col], errors="coerce")
    else:
        sub[fc_col] = 0.0

    prefix_mask = pd.Series(False, index=sub.index)
    for pref in prefixes:
        prefix_mask = prefix_mask | sub[gene_col].str.startswith(pref)

    fam_mask = sub[gene_col].isin(candidate_set) | prefix_mask
    fam = sub.loc[fam_mask].copy()
    if fam.empty:
        return fam

    is_target = fam[gene_col].isin(set(targets_up))
    fam["is_primary_target"] = is_target
    fam["paralog_response_class"] = np.where(
        is_target,
        "Primary_Perturbed_Target",
        np.where(
            fam[fc_col] >= compensation_cutoff_log2,
            "Compensatory_Paralog_Upregulation",
            np.where(
                fam[fc_col] <= -compensation_cutoff_log2,
                "Co_Depleted_or_Shared_Peptide_OffTarget",
                "Unchanged_Paralog",
            ),
        ),
    )

    fam["paralog_gene"] = fam[gene_col]
    fam["paralog_log2FC"] = fam[fc_col]
    fam["paralog_relation"] = fam["paralog_response_class"]

    if sample_cols:
        lod_map = compute_sample_lod_bounds(sub, sample_cols)
        for scol in sample_cols:
            if scol in fam.columns:
                raw_v = pd.to_numeric(fam[scol], errors="coerce")
                fam[f"{scol}_is_lod"] = raw_v.isna()
                fam[f"{scol}_plot_val"] = raw_v.fillna(lod_map.get(scol, 10.6))

    return fam.sort_values(["is_primary_target", fc_col], ascending=[False, False]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 3: Universal Archetype & Modality Solver")
    parser.add_argument("--archetype", required=True, choices=["A1", "A2", "A3", "A4", "A5", "M2"])
    parser.add_argument("--input", "-i", required=True, help="Input CSV")
    parser.add_argument("--output", "-o", required=True, help="Output CSV")
    args, unknown = parser.parse_known_args()
    print(f"[Stage 3] Loaded Archetype Solver ({args.archetype}) for {args.input}")


if __name__ == "__main__":
    main()

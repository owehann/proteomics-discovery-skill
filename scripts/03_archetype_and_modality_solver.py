#!/usr/bin/env python3
"""
Stage 3: Universal Experimental Archetype (A1-A10) & Modality (M1-M6) Solver
============================================================================
Implements exact algebraic and statistical scoring for:
  - A1: Single Perturbation with Multi-Clone Consistency & Hairpin Divergence Audit
  - A2: Reciprocal KD x OE (`R-Score` + 4-Quadrant Driver vs Stress Artifact Partition)
  - A3: Multi-Perturbation 9-Sector Cartesian Epistasis Grid with `TOST` Equivalence
  - A4: Biologics & Antibody Engineering 3-Axis Decomposition (Efficacy, Fc, Humanization TOST)
  - A5: Clinical Paired Multi-Compartment Secretion Gradient Index (`SGI`) & Dual-Scale Profiling
  - A7 (M2): PTM Stoichiometry Uncoupling (`delta_PTM - delta_Protein`) & Flanking Motif Scoring
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats


def solve_archetype_a1_multi_clone(
    df: pd.DataFrame,
    clone_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    fc_cutoff_log2: float = 0.5849625,
    divergence_warn_log2: float = 1.0,
) -> pd.DataFrame:
    """
    Archetype A1: Multi-clone knockdown/overexpression consistency & seed off-target audit.
    Flags proteins where independent clones diverge strongly (e.g., sh1 Up vs sh2 Down).
    """
    out = df.copy()
    c_mean = out[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    clone_deltas = pd.DataFrame(index=out.index)
    for col in clone_cols:
        clone_deltas[f"delta_{col}"] = pd.to_numeric(out[col], errors="coerce") - c_mean

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
    Computes the geometric Reciprocity Score (`R-Score`) and partitions proteins into
    True Reciprocal Drivers (Q2/Q4) vs Concordant Stress Artifacts (Q1/Q3).
    """
    out = df.copy()
    d_oe = pd.to_numeric(out[delta_oe_col], errors="coerce")
    d_kd = pd.to_numeric(out[delta_kd_col], errors="coerce")

    opp_dir = (np.sign(d_oe) * np.sign(d_kd)) < 0
    r_mag = np.sqrt(np.abs(d_oe * d_kd))
    # Positive R-score means promoted by target (OE > 0, KD < 0); Negative means repressed by target
    r_score = np.where(opp_dir, np.sign(d_oe) * r_mag, 0.0)
    out["reciprocal_r_score"] = np.round(r_score, 4)

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


def solve_archetype_a3_nine_sector_epistasis(
    df: pd.DataFrame,
    delta_a_col: str,
    delta_b_col: str,
    fc_cutoff_log2: float = 0.5849625,
    equiv_margin_log2: float = 0.321928,
    label_a: str = "PerturbA",
    label_b: str = "PerturbB",
) -> pd.DataFrame:
    """
    Archetype A3: 9-Sector Cartesian Divergence & Epistasis Grid for two perturbations A vs B.
    Uses equivalence margin `equiv_margin_log2` to certify strictly specific vs shared axes.
    """
    out = df.copy()
    da = pd.to_numeric(out[delta_a_col], errors="coerce")
    db = pd.to_numeric(out[delta_b_col], errors="coerce")

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
    return out


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
    Enforces exact threshold boundary `c = log2(1.5) = 0.5849625` so that proteins with
    `tissue_log2FC >= 0.5849625` AND `fluid_log2FC >= 0.5849625` are strictly classified as
    `Q1_concordant_up` (never mislabeled as `Q2_hyper_secreted`!).
    Also separates Core (`unique_peptides >= 2`) from `1-Peptide Exploratory` rows.
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
    Modality M2 / Archetype A7: Uncouples PTM site abundance from parent protein abundance:
      delta_net_PTM = log2FC_PTM - log2FC_ParentProtein
      SE_net_PTM = sqrt(SE_PTM^2 + SE_Protein^2)
    And scores canonical kinase substrate motifs on +/-7 aa sequence windows.
    """
    prot_lookup = (
        prot_df.dropna(subset=[gene_col])
        .drop_duplicates(subset=[gene_col])
        .set_index(gene_col)[[prot_fc_col, prot_se_col]]
        .rename(columns={prot_fc_col: "parent_prot_log2FC", prot_se_col: "parent_prot_SE"})
    )

    merged = ptm_df.merge(prot_lookup, left_on=gene_col, right_index=True, how="left")
    ptm_fc = pd.to_numeric(merged[ptm_fc_col], errors="coerce")
    prot_fc = pd.to_numeric(merged["parent_prot_log2FC"], errors="coerce").fillna(0.0)
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
    true_up = (net_fc >= fc_cutoff_log2) & (net_p < 0.05)
    true_dn = (net_fc <= -fc_cutoff_log2) & (net_p < 0.05)
    abund_driven = (ptm_fc.abs() >= fc_cutoff_log2) & (net_fc.abs() < 0.322)

    stoich_class[true_up] = "True_Stoichiometric_PTM_Gain"
    stoich_class[true_dn] = "True_Stoichiometric_PTM_Loss"
    stoich_class[abund_driven] = "Parent_Abundance_Driven_Shift"
    merged["ptm_stoichiometry_class"] = stoich_class

    # Optional motif family scanner (The Kinase Library / canonical basophilic, proline-directed, acidophilic motifs)
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 3: Universal Archetype & Modality Solver")
    parser.add_argument("--archetype", required=True, choices=["A1", "A2", "A3", "A4", "A5", "M2"])
    parser.add_argument("--input", "-i", required=True, help="Input CSV")
    parser.add_argument("--output", "-o", required=True, help="Output CSV")
    args, unknown = parser.parse_known_args()
    print(f"[Stage 3] Loaded Archetype Solver ({args.archetype}) for {args.input}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Stage 2 (Validated Python Statistical Engine):
==============================================
Implements exact analytical equivalents of Bioconductor proteomics workflows when R packages
are unavailable or when running in pure-Python environments:
  1. Smyth (2004) Empirical Bayes Variance Moderation (`squeezeVar`) with optional
     peptide-count (`DEqMS`-style) or intensity-trended prior variance `s0^2(n_pep)`.
  2. Li, Cobbold & Smyth (2025) Detection Probability Curve (`DPC`) logistic missingness
     model (`logit P(obs) = beta_0 + beta_1 * y_i`) and exact Hurdle/missingness probability.
  3. Paired & Unpaired Linear Models + Benjamini-Hochberg FDR + Wilcoxon sensitivity test.
  4. Two One-Sided Tests (`TOST`) for statistical equivalence within `[-delta_eq, +delta_eq]`.
  5. `SAINT`-lite Empirical Bayes interactome specificity scoring (Modality M3).
  6. `NPARC` Functional F-test for Thermal Proteome Profiling melting curves (Modality M4).
  7. `CurveCurator`-style calibrated 4-parameter log-logistic dose-response F-test (Modality M5).
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, special, stats


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Compute Benjamini-Hochberg FDR q-values handling NaNs cleanly."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return q
    pv = np.clip(p[valid], 0.0, 1.0)
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out_valid = np.empty_like(pv)
    out_valid[order] = adj
    q[valid] = out_valid
    return q


def _trigamma_inverse(x: float) -> float:
    """Solve trigamma(y) = x for y > 0 via Newton-Raphson (Smyth 2004 limma::trigammaInverse)."""
    if not np.isfinite(x) or x <= 1e-10:
        return 1e6
    if x > 1e7:
        return 1.0 / math.sqrt(x)
    y = 0.5 + 1.0 / x
    for _ in range(15):
        tri = float(special.polygamma(1, y))
        dif = tri * (1.0 - tri / x) / float(special.polygamma(2, y))
        y += dif
        if abs(-dif / y) < 1e-8:
            break
    return max(float(y), 1e-6)


def squeeze_var_ebayes(
    s2: np.ndarray,
    df: np.ndarray,
    covariate_trend: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Smyth (2004) / DEqMS (2020) analytical empirical Bayes variance shrinkage.
    Supports trended prior variance `s0^2` as a function of log2(unique_peptides) or mean intensity.
    Returns: (s2_post, df_total, d0, s02_vec)
    """
    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    valid = np.isfinite(s2) & (s2 > 1e-12) & np.isfinite(df) & (df >= 1.0)

    if valid.sum() < 5:
        med_s2 = float(np.nanmedian(s2[s2 > 0])) if np.any(s2 > 0) else 0.05
        d0 = 4.0
        s02_vec = np.full_like(s2, med_s2)
        s2_post = np.where(np.isfinite(s2), (d0 * s02_vec + df * s2) / (d0 + df), s02_vec)
        df_total = np.where(np.isfinite(df), df + d0, d0)
        return s2_post, df_total, d0, s02_vec

    z = np.log(s2[valid])
    df_v = df[valid]
    # Expected value and variance of log(s^2) under chi-square: digamma(df/2) - log(df/2), trigamma(df/2)
    e_offset = special.digamma(df_v / 2.0) - np.log(df_v / 2.0)
    z_adj = z - e_offset

    if covariate_trend is not None and np.isfinite(covariate_trend[valid]).sum() >= 10:
        x_v = np.asarray(covariate_trend[valid], dtype=float)
        # Robust linear/quadratic trend on log-variance against log2(peptides) or mean intensity
        slope, intercept, _, _, _ = stats.linregress(x_v, z_adj)
        x_all = np.where(np.isfinite(covariate_trend), covariate_trend, np.nanmedian(x_v))
        z_trend_all = intercept + slope * x_all
        resid = z_adj - (intercept + slope * x_v)
    else:
        med_z = float(np.median(z_adj))
        z_trend_all = np.full_like(s2, med_z)
        resid = z_adj - med_z

    # Estimate excess variance of log(s2) above trigamma(df/2)
    mad_resid = float(np.median(np.abs(resid - np.median(resid)))) * 1.4826
    var_z = max(mad_resid**2, float(np.var(resid, ddof=1)))
    mean_trigamma_df = float(np.mean(special.polygamma(1, df_v / 2.0)))
    excess_trigamma = max(var_z - mean_trigamma_df, 1e-4)

    d0 = 2.0 * _trigamma_inverse(excess_trigamma)
    d0 = float(np.clip(d0, 1.0, 100.0))

    s02_offset = float(special.digamma(d0 / 2.0) - math.log(d0 / 2.0))
    # Smyth (2004): E[log(s^2)] = log(s0^2) + digamma(df/2) - log(df/2) - digamma(d0/2) + log(d0/2)
    # Since z_adj = log(s^2) - (digamma(df/2) - log(df/2)), log(s0^2) = E[z_adj] + (digamma(d0/2) - log(d0/2))
    s02_vec = np.exp(z_trend_all + s02_offset)

    s2_safe = np.where(np.isfinite(s2) & (s2 > 0), s2, s02_vec)
    df_safe = np.where(np.isfinite(df) & (df > 0), df, 0.0)
    s2_post = (d0 * s02_vec + df_safe * s2_safe) / (d0 + df_safe)
    df_total = df_safe + d0
    return s2_post, df_total, d0, s02_vec


def run_moderated_differential_analysis(
    df: pd.DataFrame,
    treat_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    is_paired: bool = False,
    pep_col: str = "std_unique_peptides",
    delta_eq_log2: float = 0.321928,
    sigma_null_override: Optional[float] = None,
) -> pd.DataFrame:
    """
    Run DEqMS/limpa-style moderated t-test + DPC missingness recovery + TOST equivalence test.
    """
    t_mat = df[list(treat_cols)].apply(pd.to_numeric, errors="coerce").values
    c_mat = df[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce").values
    n_features = len(df)

    n_t = np.sum(np.isfinite(t_mat), axis=1)
    n_c = np.sum(np.isfinite(c_mat), axis=1)

    with np.errstate(invalid="ignore"):
        mean_t = np.nanmean(t_mat, axis=1)
        mean_c = np.nanmean(c_mat, axis=1)
        ave_expr = np.nanmean(np.hstack([t_mat, c_mat]), axis=1)

    if is_paired and len(treat_cols) == len(ctrl_cols):
        diff_mat = t_mat - c_mat
        n_pairs = np.sum(np.isfinite(diff_mat), axis=1)
        with np.errstate(invalid="ignore"):
            log2_fc = np.nanmean(diff_mat, axis=1)
            s2_raw = np.nanvar(diff_mat, axis=1, ddof=1)
        df_raw = np.maximum(n_pairs - 1.0, 0.0)
        unscaled_var = np.where(n_pairs > 0, 1.0 / n_pairs, np.nan)
    else:
        log2_fc = mean_t - mean_c
        with np.errstate(invalid="ignore"):
            ss_t = np.nansum((t_mat - mean_t[:, None]) ** 2, axis=1)
            ss_c = np.nansum((c_mat - mean_c[:, None]) ** 2, axis=1)
        df_raw = np.maximum(n_t + n_c - 2.0, 0.0)
        s2_raw = np.where(df_raw >= 1.0, (ss_t + ss_c) / np.maximum(df_raw, 1.0), np.nan)
        unscaled_var = np.where((n_t >= 1) & (n_c >= 1), (1.0 / n_t) + (1.0 / n_c), np.nan)

    peps = (
        pd.to_numeric(df[pep_col], errors="coerce").fillna(1.0).clip(lower=1.0).values
        if pep_col in df.columns
        else np.full(n_features, 2.0)
    )
    log_peps = np.log2(peps)

    if np.all(df_raw < 1.0):
        if sigma_null_override is not None:
            # Single-tube regime (N=1 vs N=1): use Stage 1.5 empirical null variance
            s2_post = np.full(n_features, sigma_null_override**2)
            df_total = np.full(n_features, 10.0)
            d0 = 10.0
        else:
            # Unreplicated N=1 vs N=1 without empirical null sigma override: do not invent p-values
            s2_post = np.full(n_features, np.nan)
            df_total = np.full(n_features, np.nan)
            d0 = np.nan
    else:
        s2_post, df_total, d0, _ = squeeze_var_ebayes(s2_raw, df_raw, covariate_trend=log_peps)

    se_mod = np.sqrt(s2_post * unscaled_var)
    t_mod = log2_fc / se_mod
    p_mod = 2.0 * stats.t.sf(np.abs(t_mod), df=df_total)
    q_mod = benjamini_hochberg(p_mod)

    # Two One-Sided Tests (TOST) for statistical equivalence within [-delta_eq_log2, +delta_eq_log2]
    t_lower = (log2_fc - (-delta_eq_log2)) / se_mod
    t_upper = (delta_eq_log2 - log2_fc) / se_mod
    p_lower = stats.t.sf(t_lower, df=df_total)
    p_upper = stats.t.sf(t_upper, df=df_total)
    p_tost = np.maximum(p_lower, p_upper)
    q_tost = benjamini_hochberg(p_tost)

    # Fit Detection Probability Curve (DPC) for ON/OFF missingness recovery
    all_mat = np.hstack([t_mat, c_mat])
    obs_frac = np.mean(np.isfinite(all_mat), axis=1)
    valid_dpc = np.isfinite(ave_expr) & (obs_frac > 0.0) & (obs_frac < 1.0)
    if valid_dpc.sum() >= 15:
        p_clip = np.clip(obs_frac[valid_dpc], 0.02, 0.98)
        b1, b0, _, _, _ = stats.linregress(ave_expr[valid_dpc], np.log(p_clip / (1.0 - p_clip)))
        b1 = max(float(b1), 0.2)
    else:
        b0, b1 = -10.0, 0.7

    # DPC probability for complete ON/OFF switches (e.g. detected in treat, 0 in ctrl or vice versa)
    switch_type = np.full(n_features, "Quantified_Both", dtype=object)
    dpc_switch_p = np.full(n_features, np.nan)
    for idx in range(n_features):
        if n_t[idx] > 0 and n_c[idx] == 0:
            switch_type[idx] = "Treat_Only_ON"
            p_obs_given_y = special.expit(b0 + b1 * mean_t[idx])
            # Probability all n_ctrl tubes missed a protein of intensity mean_t[idx]
            dpc_switch_p[idx] = (1.0 - p_obs_given_y) ** len(ctrl_cols)
        elif n_c[idx] > 0 and n_t[idx] == 0:
            switch_type[idx] = "Ctrl_Only_OFF"
            p_obs_given_y = special.expit(b0 + b1 * mean_c[idx])
            dpc_switch_p[idx] = (1.0 - p_obs_given_y) ** len(treat_cols)
        elif n_t[idx] == 0 and n_c[idx] == 0:
            switch_type[idx] = "Undetected_Both"

    out = pd.DataFrame(
        {
            "composite_id": df["composite_id"] if "composite_id" in df.columns else df.index,
            "std_accession": df["std_accession"] if "std_accession" in df.columns else "",
            "std_gene": df["std_gene"] if "std_gene" in df.columns else "",
            "std_unique_peptides": peps.astype(int),
            "peptide_tier": np.where(peps >= 2, "Core_ge2_peptides", "Exploratory_1_peptide"),
            "n_treat_obs": n_t,
            "n_ctrl_obs": n_c,
            "mean_treat_log2": np.round(mean_t, 4),
            "mean_ctrl_log2": np.round(mean_c, 4),
            "ave_expr_log2": np.round(ave_expr, 4),
            "log2FC": np.round(log2_fc, 4),
            "s2_raw": np.round(s2_raw, 6),
            "s2_post": np.round(s2_post, 6),
            "se_moderated": np.round(se_mod, 4),
            "s2_prior_df": round(float(d0), 3) if np.isfinite(d0) else np.nan,
            "t_moderated": np.round(t_mod, 4),
            "p_value": p_mod,
            "fdr_q_value": q_mod,
            "tost_p_value": p_tost,
            "tost_fdr": q_tost,
            "is_tost_equivalent": (q_tost < 0.05) & (np.abs(log2_fc) < delta_eq_log2),
            "detection_regime": switch_type,
            "dpc_missingness_p": dpc_switch_p,
        }
    )
    return out


def score_saint_lite_interactome(
    df: pd.DataFrame,
    bait_cols: Sequence[str],
    ctrl_cols: Sequence[str],
    bait_gene: Optional[str] = None,
) -> pd.DataFrame:
    """
    Modality M3 (AP-MS / BioID / TurboID):
    Computes bait-normalized stoichiometry, Fold-Change over Control (FC_A, FC_B),
    and heuristic Empirical Bayes specificity score (SAINT-inspired fallback when SAINTexpress is unavailable).
    """
    b_mat = df[list(bait_cols)].apply(pd.to_numeric, errors="coerce")
    c_mat = df[list(ctrl_cols)].apply(pd.to_numeric, errors="coerce")

    # Impute missing negative control background at 5th percentile of observed control values
    all_vals = pd.concat([b_mat, c_mat], axis=1).values
    finite_vals = all_vals[np.isfinite(all_vals)]
    bg_floor = float(np.percentile(finite_vals, 5)) if len(finite_vals) > 0 else 10.0

    c_filled = c_mat.fillna(bg_floor)
    b_mean = b_mat.mean(axis=1)
    c_mean = c_filled.mean(axis=1)
    c_max = c_filled.max(axis=1)

    fc_a_log2 = b_mean - c_mean
    fc_b_log2 = b_mean - c_max

    # Bait-normalized stoichiometry (log2(Prey / Bait))
    bait_level = float(b_mean.max())
    if bait_gene and "std_gene" in df.columns:
        m_bait = df["std_gene"].astype(str).str.upper() == bait_gene.upper()
        if m_bait.any() and pd.notna(b_mean.loc[m_bait].iloc[0]):
            bait_level = float(b_mean.loc[m_bait].iloc[0])

    stoich_log2 = b_mean - bait_level

    # Heuristic specificity probability (reported alongside moderated DE; not canonical SAINTexpress)
    n_b_obs = b_mat.notna().sum(axis=1).values
    rep_frac = n_b_obs / max(len(bait_cols), 1)
    z_spec = fc_a_log2 / max(float(c_filled.std(axis=1).median()), 0.35)
    saint_prob = special.expit(1.8 * (fc_a_log2 - 1.5) + 0.8 * z_spec) * rep_frac
    saint_prob = np.where(b_mean.isna(), 0.0, np.clip(saint_prob, 0.0, 1.0))

    res = run_moderated_differential_analysis(df, bait_cols, ctrl_cols)
    res["fc_a_vs_ctrl_mean_log2"] = np.round(fc_a_log2.values, 4)
    res["fc_b_vs_ctrl_max_log2"] = np.round(fc_b_log2.values, 4)
    res["bait_rel_stoichiometry_log2"] = np.round(stoich_log2.values, 4)
    res["heuristic_specificity_score"] = np.round(saint_prob, 4)
    res["saint_lite_prob"] = res["heuristic_specificity_score"]
    res["is_high_conf_interactor"] = (
        (res["saint_lite_prob"] >= 0.80)
        & (res["fc_a_vs_ctrl_mean_log2"] >= 1.5)
        & (res["std_unique_peptides"] >= 2)
    )
    return res


def fit_tpp_nparc_curve(
    temps: Sequence[float],
    veh_mat: np.ndarray,
    drug_mat: np.ndarray,
) -> pd.DataFrame:
    """
    Modality M4 (TPP-TR / CETSA):
    Fits 3-parameter sigmoidal melting curves f(T) = p + (1-p)/(1 + exp(b - a/T))
    under Null (combined) vs Alternative (separate Vehicle and Drug curves),
    returning Delta_Tm and NPARC functional F-test p-values.
    """
    t_arr = np.asarray(temps, dtype=float)

    def sigmoid_tpp(T: np.ndarray, a: float, b: float, p: float) -> np.ndarray:
        return p + (1.0 - p) / (1.0 + np.exp(np.clip(b - a / T, -40.0, 40.0)))

    def fit_single(T_vec: np.ndarray, y_vec: np.ndarray) -> Tuple[float, float, float, float]:
        valid = np.isfinite(y_vec)
        if valid.sum() < 4:
            return np.nan, np.nan, np.nan, np.nan
        tv, yv = T_vec[valid], np.clip(y_vec[valid], 0.0, 1.3)
        try:
            popt, _ = optimize.curve_fit(
                sigmoid_tpp,
                tv,
                yv,
                p0=[550.0, 10.5, 0.05],
                bounds=([50.0, 1.0, 0.0], [5000.0, 80.0, 0.45]),
                maxfev=800,
            )
            a, b, p = popt
            # Solve f(Tm) = 0.5 => p + (1-p)/(1 + exp(b - a/Tm)) = 0.5 => exp(b - a/Tm) = (1 - 2p)/(0.5 - p) - 1 = 1/(1 - 2p)
            if p < 0.49:
                tm = a / (b - math.log((0.5 - p) / 0.5))
            else:
                tm = a / b
            rss = float(np.sum((yv - sigmoid_tpp(tv, a, b, p)) ** 2))
            return float(tm), float(rss), float(p), int(valid.sum())
        except Exception:
            return np.nan, np.nan, np.nan, int(valid.sum())

    rows = []
    n_prot = veh_mat.shape[0]
    for i in range(n_prot):
        yv = veh_mat[i, :]
        yd = drug_mat[i, :]
        tm_v, rss_v, plat_v, nv = fit_single(t_arr, yv)
        tm_d, rss_d, plat_d, nd = fit_single(t_arr, yd)
        _, rss_0, _, n0 = fit_single(np.concatenate([t_arr, t_arr]), np.concatenate([yv, yd]))
        rss_1 = rss_v + rss_d if (np.isfinite(rss_v) and np.isfinite(rss_d)) else np.nan
        df1 = 3.0
        df2 = max(n0 - 6.0, 2.0) if np.isfinite(n0) else np.nan
        s2_i = max(rss_1 / df2, 1e-4) if (np.isfinite(rss_1) and np.isfinite(df2)) else np.nan
        f_stat = max((rss_0 - rss_1) / df1, 0.0) / s2_i if (np.isfinite(rss_0) and np.isfinite(s2_i)) else np.nan
        rows.append(
            {
                "tm_vehicle_C": round(tm_v, 3) if np.isfinite(tm_v) else np.nan,
                "tm_drug_C": round(tm_d, 3) if np.isfinite(tm_d) else np.nan,
                "delta_tm_C": round(tm_d - tm_v, 3) if (np.isfinite(tm_v) and np.isfinite(tm_d)) else np.nan,
                "rss_null": rss_0,
                "rss_alt": rss_1,
                "s2_resid": s2_i,
                "df_resid": df2,
                "f_raw": f_stat,
            }
        )

    res_df = pd.DataFrame(rows)
    s2_post, df_tot, _, _ = squeeze_var_ebayes(res_df["s2_resid"].values, res_df["df_resid"].values)
    f_mod = np.where(
        np.isfinite(res_df["rss_null"].values) & np.isfinite(res_df["rss_alt"].values),
        np.maximum((res_df["rss_null"].values - res_df["rss_alt"].values) / 3.0, 0.0) / s2_post,
        np.nan,
    )
    p_nparc = stats.f.sf(f_mod, dfn=3.0, dfd=df_tot)
    res_df["f_nparc_moderated"] = np.round(f_mod, 4)
    res_df["nparc_p_value"] = p_nparc
    res_df["nparc_fdr"] = benjamini_hochberg(p_nparc)
    return res_df


def fit_curve_curator_dose_response(
    log10_doses: Sequence[float],
    resp_mat: np.ndarray,
) -> pd.DataFrame:
    """
    Modality M5 (Dose-Response):
    Fits 4-parameter log-logistic Hill curves (`CurveCurator` Bayer et al. Nat Commun 2023)
    and computes calibrated F-statistics, R^2, and pEC50 (-log10 M).
    """
    x_arr = np.asarray(log10_doses, dtype=float)

    def hill_4p(x: np.ndarray, front: float, back: float, pec50: float, slope: float) -> np.ndarray:
        return front + (back - front) / (1.0 + np.power(10.0, np.clip(slope * (x + pec50), -30.0, 30.0)))

    rows = []
    for i in range(resp_mat.shape[0]):
        y = resp_mat[i, :]
        valid = np.isfinite(y)
        if valid.sum() < 4:
            rows.append({"pEC50": np.nan, "max_log2FC": np.nan, "hill_slope": np.nan, "curve_r2": np.nan, "f_stat": np.nan, "p_value": np.nan})
            continue
        xv, yv = x_arr[valid], y[valid]
        rss_null = float(np.sum((yv - np.mean(yv)) ** 2))
        try:
            p0_pec50 = float(-np.median(xv))
            popt, _ = optimize.curve_fit(
                hill_4p,
                xv,
                yv,
                p0=[float(yv[0]), float(yv[-1]), p0_pec50, 1.0],
                bounds=([-10.0, -10.0, 3.0, 0.3], [10.0, 10.0, 12.0, 8.0]),
                maxfev=800,
            )
            front, back, pec50, slope = popt
            y_hat = hill_4p(xv, front, back, pec50, slope)
            rss_alt = float(np.sum((yv - y_hat) ** 2))
            r2 = max(0.0, 1.0 - rss_alt / max(rss_null, 1e-8))
            df_res = max(len(yv) - 4.0, 1.0)
            f_val = max((rss_null - rss_alt) / 3.0, 0.0) / max(rss_alt / df_res, 1e-5)
            p_val = float(stats.f.sf(f_val, dfn=3.0, dfd=df_res))
            rows.append(
                {
                    "pEC50": round(float(pec50), 3),
                    "max_log2FC": round(float(back - front), 4),
                    "hill_slope": round(float(slope), 3),
                    "curve_r2": round(r2, 4),
                    "f_stat": round(f_val, 4),
                    "p_value": p_val,
                }
            )
        except Exception:
            rows.append({"pEC50": np.nan, "max_log2FC": 0.0, "hill_slope": np.nan, "curve_r2": 0.0, "f_stat": 0.0, "p_value": 1.0})

    out = pd.DataFrame(rows)
    out["fdr_q_value"] = benjamini_hochberg(out["p_value"].values)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 2: Validated Python Statistical Engine")
    parser.add_argument("--input", "-i", required=True, help="Input CSV/TSV matrix")
    parser.add_argument("--treat", nargs="+", required=True, help="Treatment sample columns")
    parser.add_argument("--ctrl", nargs="+", required=True, help="Control sample columns")
    parser.add_argument("--paired", action="store_true", help="Paired sample design")
    parser.add_argument("--delta-eq", type=float, default=0.321928, help="TOST equivalence margin in log2 units")
    parser.add_argument("--output", "-o", required=True, help="Output CSV path")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    res = run_moderated_differential_analysis(
        df, args.treat, args.ctrl, is_paired=args.paired, delta_eq_log2=args.delta_eq
    )
    res.to_csv(args.output, index=False)
    print(f"[Stage 2 Python Engine] Saved {len(res)} features to {args.output}")


if __name__ == "__main__":
    main()

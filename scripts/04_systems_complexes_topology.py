#!/usr/bin/env python3
"""
Stage 4: Systems Biology, CORUM Macromolecular Complexes, Biophysical Topology & STRING KDA
===========================================================================================
Implements Engines [04]-[08] end-to-end directly from standard protein tables:
  - Engine [04]: Custom Detected-Proteome Background ORA + Greedy Jaccard Redundancy De-duplication
                 (`run_detected_background_ora`, `cluster_redundant_pathways_jaccard`,
                  `run_multi_group_comparative_ora`, `load_comprehensive_gene_sets`)
  - Engine [05]: CORUM Macromolecular Complex Stoichiometry & Co-Shift Permutation Audit
                 (`audit_corum_complex_co_shifts`, `extract_complex_subunit_matrix`)
  - Engine [06]: STRING Physical/Functional PPI Network Construction, Degree-Corrected RWR,
                 1-Hop Key Driver Analysis (KDA), & Cross-Sector Lock-and-Key Bridges
                 (`fetch_string_ppi_network`, `degree_corrected_rwr`, `run_key_driver_analysis_kda`)
  - Engine [07]: Direct Protein-Table Biophysical & Membrane Topology Enrichment (UniProt/HPA
                 Signal Peptide, Secreted, ECM, 1-TM, Multi-TM, 4-TM MARVEL/Tetraspanin,
                 GPI-anchor, Lipidated Raft Fisher's Exact Tests + Ectodomain Shedding)
                 (`annotate_uniprot_membrane_topology`, `audit_membrane_and_secretory_topology`,
                  `audit_ectodomain_shedding_topology`)
  - Engine [08]: Signed Univariate Linear Model (ULM) Upstream TF (CollecTRI) & Kinase Activity
                 Inference with Modality Scope Gate (`infer_upstream_regulator_ulm`)
"""

from __future__ import annotations
import argparse
import json
import os
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats


SKILL_ROOT = Path(__file__).resolve().parent.parent
SIGNATURES_PATH = SKILL_ROOT / "resources" / "curated_proteomics_signatures.json"


def _load_curated_signatures() -> Dict[str, Any]:
    if SIGNATURES_PATH.exists():
        try:
            with open(SIGNATURES_PATH, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass
    return {}


def _locate_external_data_dir(hint_dir: Optional[str | Path] = None) -> Optional[Path]:
    """
    Search for the local `external_data/` repository containing UniProt, HPA, STRING,
    CollecTRI, CPTAC, DepMap, and MSigDB/GO gene set caches.
    """
    candidates: List[Path] = []
    env_dir = os.environ.get("PROTEOSKILL_EXTERNAL_DATA")
    if env_dir:
        candidates.append(Path(env_dir))
    if hint_dir:
        p = Path(hint_dir).resolve()
        candidates.extend([p, p / "external_data", p.parent / "external_data", p.parent.parent / "external_data"])
    cwd = Path.cwd().resolve()
    candidates.extend([
        cwd / "external_data",
        cwd.parent / "external_data",
        SKILL_ROOT.parent.parent.parent / "external_data",
    ])
    for c in candidates:
        if c.exists() and c.is_dir():
            return c
    return None


# Built-in fallback macromolecular complexes (augmented by curated_proteomics_signatures.json)
CORE_MACROMOLECULAR_COMPLEXES: Dict[str, List[str]] = {
    "26S_Proteasome_Core_and_Lid": [
        "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7",
        "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7",
        "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6",
        "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD6", "PSMD7", "PSMD8", "PSMD11", "PSMD12", "PSMD13", "PSMD14"
    ],
    "Integrin_Tetraspanin_Web": [
        "ITGA1", "ITGA2", "ITGA3", "ITGA5", "ITGA6", "ITGAV", "ITGB1", "ITGB3", "ITGB4", "ITGB5",
        "CD9", "CD63", "CD81", "CD82", "CD151", "TSPAN1", "TSPAN4", "TSPAN8", "TSPAN9", "TSPAN14", "TSPAN15", "PTGFRN", "IGSF8", "ADAM10"
    ],
    "ESCRT_MVB_Biogenesis_Machinery": [
        "HGS", "STAM", "STAM2", "TSG101", "VPS28", "VPS37A", "VPS37B", "MVB12A",
        "VPS25", "VPS36", "SNF8", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP6", "PDCD6IP", "VPS4A", "VPS4B", "SDCBP"
    ],
    "Rab_SNARE_Exocytic_Fusion_Complex": [
        "RAB3A", "RAB27A", "RAB27B", "RAB11A", "RAB11B", "RAB7A", "RAB5A", "RAB8A", "RAB10",
        "SYNGR2", "VAMP2", "VAMP3", "VAMP7", "VAMP8", "STX1A", "STX3", "STX4", "STX7", "SNAP23", "SNAP25", "STXBP1", "STXBP2"
    ],
    "MARVEL_Raft_Microdomain_Complex": [
        "SYNGR2", "MYADM", "MAL2", "PLLP", "CMTM6", "FLOT1", "FLOT2", "CAV1", "STOM", "ERLIN1", "ERLIN2"
    ],
    "ECM_Basement_Membrane_Corona": [
        "FN1", "THBS1", "LGALS3BP", "GPC3", "MFGE8", "SPP1", "TGFBI", "LAMA5", "LAMB1", "LAMC1", "COL6A1", "COL6A2", "COL6A3", "HSPG2", "NID1", "NID2", "AGRIN"
    ],
    "COPII_Anterograde_ER_Golgi_Coat": [
        "SAR1A", "SAR1B", "SEC23A", "SEC23B", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC13", "SEC31A"
    ],
    "COPI_Retrograde_Coatomer": [
        "COPA", "COPB1", "COPB2", "COPD", "COPE", "COPG1", "COPG2", "COPZ1", "ARF1"
    ],
    "Exocyst_Tethering_Complex": [
        "EXOC1", "EXOC2", "EXOC3", "EXOC4", "EXOC5", "EXOC6", "EXOC7", "EXOC8"
    ],
    "Cytoplasmic_Ribosome_80S": [
        "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL10", "RPL11", "RPL13", "RPL15",
        "RPS2", "RPS3", "RPS4X", "RPS6", "RPS8", "RPS9", "RPS14", "RPS16", "RPS18", "RPS20"
    ],
    "Ferritin_Iron_Storage_Autophagy_Complex": [
        "FTH1", "FTL", "NCOA4", "TFRC", "SLC3A2", "SLC7A11", "HMOX1", "FDFT1", "GPX4", "ACSL4", "GCLM", "NQO1", "TXNRD1"
    ],
    "Chromatin_Methylation_Helicase_Complex": [
        "HELLS", "DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "SMARCA4", "SMARCA5", "CBX3", "CBX5", "SUV39H1", "EZH2", "SUZ12", "EED"
    ],
    "MCM_CMG_Replisome_Helicase": [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC45", "GINS1", "GINS2", "GINS3", "GINS4", "PCNA", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5"
    ],
    "ISGF3_Interferon_Antiviral_Signaling": [
        "STAT1", "STAT2", "IRF9", "ISG15", "PARP9", "DTX3L", "PARP14", "OAS1", "OAS2", "OAS3", "MX1", "IFIT1", "IFIT2", "IFIT3", "BST2", "TAP1", "TAP2"
    ],
}

_SIGS = _load_curated_signatures()
if "corum_complexes" in _SIGS and isinstance(_SIGS["corum_complexes"], dict):
    for _k, _v in _SIGS["corum_complexes"].items():
        if isinstance(_v, list) and len(_v) >= 3:
            CORE_MACROMOLECULAR_COMPLEXES[_k] = [str(g).upper() for g in _v]


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return q
    pv = np.clip(p[valid], 0.0, 1.0)
    n = len(pv)
    order = np.argsort(pv)
    adj = pv[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out_v = np.empty_like(pv)
    out_v[order] = np.clip(adj, 0.0, 1.0)
    q[valid] = out_v
    return q


# =============================================================================
# ENGINE [04]: Custom Detected-Background ORA & Jaccard De-duplication
# =============================================================================

def load_comprehensive_gene_sets(
    detected_universe: Optional[Sequence[str]] = None,
    external_data_dir: Optional[str | Path] = None,
    include_go_msigdb: bool = True,
) -> Dict[str, List[str]]:
    """
    Assemble a comprehensive, multi-ontology gene set dictionary combining:
      1. `resources/curated_proteomics_signatures.json` (`pathway_gene_sets`, `corum_complexes`,
         `subcellular_compartments`, `misev2023_markers`)
      2. `external_data/genesets/gs_*.txt` (`Hallmark`, `KEGG`, `Reactome`, `GO_BP`, `GO_CC`, `GO_MF`)
         when present on disk.
    """
    u_set: Optional[Set[str]] = None
    if detected_universe is not None:
        u_set = {str(g).strip().upper() for g in detected_universe if str(g).strip()}

    gene_sets: Dict[str, List[str]] = {}
    for k, v in CORE_MACROMOLECULAR_COMPLEXES.items():
        genes = [str(g).strip().upper() for g in v if str(g).strip()]
        if u_set is not None:
            genes = [g for g in genes if g in u_set]
        if len(genes) >= 2:
            gene_sets[f"CORUM:{k}"] = genes

    sigs = _load_curated_signatures()
    for section, prefix in [
        ("pathway_gene_sets", "PATHWAY"),
        ("subcellular_compartments", "COMPARTMENT"),
        ("misev2023_markers", "MISEV2023"),
    ]:
        sec_dict = sigs.get(section, {})
        if isinstance(sec_dict, dict):
            for term, members in sec_dict.items():
                if not isinstance(members, list):
                    continue
                genes = [str(g).strip().upper() for g in members if str(g).strip()]
                if u_set is not None:
                    genes = [g for g in genes if g in u_set]
                if len(genes) >= 2:
                    gene_sets[f"{prefix}:{term}"] = genes

    if include_go_msigdb:
        ext_dir = _locate_external_data_dir(external_data_dir)
        if ext_dir is not None and (ext_dir / "genesets").exists():
            gs_dir = ext_dir / "genesets"
            for gs_file in sorted(gs_dir.glob("gs_*.txt")):
                lib_name = gs_file.stem.replace("gs_", "").upper()
                try:
                    with open(gs_file, "r", encoding="utf-8", errors="ignore") as fh:
                        for line in fh:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) < 2:
                                continue
                            term_name = parts[0].strip()
                            raw_genes = [
                                g.strip().upper()
                                for g in re.split(r"[,;\s\t]+", "\t".join(parts[1:]))
                                if g.strip()
                            ]
                            if u_set is not None:
                                raw_genes = [g for g in raw_genes if g in u_set]
                            if 4 <= len(raw_genes) <= 600:
                                full_term = term_name if ":" in term_name else f"{lib_name}:{term_name}"
                                gene_sets[full_term] = raw_genes
                except Exception:
                    continue

    return gene_sets


def run_detected_background_ora(
    hit_genes: Sequence[str],
    detected_universe: Sequence[str],
    gene_sets: Optional[Dict[str, Sequence[str]]] = None,
    min_overlap: int = 2,
    external_data_dir: Optional[str | Path] = None,
) -> pd.DataFrame:
    """
    Run Fisher's Exact Over-Representation Analysis (ORA) strictly against the
    EMPIRICALLY DETECTED proteome universe `U` (never the whole human genome!).
    """
    u_set: Set[str] = {str(g).strip().upper() for g in detected_universe if str(g).strip()}
    h_set: Set[str] = {str(g).strip().upper() for g in hit_genes if str(g).strip()} & u_set
    n_u = len(u_set)
    n_h = len(h_set)

    if gene_sets is None:
        gene_sets = load_comprehensive_gene_sets(detected_universe=list(u_set), external_data_dir=external_data_dir)
        if not gene_sets:
            gene_sets = CORE_MACROMOLECULAR_COMPLEXES

    rows = []
    for term_name, members in gene_sets.items():
        m_set = {str(g).strip().upper() for g in members} & u_set
        n_m = len(m_set)
        if n_m < min_overlap:
            continue
        overlap_genes = sorted(h_set & m_set)
        k = len(overlap_genes)
        if k < min_overlap:
            continue

        a = k
        b = max(n_h - k, 0)
        c = max(n_m - k, 0)
        d = max(n_u - n_h - n_m + k, 0)
        odds_ratio, p_val = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        expected = (n_h * n_m) / max(n_u, 1)
        fold_enrich = k / max(expected, 1e-8)

        rows.append(
            {
                "term": term_name,
                "n_detected_in_universe": n_m,
                "n_hits_in_term": k,
                "expected_hits": round(expected, 2),
                "fold_enrichment": round(fold_enrich, 3),
                "odds_ratio": round(float(odds_ratio), 3) if np.isfinite(odds_ratio) else 999.0,
                "p_value": float(p_val),
                "overlap_genes": ";".join(overlap_genes),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "term", "n_detected_in_universe", "n_hits_in_term", "expected_hits",
                "fold_enrichment", "odds_ratio", "p_value", "fdr_q_value", "overlap_genes",
            ]
        )

    out = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    out["fdr_q_value"] = benjamini_hochberg(out["p_value"].values)
    return out


def cluster_redundant_pathways_jaccard(
    ora_df: pd.DataFrame,
    genes_col: str = "overlap_genes",
    jaccard_threshold: float = 0.50,
) -> pd.DataFrame:
    """
    De-duplicate redundant enriched pathways using greedy Jaccard similarity clustering
    so that only the most statistically significant parent representative per cluster is promoted.
    """
    if ora_df.empty or genes_col not in ora_df.columns:
        return ora_df
    out = ora_df.sort_values("p_value").reset_index(drop=True).copy()
    gene_sets = [set(str(x).replace(",", ";").split(";")) - {""} for x in out[genes_col]]
    cluster_ids = [-1] * len(out)
    is_rep = [False] * len(out)
    cid = 0
    for i in range(len(out)):
        if cluster_ids[i] != -1:
            continue
        cid += 1
        cluster_ids[i] = cid
        is_rep[i] = True
        gi = gene_sets[i]
        for j in range(i + 1, len(out)):
            if cluster_ids[j] != -1:
                continue
            gj = gene_sets[j]
            union_len = len(gi | gj)
            if union_len > 0 and (len(gi & gj) / union_len) >= jaccard_threshold:
                cluster_ids[j] = cid
    out["jaccard_cluster_id"] = cluster_ids
    out["is_cluster_representative"] = is_rep
    return out


def run_multi_group_comparative_ora(
    df: pd.DataFrame,
    group_col: str = "Divergence_Group",
    gene_col: str = "std_gene",
    gene_sets: Optional[Dict[str, Sequence[str]]] = None,
    min_overlap: int = 2,
    jaccard_threshold: float = 0.35,
    external_data_dir: Optional[str | Path] = None,
) -> pd.DataFrame:
    """
    Run detected-background ORA across all biological sectors/groups in `df[group_col]`
    and apply per-group Jaccard deduplication. Ready for `plot_12_jaccard_comparative_ora_and_regulators`.
    """
    if df.empty or gene_col not in df.columns or group_col not in df.columns:
        return pd.DataFrame()
    universe = df[gene_col].dropna().astype(str).str.strip().str.upper().unique().tolist()
    if gene_sets is None:
        gene_sets = load_comprehensive_gene_sets(detected_universe=universe, external_data_dir=external_data_dir)

    all_frames: List[pd.DataFrame] = []
    for grp_name, grp_df in df.groupby(group_col):
        g_label = str(grp_name)
        if g_label in ("G0_Background", "S0_Background", "Unchanged_Background", "nan", "None"):
            continue
        hits = grp_df[gene_col].dropna().astype(str).str.strip().str.upper().unique().tolist()
        if len(hits) < min_overlap:
            continue
        ora_res = run_detected_background_ora(
            hit_genes=hits,
            detected_universe=universe,
            gene_sets=gene_sets,
            min_overlap=min_overlap,
        )
        if ora_res.empty:
            continue
        ora_res = cluster_redundant_pathways_jaccard(ora_res, jaccard_threshold=jaccard_threshold)
        ora_res.insert(0, "Group", g_label)
        ora_res.insert(1, "group_size", len(hits))
        all_frames.append(ora_res)

    if not all_frames:
        return pd.DataFrame()
    combined = pd.concat(all_frames, ignore_index=True)
    return combined.sort_values(["Group", "p_value"]).reset_index(drop=True)


# =============================================================================
# ENGINE [05]: CORUM Macromolecular Complex Stoichiometry & Subunit Matrix
# =============================================================================

def audit_corum_complex_co_shifts(
    df: pd.DataFrame,
    fc_col: str = "log2FC",
    gene_col: str = "std_gene",
    complex_dict: Optional[Dict[str, Sequence[str]]] = None,
    min_detected_subunits: int = 3,
    n_perm: int = 5000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Evaluate macromolecular complex coordinate shifts, subunit directional coherence,
    and stoichiometric uncoupling (orphan subunit shifts) against size-matched permutations
    drawn from the detected proteome.
    """
    if complex_dict is None:
        complex_dict = CORE_MACROMOLECULAR_COMPLEXES

    if fc_col not in df.columns:
        cand_cols = [c for c in df.columns if "log2fc" in c.lower() or "eff_fc" in c.lower()]
        if cand_cols:
            fc_col = cand_cols[0]
        else:
            return pd.DataFrame()

    sub = df[[gene_col, fc_col]].dropna().copy()
    sub[gene_col] = sub[gene_col].astype(str).str.strip().str.upper()
    sub[fc_col] = pd.to_numeric(sub[fc_col], errors="coerce")
    sub = sub.dropna().groupby(gene_col, as_index=False)[fc_col].mean()

    gene_to_fc = dict(zip(sub[gene_col], sub[fc_col]))
    all_fcs = sub[fc_col].values
    rng = np.random.default_rng(seed)

    rows = []
    for cx_name, members in complex_dict.items():
        det_genes = [g.upper() for g in members if g.upper() in gene_to_fc]
        m = len(det_genes)
        if m < min_detected_subunits or m >= len(all_fcs):
            continue

        vals = np.array([gene_to_fc[g] for g in det_genes], dtype=float)
        mean_shift = float(np.mean(vals))
        med_shift = float(np.median(vals))
        iqr_shift = float(np.percentile(vals, 75) - np.percentile(vals, 25))
        dir_frac = float(max(np.mean(vals > 0), np.mean(vals < 0)))

        perm_means = np.mean(rng.choice(all_fcs, size=(n_perm, m), replace=True), axis=1)
        p_perm = float((np.sum(np.abs(perm_means) >= abs(mean_shift)) + 1.0) / (n_perm + 1.0))

        max_idx = int(np.argmax(np.abs(vals)))
        orphan_gene = det_genes[max_idx]
        rest_vals = np.delete(vals, max_idx)
        is_orphan_uncoupled = bool(abs(vals[max_idx]) >= 1.0 and abs(np.median(rest_vals)) < 0.25)

        if p_perm < 0.05 and dir_frac >= 0.70 and abs(mean_shift) >= 0.25:
            regime = "Coordinated_Complex_CoShift"
        elif is_orphan_uncoupled:
            regime = f"Orphan_Subunit_Uncoupling ({orphan_gene})"
        elif abs(mean_shift) < 0.20 and iqr_shift < 0.35:
            regime = "Stoichiometrically_Preserved_Complex"
        else:
            regime = "Heterogeneous_or_Background"

        rows.append(
            {
                "complex_name": cx_name,
                "n_subunits_detected": m,
                "n_subunits_annotated": len(members),
                "mean_log2FC": round(mean_shift, 4),
                "median_log2FC": round(med_shift, 4),
                "subunit_iqr_log2FC": round(iqr_shift, 4),
                "directional_agreement_frac": round(dir_frac, 3),
                "perm_p_value": round(p_perm, 5),
                "complex_regime": regime,
                "subunit_details": "; ".join(f"{g}:{gene_to_fc[g]:+.2f}" for g in det_genes),
            }
        )

    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values("perm_p_value").reset_index(drop=True)
    out["perm_fdr"] = benjamini_hochberg(out["perm_p_value"].values)
    return out


def extract_complex_subunit_matrix(
    df: pd.DataFrame,
    fc_cols: Optional[Sequence[str]] = None,
    gene_col: str = "std_gene",
    complex_dict: Optional[Dict[str, Sequence[str]]] = None,
    min_detected_subunits: int = 3,
) -> pd.DataFrame:
    """
    Extract a tidy subunit-by-contrast fold-change table for CORUM / subpopulation complexes
    to power `plot_11_corum_complex_and_subpopulation_heatmap()`.
    """
    if complex_dict is None:
        complex_dict = CORE_MACROMOLECULAR_COMPLEXES
    if fc_cols is None:
        fc_cols = [c for c in df.columns if ("log2fc" in c.lower() or "eff_fc" in c.lower() or "eff_di" in c.lower())]
    if not fc_cols or gene_col not in df.columns:
        return pd.DataFrame()

    work = df[[gene_col] + list(fc_cols)].copy()
    work[gene_col] = work[gene_col].astype(str).str.strip().str.upper()
    for c in fc_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce")
    work = work.groupby(gene_col, as_index=False)[list(fc_cols)].mean()
    gene_Lookup = work.set_index(gene_col)

    rows = []
    for cx_name, members in complex_dict.items():
        det = [g.upper() for g in members if g.upper() in gene_Lookup.index]
        if len(det) < min_detected_subunits:
            continue
        for g in det:
            rec: Dict[str, Any] = {"complex_name": cx_name, "std_gene": g}
            for c in fc_cols:
                rec[c] = float(gene_Lookup.loc[g, c])
            rows.append(rec)
    return pd.DataFrame(rows)


# =============================================================================
# ENGINE [06]: STRING Physical/Functional PPI, Degree-Corrected RWR & KDA
# =============================================================================

# Built-in high-confidence physical & functional PPI interactions for zero-network fallback
BUILTIN_CURATED_STRING_EDGES: List[Tuple[str, str, float, float]] = [
    # (gene_a, gene_b, combined_score_0_to_1, physical_exp_score_0_to_1)
    # Integrin-Tetraspanin-ADAM-MARVEL web & ECM corona
    ("ITGB1", "ITGA2", 0.999, 0.980), ("ITGB1", "ITGA3", 0.999, 0.985), ("ITGB1", "ITGA5", 0.999, 0.990),
    ("ITGB1", "ITGA6", 0.999, 0.985), ("ITGB1", "ITGAV", 0.998, 0.960), ("ITGAV", "ITGB3", 0.999, 0.990),
    ("ITGAV", "ITGB5", 0.999, 0.980), ("ITGA3", "CD151", 0.995, 0.950), ("ITGA6", "CD151", 0.996, 0.960),
    ("ITGB1", "CD151", 0.992, 0.920), ("ITGB1", "CD81", 0.985, 0.890), ("ITGB1", "CD9", 0.988, 0.910),
    ("ITGB1", "CD63", 0.975, 0.850), ("CD151", "CD81", 0.970, 0.870), ("CD151", "CD9", 0.965, 0.860),
    ("CD81", "CD9", 0.990, 0.930), ("CD81", "IGSF8", 0.995, 0.960), ("CD9", "IGSF8", 0.992, 0.950),
    ("CD81", "PTGFRN", 0.990, 0.940), ("CD9", "PTGFRN", 0.992, 0.950), ("CD81", "TSPAN4", 0.910, 0.720),
    ("CD151", "TSPAN9", 0.880, 0.680), ("CD81", "ADAM10", 0.960, 0.880), ("CD9", "ADAM10", 0.955, 0.870),
    ("TSPAN14", "ADAM10", 0.992, 0.960), ("TSPAN15", "ADAM10", 0.990, 0.950), ("ITGB1", "ADAM9", 0.965, 0.880),
    ("ITGAV", "ADAM9", 0.950, 0.840), ("ITGB1", "FN1", 0.999, 0.990), ("ITGA5", "FN1", 0.999, 0.990),
    ("ITGAV", "FN1", 0.998, 0.970), ("ITGAV", "THBS1", 0.995, 0.940), ("ITGB1", "THBS1", 0.990, 0.920),
    ("FN1", "THBS1", 0.992, 0.910), ("FN1", "LGALS3BP", 0.965, 0.850), ("ITGB1", "LGALS3BP", 0.955, 0.840),
    ("FN1", "GPC3", 0.890, 0.680), ("THBS1", "GPC3", 0.875, 0.650), ("ITGAV", "MFGE8", 0.995, 0.960),
    ("ITGB3", "MFGE8", 0.996, 0.960), ("ITGB5", "MFGE8", 0.994, 0.950), ("FN1", "SPP1", 0.985, 0.890),
    ("ITGAV", "SPP1", 0.995, 0.950), ("ITGB1", "SPP1", 0.992, 0.930), ("FN1", "TGFBI", 0.960, 0.830),
    ("ITGB1", "LAMB1", 0.990, 0.920), ("ITGB1", "LAMC1", 0.990, 0.920), ("LAMB1", "LAMC1", 0.999, 0.995),
    ("FN1", "COL6A1", 0.975, 0.880), ("COL6A1", "COL6A2", 0.999, 0.990), ("COL6A1", "COL6A3", 0.999, 0.990),
    # Vesicle / SNARE / MARVEL trafficking axis
    ("SYNGR2", "VAMP3", 0.925, 0.810), ("SYNGR2", "VAMP2", 0.910, 0.790), ("SYNGR2", "SNAP23", 0.890, 0.760),
    ("RAB27A", "RAB27B", 0.940, 0.820), ("RAB27A", "VAMP2", 0.935, 0.840), ("RAB27A", "VAMP3", 0.920, 0.810),
    ("RAB27A", "SNAP23", 0.915, 0.800), ("RAB27A", "EXOC4", 0.885, 0.740),
    ("RAB27A", "SYTL2", 0.995, 0.970), ("RAB27A", "STX4", 0.920, 0.810), ("VAMP3", "SNAP23", 0.995, 0.970),
    ("STX4", "SNAP23", 0.998, 0.985), ("VAMP3", "STX4", 0.992, 0.960), ("VAMP7", "SNAP23", 0.980, 0.920),
    ("VAMP7", "CD63", 0.910, 0.790), ("RAB11A", "EXOC6", 0.965, 0.890), ("RAB7A", "CD63", 0.920, 0.810),
    # ESCRT MVB Biogenesis
    ("TSG101", "PDCD6IP", 0.999, 0.995), ("TSG101", "HGS", 0.998, 0.985), ("TSG101", "VPS28", 0.999, 0.995),
    ("PDCD6IP", "CHMP4B", 0.999, 0.995), ("PDCD6IP", "SDCBP", 0.999, 0.995), ("CHMP4B", "VPS4A", 0.998, 0.985),
    ("CHMP4B", "VPS4B", 0.998, 0.985), ("VPS4A", "VPS4B", 0.995, 0.960), ("SDCBP", "CD63", 0.985, 0.940),
    # Redox / Glutathione / Ferroptosis / Iron Storage axis
    ("GCLM", "GSTP1", 0.920, 0.810), ("GSTP1", "PRDX6", 0.915, 0.780),
    ("GCLM", "GCLC", 0.999, 0.995), ("GCLM", "NQO1", 0.955, 0.820),
    ("GCLM", "TXNRD1", 0.960, 0.840), ("GCLM", "HMOX1", 0.950, 0.810), ("NQO1", "TXNRD1", 0.965, 0.850),
    ("HMOX1", "FTH1", 0.975, 0.880), ("FTH1", "FTL", 0.999, 0.998), ("FTH1", "NCOA4", 0.998, 0.990),
    ("FTL", "NCOA4", 0.995, 0.970), ("FTH1", "TFRC", 0.980, 0.910), ("GPX4", "SLC7A11", 0.975, 0.860),
    # Chromatin / Epigenetic / IFN-STAT1 / Replisome axis
    ("HELLS", "DNMT1", 0.985, 0.940), ("HELLS", "DNMT3B", 0.975, 0.910), ("DNMT1", "UHRF1", 0.999, 0.995),
    ("HELLS", "UHRF1", 0.960, 0.870), ("DNMT1", "PCNA", 0.998, 0.990), ("PCNA", "MCM2", 0.990, 0.940),
    ("MCM2", "MCM4", 0.999, 0.998), ("MCM2", "MCM6", 0.999, 0.998), ("MCM4", "MCM6", 0.999, 0.998),
    ("STAT1", "STAT2", 0.999, 0.995), ("STAT1", "ISG15", 0.995, 0.950), ("STAT1", "PARP9", 0.990, 0.940),
    ("PARP9", "DTX3L", 0.999, 0.995), ("STAT1", "IFIT1", 0.985, 0.900), ("STAT1", "MX1", 0.985, 0.900),
]


def fetch_string_ppi_network(
    genes: Sequence[str],
    external_data_dir: Optional[str | Path] = None,
    min_score: int = 450,
    query_online: bool = False,
) -> pd.DataFrame:
    """
    Build a high-confidence protein-protein interaction (PPI) network (`source`, `target`,
    `weight`, `physical_score`, `is_physical`, `edge_source`) for `genes` using a 3-tier resolution engine:
      1. Local `external_data/string/string_map.txt` + `string_edges.txt`
      2. Built-in `BUILTIN_CURATED_STRING_EDGES` + intra-complex CORUM physical edges
      3. Optional live STRING REST API (`https://string-db.org/api/json/network`)
    """
    gene_set: Set[str] = {str(g).strip().upper() for g in genes if str(g).strip()}
    if len(gene_set) < 2:
        return pd.DataFrame(columns=["source", "target", "weight", "physical_score", "is_physical", "edge_source"])

    edge_map: Dict[Tuple[str, str], Tuple[float, float, str]] = {}

    def _add_edge(u: str, v: str, w: float, phys: float, src_tag: str = "curated_ppi") -> None:
        if u == v or u not in gene_set or v not in gene_set:
            return
        pair = (u, v) if u < v else (v, u)
        old = edge_map.get(pair)
        if old is None or w > old[0]:
            edge_map[pair] = (float(w), float(phys), src_tag)

    # Tier 1: Local `external_data/string/` cache (supports both headerless integer-index format and headered TSV)
    ext_dir = _locate_external_data_dir(external_data_dir)
    if ext_dir is not None and (ext_dir / "string" / "string_map.txt").exists() and (ext_dir / "string" / "string_edges.txt").exists():
        try:
            smap = pd.read_csv(
                ext_dir / "string" / "string_map.txt",
                sep="\t", header=None, names=["gene", "sid", "pref"]
            )
            sids = list(dict.fromkeys(smap["sid"]))
            first_gene = smap.drop_duplicates("sid").set_index("sid")["gene"]
            names = [str(first_gene[s]).strip().upper() for s in sids]
            in_set = np.array([g in gene_set for g in names], dtype=bool)
            if np.any(in_set):
                e_df = pd.read_csv(
                    ext_dir / "string" / "string_edges.txt",
                    sep=r"\s+", header=None, names=["i", "j", "score", "n", "f", "p", "a", "e", "d", "t"]
                )
                i_arr = pd.to_numeric(e_df["i"], errors="coerce").fillna(-1).astype(int).values
                j_arr = pd.to_numeric(e_df["j"], errors="coerce").fillna(-1).astype(int).values
                s_arr = pd.to_numeric(e_df["score"], errors="coerce").fillna(0.0).values
                exp_arr = pd.to_numeric(e_df["e"], errors="coerce").fillna(0.0).values
                valid_mask = (
                    (i_arr >= 0) & (i_arr < len(names))
                    & (j_arr >= 0) & (j_arr < len(names))
                    & (s_arr >= min_score)
                )
                for idx_row in np.where(valid_mask)[0]:
                    ia, jb = int(i_arr[idx_row]), int(j_arr[idx_row])
                    if in_set[ia] and in_set[jb]:
                        _add_edge(names[ia], names[jb], float(s_arr[idx_row]) / 1000.0, float(exp_arr[idx_row]) / 1000.0, "local_string_cache")
        except Exception:
            pass

    # Tier 2: Built-in curated STRING physical & functional interactions + CORUM co-complex edges
    thresh = min_score / 1000.0
    for u, v, w, phys in BUILTIN_CURATED_STRING_EDGES:
        if w >= thresh:
            _add_edge(u.upper(), v.upper(), w, phys, "curated_ppi")

    for cx_name, members in CORE_MACROMOLECULAR_COMPLEXES.items():
        present = [g.upper() for g in members if g.upper() in gene_set]
        if 2 <= len(present) <= 25:
            for i in range(len(present)):
                for j in range(i + 1, len(present)):
                    _add_edge(present[i], present[j], 0.85, 0.80, "curated_complex")

    # Tier 3: Optional live STRING API query if few edges found and query_online=True
    if query_online and len(edge_map) < 5 and len(gene_set) <= 300:
        try:
            identifiers = "%0d".join(sorted(list(gene_set))[:200])
            url = (
                f"https://string-db.org/api/json/network?identifiers={identifiers}"
                f"&species=9606&required_score={int(min_score)}"
            )
            req = urllib.request.Request(url, headers={"User-Agent": "ProteoSkill/4.1"})
            with urllib.request.urlopen(req, timeout=8) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                for item in data:
                    u = str(item.get("preferredName_A", "")).upper()
                    v = str(item.get("preferredName_B", "")).upper()
                    w = float(item.get("score", 0.0))
                    phys = float(item.get("escore", 0.0))
                    _add_edge(u, v, w, phys, "string_api")
        except Exception:
            pass

    if not edge_map:
        return pd.DataFrame(columns=["source", "target", "weight", "physical_score", "is_physical", "edge_source"])

    rows = [
        {
            "source": pair[0],
            "target": pair[1],
            "weight": round(vals[0], 4),
            "physical_score": round(vals[1], 4),
            "is_physical": bool(vals[1] >= 0.30),
            "edge_source": vals[2],
        }
        for pair, vals in sorted(edge_map.items())
    ]
    return pd.DataFrame(rows)


def degree_corrected_rwr(
    edges_df: pd.DataFrame,
    seed_genes: Sequence[str],
    restart_prob: float = 0.30,
    n_degree_perms: int = 500,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Run Random Walk with Restart (RWR) on a PPI network (`source`, `target`, `weight`)
    and calibrate every node's stationary probability against a DEGREE-PRESERVING
    random seed null distribution to eliminate high-degree hub bias.
    """
    if edges_df.empty:
        return pd.DataFrame(columns=["gene", "degree", "rwr_score", "degree_null_z", "degree_null_p", "is_true_bridge"])

    u_nodes = sorted(
        set(edges_df["source"].astype(str).str.upper()).union(set(edges_df["target"].astype(str).str.upper()))
    )
    node_idx = {g: i for i, g in enumerate(u_nodes)}
    n = len(u_nodes)

    adj = np.zeros((n, n), dtype=float)
    for _, row in edges_df.iterrows():
        u = str(row["source"]).upper()
        v = str(row["target"]).upper()
        w = float(row.get("weight", 1.0))
        if u in node_idx and v in node_idx and u != v:
            i, j = node_idx[u], node_idx[v]
            adj[i, j] = max(adj[i, j], w)
            adj[j, i] = max(adj[j, i], w)

    deg = np.sum(adj > 0, axis=1)
    col_sums = np.sum(adj, axis=0)
    col_sums_safe = np.where(col_sums > 0, col_sums, 1.0)
    w_mat = adj / col_sums_safe[None, :]

    inv_op = restart_prob * np.linalg.inv(np.eye(n) - (1.0 - restart_prob) * w_mat)

    valid_seeds = [s.upper() for s in seed_genes if s.upper() in node_idx]
    if not valid_seeds:
        return pd.DataFrame(columns=["gene", "degree", "rwr_score", "degree_null_z", "degree_null_p", "is_true_bridge"])

    p0 = np.zeros(n, dtype=float)
    for s in valid_seeds:
        p0[node_idx[s]] = 1.0 / len(valid_seeds)
    rwr_obs = inv_op @ p0

    rng = np.random.default_rng(random_seed)
    deg_bins = pd.qcut(pd.Series(deg).rank(method="first"), q=min(5, max(1, n // 5)), labels=False).values
    seed_indices = [node_idx[s] for s in valid_seeds]
    seed_bins = [deg_bins[idx] for idx in seed_indices]

    bin_to_candidates = {b: np.where(deg_bins == b)[0] for b in set(deg_bins)}
    perm_scores = np.zeros((n_degree_perms, n), dtype=float)
    for b_idx in range(n_degree_perms):
        p0_perm = np.zeros(n, dtype=float)
        for sb in seed_bins:
            chosen = rng.choice(bin_to_candidates[sb])
            p0_perm[chosen] += 1.0 / len(seed_bins)
        perm_scores[b_idx, :] = inv_op @ p0_perm

    null_mean = np.mean(perm_scores, axis=0)
    null_sd = np.maximum(np.std(perm_scores, axis=0, ddof=1), 1e-8)
    z_deg = (rwr_obs - null_mean) / null_sd
    p_deg = (np.sum(perm_scores >= rwr_obs[None, :], axis=0) + 1.0) / (n_degree_perms + 1.0)

    out = pd.DataFrame(
        {
            "gene": u_nodes,
            "is_seed": [g in set(valid_seeds) for g in u_nodes],
            "degree": deg.astype(int),
            "rwr_score": np.round(rwr_obs, 6),
            "degree_null_z": np.round(z_deg, 3),
            "degree_null_p": np.round(p_deg, 5),
        }
    )
    out["is_true_bridge"] = (~out["is_seed"]) & (out["degree_null_p"] < 0.05) & (out["degree_null_z"] >= 2.0)
    return out.sort_values("degree_null_z", ascending=False).reset_index(drop=True)


def run_key_driver_analysis_kda(
    edges_df: pd.DataFrame,
    seed_genes: Sequence[str],
    node_group_map: Optional[Dict[str, str]] = None,
    restart_prob: float = 0.30,
    n_degree_perms: int = 500,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Combine Degree-Corrected RWR with 1-Hop Hypergeometric Key Driver Analysis (KDA)
    and annotate Cross-Sector "Lock-and-Key" physical/functional PPI bridges.
    Returns `(kda_nodes_df, annotated_edges_df)`.
    """
    if edges_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    rwr_df = degree_corrected_rwr(
        edges_df=edges_df,
        seed_genes=seed_genes,
        restart_prob=restart_prob,
        n_degree_perms=n_degree_perms,
    )
    if rwr_df.empty:
        return pd.DataFrame(), edges_df.copy()

    # Build 1-hop adjacency list
    neighbors: Dict[str, Set[str]] = {g: set() for g in rwr_df["gene"]}
    phys_neighbors: Dict[str, Set[str]] = {g: set() for g in rwr_df["gene"]}
    for _, r in edges_df.iterrows():
        u = str(r["source"]).upper()
        v = str(r["target"]).upper()
        is_p = bool(r.get("is_physical", float(r.get("physical_score", 0.0)) >= 0.30))
        if u in neighbors and v in neighbors and u != v:
            neighbors[u].add(v)
            neighbors[v].add(u)
            if is_p:
                phys_neighbors[u].add(v)
                phys_neighbors[v].add(u)

    seed_set = {str(s).strip().upper() for s in seed_genes} & set(neighbors.keys())
    n_total = len(neighbors)
    n_seeds = len(seed_set)

    group_lookup: Dict[str, str] = {}
    if node_group_map:
        group_lookup = {str(k).strip().upper(): str(v) for k, v in node_group_map.items()}

    kda_hits_list = []
    kda_fold_list = []
    kda_p_list = []
    phys_deg_list = []
    cross_sector_deg_list = []
    node_group_list = []

    for g in rwr_df["gene"]:
        nb = neighbors.get(g, set())
        deg_i = len(nb)
        hits_i = len(nb & seed_set)
        expected_i = (deg_i * n_seeds) / max(n_total, 1)
        fold_i = hits_i / max(expected_i, 1e-8)
        if deg_i >= 2 and hits_i >= 1:
            p_kda = float(stats.hypergeom.sf(hits_i - 1, n_total, n_seeds, deg_i))
        else:
            p_kda = 1.0

        g_grp = group_lookup.get(g, "Seed" if g in seed_set else "Bridge")
        cross_cnt = sum(
            1 for n_gene in nb
            if group_lookup.get(n_gene, "Background") not in (g_grp, "Background", "G0_Background", "S0_Background")
        )
        kda_hits_list.append(hits_i)
        kda_fold_list.append(round(fold_i, 3))
        kda_p_list.append(p_kda)
        phys_deg_list.append(len(phys_neighbors.get(g, set())))
        cross_sector_deg_list.append(cross_cnt)
        node_group_list.append(g_grp)

    rwr_df["node_group"] = node_group_list
    rwr_df["physical_degree"] = phys_deg_list
    rwr_df["kda_seed_hits"] = kda_hits_list
    rwr_df["kda_fold_enrichment"] = kda_fold_list
    rwr_df["kda_p_value"] = np.round(kda_p_list, 6)
    rwr_df["kda_fdr"] = np.round(benjamini_hochberg(np.array(kda_p_list)), 6)
    rwr_df["cross_sector_edges"] = cross_sector_deg_list
    rwr_df["is_kda_driver"] = (
        ((rwr_df["kda_p_value"] < 0.05) & (rwr_df["degree"] >= 3))
        | rwr_df["is_true_bridge"]
        | (rwr_df["cross_sector_edges"] >= 2)
    )

    # Annotate edges with intra-sector vs cross-sector Lock-and-Key status
    ann_edges = edges_df.copy()
    ann_edges["source_group"] = ann_edges["source"].astype(str).str.upper().map(lambda x: group_lookup.get(x, "Background"))
    ann_edges["target_group"] = ann_edges["target"].astype(str).str.upper().map(lambda x: group_lookup.get(x, "Background"))
    ann_edges["edge_type"] = np.where(
        (ann_edges["source_group"] != ann_edges["target_group"])
        & (~ann_edges["source_group"].isin(["Background", "G0_Background", "S0_Background"]))
        & (~ann_edges["target_group"].isin(["Background", "G0_Background", "S0_Background"])),
        "Cross_Sector_Lock_and_Key",
        "Intra_Sector_CoModule",
    )

    rwr_df = rwr_df.sort_values(["is_kda_driver", "rwr_score", "degree"], ascending=[False, False, False]).reset_index(drop=True)
    return rwr_df, ann_edges


# =============================================================================
# ENGINE [07]: Direct Protein-Table Biophysical & Membrane Topology Enrichment
# =============================================================================

# Curated fallback structural annotations for offline/standalone runs
_CURATED_SIGNAL_PEPTIDE_SECRETED: Set[str] = {
    "FN1", "THBS1", "LGALS3BP", "GPC3", "GPC1", "MFGE8", "SPP1", "TGFBI", "POSTN", "CCN1", "CCN2",
    "LAMA1", "LAMA5", "LAMB1", "LAMC1", "COL1A1", "COL1A2", "COL4A1", "COL4A2", "COL6A1", "COL6A2", "COL6A3", "COL18A1",
    "HSPG2", "NID1", "NID2", "AGRIN", "VTN", "TNC", "FBN1", "EMILIN1", "TIMP1", "TIMP2", "MMP2", "MMP9",
    "SERPINA1", "SERPINA3", "SERPINC1", "SERPINE1", "SERPING1", "A2M", "ALB", "TF", "APOA1", "APOA2", "APOB", "APOE", "APOH",
    "FGA", "FGB", "FGG", "F2", "PLG", "C3", "C4A", "C4B", "CFB", "CFH", "CLU", "IGFBP2", "IGFBP3", "GDF15",
    "CXCL8", "CCL2", "IL6", "VEGFA", "HDGF", "MDK", "PTN", "GALNT2", "PCSK9", "AFP", "GPC4", "SPARC", "LUM", "DCN", "BGN",
}

_CURATED_ECM_MATRISOME: Set[str] = {
    "FN1", "THBS1", "LGALS3BP", "GPC3", "GPC1", "MFGE8", "SPP1", "TGFBI", "POSTN", "CCN1", "CCN2",
    "LAMA1", "LAMA5", "LAMB1", "LAMC1", "COL1A1", "COL1A2", "COL4A1", "COL4A2", "COL6A1", "COL6A2", "COL6A3", "COL18A1",
    "HSPG2", "NID1", "NID2", "AGRIN", "VTN", "TNC", "FBN1", "EMILIN1", "TIMP1", "TIMP2", "MMP2", "MMP9", "SPARC", "LUM", "DCN", "BGN",
}

_CURATED_SINGLE_PASS_TM: Set[str] = {
    "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB2", "ITGB6", "ITGB8",
    "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGAV",
    "ADAM10", "ADAM9", "ADAM17", "ADAM15", "IGSF8", "PTGFRN", "BSG", "CD44", "ICAM1", "VCAM1", "EPCAM", "CD47",
    "EGFR", "ERBB2", "MET", "FGFR1", "FGFR2", "IGF1R", "INSR", "AXL", "LDLR", "LRP1", "TFRC", "LAMP1", "LAMP2",
    "STX1A", "STX3", "STX4", "STX7", "VAMP2", "VAMP3", "VAMP7", "VAMP8", "HLA-A", "HLA-B", "HLA-C", "B2M",
    "SLC3A2", "CD274", "PDCD1LG2", "NOTCH1", "NOTCH2", "APP", "SDC1", "SDC4", "CADM1", "CDH1", "CDH2",
}

_CURATED_4TM_MARVEL_TSPAN: Set[str] = {
    "SYNGR1", "SYNGR2", "SYNGR3", "SYNGR4", "SYP", "MYADM", "MAL", "MAL2", "PLLP", "MARVELD2", "MARVELD3", "OCLN",
    "CMTM3", "CMTM4", "CMTM6", "CMTM7", "CMTM8",
    "CD9", "CD63", "CD81", "CD82", "CD151", "CD37", "CD53",
    "TSPAN1", "TSPAN2", "TSPAN3", "TSPAN4", "TSPAN5", "TSPAN6", "TSPAN7", "TSPAN8", "TSPAN9", "TSPAN12", "TSPAN13", "TSPAN14", "TSPAN15", "TSPAN18", "TSPAN31", "TSPAN33",
}

_CURATED_GPI_ANCHORED: Set[str] = {
    "GPC3", "GPC1", "GPC4", "GPC6", "CD55", "CD59", "NT5E", "CD109", "ALPL", "ALPP", "CEACAM5", "CEACAM6",
    "PLAUR", "THY1", "LY6E", "LY6K", "PRNP", "BST2", "FOLR1", "EFNA1", "EFNA5", "RECK", "GAS1", "MSLN", "ULBP2",
}

_CURATED_LIPIDATED_RAFT: Set[str] = {
    "RAB3A", "RAB27A", "RAB27B", "RAB11A", "RAB11B", "RAB7A", "RAB5A", "RAB8A", "RAB10", "RAB1A", "RAB1B", "RAB2A", "RAB6A", "RAB14", "RAB18", "RAB35",
    "RHOA", "RHOB", "RHOC", "RAC1", "CDC42", "KRAS", "NRAS", "HRAS", "RAP1A", "RAP1B", "ARF1", "ARF6",
    "SRC", "FYN", "YES1", "LYN", "LCK", "GNAI1", "GNAI2", "GNAI3", "GNAQ", "GNA11", "GNAS", "GNA12", "GNA13",
    "CAV1", "CAV2", "FLOT1", "FLOT2", "STOM", "ERLIN1", "ERLIN2", "SNAP23", "SNAP25", "CD9", "CD63", "CD81", "CD82", "CD151", "SYNGR2", "MYADM",
}


def annotate_uniprot_membrane_topology(
    df: pd.DataFrame,
    gene_col: str = "std_gene",
    external_data_dir: Optional[str | Path] = None,
) -> pd.DataFrame:
    """
    Annotate every protein in `df` with 10 structural & membrane topology features:
      - `has_signal_peptide` (UniProt FT SIGNAL)
      - `is_secreted` (UniProt KW-0964 / HPA Secretome)
      - `is_ecm` (UniProt KW-0272 Extracellular Matrix / Matrisome)
      - `tm_count` (Number of UniProt TRANSMEM helices)
      - `is_single_tm` (`tm_count == 1`)
      - `is_multi_tm` (`tm_count >= 2`)
      - `is_4tm_marvel_tspan` (`tm_count == 4` or MARVEL/Tetraspanin family)
      - `is_gpi_anchored` (UniProt KW-0336 GPI-anchor)
      - `is_lipidated_raft` (UniProt KW-0449 Lipoprotein / Palmitoylated / Prenylated / Raft)
      - `is_classical_cargo_sp_or_1tm` (`has_signal_peptide | is_single_tm | is_secreted | is_ecm`)
      - `is_raft_multipass_lipidated_gpi` (`is_multi_tm | is_4tm_marvel_tspan | is_gpi_anchored | is_lipidated_raft`)
    Automatically uses `external_data/annotation/uniprot_human_reviewed.tsv` & `hpa_annotation.tsv`
    when available, with deterministic fallback to curated structural dictionaries & family rules.
    """
    out = df.copy()
    if gene_col not in out.columns:
        return out

    genes_upper = out[gene_col].astype(str).str.strip().str.upper()

    # Load UniProt & HPA caches if present (supports headerless TSV format: acc, gene, kw, signal, tm)
    up_lookup: Dict[str, Dict[str, Any]] = {}
    hpa_sec_genes: Set[str] = set()
    ext_dir = _locate_external_data_dir(external_data_dir)
    if ext_dir is not None:
        up_file = ext_dir / "annotation" / "uniprot_human_reviewed.tsv"
        if up_file.exists():
            try:
                up_df = pd.read_csv(up_file, sep="\t", header=None, names=["acc", "gene", "kw", "signal", "tm"])
                for _, r in up_df.iterrows():
                    raw_g = str(r.get("gene", "")).strip()
                    if not raw_g or raw_g == "nan" or raw_g.lower() == "gene":
                        continue
                    primary_g = raw_g.split()[0].upper()
                    kw_str = str(r.get("kw", ""))
                    sig_val = r.get("signal", 0)
                    sig_str = str(sig_val).strip()
                    has_sp_flag = (sig_str in ("1", "1.0", "True", "true")) or ("SIGNAL" in sig_str) or ("KW-0732" in kw_str)
                    tm_val = r.get("tm", 0)
                    tm_str = str(tm_val).strip()
                    if "TRANSMEM" in tm_str:
                        tm_cnt = tm_str.count("TRANSMEM")
                    else:
                        try:
                            tm_cnt = int(float(tm_str))
                        except Exception:
                            tm_cnt = 0
                    up_lookup[primary_g] = {
                        "has_sp": has_sp_flag,
                        "is_sec": "KW-0964" in kw_str,
                        "is_ecm": ("KW-0272" in kw_str) or ("KW-0084" in kw_str),
                        "is_gpi": "KW-0336" in kw_str,
                        "is_lipid": "KW-0449" in kw_str,
                        "tm_cnt": tm_cnt,
                    }
            except Exception:
                pass

        hpa_file = ext_dir / "annotation" / "hpa_annotation.tsv"
        if hpa_file.exists():
            try:
                hpa_df = pd.read_csv(hpa_file, sep="\t", header=None)
                for _, r in hpa_df.iterrows():
                    g = str(r.iloc[0]).strip().upper()
                    sec_flag = str(r.iloc[1]).strip() if len(r) > 1 else "0"
                    if g and g != "GENE" and (sec_flag in ("1", "1.0", "True") or "Secreted" in sec_flag):
                        hpa_sec_genes.add(g)
            except Exception:
                pass

    has_sp_arr = []
    is_sec_arr = []
    is_ecm_arr = []
    tm_cnt_arr = []
    is_1tm_arr = []
    is_mtm_arr = []
    is_4tm_arr = []
    is_gpi_arr = []
    is_lip_arr = []
    topo_class_arr = []

    for g in genes_upper:
        rec = up_lookup.get(g)
        if rec is not None:
            has_sp = bool(rec["has_sp"] or g in _CURATED_SIGNAL_PEPTIDE_SECRETED)
            is_sec = bool(rec["is_sec"] or g in hpa_sec_genes or g in _CURATED_SIGNAL_PEPTIDE_SECRETED)
            is_ecm = bool(rec["is_ecm"] or g in _CURATED_ECM_MATRISOME)
            tm_cnt = int(rec["tm_cnt"])
            if tm_cnt == 0 and g in _CURATED_4TM_MARVEL_TSPAN:
                tm_cnt = 4
            elif tm_cnt == 0 and g in _CURATED_SINGLE_PASS_TM:
                tm_cnt = 1
            is_1tm = bool(tm_cnt == 1)
            is_mtm = bool(tm_cnt >= 2)
            is_4tm = bool(tm_cnt == 4 or g in _CURATED_4TM_MARVEL_TSPAN)
            is_gpi = bool(rec["is_gpi"] or g in _CURATED_GPI_ANCHORED)
            is_lip = bool(rec["is_lipid"] or g in _CURATED_LIPIDATED_RAFT)
        else:
            # Deterministic curated & gene-family rule inference
            has_sp = bool(
                g in _CURATED_SIGNAL_PEPTIDE_SECRETED
                or (g.startswith(("COL", "LAMA", "LAMB", "LAMC", "SERPINA", "SERPINC", "SERPINE", "SERPING", "APOA", "APOB", "APOC", "APOE", "APOH", "IGH", "IGK", "IGL", "MMP", "TIMP", "ADAM", "ITGA", "ITGB", "GPC")))
            )
            is_ecm = bool(
                g in _CURATED_ECM_MATRISOME
                or g.startswith(("COL", "LAMA", "LAMB", "LAMC", "NID", "MMP", "TIMP", "CCN"))
            )
            is_sec = bool(has_sp or is_ecm or g in hpa_sec_genes)
            is_4tm = bool(
                g in _CURATED_4TM_MARVEL_TSPAN
                or g.startswith(("TSPAN", "SYNGR", "MARVELD", "CMTM"))
            )
            is_multipass_family = bool(
                is_4tm or (g.startswith(("SLC", "ABCA", "ABCB", "ABCC", "ABCD", "ABCG", "ATP1", "ATP2", "CLDN", "GJB", "FZD", "GPR", "TMEM")))
            )
            is_1tm = bool(
                not is_multipass_family
                and (
                    g in _CURATED_SINGLE_PASS_TM
                    or g.startswith(("ITGA", "ITGB", "ADAM", "HLA-", "CDH", "SDC", "LAMP", "ICAM", "VCAM"))
                )
            )
            tm_cnt = 4 if is_4tm else (6 if is_multipass_family else (1 if is_1tm else 0))
            is_mtm = bool(tm_cnt >= 2)
            is_gpi = bool(g in _CURATED_GPI_ANCHORED or g.startswith("GPC"))
            is_lip = bool(
                g in _CURATED_LIPIDATED_RAFT
                or g.startswith(("RAB", "RHOA", "RHOB", "RHOC", "RAC1", "RAC2", "GNA", "FLOT", "CAV"))
            )

        if is_4tm:
            t_cls = "4-TM MARVEL / Tetraspanin"
        elif is_gpi:
            t_cls = "GPI-Anchored Surface"
        elif is_mtm:
            t_cls = f"Multi-Pass TM ({tm_cnt}-TM)"
        elif is_1tm:
            t_cls = "Single-Pass 1-TM Receptor"
        elif is_ecm:
            t_cls = "Secreted ECM / Matrisome"
        elif is_sec or has_sp:
            t_cls = "Secreted Signal-Peptide"
        elif is_lip:
            t_cls = "Lipidated Membrane Raft"
        else:
            t_cls = "Intracellular / Cytosolic"

        has_sp_arr.append(has_sp)
        is_sec_arr.append(is_sec)
        is_ecm_arr.append(is_ecm)
        tm_cnt_arr.append(tm_cnt)
        is_1tm_arr.append(is_1tm)
        is_mtm_arr.append(is_mtm)
        is_4tm_arr.append(is_4tm)
        is_gpi_arr.append(is_gpi)
        is_lip_arr.append(is_lip)
        topo_class_arr.append(t_cls)

    out["has_signal_peptide"] = has_sp_arr
    out["is_secreted"] = is_sec_arr
    out["is_ecm"] = is_ecm_arr
    out["tm_count"] = tm_cnt_arr
    out["is_single_tm"] = is_1tm_arr
    out["is_multi_tm"] = is_mtm_arr
    out["is_4tm_marvel_tspan"] = is_4tm_arr
    out["is_gpi_anchored"] = is_gpi_arr
    out["is_lipidated_raft"] = is_lip_arr
    out["topology_class"] = topo_class_arr
    out["is_classical_cargo_sp_or_1tm"] = (
        out["has_signal_peptide"] | out["is_single_tm"] | out["is_secreted"] | out["is_ecm"]
    )
    out["is_raft_multipass_lipidated_gpi"] = (
        out["is_multi_tm"] | out["is_4tm_marvel_tspan"] | out["is_gpi_anchored"] | out["is_lipidated_raft"]
    )
    return out


def audit_membrane_and_secretory_topology(
    df: pd.DataFrame,
    group_col: str = "Divergence_Group",
    gene_col: str = "std_gene",
    fc_col: str = "log2FC",
    external_data_dir: Optional[str | Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Engine [07]: Directly annotate `df` with UniProt/HPA membrane & secretory structural
    topology features and compute 2-sided Fisher's Exact Tests (`Odds_Ratio`, `OR_Haldane_Anscombe`,
    `log2_OR`, `CI95_low`, `CI95_high`, `Fold_Enrichment`, `p_value`, `FDR_BH`) across all groups/sectors
    plus head-to-head antagonistic contrasts.
    Returns `(annotated_df, fisher_enrichment_df)`.
    """
    ann_df = annotate_uniprot_membrane_topology(df, gene_col=gene_col, external_data_dir=external_data_dir)
    if ann_df.empty:
        return ann_df, pd.DataFrame()

    # Ensure a valid grouping column exists; if not, derive from fold change
    work_group_col = group_col
    if work_group_col not in ann_df.columns:
        cand_fc = fc_col if fc_col in ann_df.columns else next(
            (c for c in ann_df.columns if "log2fc" in c.lower() or "eff_fc" in c.lower()), None
        )
        if cand_fc is not None:
            fcs = pd.to_numeric(ann_df[cand_fc], errors="coerce").fillna(0.0)
            ann_df["_auto_topology_group"] = np.where(
                fcs >= 0.585, "Up_Regulated", np.where(fcs <= -0.585, "Down_Regulated", "G0_Background")
            )
            work_group_col = "_auto_topology_group"
        else:
            return ann_df, pd.DataFrame()

    features = [
        ("has_signal_peptide", "Signal Peptide (ER/Golgi Secretory)"),
        ("is_secreted", "Secreted / Extracellular (UniProt/HPA)"),
        ("is_ecm", "Extracellular Matrix (ECM / Matrisome)"),
        ("is_single_tm", "Single-Pass Transmembrane (1-TM Receptor/Integrin)"),
        ("is_multi_tm", "Multi-Pass Transmembrane (>=2 TM Helices)"),
        ("is_4tm_marvel_tspan", "4-TM MARVEL / Tetraspanin Web (CD63/CD81/CD151)"),
        ("is_gpi_anchored", "GPI-Anchored Surface Protein (KW-0336)"),
        ("is_lipidated_raft", "Lipidated / Palmitoylated / Prenylated Raft (KW-0449)"),
        ("is_classical_cargo_sp_or_1tm", "Classical Secretory / 1-TM / ECM Cargo"),
        ("is_raft_multipass_lipidated_gpi", "Non-Classical Multi-TM / Raft / Lipidated / GPI Cargo"),
    ]

    dedup = ann_df.drop_duplicates(subset=[gene_col]).copy()
    n_total = len(dedup)
    rows: List[Dict[str, Any]] = []

    def _run_fisher_2x2(a: int, b: int, c: int, d: int) -> Tuple[float, float, float, float, float, float]:
        odds_r, p_v = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        # Haldane-Anscombe correction for finite log2(OR) and 95% Wald CI
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        or_ha = float((aa * dd) / (bb * cc))
        log2_or = float(np.log2(or_ha))
        se_ln = float(np.sqrt(1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd))
        ci_low = float(np.log2(np.exp(np.log(or_ha) - 1.96 * se_ln)))
        ci_high = float(np.log2(np.exp(np.log(or_ha) + 1.96 * se_ln)))
        return (float(odds_r) if np.isfinite(odds_r) else or_ha, or_ha, log2_or, ci_low, ci_high, float(p_v))

    for grp_name, grp_df in dedup.groupby(work_group_col):
        g_str = str(grp_name)
        if g_str in ("G0_Background", "S0_Background", "Unchanged_Background", "nan", "None"):
            continue
        n_grp = len(grp_df)
        if n_grp < 3:
            continue
        rest_df = dedup[dedup[work_group_col] != grp_name]
        for feat_col, feat_label in features:
            a = int(grp_df[feat_col].sum())
            b = int(n_grp - a)
            c = int(rest_df[feat_col].sum())
            d = int(len(rest_df) - c)
            tot_feat = a + c
            if tot_feat < 2:
                continue
            or_val, or_ha, log2_or, ci_l, ci_h, p_val = _run_fisher_2x2(a, b, c, d)
            bg_pct = 100.0 * tot_feat / max(n_total, 1)
            grp_pct = 100.0 * a / max(n_grp, 1)
            fold_enr = grp_pct / max(bg_pct, 1e-6)
            hit_genes = sorted(grp_df[grp_df[feat_col]][gene_col].astype(str).tolist())
            rows.append(
                {
                    "Comparison_Type": "Group_vs_Background",
                    "Group": g_str,
                    "Feature_Col": feat_col,
                    "Feature": feat_label,
                    "Group_Hits": a,
                    "Group_Total": n_grp,
                    "Group_Pct": round(grp_pct, 2),
                    "Bg_Hits": tot_feat,
                    "Bg_Total": n_total,
                    "Bg_Pct": round(bg_pct, 2),
                    "Fold_Enrichment": round(fold_enr, 3),
                    "Odds_Ratio": round(or_val, 3),
                    "OR_Haldane_Anscombe": round(or_ha, 3),
                    "log2_OR": round(log2_or, 3),
                    "CI95_log2_low": round(ci_l, 3),
                    "CI95_log2_high": round(ci_h, 3),
                    "p_value": p_val,
                    "Hit_Genes": ";".join(hit_genes[:25]),
                }
            )

    # Head-to-head antagonistic contrast if G7 sub-directions or Up vs Down exist
    head_pairs = [
        ("G7_A_up_B_down", "G7_A_down_B_up"),
        ("S7_Seesaw_A_Up_B_Down", "S7_Seesaw_A_Down_B_Up"),
        ("G5_A_specific", "G6_B_specific"),
        ("Up_Regulated", "Down_Regulated"),
    ]
    present_groups = set(dedup[work_group_col].astype(str).unique())
    for g_pos, g_neg in head_pairs:
        if g_pos in present_groups and g_neg in present_groups:
            df_pos = dedup[dedup[work_group_col].astype(str) == g_pos]
            df_neg = dedup[dedup[work_group_col].astype(str) == g_neg]
            if len(df_pos) >= 3 and len(df_neg) >= 3:
                for feat_col, feat_label in features:
                    a = int(df_pos[feat_col].sum())
                    b = int(len(df_pos) - a)
                    c = int(df_neg[feat_col].sum())
                    d = int(len(df_neg) - c)
                    if (a + c) < 2:
                        continue
                    or_val, or_ha, log2_or, ci_l, ci_h, p_val = _run_fisher_2x2(a, b, c, d)
                    pos_pct = 100.0 * a / max(len(df_pos), 1)
                    neg_pct = 100.0 * c / max(len(df_neg), 1)
                    rows.append(
                        {
                            "Comparison_Type": "HeadToHead",
                            "Group": f"{g_pos}_vs_{g_neg}",
                            "Feature_Col": feat_col,
                            "Feature": feat_label,
                            "Group_Hits": a,
                            "Group_Total": len(df_pos),
                            "Group_Pct": round(pos_pct, 2),
                            "Bg_Hits": c,
                            "Bg_Total": len(df_neg),
                            "Bg_Pct": round(neg_pct, 2),
                            "Fold_Enrichment": round((pos_pct + 0.5) / (neg_pct + 0.5), 3),
                            "Odds_Ratio": round(or_val, 3),
                            "OR_Haldane_Anscombe": round(or_ha, 3),
                            "log2_OR": round(log2_or, 3),
                            "CI95_log2_low": round(ci_l, 3),
                            "CI95_log2_high": round(ci_h, 3),
                            "p_value": p_val,
                            "Hit_Genes": ";".join(sorted(df_pos[df_pos[feat_col]][gene_col].astype(str).tolist())[:25]),
                        }
                    )

    if not rows:
        return ann_df, pd.DataFrame()
    fisher_df = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    fisher_df["FDR_BH"] = benjamini_hochberg(fisher_df["p_value"].values)
    return ann_df, fisher_df


def audit_ectodomain_shedding_topology(
    peptide_df: pd.DataFrame,
    gene_col: str = "std_gene",
    pep_start_col: str = "start_pos",
    pep_end_col: str = "end_pos",
    tm_start_col: str = "tm_helix_start",
    tm_end_col: str = "tm_helix_end",
    ectodomain_side: str = "N_terminal",
) -> pd.DataFrame:
    """
    Map detected peptides against UniProt transmembrane helix coordinates (when peptide-level
    coordinates are available) OR fall back to protein-level membrane shedding annotation
    when given a standard protein table without peptide coordinates.
    """
    if peptide_df.empty:
        return pd.DataFrame()
    if tm_start_col not in peptide_df.columns or pep_start_col not in peptide_df.columns:
        # Protein-table mode: run protein-level membrane & secretory topology audit
        _, fisher_df = audit_membrane_and_secretory_topology(peptide_df, gene_col=gene_col)
        return fisher_df

    rows = []
    for gene, grp in peptide_df.groupby(gene_col):
        s_clean = pd.to_numeric(grp[tm_start_col], errors="coerce").dropna() if tm_start_col in grp else pd.Series(dtype=float)
        e_clean = pd.to_numeric(grp[tm_end_col], errors="coerce").dropna() if tm_end_col in grp else pd.Series(dtype=float)
        tm_s = float(s_clean.iloc[0]) if not s_clean.empty else np.nan
        tm_e = float(e_clean.iloc[0]) if not e_clean.empty else np.nan
        if not np.isfinite(tm_s) or not np.isfinite(tm_e):
            continue
        p_starts = pd.to_numeric(grp[pep_start_col], errors="coerce")
        p_ends = pd.to_numeric(grp[pep_end_col], errors="coerce")
        if ectodomain_side == "N_terminal":
            n_ecto = int((p_ends <= tm_s).sum())
            n_cyto = int((p_starts >= tm_e).sum())
        else:
            n_ecto = int((p_starts >= tm_e).sum())
            n_cyto = int((p_ends <= tm_s).sum())
        mechanism = (
            "Confirmed_Ectodomain_Shedding (100% Extracellular Peptides)"
            if (n_ecto >= 2 and n_cyto == 0)
            else ("Full_Length_or_Vesicular_Release" if n_cyto > 0 else "Low_Coverage_Ambiguous")
        )
        rows.append(
            {
                "std_gene": str(gene),
                "n_ectodomain_peptides": n_ecto,
                "n_cytoplasmic_peptides": n_cyto,
                "topology_mechanism": mechanism,
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# ENGINE [08]: Signed ULM Upstream TF (CollecTRI) & Kinase Activity Inference
# =============================================================================

def load_default_regulon_edges(
    external_data_dir: Optional[str | Path] = None,
    modality: str = "WCL",
) -> pd.DataFrame:
    """
    Load signed regulator-target edges (`regulator`, `target`, `mor`, `regulator_class`)
    from `external_data/tf_regulons/collectri.tsv` and `curated_proteomics_signatures.json`
    (`collectri_tf_regulons` + `kinase_substrate_regulons`).
    Enforces Modality Scope Gate: for `EV` / `Secretome` / `Plasma`, nuclear TFs are flagged
    or filtered unless explicitly requested, prioritizing trafficking/kinase regulators.
    """
    rows: List[Dict[str, Any]] = []
    sigs = _load_curated_signatures()

    def _parse_regulon_section(sec: Any, reg_class: str) -> None:
        if not isinstance(sec, dict):
            return
        for reg, spec in sec.items():
            reg_u = str(reg).strip().upper()
            if isinstance(spec, list):
                for item in spec:
                    if isinstance(item, (list, tuple)) and len(item) >= 2:
                        rows.append({
                            "regulator": reg_u,
                            "target": str(item[0]).strip().upper(),
                            "mor": float(item[1]),
                            "regulator_class": reg_class,
                        })
            elif isinstance(spec, dict):
                for t in spec.get("positive_targets", []):
                    rows.append({"regulator": reg_u, "target": str(t).upper(), "mor": 1.0, "regulator_class": reg_class})
                for t in spec.get("negative_targets", []):
                    rows.append({"regulator": reg_u, "target": str(t).upper(), "mor": -1.0, "regulator_class": reg_class})

    # 1. Built-in Kinase/Trafficking regulons (valid across all modalities)
    _parse_regulon_section(sigs.get("kinase_substrate_regulons", {}), "Kinase_or_Trafficking")

    # 2. Built-in CollecTRI TF regulons
    _parse_regulon_section(sigs.get("collectri_tf_regulons", {}), "Nuclear_TF")

    # 3. Local `external_data/tf_regulons/collectri.tsv` (43,156 signed interactions)
    ext_dir = _locate_external_data_dir(external_data_dir)
    if ext_dir is not None and (ext_dir / "tf_regulons" / "collectri.tsv").exists():
        try:
            ct_df = pd.read_csv(ext_dir / "tf_regulons" / "collectri.tsv", sep="\t")
            src_col = "source" if "source" in ct_df.columns else ct_df.columns[0]
            tgt_col = "target" if "target" in ct_df.columns else ct_df.columns[1]
            mor_col = "weight" if "weight" in ct_df.columns else ("mor" if "mor" in ct_df.columns else ct_df.columns[2])
            for _, r in ct_df.iterrows():
                rows.append(
                    {
                        "regulator": str(r[src_col]).strip().upper(),
                        "target": str(r[tgt_col]).strip().upper(),
                        "mor": float(r[mor_col]) if pd.notna(r[mor_col]) else 1.0,
                        "regulator_class": "Nuclear_TF",
                    }
                )
        except Exception:
            pass

    if not rows:
        return pd.DataFrame(columns=["regulator", "target", "mor", "regulator_class"])

    reg_df = pd.DataFrame(rows).drop_duplicates(subset=["regulator", "target"], keep="first")
    if str(modality).upper() in ("EV", "EXOSOME", "SECRETOME", "PLASMA"):
        # Modality Scope Gate: in cell-free EVs, prioritize Kinase/Trafficking regulators
        # while tagging Nuclear_TFs as indirect donor-cell reflections
        reg_df["modality_scope_note"] = np.where(
            reg_df["regulator_class"] == "Nuclear_TF",
            "Donor_Cell_Reflection_Only (EV_Scope_Gate)",
            "Direct_Vesicular_or_Signaling_Regulator",
        )
    else:
        reg_df["modality_scope_note"] = "Direct_Cellular_Regulator"
    return reg_df


def infer_upstream_regulator_ulm(
    df: pd.DataFrame,
    regulon_edges: Optional[pd.DataFrame] = None,
    fc_col: str = "log2FC",
    gene_col: str = "std_gene",
    min_targets: int = 4,
    modality: str = "WCL",
    external_data_dir: Optional[str | Path] = None,
) -> pd.DataFrame:
    """
    Engine [08]: Signed Univariate Linear Model (ULM) upstream regulator activity inference
    for Kinases, Trafficking Hubs, and Transcription Factors (CollecTRI/OmniPath).
    If `regulon_edges` is None, automatically loads `load_default_regulon_edges()`.
    Also cross-checks whether the regulator protein itself is detected in `df`.
    """
    if regulon_edges is None or regulon_edges.empty:
        regulon_edges = load_default_regulon_edges(external_data_dir=external_data_dir, modality=modality)

    if fc_col not in df.columns:
        cand_cols = [c for c in df.columns if "log2fc" in c.lower() or "eff_fc" in c.lower()]
        if cand_cols:
            fc_col = cand_cols[0]
        else:
            return pd.DataFrame(columns=["regulator", "n_detected_targets", "activity_t_score", "p_value", "fdr_q_value", "status"])

    sub = df[[gene_col, fc_col]].dropna().copy()
    sub[gene_col] = sub[gene_col].astype(str).str.strip().str.upper()
    sub[fc_col] = pd.to_numeric(sub[fc_col], errors="coerce")
    sub = sub.dropna().groupby(gene_col, as_index=False)[fc_col].mean()
    y_vec = sub.set_index(gene_col)[fc_col]
    n_total = len(y_vec)
    if n_total < 8 or regulon_edges.empty:
        return pd.DataFrame(columns=["regulator", "n_detected_targets", "activity_t_score", "p_value", "fdr_q_value", "status"])

    has_class = "regulator_class" in regulon_edges.columns
    rows = []
    for reg, grp in regulon_edges.groupby("regulator"):
        reg_u = str(reg).strip().upper()
        t_map = dict(zip(grp["target"].astype(str).str.upper(), pd.to_numeric(grp["mor"], errors="coerce").fillna(1.0)))
        shared = [g for g in y_vec.index if g in t_map]
        if len(shared) < min_targets:
            continue
        x = np.array([t_map.get(g, 0.0) for g in y_vec.index], dtype=float)
        y = y_vec.values
        if np.std(x) < 1e-8 or np.std(y) < 1e-8:
            continue
        r = float(np.corrcoef(x, y)[0, 1])
        df_deg = max(n_total - 2, 1)
        t_stat = r * np.sqrt(df_deg / max(1.0 - r * r, 1e-12))
        p_val = float(2.0 * stats.t.sf(abs(t_stat), df=df_deg))
        reg_class = str(grp["regulator_class"].iloc[0]) if has_class else "Regulator"
        self_fc = float(y_vec.loc[reg_u]) if reg_u in y_vec.index else np.nan
        # Sort top responsive targets by |fc * mor|
        scored_targets = sorted(shared, key=lambda g: abs(float(y_vec.loc[g]) * t_map[g]), reverse=True)
        rows.append(
            {
                "regulator": str(reg),
                "regulator_class": reg_class,
                "n_detected_targets": len(shared),
                "activity_t_score": round(float(t_stat), 4),
                "pearson_r": round(r, 4),
                "p_value": p_val,
                "regulator_self_log2FC": round(self_fc, 3) if np.isfinite(self_fc) else np.nan,
                "status": "Activated" if t_stat >= 2.0 else ("Inhibited" if t_stat <= -2.0 else "Neutral"),
                "top_targets": ";".join(scored_targets[:10]),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["regulator", "n_detected_targets", "activity_t_score", "p_value", "fdr_q_value", "status"])
    out = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    out["fdr_q_value"] = benjamini_hochberg(out["p_value"].values)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 4: Systems Biology, CORUM Complexes, Topology & STRING KDA")
    parser.add_argument("--input", "-i", required=True, help="Differential results CSV")
    parser.add_argument("--output-prefix", "-o", default="stage4_systems", help="Output file prefix")
    parser.add_argument("--modality", "-m", default="WCL", help="Modality (WCL, EV, Secretome, etc.)")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    cx_df = audit_corum_complex_co_shifts(df)
    if not cx_df.empty:
        cx_df.to_csv(f"{args.output_prefix}_complexes.csv", index=False)

    ann_df, topo_df = audit_membrane_and_secretory_topology(df)
    if not topo_df.empty:
        topo_df.to_csv(f"{args.output_prefix}_topology_fisher.csv", index=False)

    ulm_df = infer_upstream_regulator_ulm(df, modality=args.modality)
    if not ulm_df.empty:
        ulm_df.to_csv(f"{args.output_prefix}_upstream_ulm.csv", index=False)

    print(
        f"[Stage 4] Evaluated {len(cx_df)} CORUM complexes, "
        f"{len(topo_df)} biophysical topology tests, and {len(ulm_df)} upstream regulators."
    )


if __name__ == "__main__":
    main()

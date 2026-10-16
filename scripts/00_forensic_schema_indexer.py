#!/usr/bin/env python3
"""
Stage 0: Multi-Format Schema Indexer & Experimental Design Inferencer
=====================================================================
Inspects raw or processed proteomics tables (.xlsx, .tsv, .csv, .txt, .parquet)
from DIA-NN, MaxQuant, FragPipe, Proteome Discoverer, Spectronaut, or Olink/SomaScan.
Automatically:
  1. Sanitizes HGNC gene symbols corrupted into Excel dates (e.g., 1-Mar -> MARCHF1).
  2. Constructs composite primary keys (Accession|Gene or Accession|Gene|Site).
  3. Infers Proteomics Modality (M1-M6), Experimental Archetype (A1-A10),
     sample groups, replicate counts, peptide support columns, and intensity scale.
  4. Emits `design.yaml` and `analysis_spec.yaml` for Gate 1 confirmation.
"""

from __future__ import annotations
import argparse
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:
    yaml = None


SKILL_ROOT = Path(__file__).resolve().parent.parent
SIGNATURES_PATH = SKILL_ROOT / "resources" / "curated_proteomics_signatures.json"


def load_curated_signatures() -> Dict[str, Any]:
    if SIGNATURES_PATH.exists():
        with open(SIGNATURES_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
            if "excel_date_gene_fixes" in data and "excel_date_gene_map" not in data:
                data["excel_date_gene_map"] = data["excel_date_gene_fixes"]
            return data
    return {"excel_date_gene_map": {}}


def sanitize_gene_symbol(val: Any, date_map: Optional[Dict[str, str]] = None) -> str:
    """Recover HGNC symbols corrupted by Excel date auto-formatting and strip isoform suffixes."""
    if pd.isna(val):
        return ""
    s = str(val).strip()
    if not s or s.lower() in {"nan", "none", "null", "---"}:
        return ""
    if date_map is None:
        date_map = load_curated_signatures().get("excel_date_gene_map", {})
    if s in date_map:
        return date_map[s]
    # Also check case-insensitive match for Excel date patterns like '01-Mar', 'Mar-01'
    m1 = re.match(r"^0?(\d{1,2})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$", s, re.IGNORECASE)
    if m1:
        key = f"{int(m1.group(1))}-{m1.group(2).capitalize()}"
        if key in date_map:
            return date_map[key]
    m2 = re.match(r"^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-0?(\d{1,2})$", s, re.IGNORECASE)
    if m2:
        key = f"{m2.group(2)}-{m2.group(1).capitalize()}"
        if key in date_map:
            return date_map[key]
        key_rev = f"{m2.group(1).capitalize()}-{int(m2.group(2))}"
        if key_rev in date_map:
            return date_map[key_rev]
    # If semicolon-delimited gene group, take primary gene symbol
    if ";" in s:
        s = s.split(";")[0].strip()
    return s


def extract_gene_from_description(desc: Any) -> str:
    """Extract GN=GENE from UniProt FASTA description header."""
    if pd.isna(desc):
        return ""
    m = re.search(r"GN=([^\s;]+)", str(desc))
    return m.group(1).strip() if m else ""


def detect_Header_and_load_sheet(
    file_path: Path, sheet_name: Optional[str] = None
) -> Dict[str, pd.DataFrame]:
    """Load one or all sheets from an .xlsx, .tsv, .csv, or .parquet file, auto-detecting header row."""
    suffix = file_path.suffix.lower()
    sheets: Dict[str, pd.DataFrame] = {}

    if suffix in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(file_path)
        target_sheets = [sheet_name] if sheet_name else xl.sheet_names
        for sh in target_sheets:
            raw_preview = xl.parse(sh, nrows=15, header=None)
            header_idx = 0
            best_score = -1
            for r_idx in range(min(10, len(raw_preview))):
                row_vals = [str(v).strip().lower() for v in raw_preview.iloc[r_idx].values if pd.notna(v)]
                score = sum(
                    1
                    for v in row_vals
                    if any(
                        k in v
                        for k in [
                            "accession",
                            "gene",
                            "protein",
                            "description",
                            "peptide",
                            "abundance",
                            "intensity",
                            "coverage",
                            "psm",
                        ]
                    )
                )
                if score > best_score and score >= 2:
                    best_score = score
                    header_idx = r_idx
            df = xl.parse(sh, header=header_idx)
            sheets[str(sh)] = df
    elif suffix == ".parquet":
        sheets["main"] = pd.read_parquet(file_path)
    elif suffix in {".tsv", ".txt"}:
        sheets["main"] = pd.read_csv(file_path, sep="\t", low_memory=False)
    else:
        sheets["main"] = pd.read_csv(file_path, low_memory=False)

    return sheets


def infer_columns_and_modality(df: pd.DataFrame, sheet_label: str = "main") -> Dict[str, Any]:
    """Classify metadata vs quantitative sample columns and infer modality/archetype."""
    date_map = load_curated_signatures().get("excel_date_gene_map", {})
    cols = [str(c) for c in df.columns]

    acc_col = next(
        (
            c
            for c in cols
            if any(
                k in c.lower()
                for k in ["accession", "protein.group", "protein.ids", "majority protein ids", "pg.proteinaccessions", "uniprot"]
            )
        ),
        None,
    )
    gene_col = next(
        (
            c
            for c in cols
            if any(k in c.lower() for k in ["gene symbol", "gene_symbol", "genes", "pg.genes", "gene name", "gene"])
            and "description" not in c.lower()
        ),
        None,
    )
    desc_col = next((c for c in cols if "description" in c.lower() or "fasta" in c.lower() or "protein name" in c.lower()), None)
    pep_col = next(
        (
            c
            for c in cols
            if any(
                k in c.lower()
                for k in [
                    "unique peptide count",
                    "unique peptides",
                    "unique_peptides",
                    "unique peptide",
                    "razor + unique peptides",
                    "peptides",
                ]
            )
            and "sequence" not in c.lower()
        ),
        None,
    )
    psm_col = next(
        (
            c
            for c in cols
            if any(k in c.lower() for k in ["psm", "ms/ms count", "precursor", "total peptide count", "total peptides"])
        ),
        None,
    )
    ptm_loc_col = next(
        (
            c
            for c in cols
            if any(k in c.lower() for k in ["localization prob", "ptm.site.confidence", "assayprobability", "phospho", "glygly"])
        ),
        None,
    )

    # Identify quantitative sample columns
    metadata_keywords = {
        "row_id",
        "unnamed",
        "index",
        "rank",
        "length",
        "accession",
        "description",
        "gene",
        "coverage",
        "peptide",
        "psm",
        "unique",
        "razor",
        "mw",
        "kda",
        "calc. pi",
        "score",
        "fdr",
        "q-value",
        "q.value",
        "contaminant",
        "reverse",
        "master",
        "marked as",
        "number of",
        "sequence",
        "modifications",
        "positions",
        "pfam",
        "go ",
        "kegg",
        "reactome",
        "ensembl",
        "entrez",
    }

    precomputed_stat_patterns = (
        "log2 fc",
        "log 2fc",
        "log2fc",
        "logfc",
        "log2 (",
        "log2(",
        "fold change",
        "fold_change",
        "p-value",
        "pvalue",
        "p_value",
        "q-value",
        "qvalue",
        "adj.p",
    )

    sample_cols: List[str] = []
    precomputed_stat_cols: List[str] = []
    for c in cols:
        cl = c.lower().strip()
        if c in {acc_col, gene_col, desc_col, pep_col, psm_col, ptm_loc_col}:
            continue
        if pd.api.types.is_bool_dtype(df[c]) or cl.startswith("is_"):
            continue
        if any(pk in cl for pk in precomputed_stat_patterns):
            precomputed_stat_cols.append(c)
            continue
        if any(mk in cl for mk in metadata_keywords) and not any(
            ak in cl for ak in ["abundance", "intensity", "lfq", "ibaq", "npx", "ratio"]
        ):
            continue
        # Check if column is predominantly numeric or convertible to numeric
        s_num = pd.to_numeric(df[c], errors="coerce")
        non_na_frac = s_num.notna().mean()
        if non_na_frac >= 0.08:
            sample_cols.append(c)

    # Check for Excel date-corrupted genes
    corrupted_examples: List[str] = []
    if gene_col and gene_col in df.columns:
        raw_genes = df[gene_col].dropna().astype(str).tolist()
        for g in raw_genes:
            cleaned = sanitize_gene_symbol(g, date_map)
            if cleaned != g.split(";")[0].strip() and cleaned:
                corrupted_examples.append(f"{g} -> {cleaned}")

    # Inspect numeric scale (raw linear vs log2)
    sample_matrix = (
        df[sample_cols].apply(pd.to_numeric, errors="coerce").astype(float)
        if sample_cols
        else pd.DataFrame()
    )
    if not sample_matrix.empty:
        arr = sample_matrix.to_numpy(dtype=float)
        pos_vals = arr[np.isfinite(arr) & (arr > 0)]
        med_val = float(np.median(pos_vals)) if len(pos_vals) > 0 else 0.0
        max_val = float(np.max(pos_vals)) if len(pos_vals) > 0 else 0.0
        is_log2_scale = bool(med_val < 35.0 and max_val < 65.0)
        missing_frac = float(sample_matrix.isna().mean().mean())
    else:
        med_val, max_val, is_log2_scale, missing_frac = 0.0, 0.0, False, 1.0

    # Infer Modality (M1-M6)
    joined_names = " ".join(sample_cols + [sheet_label] + ([ptm_loc_col] if ptm_loc_col else [])).lower()
    if ptm_loc_col or any(k in joined_names for k in ["phospho", "sty", "ubiquitin", "glygly", "acetyl"]):
        modality = "M2_PTM_Signaling"
    elif any(k in joined_names for k in ["ap-ms", "co-ip", "bioid", "turboid", "apex", "bait", "prey"]):
        modality = "M3_Interactome_Proximity"
    elif any(k in joined_names for k in ["tpp", "cetsa", "pisa", "lip-ms", "37c", "52c", "temp"]):
        modality = "M4_Thermal_Chemoproteomics"
    elif any(k in joined_names for k in ["dose", "ec50", "ic50", "timecourse", "0h", "24h", "48h", "uM", "nM"]):
        modality = "M5_Dose_Time_Course"
    elif any(k in joined_names for k in ["npx", "olink", "somascan", "prm", "mrm"]):
        modality = "M6_Targeted_Affinity"
    else:
        modality = "M1_Global_Abundance"

    # Infer Experimental Archetype (A1-A10)
    if any(k in joined_names for k in ["biofluid", "fluid", "plasma", "serum"]) and any(
        k in joined_names for k in ["tissue", "tumor", "t1", "n1"]
    ):
        archetype = "A5_Clinical_Paired_MultiCompartment_SGI"
    elif any(k in joined_names for k in ["igg", "humanized", "parental", "chimeric", "mab", "antibody"]):
        archetype = "A4_Biologics_Antibody_Engineering"
    elif any(k in joined_names for k in ["sh", "kd", "si"]) and any(k in joined_names for k in ["oe", "overexpress"]):
        archetype = "A2_Reciprocal_KD_x_OE"
    elif len(set(re.findall(r"sh[a-z0-9]+", joined_names))) >= 2:
        archetype = "A3_Multi_Perturbation_9Sector_Epistasis"
    elif modality == "M2_PTM_Signaling":
        archetype = "A7_PTM_Stoichiometry_Kinase_E3"
    elif modality == "M3_Interactome_Proximity":
        archetype = "A8_Interactome_SAINT_CompPASS"
    elif modality == "M4_Thermal_Chemoproteomics":
        archetype = "A9_Thermal_Stability_NPARC_PISA"
    elif modality == "M5_Dose_Time_Course":
        archetype = "A10_Dose_Response_CurveCurator_Spline"
    else:
        archetype = "A1_Single_Perturbation_or_Cohort_DE"

    # Group sample columns by stripping replicate suffix (e.g., -1, _rep1, .1)
    inferred_groups: Dict[str, List[str]] = {}
    for sc in sample_cols:
        clean_sc = re.sub(r"^abundance[s]?(\s*\([^)]+\))?\s*:\s*", "", sc, flags=re.IGNORECASE).strip()
        grp = re.sub(r"([_-]rep\d+|[_-]\d+|\s+\d+)$", "", clean_sc, flags=re.IGNORECASE).strip()
        inferred_groups.setdefault(grp, []).append(sc)

    max_reps = max((len(v) for v in inferred_groups.values()), default=0)
    if max_reps >= 5:
        stat_regime = "Large_Cohort_limpa_DEqMS_Hurdle"
    elif max_reps >= 3:
        stat_regime = "Standard_Replicate_limpa_DEqMS"
    elif max_reps == 2:
        stat_regime = "Small_Replicate_eBayes_SharedControl_Null"
    else:
        stat_regime = "Single_Tube_Empirical_Null_Calibrated_E3"

    return {
        "sheet_name": sheet_label,
        "n_rows": int(len(df)),
        "n_sample_cols": len(sample_cols),
        "accession_col": acc_col,
        "gene_col": gene_col,
        "description_col": desc_col,
        "unique_peptides_col": pep_col,
        "peptide_support": "present" if (pep_col and pep_col in df.columns) else "unknown",
        "psm_col": psm_col,
        "ptm_localization_col": ptm_loc_col,
        "sample_cols": sample_cols,
        "precomputed_stat_cols": precomputed_stat_cols,
        "inferred_groups": inferred_groups,
        "max_replicates_per_group": max_reps,
        "is_log2_scale": is_log2_scale,
        "median_positive_value": round(med_val, 3),
        "overall_missing_fraction": round(missing_frac, 4),
        "modality": modality,
        "archetype": archetype,
        "recommended_statistical_engine": stat_regime,
        "excel_date_corruptions_detected": corrupted_examples[:10],
    }


def build_composite_Feature_table(df: pd.DataFrame, schema: Dict[str, Any]) -> pd.DataFrame:
    """Standardize accession, sanitized gene symbol, composite key, and peptide counts."""
    date_map = load_curated_signatures().get("excel_date_gene_map", {})
    out = df.copy()
    acc_col = schema.get("accession_col")
    gene_col = schema.get("gene_col")
    desc_col = schema.get("description_col")
    pep_col = schema.get("unique_peptides_col")

    if acc_col and acc_col in out.columns:
        out["std_accession"] = out[acc_col].astype(str).str.strip()
    else:
        out["std_accession"] = [f"PROT_{i:05d}" for i in range(len(out))]

    if gene_col and gene_col in out.columns:
        out["std_gene"] = [sanitize_gene_symbol(g, date_map) for g in out[gene_col]]
    elif desc_col and desc_col in out.columns:
        out["std_gene"] = [ sanitize_gene_symbol(extract_gene_from_description(d), date_map) for d in out[desc_col] ]
    else:
        out["std_gene"] = out["std_accession"]

    # Fallback to accession when gene is empty
    out["std_gene"] = [
        g if g else str(a).split(";")[0].strip()
        for g, a in zip(out["std_gene"], out["std_accession"])
    ]
    out["composite_id"] = out["std_accession"] + "|" + out["std_gene"]

    if pep_col and pep_col in out.columns:
        out["std_unique_peptides"] = pd.to_numeric(out[pep_col], errors="coerce").fillna(1).astype(int)
        out["peptide_support_known"] = True
        out["peptide_tier"] = np.where(out["std_unique_peptides"] >= 2, "Core_ge2_peptides", "Exploratory_1_peptide")
    else:
        out["std_unique_peptides"] = 1
        out["peptide_support_known"] = False
        out["peptide_tier"] = "Unverified_Peptide_Count"

    return out


def index_proteomics_file(
    input_path: str | Path,
    output_dir: Optional[str | Path] = None,
    organ_preset: Optional[str] = None,
    report_language: str = "bilingual",
) -> Dict[str, Any]:
    """Run Stage 0 schema indexing on an input file and emit design.yaml + analysis_spec.yaml."""
    in_path = Path(input_path).resolve()
    sheets = detect_Header_and_load_sheet(in_path)
    sheet_summaries = {}
    for sh_name, df in sheets.items():
        if len(df) > 0 and len(df.columns) >= 3:
            sheet_summaries[sh_name] = infer_columns_and_modality(df, sh_name)

    summary_payload = {
        "input_file": str(in_path),
        "n_sheets": len(sheet_summaries),
        "sheets": sheet_summaries,
        "gate_1_confirmation_required": True,
        "organ_context_preset": organ_preset,
        "deliverable_language": report_language,
    }

    if output_dir:
        out_dir = Path(output_dir).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        json_out = out_dir / "stage0_schema_report.json"
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump(summary_payload, f, indent=2)

        design_spec = {
            "project_name": in_path.stem,
            "input_file": str(in_path),
            "organ_preset": organ_preset,
            "deliverable_language": report_language,
            "qc_rules": {
                "min_unique_peptides_core": 2,
                "fc_cutoff_log2": 0.5849625,
                "fdr_cutoff": 0.05,
                "tost_equivalence_margin_log2": 0.321928,
                "max_median_shift_for_naive_norm_log2": 0.15,
                "empirical_null_fdr_e2_ceiling": 0.25,
            },
            "sheets": sheet_summaries,
        }
        if yaml is not None:
            for yname in ("analysis_spec.yaml", "design.yaml"):
                with open(out_dir / yname, "w", encoding="utf-8") as f:
                    yaml.safe_dump(design_spec, f, sort_keys=False)
        else:
            for jname in ("analysis_spec.json", "design.json"):
                with open(out_dir / jname, "w", encoding="utf-8") as f:
                    json.dump(design_spec, f, indent=2)

    return summary_payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 0: Proteomics Schema Indexer & Design Inferencer")
    parser.add_argument("--input", "-i", required=True, help="Path to input .xlsx/.tsv/.csv/.parquet file")
    parser.add_argument("--outdir", "-o", default=".", help="Directory to write design.yaml & stage0_schema_report.json")
    parser.add_argument("--organ-preset", default=None, help="Optional organ/disease context preset")
    parser.add_argument("--lang", default="bilingual", choices=["en", "zh", "bilingual"], help="Deliverable language")
    args = parser.parse_args()

    res = index_proteomics_file(args.input, args.outdir, args.organ_preset, args.lang)
    print(f"[Stage 0] Indexed {res['input_file']} ({res['n_sheets']} sheet(s)):")
    for sh, info in res["sheets"].items():
        print(
            f"  - Sheet '{sh}': {info['n_rows']} rows, {info['n_sample_cols']} sample cols, "
            f"Modality={info['modality']}, Archetype={info['archetype']}, "
            f"Engine={info['recommended_statistical_engine']}"
        )


if __name__ == "__main__":
    main()

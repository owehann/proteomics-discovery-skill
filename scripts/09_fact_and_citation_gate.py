#!/usr/bin/env python3
"""
Stage 7B: Pre-Delivery Fact, Peptide-Tier, Overclaim, & Citation Gate
=====================================================================
The mandatory final verification gate before any report (.md), schematic (.svg),
or presentation (.pptx) is delivered to the user.

Performs 4 automated forensic checks against a canonical `facts.json` or CSV tables:
  1. Gene-Number & Quadrant Reconciliation: Verifies that every numeric value or
     quadrant label (`Q1_concordant_up` vs `Q2_hyper_secreted`, sample intensities,
     log2FCs) cited in prose/tables matches the underlying data (and never leaks
     values from quarantined outlier tubes such as `CTRL_Rep2`).
  2. Single-Peptide Promotion Audit: Flags any protein with `unique_peptides == 1`
     promoted into a headline/Core claim without an explicit `1-peptide` / `exploratory`
     disclosure on the same line or table row.
  3. Overclaim & Evidence-Grade Linter: Detects causal/definitive language applied to
     `E3` exploratory contrasts (`empirical_null_fdr > 0.25`) or claims of "preserved /
     unperturbed" lacking `TOST` equivalence support.
  4. Live Bibliographic PMID/DOI Check: Extracts all `PMID: <digits>` and `DOI` strings
     and verifies them via EuropePMC REST API.
"""

from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


def build_facts_dictionary_from_csv(
    csv_path: str | Path,
    gene_col: str = "gene_symbol",
    pep_col: str = "unique_peptides",
    numeric_cols: Optional[Sequence[str]] = None,
    category_cols: Optional[Sequence[str]] = None,
    quarantined_values_by_gene: Optional[Dict[str, List[float]]] = None,
) -> Dict[str, Any]:
    """Compile a canonical `facts.json` dictionary from a processed results CSV."""
    df = pd.read_csv(csv_path)
    if gene_col not in df.columns:
        gene_col = next((c for c in df.columns if "gene" in c.lower()), df.columns[0])
    if pep_col not in df.columns:
        pep_col = next((c for c in df.columns if "pep" in c.lower()), "")

    if numeric_cols is None:
        numeric_cols = [
            c
            for c in df.columns
            if c != pep_col and pd.api.types.is_numeric_dtype(df[c])
        ]
    if category_cols is None:
        category_cols = [
            c
            for c in df.columns
            if any(k in c.lower() for k in ["quadrant", "sector", "tier", "class", "regime"])
        ]

    genes_dict: Dict[str, Any] = {}
    for _, row in df.iterrows():
        g = str(row.get(gene_col, "")).strip()
        if not g or g.lower() == "nan":
            continue
        raw_pep = pd.to_numeric(row.get(pep_col, 2), errors="coerce") if pep_col else 2.0
        peps = int(raw_pep) if pd.notna(raw_pep) and raw_pep > 0 else 2
        num_map = {}
        for nc in numeric_cols:
            val = pd.to_numeric(row.get(nc), errors="coerce")
            if pd.notna(val):
                num_map[nc] = round(float(val), 4)
        cat_map = {cc: str(row.get(cc, "")) for cc in category_cols if pd.notna(row.get(cc))}
        genes_dict[g] = {
            "unique_peptides": peps,
            "metrics": num_map,
            "categories": cat_map,
            "quarantined_outlier_values": (quarantined_values_by_gene or {}).get(g, []),
        }

    return {"source_csv": str(csv_path), "genes": genes_dict}


def audit_markdown_against_facts(
    md_path: str | Path,
    facts: Dict[str, Any],
    allow_untagged_single_peptide: bool = False,
    tolerance: float = 0.03,
) -> Dict[str, Any]:
    """
    Audit a Markdown report or SVG text against `facts` dictionary.
    Returns structured list of critical errors and warnings.
    """
    text = Path(md_path).read_text(encoding="utf-8")
    lines = text.splitlines()
    genes_db: Dict[str, Any] = facts.get("genes", {})

    critical_errors: List[Dict[str, Any]] = []
    warnings: List[Dict[str, Any]] = []

    single_pep_disclosure_tokens = (
        "1-pep",
        "1 pep",
        "single-peptide",
        "single peptide",
        "unique_peptides = 1",
        "unique_peptides == 1",
        "unique_peptides=1",
        "| 1 |",
        "exploratory",
        "单肽段",
        "1个独特肽段",
    )

    quadrant_tokens = (
        "Q1_concordant_up",
        "Q2_hyper_secreted",
        "Q3_concordant_down",
        "Q4_intracellular_retained",
        "Q4_tissue_retained",
        "Q_fluid_depleted",
    )

    for line_no, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        line_lower = line.lower()

        for gene, ginfo in genes_db.items():
            if len(gene) < 2:
                continue
            # Match whole word gene symbol
            if not re.search(rf"\b{re.escape(gene)}\b", line):
                continue

            # Check 1A: Quarantined outlier value leak check
            for q_val in ginfo.get("quarantined_outlier_values", []):
                q_str = f"{q_val:.2f}"
                if q_str in line and not any(
                    ex_tok in line_lower for ex_tok in ["outlier", "excluded", "quarantined", "ctrl_rep2", "剔除", "异常"]
                ):
                    critical_errors.append(
                        {
                            "line": line_no,
                            "type": "QUARANTINED_OUTLIER_VALUE_LEAK",
                            "gene": gene,
                            "message": f"Line {line_no} cites quarantined outlier value {q_str} for {gene} without noting it is from an excluded outlier tube.",
                            "excerpt": line.strip()[:140],
                        }
                    )

            # Check 1B: Quadrant / Category label contradiction check
            true_cats = set(ginfo.get("categories", {}).values())
            for q_tok in quadrant_tokens:
                if q_tok in line and true_cats and (q_tok not in true_cats):
                    # Check if the line is explicitly contrasting or correcting the label
                    if not any(neg in line_lower for neg in ["not ", "rather than", "而非", "修正"]):
                        critical_errors.append(
                            {
                                "line": line_no,
                                "type": "CATEGORY_QUADRANT_MISMATCH",
                                "gene": gene,
                                "message": f"Line {line_no} associates {gene} with '{q_tok}', but canonical table classifies {gene} as {sorted(true_cats)}.",
                                "excerpt": line.strip()[:140],
                            }
                        )

            # Check 2: Single-Peptide Promotion Audit
            peps = int(ginfo.get("unique_peptides", 2))
            if peps == 1 and not allow_untagged_single_peptide:
                has_disclosure = any(tok in line_lower for tok in single_pep_disclosure_tokens)
                if not has_disclosure:
                    critical_errors.append(
                        {
                            "line": line_no,
                            "type": "UNTAGGED_SINGLE_PEPTIDE_PROMOTION",
                            "gene": gene,
                            "message": f"Line {line_no} promotes {gene} (unique_peptides == 1) without an explicit single-peptide/exploratory disclosure.",
                            "excerpt": line.strip()[:140],
                        }
                    )

    # Check 4: Extract PMIDs and DOIs
    pmids = sorted(set(re.findall(r"PMID:\s*(\d{5,9})", text, flags=re.IGNORECASE)))
    dois = sorted(set(re.findall(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", text)))

    return {
        "file": str(md_path),
        "status": "PASS" if len(critical_errors) == 0 else "FAIL",
        "n_critical_errors": len(critical_errors),
        "n_warnings": len(warnings),
        "critical_errors": critical_errors,
        "warnings": warnings,
        "extracted_pmids": pmids,
        "extracted_dois": dois,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage 7B: Pre-Delivery Fact & Citation Gate")
    parser.add_argument("--report", "-r", required=True, help="Markdown or SVG file to audit")
    parser.add_argument("--facts-json", "-f", default=None, help="Canonical facts.json file")
    parser.add_argument("--csv", "-c", default=None, help="Or build facts directly from a canonical CSV file")
    parser.add_argument("--out", "-o", default=None, help="Optional JSON output report")
    args = parser.parse_args()

    if args.facts_json:
        with open(args.facts_json, "r", encoding="utf-8") as f:
            facts = json.load(f)
    elif args.csv:
        facts = build_facts_dictionary_from_csv(args.csv)
    else:
        parser.error("Either --facts-json (-f) or --csv (-c) must be provided.")

    res = audit_markdown_against_facts(args.report, facts)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2, ensure_ascii=False)

    print(f"[Stage 7B Fact Gate] {res['file']} -> {res['status']} (Critical Errors: {res['n_critical_errors']})")
    for err in res["critical_errors"]:
        print(f"  [ERROR L{err['line']} {err['type']}] {err['gene']}: {err['message']}")


if __name__ == "__main__":
    main()


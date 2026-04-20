"""
Audit the manuscript's §2.2 claim: "2,145 of 2,209 Exparel primary-suspect
records (97%) would have been misclassified without drugname-first precedence."

Compares the naive classifier (prod_ai.fillna(drugname) → classify_drug) against
the correct classifier (drugname first, prod_ai as fallback for pre-2014 Q3).
Prints both counts and the misclassification rate.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from ascii_pipeline import LIPOSOMAL_PATTERN, PLAIN_BUPI_PATTERN, ROPI_PATTERN


def classify_naive(drugname, prod_ai) -> str | None:
    """What the pipeline would produce with prod_ai.fillna(drugname)."""
    s = prod_ai if isinstance(prod_ai, str) and prod_ai.strip() else drugname
    if not isinstance(s, str):
        return None
    if LIPOSOMAL_PATTERN.search(s):
        return "exparel"
    if ROPI_PATTERN.search(s):
        return "ropivacaine"
    if PLAIN_BUPI_PATTERN.search(s):
        return "bupivacaine_plain"
    return None


def classify_correct(drugname, prod_ai) -> str | None:
    """drugname-first, most-specific-wins."""
    for field in (drugname, prod_ai):
        if isinstance(field, str) and LIPOSOMAL_PATTERN.search(field):
            return "exparel"
    for field in (drugname, prod_ai):
        if isinstance(field, str) and ROPI_PATTERN.search(field):
            return "ropivacaine"
    for field in (drugname, prod_ai):
        if isinstance(field, str) and PLAIN_BUPI_PATTERN.search(field):
            return "bupivacaine_plain"
    return None


def _load_drug_from_quarter(q: Path) -> pd.DataFrame | None:
    for sub in ("ASCII", "ascii"):
        if (q / sub).exists():
            ascii_dir = q / sub
            break
    else:
        ascii_dir = q
    drug_f = [p for p in (list(ascii_dir.glob("*.txt")) + list(ascii_dir.glob("*.TXT")))
              if p.stem.upper().startswith("DRUG")]
    if not drug_f:
        return None
    df = pd.read_csv(drug_f[0], sep="$", encoding="latin-1",
                     low_memory=False, on_bad_lines="skip")
    df.columns = [c.lstrip("\ufeff").lstrip("ï»¿").lower() for c in df.columns]
    keep = [c for c in ("primaryid", "role_cod", "drugname", "prod_ai") if c in df.columns]
    return df[keep]


def main() -> None:
    data = Path("data/raw")
    quarter_dirs = sorted(
        [p for p in data.iterdir() if p.is_dir() and "faers" in p.name.lower()]
    )
    frames = []
    for q in quarter_dirs:
        d = _load_drug_from_quarter(q)
        if d is not None:
            frames.append(d)
    drug = pd.concat(frames, ignore_index=True)
    print(f"loaded DRUG: {len(drug):,} rows")

    # Restrict to PS role
    drug_ps = drug[drug["role_cod"] == "PS"].copy()
    if "prod_ai" not in drug_ps.columns:
        drug_ps["prod_ai"] = None
    print(f"PS-role records: {len(drug_ps):,}")

    # Correct classifier: find all Exparel records
    drug_ps["correct"] = [
        classify_correct(dn, pa)
        for dn, pa in zip(drug_ps["drugname"], drug_ps["prod_ai"])
    ]
    drug_ps["naive"] = [
        classify_naive(dn, pa)
        for dn, pa in zip(drug_ps["drugname"], drug_ps["prod_ai"])
    ]

    exparel_correct = drug_ps[drug_ps["correct"] == "exparel"]
    misclassified = exparel_correct[exparel_correct["naive"] != "exparel"]

    n_total_ps = len(exparel_correct)
    n_mis = len(misclassified)
    print(f"\n=== Record-level (raw PS records) ===")
    print(f"Exparel under correct classifier:  {n_total_ps:,}")
    print(f"Misclassified under naive:         {n_mis:,}")
    print(f"Misclassification rate:            {100 * n_mis / n_total_ps:.1f}%")

    # What were they misclassified as?
    print(f"\nNaive-classifier labels for true Exparel records:")
    print(misclassified["naive"].fillna("(None/unclassified)").value_counts().to_string())

    # Case-level (after caseid-level dedup via primaryid uniqueness)
    exparel_pids = set(exparel_correct["primaryid"].unique())
    # Use analysis_dataset.parquet if available to get deduplicated case count
    try:
        ds = pd.read_parquet("out/analysis_dataset.parquet")
        ex_case_pids = set(ds.loc[ds["drug_group"] == "exparel", "primaryid"].unique())
        print(f"\n=== Case-level (after caseid dedup) ===")
        print(f"Exparel unique primaryids (raw):   {len(exparel_pids):,}")
        print(f"Exparel cases in analysis_dataset: {len(ex_case_pids):,}")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()

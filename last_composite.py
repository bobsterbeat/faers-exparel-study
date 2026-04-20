"""
LAST-composite endpoint analysis
================================
Aggregates all PTs in the Local Anesthetic Systemic Toxicity (LAST) spectrum
into a single composite case-level endpoint, then runs disproportionality
vs both bupivacaine and ropivacaine comparators — and stratified by era.

The case-level narrative review of the 64 early Exparel reports showed that
individual LAST PTs were each at ~ROR 1 vs bupivacaine but the *pattern*
(cardiac + CNS toxicity together) was unmistakable. Composite endpoints
capture this where PT-level analysis misses it.

Protocol note: the LAST PT list is pre-specified (locked before analysis),
not data-driven. Disclose in methods per READUS-PV.

Usage:
    python last_composite.py \\
        --dataset out/analysis_dataset.parquet \\
        --cases   out/exparel_cases.csv \\
        --out     out/last_composite.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# PTs in the LAST spectrum. Extended from the protocol's minimal list based
# on what the case-narrative review surfaced (Depressed level of consciousness,
# Cardio-respiratory arrest, Hypoxia, etc.) — still pre-specified.
LAST_PTS = {
    # Cardiac
    "cardiac arrest",
    "cardio-respiratory arrest",
    "cardiac disorder",
    "ventricular tachycardia",
    "ventricular fibrillation",
    "ventricular arrhythmia",
    "arrhythmia",
    "bradycardia",
    "hypotension",
    "blood pressure decreased",
    "cardiovascular insufficiency",
    # CNS
    "seizure",
    "convulsion",
    "generalised tonic-clonic seizure",
    "loss of consciousness",
    "depressed level of consciousness",
    "unresponsive to stimuli",
    # Respiratory (often a manifestation of CNS-depression-induced hypoventilation)
    "hypoxia",
    "dyspnoea",
    "hypercapnia",
    "respiratory depression",
    "respiratory failure",
    "pulmonary oedema",
    "acute pulmonary oedema",
    # Metabolic marker of end-organ damage
    "acidosis",
    # Nonspecific toxicity marker used by FAERS reporters
    "toxicity to various agents",
}


def case_level_contingency(
    ds: pd.DataFrame, pts: set[str], target: str, comparator: str,
) -> dict:
    """
    a = unique target cases with ≥1 PT in `pts`
    b = unique comparator cases with ≥1 PT in `pts`
    c = target cases without any PT in `pts`
    d = comparator cases without any PT in `pts`
    """
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    ds["hit"] = ds["pt_lower"].isin(pts)

    target_cases = ds.loc[ds["drug_group"] == target, "primaryid"].nunique()
    comp_cases = ds.loc[ds["drug_group"] == comparator, "primaryid"].nunique()

    a = ds.loc[(ds["drug_group"] == target) & ds["hit"], "primaryid"].nunique()
    b = ds.loc[(ds["drug_group"] == comparator) & ds["hit"], "primaryid"].nunique()
    c = target_cases - a
    d = comp_cases - b
    return {
        "target": target, "comparator": comparator,
        "a": a, "b": b, "c": c, "d": d,
        "target_total": target_cases, "comparator_total": comp_cases,
    }


def compute_metrics(row: dict) -> dict:
    a, b, c, d = row["a"], row["b"], row["c"], row["d"]
    aa, bb, cc, dd = [x + 0.5 if x == 0 else x for x in (a, b, c, d)]
    ror = (aa / cc) / (bb / dd)
    se = np.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    ln_ror = np.log(ror)
    prr = (aa / (aa + cc)) / (bb / (bb + dd))
    N = a + b + c + d
    expected_a = (a + b) * (a + c) / N if N else 0
    ic = np.log2((a + 0.5) / (expected_a + 0.5)) if expected_a else np.nan
    var_ic = (1 / np.log(2) ** 2) * (
        (N - a + 0.5) / ((a + 0.5) * (1 + N))
        + (N - (a + b) + 0.5) / (((a + b) + 0.5) * (1 + N))
        + (N - (a + c) + 0.5) / (((a + c) + 0.5) * (1 + N))
    ) if N else np.nan
    ic025 = ic - 1.96 * np.sqrt(var_ic) if not np.isnan(ic) else np.nan
    chi2 = ((a - expected_a) ** 2 / expected_a) if expected_a else 0
    signal = (a >= 3) and (np.exp(ln_ror - 1.96 * se) > 1) and (prr >= 2) and (chi2 >= 4) and (ic025 > 0)
    return {
        **row,
        "ROR": ror,
        "ROR_lower": float(np.exp(ln_ror - 1.96 * se)),
        "ROR_upper": float(np.exp(ln_ror + 1.96 * se)),
        "PRR": prr,
        "chi2": chi2,
        "IC": ic,
        "IC025": ic025,
        "signal": bool(signal),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, type=Path,
                    help="analysis_dataset.parquet from ascii_pipeline.py")
    ap.add_argument("--cases", required=True, type=Path,
                    help="exparel_cases.csv (for era stratification)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ds = pd.read_parquet(args.dataset)
    cases = pd.read_csv(args.cases)

    rows = []

    # Overall LAST-composite vs each comparator
    for comparator, tier in [("bupivacaine_plain", "primary"),
                             ("ropivacaine", "secondary")]:
        r = case_level_contingency(ds, LAST_PTS, "exparel", comparator)
        r = compute_metrics(r)
        r["tier"] = tier
        r["era"] = "all"
        rows.append(r)

    # Era stratification with ERA-MATCHED COMPARATOR POOLS. Earlier version
    # restricted only Exparel to the era and used the full-period comparator
    # pool as background — which is methodologically fragile (bupivacaine's
    # decade-long pre-Exparel baseline dominated). Now both sides share the
    # same era, using out/case_meta.parquet (primaryid → event_year for all
    # 17M FAERS cases).
    #   pre-2018 : before FDA approved interscalene block (ISB) indication
    #   2018-2020: post-ISB, pre-Anesthesiology-2021-publication
    #   2021+    : post-Anesthesiology ruling and Pacira v. ASA litigation
    era_bins = {
        "pre_2018":   (0, 2017),
        "2018_2020":  (2018, 2020),
        "2021_plus":  (2021, 9999),
    }
    case_meta = pd.read_parquet("out/case_meta.parquet")
    year_map = case_meta.set_index("primaryid")["event_year"]

    for era, (lo, hi) in era_bins.items():
        era_pids = set(year_map[(year_map >= lo) & (year_map <= hi)].index)
        ds_era = ds[ds["primaryid"].isin(era_pids)]
        for comparator, tier in [("bupivacaine_plain", "primary"),
                                 ("ropivacaine", "secondary")]:
            r = case_level_contingency(ds_era, LAST_PTS, "exparel", comparator)
            r = compute_metrics(r)
            r["tier"] = tier
            r["era"] = era
            rows.append(r)

    out = pd.DataFrame(rows)
    cols = ["era", "tier", "comparator", "target_total", "comparator_total",
            "a", "b", "c", "d",
            "ROR", "ROR_lower", "ROR_upper", "PRR", "chi2",
            "IC", "IC025", "signal"]
    out = out[cols]
    out.to_csv(args.out, index=False)
    print(out.to_string(index=False))
    print(f"\nsaved → {args.out}")


if __name__ == "__main__":
    main()

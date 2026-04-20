"""
Duration and permanence stratification for the prolonged-block composite.

FAERS does not record block duration in hours. This script produces three
indirect analyses that together address whether the signal reflects abnormal
or permanent deficit rather than expected Exparel pharmacology:

  1. Pathology-implying vs expected-pharmacology PT split (composite ROR)
  2. Disability (DS) and other serious outcomes among flagged cases
  3. Event-to-report lag distribution (>30d, >90d, >180d)

Usage:
    python duration_analysis.py \\
        --dataset out/analysis_dataset.parquet \\
        --cases   out/exparel_cases.csv \\
        --out-dir out/duration/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# PTs whose MedDRA definition implies abnormal duration or permanent deficit.
# These cannot be expected <72h pharmacology.
PATHOLOGY_PTS = {
    "peroneal nerve palsy", "peroneal nerve injury",
    "neuromuscular block prolonged",
    "nerve injury", "neuropathy peripheral", "peripheral nerve injury",
    "peripheral motor neuropathy", "peripheral sensory neuropathy",
    "brachial plexopathy",
    "paresis", "hemiparesis", "monoparesis", "paralysis",
    "neuralgia",
    "motor dysfunction",
}

# PTs that could reflect expected Exparel pharmacology at <72h exposure.
# Disproportional reporting vs bupivacaine is still interesting (it suggests
# greater-than-expected rates) but cannot distinguish normal from prolonged.
# NOTE: "anaesthesia" is intentionally excluded from this set — the PT is
# terminologically ambiguous (the procedural state vs. an unexpected
# reported reaction) and cannot be cleanly classified.
EXPECTED_PTS = {
    "hypoaesthesia", "paraesthesia",
    "muscular weakness", "sensory loss",
}

# PTs that suggest systemic absorption rather than local/perineural block —
# e.g., oral hypoaesthesia after a non-oral injection implies circulating
# drug. Analyzed separately since the mechanism and clinical implication
# differ from the local prolonged-block hypothesis.
SYSTEMIC_ABSORPTION_PTS = {
    "hypoaesthesia oral",
}


def _compute_2x2(a: int, b: int, c: int, d: int) -> dict:
    aa, bb, cc, dd = [x + 0.5 if x == 0 else x for x in (a, b, c, d)]
    ror = (aa / cc) / (bb / dd)
    se = np.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    return {
        "a": a, "b": b, "c": c, "d": d,
        "ROR": ror,
        "ROR_lower": float(np.exp(np.log(ror) - 1.96 * se)),
        "ROR_upper": float(np.exp(np.log(ror) + 1.96 * se)),
    }


def _composite(ds: pd.DataFrame, pts: set[str], comparator: str) -> dict:
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    ex_total = ds.loc[ds["drug_group"] == "exparel", "primaryid"].nunique()
    comp_total = ds.loc[ds["drug_group"] == comparator, "primaryid"].nunique()
    a = ds.loc[(ds["drug_group"] == "exparel") & ds["pt_lower"].isin(pts), "primaryid"].nunique()
    b = ds.loc[(ds["drug_group"] == comparator) & ds["pt_lower"].isin(pts), "primaryid"].nunique()
    c = ex_total - a
    d = comp_total - b
    return {
        "exparel_total": ex_total,
        "comparator_total": comp_total,
        **_compute_2x2(a, b, c, d),
        "exparel_pct": 100 * a / ex_total if ex_total else np.nan,
        "comparator_pct": 100 * b / comp_total if comp_total else np.nan,
    }


def analysis_1_pt_split(ds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, pts in [("pathology_implying", PATHOLOGY_PTS),
                       ("expected_pharmacology", EXPECTED_PTS),
                       ("systemic_absorption", SYSTEMIC_ABSORPTION_PTS)]:
        for comp in ["bupivacaine_plain", "ropivacaine"]:
            r = _composite(ds, pts, comp)
            r["stratum"] = label
            r["comparator"] = comp
            rows.append(r)
    return pd.DataFrame(rows)


def analysis_2_outcomes(ds: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()

    rows = []
    for label, pts in [("all_prolonged_block", PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS),
                       ("pathology_implying", PATHOLOGY_PTS),
                       ("expected_pharmacology", EXPECTED_PTS),
                       ("systemic_absorption", SYSTEMIC_ABSORPTION_PTS)]:
        pids = set(ds.loc[
            (ds["drug_group"] == "exparel") & ds["pt_lower"].isin(pts),
            "primaryid"
        ])
        sub = cases[cases["primaryid"].isin(pids)]
        n = len(sub)
        outcome_counts = {
            code: int(sub["outcomes"].str.contains(code, na=False).sum())
            for code in ("DE", "LT", "DS", "HO", "RI", "CA")
        }
        rows.append({
            "stratum": label,
            "n_cases": n,
            **{f"n_{code}": v for code, v in outcome_counts.items()},
            **{f"pct_{code}": round(100*v/n, 2) if n else np.nan
               for code, v in outcome_counts.items()},
        })
    return pd.DataFrame(rows)


def analysis_3_lag(ds: pd.DataFrame, cases: pd.DataFrame) -> pd.DataFrame:
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    pids = set(ds.loc[
        (ds["drug_group"] == "exparel")
        & ds["pt_lower"].isin(PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS),
        "primaryid"
    ])
    sub = cases[cases["primaryid"].isin(pids)].copy()
    sub["event_dt_p"] = pd.to_datetime(sub["event_dt"], format="%Y%m%d", errors="coerce")
    sub["init_fda_dt_p"] = pd.to_datetime(sub["init_fda_dt"], format="%Y%m%d", errors="coerce")
    sub["lag_days"] = (sub["init_fda_dt_p"] - sub["event_dt_p"]).dt.days
    lag = sub["lag_days"].dropna()
    return pd.DataFrame([{
        "n_cases_with_event_dt": int(len(lag)),
        "median_lag_days": float(lag.median()),
        "q25_lag_days": float(lag.quantile(0.25)),
        "q75_lag_days": float(lag.quantile(0.75)),
        "n_lag_gt_30d":  int((lag > 30).sum()),
        "n_lag_gt_90d":  int((lag > 90).sum()),
        "n_lag_gt_180d": int((lag > 180).sum()),
        "pct_lag_gt_30d":  round(100 * (lag > 30).sum() / len(lag), 1),
        "pct_lag_gt_90d":  round(100 * (lag > 90).sum() / len(lag), 1),
        "pct_lag_gt_180d": round(100 * (lag > 180).sum() / len(lag), 1),
    }])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, type=Path)
    ap.add_argument("--cases", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(exist_ok=True, parents=True)

    ds = pd.read_parquet(args.dataset)
    cases = pd.read_csv(args.cases)

    a1 = analysis_1_pt_split(ds)
    a1.to_csv(args.out_dir / "pt_split_pathology_vs_expected.csv", index=False)
    print("=== 1. Pathology-implying vs expected-pharmacology composite ===")
    print(a1.to_string(index=False))

    a2 = analysis_2_outcomes(ds, cases)
    a2.to_csv(args.out_dir / "outcomes_by_stratum.csv", index=False)
    print("\n=== 2. Serious-outcome distribution by stratum ===")
    print(a2.to_string(index=False))

    a3 = analysis_3_lag(ds, cases)
    a3.to_csv(args.out_dir / "event_to_report_lag.csv", index=False)
    print("\n=== 3. Event-to-report lag (prolonged-block cases) ===")
    print(a3.to_string(index=False))

    print(f"\nAll outputs → {args.out_dir}/")


if __name__ == "__main__":
    main()

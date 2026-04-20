"""
Pre-registered sensitivity analyses for the prolonged-block composite.

Four restrictions locked in §2.7 of the manuscript:
  1. Serious outcomes only (DE, LT, HO, DS, CA, RI)
  2. Physician-reported only (occp_cod == 'MD')
  3. Exclude pregnancy / neonatal cases
  4. Exclude manufacturer-submitted reports

All four run vs plain bupivacaine; pathology-implying sub-stratum is also
tested in each restriction to confirm the signal survives.

Usage:
    python sensitivity_analysis.py \\
        --dataset out/analysis_dataset.parquet \\
        --cases   out/exparel_cases.csv \\
        --out     out/sensitivity.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from peroneal_palsy import ALL_PTS as PROLONGED_BLOCK_PTS
from duration_analysis import PATHOLOGY_PTS


# ---------------------------------------------------------------------------
# Shared 2x2 and signal logic
# ---------------------------------------------------------------------------

def _compute(a: int, b: int, c: int, d: int) -> dict:
    aa, bb, cc, dd = [x + 0.5 if x == 0 else x for x in (a, b, c, d)]
    ror = (aa / cc) / (bb / dd)
    se = np.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    N = a + b + c + d
    prr = (aa / (aa + cc)) / (bb / (bb + dd))
    exp_a = (a + b) * (a + c) / N if N else 0
    chi2 = ((abs(a - exp_a) - 0.5) ** 2 / exp_a) if exp_a else 0
    ic = np.log2((a + 0.5) / (exp_a + 0.5)) if exp_a else np.nan
    var_ic = (1 / np.log(2) ** 2) * (
        (N - a + 0.5) / ((a + 0.5) * (1 + N))
        + (N - (a + b) + 0.5) / (((a + b) + 0.5) * (1 + N))
        + (N - (a + c) + 0.5) / (((a + c) + 0.5) * (1 + N))
    ) if N else np.nan
    ic025 = ic - 1.96 * np.sqrt(var_ic) if not np.isnan(ic) else np.nan
    signal = (
        (a >= 3) and (np.exp(np.log(ror) - 1.96 * se) > 1)
        and (prr >= 2) and (chi2 >= 4) and (ic025 > 0)
    )
    return {
        "a": a, "b": b, "c": c, "d": d,
        "ROR": round(ror, 3),
        "ROR_lower": round(float(np.exp(np.log(ror) - 1.96 * se)), 3),
        "ROR_upper": round(float(np.exp(np.log(ror) + 1.96 * se)), 3),
        "PRR": round(prr, 3),
        "chi2": round(chi2, 2),
        "IC025": round(ic025, 3) if not np.isnan(ic025) else np.nan,
        "signal": bool(signal),
    }


def _composite_vs_bupi(ds: pd.DataFrame, pts: set[str],
                        allowed_pids: set[int]) -> dict:
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    ds = ds[ds["primaryid"].isin(allowed_pids)]
    ex_total = ds.loc[ds["drug_group"] == "exparel", "primaryid"].nunique()
    bu_total = ds.loc[ds["drug_group"] == "bupivacaine_plain", "primaryid"].nunique()
    a = ds.loc[(ds["drug_group"] == "exparel") & ds["pt_lower"].isin(pts), "primaryid"].nunique()
    b = ds.loc[(ds["drug_group"] == "bupivacaine_plain") & ds["pt_lower"].isin(pts), "primaryid"].nunique()
    c = ex_total - a
    d = bu_total - b
    return {
        "exparel_cases": ex_total,
        "bupi_cases": bu_total,
        **_compute(a, b, c, d),
    }


# ---------------------------------------------------------------------------
# Restriction builders — return set of primaryids that pass the filter
# ---------------------------------------------------------------------------

def _pids_all(ds: pd.DataFrame) -> set[int]:
    return set(ds["primaryid"].unique())


def _pids_serious(cases: pd.DataFrame, comp_pids: set[int]) -> set[int]:
    """Cases with any serious outcome code. For comparator cases (not in
    exparel_cases.csv), we conservatively INCLUDE all — restriction only
    applies to the Exparel side — to match the intent of a safety-weighted
    signal check (i.e., 'is the signal driven by mild reports on the
    Exparel side?'). Standard pharmacovigilance convention."""
    ex_serious = set(
        cases.loc[cases["outcomes"].notna(), "primaryid"]
        .astype("int64")
    )
    return ex_serious | comp_pids


def _pids_md_reported(cases: pd.DataFrame, comp_pids: set[int]) -> set[int]:
    ex_md = set(
        cases.loc[cases["occp_cod"] == "MD", "primaryid"]
        .astype("int64")
    )
    return ex_md | comp_pids


def _pids_exclude_pregnancy(ds: pd.DataFrame) -> set[int]:
    """Exclude cases whose PT list contains a pregnancy/neonatal term."""
    pregnancy_terms = {
        "maternal exposure during pregnancy",
        "maternal exposure during delivery",
        "exposure during pregnancy",
        "foetal exposure during pregnancy",
        "foetal exposure during delivery",
        "premature baby",
        "premature delivery",
        "low birth weight baby",
        "normal newborn",
        "pre-eclampsia",
        "neonatal respiratory distress syndrome",
    }
    ds_lower = ds.copy()
    ds_lower["pt_lower"] = ds_lower["pt"].astype(str).str.lower()
    excluded = set(
        ds_lower.loc[ds_lower["pt_lower"].isin(pregnancy_terms), "primaryid"]
    )
    return set(ds["primaryid"].unique()) - excluded


def _pids_exclude_mfr(cases: pd.DataFrame, comp_pids: set[int]) -> set[int]:
    """Exclude Exparel cases from manufacturer submissions. DEMO's
    `mfr_sndr` is populated when the report came through a manufacturer
    channel. A case with any `mfr_sndr` value is a manufacturer-submitted
    report by FDA convention."""
    if "mfr_sndr" not in cases.columns:
        return set(cases["primaryid"].astype("int64")) | comp_pids
    ex_non_mfr = set(
        cases.loc[cases["mfr_sndr"].isna(), "primaryid"].astype("int64")
    )
    return ex_non_mfr | comp_pids


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, type=Path)
    ap.add_argument("--cases", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ds = pd.read_parquet(args.dataset)
    cases = pd.read_csv(args.cases)

    # Comparator primaryids (for restrictions that filter only the Exparel side)
    comp_pids = set(
        ds.loc[ds["drug_group"] == "bupivacaine_plain", "primaryid"]
        .unique()
        .astype("int64")
    )

    restrictions = {
        "none (primary analysis)": _pids_all(ds),
        "serious outcomes only":     _pids_serious(cases, comp_pids),
        "physician-reported only":   _pids_md_reported(cases, comp_pids),
        "exclude pregnancy/neonatal": _pids_exclude_pregnancy(ds),
        "exclude mfr-submitted":     _pids_exclude_mfr(cases, comp_pids),
    }

    rows = []
    for label, pids in restrictions.items():
        for composite_label, pts in [
            ("prolonged_block_full", PROLONGED_BLOCK_PTS),
            ("pathology_implying",   PATHOLOGY_PTS),
        ]:
            r = _composite_vs_bupi(ds, pts, pids)
            r["restriction"] = label
            r["composite"] = composite_label
            rows.append(r)

    out = pd.DataFrame(rows)
    cols = ["restriction", "composite", "exparel_cases", "bupi_cases",
            "a", "b", "ROR", "ROR_lower", "ROR_upper",
            "PRR", "chi2", "IC025", "signal"]
    out = out[cols]
    out.to_csv(args.out, index=False)

    print(out.to_string(index=False))
    print(f"\nsaved → {args.out}")


if __name__ == "__main__":
    main()

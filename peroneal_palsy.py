"""
Peroneal palsy / prolonged-block subgroup deep-dive
===================================================
The pilot flagged PERONEAL NERVE PALSY (a=40, b=29, ROR 2.4 vs bupi). The
case narrative from the 2012-2014 partial data found a single "Neuromuscular
block prolonged + Off label use" case that matches the phenotype. This
script takes the full ASCII pipeline output and produces the figure-1
table for the manuscript:

  - PT-level table: peroneal / prolonged-block / neurotoxicity PTs
    with a, b, ROR, ROR_lower, PRR, IC025 vs both comparators
  - Era stratification (pre-2018 / 2018-2020 / 2021+) to test whether
    the signal strengthens after the FDA's interscalene block indication
  - Case-level listing of every Exparel case with a matching PT, with
    demographics and indication for narrative review

Usage:
    python peroneal_palsy.py \\
        --dataset out/analysis_dataset.parquet \\
        --cases   out/exparel_cases.csv \\
        --out-dir out/peroneal/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Pre-specified PTs for this subgroup. Grouped by sub-phenotype so the
# manuscript table can show which part of the spectrum is driving signal.
PT_GROUPS = {
    "peroneal_specific": {
        "peroneal nerve palsy",
        "peroneal nerve injury",
    },
    "prolonged_block": {
        "neuromuscular block prolonged",
        "hypoaesthesia",
        "paraesthesia",
        "muscular weakness",
        "hypoaesthesia oral",
        "hemiparesis",
        "paresis",
        "motor dysfunction",
        "sensory loss",
        "anaesthesia",
    },
    "peripheral_neuropathy": {
        "neuropathy peripheral",
        "peripheral motor neuropathy",
        "peripheral sensory neuropathy",
        "peripheral nerve injury",
        "nerve injury",
        "neuralgia",
        "brachial plexopathy",
    },
}

ALL_PTS = set().union(*PT_GROUPS.values())


def compute_2x2(a: int, b: int, c: int, d: int) -> dict:
    aa, bb, cc, dd = [x + 0.5 if x == 0 else x for x in (a, b, c, d)]
    ror = (aa / cc) / (bb / dd)
    se = np.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    prr = (aa / (aa + cc)) / (bb / (bb + dd))
    N = a + b + c + d
    exp_a = (a + b) * (a + c) / N if N else 0
    ic = np.log2((a + 0.5) / (exp_a + 0.5)) if exp_a else np.nan
    var_ic = (1 / np.log(2) ** 2) * (
        (N - a + 0.5) / ((a + 0.5) * (1 + N))
        + (N - (a + b) + 0.5) / (((a + b) + 0.5) * (1 + N))
        + (N - (a + c) + 0.5) / (((a + c) + 0.5) * (1 + N))
    ) if N else np.nan
    return {
        "a": a, "b": b, "c": c, "d": d,
        "ROR": ror,
        "ROR_lower": float(np.exp(np.log(ror) - 1.96 * se)),
        "ROR_upper": float(np.exp(np.log(ror) + 1.96 * se)),
        "PRR": prr,
        "IC025": ic - 1.96 * np.sqrt(var_ic) if not np.isnan(ic) else np.nan,
    }


def pt_level_table(ds: pd.DataFrame, pts: set[str], comparator: str) -> pd.DataFrame:
    """One row per PT in `pts`."""
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    target_cases = ds.loc[ds["drug_group"] == "exparel", "primaryid"].nunique()
    comp_cases = ds.loc[ds["drug_group"] == comparator, "primaryid"].nunique()

    rows = []
    for pt in sorted(pts):
        a = ds.loc[(ds["drug_group"] == "exparel") & (ds["pt_lower"] == pt),
                   "primaryid"].nunique()
        b = ds.loc[(ds["drug_group"] == comparator) & (ds["pt_lower"] == pt),
                   "primaryid"].nunique()
        if a + b == 0:
            continue  # PT never reported for either drug
        c = target_cases - a
        d = comp_cases - b
        rows.append({"pt": pt, **compute_2x2(a, b, c, d)})

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # Attach the PT group label
    pt_to_group = {pt: g for g, pts_set in PT_GROUPS.items() for pt in pts_set}
    df["group"] = df["pt"].map(pt_to_group)
    return df.sort_values(["group", "ROR"], ascending=[True, False])


def composite_by_era(
    ds: pd.DataFrame,
    cases: pd.DataFrame,
    pts: set[str],
    comparator: str,
) -> pd.DataFrame:
    """Composite subgroup signal vs comparator, stratified by era with
    ERA-MATCHED COMPARATOR POOLS. Earlier version left comparators
    unfiltered (full-period background); this version filters both sides
    using out/case_meta.parquet (17M FAERS primaryid → event_year map)."""
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    ds["hit"] = ds["pt_lower"].isin(pts)

    case_meta = pd.read_parquet("out/case_meta.parquet")
    year_map = case_meta.set_index("primaryid")["event_year"]

    era_bins = {
        "all":        (0, 9999),
        "pre_2018":   (0, 2017),
        "2018_2020":  (2018, 2020),
        "2021_plus":  (2021, 9999),
    }
    rows = []
    for era, (lo, hi) in era_bins.items():
        if era == "all":
            ds_era = ds
        else:
            era_pids = set(year_map[(year_map >= lo) & (year_map <= hi)].index)
            ds_era = ds[ds["primaryid"].isin(era_pids)]
        target_cases = ds_era.loc[ds_era["drug_group"] == "exparel", "primaryid"].nunique()
        comp_cases = ds_era.loc[ds_era["drug_group"] == comparator, "primaryid"].nunique()
        a = ds_era.loc[(ds_era["drug_group"] == "exparel") & ds_era["hit"], "primaryid"].nunique()
        b = ds_era.loc[(ds_era["drug_group"] == comparator) & ds_era["hit"], "primaryid"].nunique()
        c = target_cases - a
        d = comp_cases - b
        rows.append({
            "era": era, "comparator": comparator,
            "target_cases": target_cases, "comparator_cases": comp_cases,
            **compute_2x2(a, b, c, d),
        })
    return pd.DataFrame(rows)


def matching_exparel_cases(cases: pd.DataFrame, ds: pd.DataFrame, pts: set[str]) -> pd.DataFrame:
    """Case-level listing of Exparel reports with ≥1 PT in the subgroup."""
    ds = ds.copy()
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    hit_pids = ds.loc[
        (ds["drug_group"] == "exparel") & ds["pt_lower"].isin(pts),
        "primaryid",
    ].unique()
    return cases[cases["primaryid"].isin(hit_pids)].copy()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, type=Path)
    ap.add_argument("--cases", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(exist_ok=True, parents=True)

    ds = pd.read_parquet(args.dataset)
    cases = pd.read_csv(args.cases)

    # PT-level tables vs each comparator
    for comp in ["bupivacaine_plain", "ropivacaine"]:
        tbl = pt_level_table(ds, ALL_PTS, comp)
        out_path = args.out_dir / f"pt_level_vs_{comp}.csv"
        tbl.to_csv(out_path, index=False)
        print(f"\n=== PT-level vs {comp} ===")
        if tbl.empty:
            print("  (no matching PTs)")
        else:
            print(tbl[["group", "pt", "a", "b", "ROR", "ROR_lower", "PRR", "IC025"]].to_string(index=False))

    # Composite by era, per comparator
    era_rows = []
    for comp in ["bupivacaine_plain", "ropivacaine"]:
        era_rows.append(composite_by_era(ds, cases, ALL_PTS, comp))
    era_all = pd.concat(era_rows, ignore_index=True)
    era_all.to_csv(args.out_dir / "composite_by_era.csv", index=False)
    print("\n=== Composite by era ===")
    print(era_all[["era", "comparator", "target_cases", "a", "b", "ROR", "ROR_lower", "PRR", "IC025"]].to_string(index=False))

    # Per-subgroup (peroneal_specific / prolonged_block / peripheral_neuropathy)
    sub_rows = []
    for group, pts in PT_GROUPS.items():
        for comp in ["bupivacaine_plain", "ropivacaine"]:
            era_df = composite_by_era(ds, cases, pts, comp)
            era_df["subgroup"] = group
            sub_rows.append(era_df)
    sub_all = pd.concat(sub_rows, ignore_index=True)
    sub_all.to_csv(args.out_dir / "composite_by_subgroup_era.csv", index=False)

    # Matching Exparel cases (for narrative review)
    hits = matching_exparel_cases(cases, ds, ALL_PTS)
    hits.to_csv(args.out_dir / "matching_exparel_cases.csv", index=False)
    print(f"\n=== Matching Exparel cases for narrative review: {len(hits)} ===")
    cols = [c for c in ["primaryid", "event_year", "age", "gndr_cod",
                        "occp_cod", "indications", "reactions"] if c in hits.columns]
    print(hits[cols].head(20).to_string(index=False))

    print(f"\nAll outputs → {args.out_dir}/")


if __name__ == "__main__":
    main()

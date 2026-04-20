"""
FAERS Exparel Disproportionality Analysis — Production Pipeline
================================================================

This is the PRIMARY analysis script. It ingests FAERS quarterly ASCII
files directly from FDA, enabling full control over deduplication,
drug-name harmonization, and sensitivity analyses that the openFDA
API cannot provide.

Download FAERS quarterly ASCII files from:
    https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html

Expected directory layout:
    data/raw/
        faers_ascii_2012q2/ASCII/DEMO12Q2.txt
        faers_ascii_2012q2/ASCII/DRUG12Q2.txt
        faers_ascii_2012q2/ASCII/REAC12Q2.txt
        ... (for every quarter 2012Q2 through most recent)

Usage:
    python ascii_pipeline.py --data-dir data/raw --out out/
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Drug-name harmonization (regex patterns)
# ---------------------------------------------------------------------------

LIPOSOMAL_PATTERN = re.compile(
    r"\b(exparel|"
    r"bupivacaine\s+liposom(e|al)(\s+injectable\s+suspension)?|"
    r"liposomal\s+bupivacaine|"
    r"bupivacaine\s+liposom)\b",
    re.IGNORECASE,
)

PLAIN_BUPI_PATTERN = re.compile(
    r"\b(bupivacaine(\s+hydrochloride|\s+hcl)?|marcaine|sensorcaine)\b",
    re.IGNORECASE,
)

ROPI_PATTERN = re.compile(
    r"\b(ropivacaine(\s+hydrochloride|\s+hcl)?|naropin)\b",
    re.IGNORECASE,
)


def classify_drug(name: str) -> str | None:
    """Classify a DRUGNAME / PROD_AI string into one of our drug groups."""
    if not isinstance(name, str):
        return None
    # Liposomal first — if matched, this takes precedence
    if LIPOSOMAL_PATTERN.search(name):
        return "exparel"
    if PLAIN_BUPI_PATTERN.search(name):
        return "bupivacaine_plain"
    if ROPI_PATTERN.search(name):
        return "ropivacaine"
    return None


# ---------------------------------------------------------------------------
# FAERS ASCII ingestion
# ---------------------------------------------------------------------------

def load_faers_quarter(quarter_dir: Path) -> dict[str, pd.DataFrame]:
    """
    Load the seven ASCII tables from a single FAERS quarter directory.
    FAERS uses '$' as the delimiter and Latin-1 encoding.
    """
    # FDA alternates between ASCII/ (older) and ascii/ (newer) nesting,
    # and a few older quarters put files at the top level.
    for sub in ("ASCII", "ascii"):
        if (quarter_dir / sub).exists():
            ascii_dir = quarter_dir / sub
            break
    else:
        ascii_dir = quarter_dir
    tables = {}
    all_txt = list(ascii_dir.glob("*.txt")) + list(ascii_dir.glob("*.TXT"))
    for prefix in ["DEMO", "DRUG", "REAC", "OUTC", "RPSR", "THER", "INDI"]:
        # Case-insensitive match — FDA alternates between DEMO13Q1.txt
        # (2013+) and demo12q4.txt (2012) without a consistent pattern.
        candidates = [p for p in all_txt if p.stem.upper().startswith(prefix)]
        if not candidates:
            print(f"  WARN: no {prefix} file in {ascii_dir}")
            continue
        tables[prefix] = pd.read_csv(
            candidates[0],
            sep="$",
            encoding="latin-1",
            low_memory=False,
            on_bad_lines="skip",
        )
        # Some older quarters (e.g. 2012 Q4) are UTF-8-with-BOM; reading
        # as Latin-1 leaves the BOM bytes as "ï»¿" at the start of the
        # first column name. Strip BOM artefacts before lowercasing.
        tables[prefix].columns = [
            c.lstrip("\ufeff").lstrip("ï»¿").lower()
            for c in tables[prefix].columns
        ]
    return tables


def load_all_quarters(data_dir: Path) -> dict[str, pd.DataFrame]:
    """Concatenate all quarters into master FAERS tables."""
    # INDI kept for indication/procedure stratification.
    # THER kept for drug start/stop dates (time-to-onset).
    # RPSR kept for report-source analysis (e.g. excluding lit reports).
    all_tables: dict[str, list[pd.DataFrame]] = {
        k: [] for k in ["DEMO", "DRUG", "REAC", "OUTC", "INDI", "THER", "RPSR"]
    }
    quarter_dirs = sorted(
        [p for p in data_dir.iterdir() if p.is_dir() and "faers" in p.name.lower()]
    )
    print(f"Found {len(quarter_dirs)} quarter directories")
    for qdir in quarter_dirs:
        print(f"  loading {qdir.name}…")
        tables = load_faers_quarter(qdir)
        for k in all_tables:
            if k in tables:
                all_tables[k].append(tables[k])
    return {
        k: pd.concat(v, ignore_index=True) if v else pd.DataFrame()
        for k, v in all_tables.items()
    }


# ---------------------------------------------------------------------------
# Deduplication per FDA protocol
# ---------------------------------------------------------------------------

def deduplicate_demo(demo: pd.DataFrame) -> pd.DataFrame:
    """
    Per FDA guidance:
      1. For identical caseids, keep the record with most recent fda_dt.
      2. Tie-break on highest primaryid.
    """
    before = len(demo)
    demo = demo.sort_values(["caseid", "fda_dt", "primaryid"], ascending=[True, False, False])
    demo = demo.drop_duplicates(subset=["caseid"], keep="first")
    print(f"  dedup DEMO: {before:,} → {len(demo):,}")
    return demo


# ---------------------------------------------------------------------------
# Build analysis dataset
# ---------------------------------------------------------------------------

def build_analysis_dataset(all_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Build a case-level dataset where each row = (primaryid, drug_group, PT).
    Restrict to drug_group in {exparel, bupivacaine_plain, ropivacaine}
    AND role_cod == 'PS' (primary suspect).
    """
    demo = deduplicate_demo(all_tables["DEMO"])
    drug = all_tables["DRUG"]
    reac = all_tables["REAC"]

    # Restrict DRUG to PS role
    drug_ps = drug[drug["role_cod"] == "PS"].copy()

    # Classify drugs. CRITICAL: we must check drugname FIRST for the
    # liposomal pattern. FDA's prod_ai field normalizes Exparel to just
    # "BUPIVACAINE" (active molecule, not formulation), so using prod_ai
    # as primary would silently re-label every Exparel report as plain
    # bupivacaine from 2014 Q3 onward — and contaminate the comparator.
    # drugname preserves the brand ("EXPAREL", "BUPIVACAINE LIPOSOME"),
    # which is the only place formulation-specific info lives.
    def classify_row(drugname, prod_ai):
        # Most-specific label wins: exparel > ropivacaine > bupivacaine_plain.
        for field in (drugname, prod_ai):
            if LIPOSOMAL_PATTERN.search(str(field)):
                return "exparel"
        for field in (drugname, prod_ai):
            if ROPI_PATTERN.search(str(field)):
                return "ropivacaine"
        for field in (drugname, prod_ai):
            if PLAIN_BUPI_PATTERN.search(str(field)):
                return "bupivacaine_plain"
        return None

    if "prod_ai" in drug_ps.columns:
        drug_ps["drug_group"] = [
            classify_row(dn, pa)
            for dn, pa in zip(drug_ps["drugname"], drug_ps["prod_ai"])
        ]
    else:
        drug_ps["drug_group"] = drug_ps["drugname"].apply(classify_drug)
    drug_ps["drug_string"] = drug_ps["drugname"]  # kept for exparel_cases export
    drug_ps = drug_ps.dropna(subset=["drug_group"])

    print(f"  drugs (PS, in scope): {len(drug_ps):,}")
    print(drug_ps["drug_group"].value_counts())

    # Merge to demo-dedup
    drug_ps = drug_ps.merge(demo[["primaryid", "caseid"]], on="primaryid", how="inner")

    # Merge to reactions
    ds = drug_ps[["primaryid", "drug_group"]].merge(
        reac[["primaryid", "pt"]], on="primaryid", how="inner"
    ).drop_duplicates()

    # A case may list the same PT in multiple DRUG lines — collapse to unique
    # (primaryid, drug_group, pt) triples
    ds = ds.drop_duplicates()
    print(f"  final analysis rows: {len(ds):,}")
    return ds


# ---------------------------------------------------------------------------
# Case-level export (for narrative review of Exparel reports)
# ---------------------------------------------------------------------------

def _agg_unique(series: pd.Series) -> str:
    vals = sorted({str(v) for v in series.dropna() if str(v).strip()})
    return "; ".join(vals)


def export_exparel_cases(
    all_tables: dict[str, pd.DataFrame],
    ds: pd.DataFrame,
    out_path: Path,
) -> pd.DataFrame:
    """
    One row per Exparel case with deduplicated demographics, PTs, indications,
    outcomes, and concomitant drugs. Used for (a) clinical narrative review,
    (b) subgroup filtering (e.g. nerve-block cases only), (c) reviewer
    requests to "show me the actual reports."
    """
    demo = deduplicate_demo(all_tables["DEMO"])
    ex_pids = set(ds.loc[ds["drug_group"] == "exparel", "primaryid"])
    ex_demo = demo[demo["primaryid"].isin(ex_pids)].copy()

    reac = all_tables["REAC"]
    reac_by_case = (
        reac[reac["primaryid"].isin(ex_pids)]
        .groupby("primaryid")["pt"].agg(_agg_unique).rename("reactions")
    )

    indi = all_tables.get("INDI", pd.DataFrame())
    if not indi.empty and "indi_pt" in indi.columns:
        indi_by_case = (
            indi[indi["primaryid"].isin(ex_pids)]
            .groupby("primaryid")["indi_pt"].agg(_agg_unique).rename("indications")
        )
    else:
        indi_by_case = pd.Series(dtype=str, name="indications")

    outc = all_tables.get("OUTC", pd.DataFrame())
    if not outc.empty and "outc_cod" in outc.columns:
        outc_by_case = (
            outc[outc["primaryid"].isin(ex_pids)]
            .groupby("primaryid")["outc_cod"].agg(_agg_unique).rename("outcomes")
        )
    else:
        outc_by_case = pd.Series(dtype=str, name="outcomes")

    drug = all_tables["DRUG"]
    ex_drug = drug[drug["primaryid"].isin(ex_pids)].copy()
    # is_exparel must match the LIPOSOMAL_PATTERN specifically — prod_ai
    # strips the liposomal formulation tag and collapses to BUPIVACAINE,
    # so we check drugname first (same precedence as the main classifier).
    def _is_exparel_row(drugname, prod_ai):
        for f in (drugname, prod_ai) if "prod_ai" in ex_drug.columns else (drugname,):
            if LIPOSOMAL_PATTERN.search(str(f)):
                return True
        return False
    if "prod_ai" in ex_drug.columns:
        ex_drug["is_exparel"] = [
            _is_exparel_row(dn, pa)
            for dn, pa in zip(ex_drug["drugname"], ex_drug["prod_ai"])
        ]
    else:
        ex_drug["is_exparel"] = ex_drug["drugname"].apply(
            lambda x: bool(LIPOSOMAL_PATTERN.search(str(x)))
        )

    exparel_row = (
        ex_drug[ex_drug["is_exparel"]]
        .sort_values("primaryid")
        .drop_duplicates("primaryid")
    )
    keep_drug_cols = [c for c in
        ["primaryid", "drugname", "role_cod", "route", "dose_vbm"]
        if c in exparel_row.columns]
    exparel_row = exparel_row[keep_drug_cols]

    concomitant = (
        ex_drug[~ex_drug["is_exparel"]]
        .groupby("primaryid")["drugname"].agg(_agg_unique).rename("concomitant_drugs")
    )

    out = (
        ex_demo
        .merge(exparel_row, on="primaryid", how="left", suffixes=("", "_drug"))
        .merge(reac_by_case, on="primaryid", how="left")
        .merge(indi_by_case, on="primaryid", how="left")
        .merge(outc_by_case, on="primaryid", how="left")
        .merge(concomitant, on="primaryid", how="left")
    )

    # Add year for era stratification. event_dt is only ~28% populated in
    # FAERS; fall back to init_fda_dt (date FDA first received the report),
    # which is 100% populated and is a faithful proxy for the era even if
    # not the exact event date. Final fallback to fda_dt.
    def _year_with_fallback(row):
        for col in ("event_dt", "init_fda_dt", "fda_dt"):
            v = row.get(col)
            if pd.notna(v):
                y = pd.to_datetime(v, format="%Y%m%d", errors="coerce")
                if pd.notna(y):
                    return y.year
        return pd.NA
    out["event_year"] = out.apply(_year_with_fallback, axis=1)

    out.to_csv(out_path, index=False)
    print(f"  exported {len(out):,} Exparel cases → {out_path.name}")
    return out


# ---------------------------------------------------------------------------
# Disproportionality (same math as pilot, vectorized)
# ---------------------------------------------------------------------------

@dataclass
class ContingencyFrame:
    """a/b/c/d for every PT in the dataset, vectorized."""
    df: pd.DataFrame  # columns: pt, a, b, c, d


def build_contingency(
    ds: pd.DataFrame,
    target: str = "exparel",
    comparator: str = "bupivacaine_plain",
) -> ContingencyFrame:
    """
    For each PT: count primaryids reporting it under each drug group.
    a = target+PT, b = comparator+PT, c = target+not-PT, d = comparator+not-PT.
    """
    target_cases = ds.loc[ds["drug_group"] == target, "primaryid"].nunique()
    comp_cases = ds.loc[ds["drug_group"] == comparator, "primaryid"].nunique()

    # Count cases per (drug_group, pt)
    counts = (
        ds[ds["drug_group"].isin([target, comparator])]
        .groupby(["drug_group", "pt"])["primaryid"]
        .nunique()
        .unstack("drug_group", fill_value=0)
        .reset_index()
    )
    counts["a"] = counts.get(target, 0)
    counts["b"] = counts.get(comparator, 0)
    counts["c"] = target_cases - counts["a"]
    counts["d"] = comp_cases - counts["b"]
    return ContingencyFrame(df=counts[["pt", "a", "b", "c", "d"]])


def compute_signals_vectorized(cf: ContingencyFrame) -> pd.DataFrame:
    df = cf.df.copy()
    # Haldane-Anscombe: add 0.5 to zero cells
    for col in ["a", "b", "c", "d"]:
        df[f"{col}_adj"] = np.where(df[col] == 0, 0.5, df[col])

    a, b, c, d = df["a_adj"], df["b_adj"], df["c_adj"], df["d_adj"]

    # ROR
    df["ROR"] = (a / c) / (b / d)
    se_ln_ror = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    df["ROR_lower"] = np.exp(np.log(df["ROR"]) - 1.96 * se_ln_ror)
    df["ROR_upper"] = np.exp(np.log(df["ROR"]) + 1.96 * se_ln_ror)

    # PRR
    df["PRR"] = (a / (a + c)) / (b / (b + d))

    # Chi-sq (Yates)
    N = a + b + c + d
    expected_a = (a + b) * (a + c) / N
    expected_b = (a + b) * (b + d) / N
    expected_c = (c + d) * (a + c) / N
    expected_d = (c + d) * (b + d) / N
    # Yates continuity correction
    df["PRR_chi2"] = (
        (np.abs(a - expected_a) - 0.5) ** 2 / expected_a
        + (np.abs(b - expected_b) - 0.5) ** 2 / expected_b
        + (np.abs(c - expected_c) - 0.5) ** 2 / expected_c
        + (np.abs(d - expected_d) - 0.5) ** 2 / expected_d
    )

    # IC (BCPNN)
    df["IC"] = np.log2((a + 0.5) / (expected_a + 0.5))
    var_ic = (1 / np.log(2) ** 2) * (
        (N - a + 0.5) / ((a + 0.5) * (1 + N))
        + (N - (a + b) + 0.5) / (((a + b) + 0.5) * (1 + N))
        + (N - (a + c) + 0.5) / (((a + c) + 0.5) * (1 + N))
    )
    df["IC025"] = df["IC"] - 1.96 * np.sqrt(var_ic)

    # Signal flag: all four criteria met
    df["signal"] = (
        (df["a"] >= 3)
        & (df["ROR_lower"] > 1)
        & (df["PRR"] >= 2)
        & (df["PRR_chi2"] >= 4)
        & (df["IC025"] > 0)
    )
    return df.drop(columns=[c for c in df.columns if c.endswith("_adj")])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--out", default=Path("out"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)

    print("=" * 70)
    print("FAERS Exparel Disproportionality Analysis — Production Pipeline")
    print("=" * 70)

    print("\n[1/5] Loading FAERS quarterly ASCII files…")
    all_tables = load_all_quarters(args.data_dir)
    for k, v in all_tables.items():
        print(f"  {k}: {len(v):,} rows")

    print("\n[2/5] Building case-level analysis dataset…")
    ds = build_analysis_dataset(all_tables)
    ds.to_parquet(args.out / "analysis_dataset.parquet")

    print("\n[3/5] Exporting Exparel case-level extract…")
    export_exparel_cases(all_tables, ds, args.out / "exparel_cases.csv")

    print("\n[4/5] Primary analysis: Exparel vs plain bupivacaine…")
    cf = build_contingency(ds, "exparel", "bupivacaine_plain")
    signals = compute_signals_vectorized(cf)
    signals = signals.sort_values(
        ["signal", "ROR"], ascending=[False, False]
    )
    signals.to_csv(args.out / "primary_signals.csv", index=False)
    print(f"  total PTs analyzed: {len(signals):,}")
    print(f"  signals meeting all 4 criteria: {signals['signal'].sum():,}")

    print("\n[5/5] Secondary analysis: Exparel vs ropivacaine…")
    cf2 = build_contingency(ds, "exparel", "ropivacaine")
    signals2 = compute_signals_vectorized(cf2)
    signals2 = signals2.sort_values(
        ["signal", "ROR"], ascending=[False, False]
    )
    signals2.to_csv(args.out / "secondary_ropivacaine_signals.csv", index=False)

    print(f"\nDone. Outputs in {args.out}/")


if __name__ == "__main__":
    main()

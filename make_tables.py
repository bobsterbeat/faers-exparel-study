"""
Generate publication-ready tables as CSV + rendered PNG.

Tables:
  Table 1. Cohort demographics by drug group
  Table 2. PT-level signals (prolonged-block + peripheral-neuropathy spectrum)
  Table 3. Composite endpoint by era, vs both comparators

Usage:
    python make_tables.py --out out/tables/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _render_table(df: pd.DataFrame, out_path: Path, title: str,
                  col_widths: list[float] | None = None,
                  row_height: float = 0.35) -> None:
    """Render a dataframe as a styled PNG table."""
    rows, cols = df.shape
    fig_h = row_height * (rows + 2) + 0.6
    fig_w = max(8, cols * 1.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=12, pad=14, fontweight="bold")

    cell_text = df.astype(str).values.tolist()
    tbl = ax.table(
        cellText=cell_text,
        colLabels=df.columns.tolist(),
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=col_widths,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.35)

    # Style header
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", fontweight="bold")
            cell.set_edgecolor("#2c3e50")
        else:
            cell.set_edgecolor("#bbb")
            if r % 2 == 0:
                cell.set_facecolor("#f8f8f8")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Table 1: Cohort demographics
# ---------------------------------------------------------------------------

def table_1_cohort(out_dir: Path) -> pd.DataFrame:
    cases = pd.read_csv("out/exparel_cases.csv")
    ds = pd.read_parquet("out/analysis_dataset.parquet")

    def _group_stats(pids: set, label: str) -> dict:
        sub = cases[cases["primaryid"].isin(pids)] if label == "Exparel" else None
        return {"n": len(pids), "label": label, "cases": sub}

    totals = {
        "Exparel": ds.loc[ds["drug_group"] == "exparel", "primaryid"].nunique(),
        "Plain bupivacaine": ds.loc[ds["drug_group"] == "bupivacaine_plain", "primaryid"].nunique(),
        "Ropivacaine": ds.loc[ds["drug_group"] == "ropivacaine", "primaryid"].nunique(),
    }

    # Exparel-specific demographics (we only have case-level export for Exparel)
    c = cases.copy()
    age_med = c["age"].median()
    age_n = int(c["age"].notna().sum())
    # FAERS renamed gndr_cod -> sex around 2014 Q3. Coalesce both columns so
    # we don't silently drop >90% of cases.
    if "sex" in c.columns:
        c["sex_coalesced"] = c["gndr_cod"].fillna(c["sex"])
    else:
        c["sex_coalesced"] = c["gndr_cod"]
    sex_counts = c["sex_coalesced"].value_counts()
    rpt_counts = c["occp_cod"].value_counts()

    outc = c["outcomes"].fillna("").str.split("; ").explode()
    any_serious = c["outcomes"].notna().sum()
    death = c["outcomes"].str.contains("DE", na=False).sum()
    lt = c["outcomes"].str.contains("LT", na=False).sum()
    ho = c["outcomes"].str.contains("HO", na=False).sum()

    # Era breakdown for Exparel
    era_n = {
        "Pre-2018 (2012 Q4 – 2017)": int(((c["event_year"] >= 2012) & (c["event_year"] <= 2017)).sum()),
        "2018 – 2020": int(((c["event_year"] >= 2018) & (c["event_year"] <= 2020)).sum()),
        "2021 – 2025": int(((c["event_year"] >= 2021) & (c["event_year"] <= 2025)).sum()),
    }

    rows = [
        ("Total primary-suspect cases", f"{totals['Exparel']:,}", f"{totals['Plain bupivacaine']:,}", f"{totals['Ropivacaine']:,}"),
        ("", "", "", ""),
        ("— Exparel cohort detail —", "", "", ""),
        ("Age, years, median (n with age)", f"{age_med:.0f} (n={age_n})", "", ""),
        ("Sex, F / M / NS", f"{sex_counts.get('F',0)} / {sex_counts.get('M',0)} / {sex_counts.get('NS',0)}", "", ""),
        ("Reporter, MD / PH / RN / Other", f"{rpt_counts.get('MD',0)} / {rpt_counts.get('PH',0)} / {rpt_counts.get('RN',0)} / {rpt_counts.drop(['MD','PH','RN'], errors='ignore').sum()}", "", ""),
        ("Any serious outcome", f"{any_serious:,} ({100*any_serious/len(c):.1f}%)", "", ""),
        ("  Death (DE)", f"{int(death)}", "", ""),
        ("  Life-threatening (LT)", f"{int(lt)}", "", ""),
        ("  Hospitalization (HO)", f"{int(ho)}", "", ""),
        ("Era: Pre-2018", f"{era_n['Pre-2018 (2012 Q4 – 2017)']}", "", ""),
        ("Era: 2018 – 2020", f"{era_n['2018 – 2020']}", "", ""),
        ("Era: 2021 – 2025", f"{era_n['2021 – 2025']}", "", ""),
    ]
    df = pd.DataFrame(rows, columns=["Characteristic", "Exparel", "Plain bupivacaine", "Ropivacaine"])

    df.to_csv(out_dir / "table_1_cohort.csv", index=False)
    _render_table(df, out_dir / "table_1_cohort.png",
                  "Table 1. Cohort demographics and outcome distribution")
    return df


# ---------------------------------------------------------------------------
# Table 2: PT-level signals
# ---------------------------------------------------------------------------

def table_2_pt_signals(out_dir: Path) -> pd.DataFrame:
    tbl_bupi = pd.read_csv("out/peroneal/pt_level_vs_bupivacaine_plain.csv")
    tbl_ropi = pd.read_csv("out/peroneal/pt_level_vs_ropivacaine.csv")

    tbl_bupi = tbl_bupi[(tbl_bupi["a"] >= 3)].copy()
    m = tbl_bupi.merge(
        tbl_ropi[["pt", "ROR", "ROR_lower", "ROR_upper", "IC025"]],
        on="pt", how="left", suffixes=("_bupi", "_ropi"),
    )
    m["signal_bupi"] = (
        (m["a"] >= 3) & (m["ROR_lower_bupi"] > 1)
        & (m["PRR"] >= 2) & (m["IC025_bupi"] > 0)
    )
    m = m.sort_values(["signal_bupi", "ROR_bupi"], ascending=[False, False])

    rows = []
    for _, r in m.iterrows():
        if pd.notna(r["ROR_ropi"]):
            ropi_str = f"{r['ROR_ropi']:.2f} ({r['ROR_lower_ropi']:.2f}–{r['ROR_upper_ropi']:.2f})"
        else:
            ropi_str = "—"
        rows.append({
            "PT": r["pt"].title(),
            "Group": r["group"].replace("_", " "),
            "a / b": f"{int(r['a'])} / {int(r['b'])}",
            "ROR vs bupivacaine (95% CI)":
                f"{r['ROR_bupi']:.2f} ({r['ROR_lower_bupi']:.2f}–{r['ROR_upper_bupi']:.2f})",
            "IC025 bupi": f"{r['IC025_bupi']:.2f}",
            "ROR vs ropivacaine (95% CI)": ropi_str,
            "Signal (bupi)": "✓" if r["signal_bupi"] else "",
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "table_2_pt_signals.csv", index=False)
    _render_table(df, out_dir / "table_2_pt_signals.png",
                  "Table 2. PT-level signals for prolonged-block and peripheral-neuropathy spectrum")
    return df


# ---------------------------------------------------------------------------
# Table 3: Composite endpoint by era
# ---------------------------------------------------------------------------

def table_3_composite(out_dir: Path) -> pd.DataFrame:
    per = pd.read_csv("out/peroneal/composite_by_era.csv")
    lst = pd.read_csv("out/last_composite.csv")
    cases = pd.read_csv("out/exparel_cases.csv")
    ds = pd.read_parquet("out/analysis_dataset.parquet")
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()

    from duration_analysis import PATHOLOGY_PTS, EXPECTED_PTS, SYSTEMIC_ABSORPTION_PTS
    from last_composite import LAST_PTS
    PROLONGED_BLOCK_PTS = PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS

    era_order = ["all", "pre_2018", "2018_2020", "2021_plus"]
    era_labels = {
        "all": "All time",
        "pre_2018": "Pre-2018",
        "2018_2020": "2018–2020",
        "2021_plus": "2021+",
    }
    era_bins = {
        "all": (0, 9999),
        "pre_2018": (0, 2017),
        "2018_2020": (2018, 2020),
        "2021_plus": (2021, 9999),
    }

    def _serious_rate(pts: set[str], era: str) -> str:
        """% of flagged Exparel cases in era with any serious outcome."""
        lo, hi = era_bins[era]
        hit_pids = ds.loc[
            (ds["drug_group"] == "exparel") & ds["pt_lower"].isin(pts),
            "primaryid",
        ].unique()
        sub = cases[cases["primaryid"].isin(hit_pids)]
        if era != "all":
            sub = sub[(sub["event_year"] >= lo) & (sub["event_year"] <= hi)]
        n = len(sub)
        if n == 0:
            return "—"
        n_ser = int(sub["outcomes"].notna().sum())
        return f"{n_ser}/{n} ({100*n_ser/n:.0f}%)"

    def _inverse_marker(ror: float, ror_lower: float, ror_upper: float) -> str:
        """Mark statistically below 1 (inverse) distinct from null."""
        if ror_upper < 1:
            return " ↓"  # entire CI below 1 — inverse signal
        return ""

    def _rows(df, composite_name, pts):
        out = []
        for era in era_order:
            rows_era = df[df["era"] == era]
            if rows_era.empty:
                continue
            bupi = rows_era[rows_era["comparator"] == "bupivacaine_plain"].iloc[0]
            ropi_row = rows_era[rows_era["comparator"] == "ropivacaine"]
            ropi = ropi_row.iloc[0] if len(ropi_row) else None
            signal_bupi = (
                (bupi["a"] >= 3) and (bupi["ROR_lower"] > 1)
                and (bupi["PRR"] >= 2) and (bupi["IC025"] > 0)
            )
            ror_bupi_str = (
                f"{bupi['ROR']:.2f} ({bupi['ROR_lower']:.2f}–{bupi['ROR_upper']:.2f})"
                + _inverse_marker(bupi["ROR"], bupi["ROR_lower"], bupi["ROR_upper"])
            )
            ror_ropi_str = (
                f"{ropi['ROR']:.2f} ({ropi['ROR_lower']:.2f}–{ropi['ROR_upper']:.2f})"
                + _inverse_marker(ropi["ROR"], ropi["ROR_lower"], ropi["ROR_upper"])
                if ropi is not None else "—"
            )
            out.append({
                "Composite": composite_name,
                "Era": era_labels[era],
                "Exparel cases": int(bupi.get("target_cases", bupi.get("target_total"))),
                "a / b": f"{int(bupi['a'])} / {int(bupi['b'])}",
                "ROR vs bupi (95% CI)": ror_bupi_str,
                "Signal (bupi)": "✓" if signal_bupi else "",
                "ROR vs ropi (95% CI)": ror_ropi_str,
                "Serious outcome (Exparel)": _serious_rate(pts, era),
            })
        return out

    all_rows = _rows(per, "Prolonged sensory/motor block", PROLONGED_BLOCK_PTS)
    all_rows += _rows(lst, "LAST spectrum", LAST_PTS)
    df = pd.DataFrame(all_rows)
    df.to_csv(out_dir / "table_3_composite_era.csv", index=False)
    _render_table(df, out_dir / "table_3_composite_era.png",
                  "Table 3. Composite-endpoint disproportionality by reporting era "
                  "(↓ = entire 95% CI below 1, i.e. inverse signal)")
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=Path("out/tables"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)

    t1 = table_1_cohort(args.out)
    print(f"\nTable 1 — Cohort ({len(t1)} rows)")
    print(t1.to_string(index=False))

    t2 = table_2_pt_signals(args.out)
    print(f"\nTable 2 — PT signals ({len(t2)} rows)")
    print(t2.head(10).to_string(index=False))

    t3 = table_3_composite(args.out)
    print(f"\nTable 3 — Composite by era ({len(t3)} rows)")
    print(t3.to_string(index=False))

    print(f"\nAll tables → {args.out}/")


if __name__ == "__main__":
    main()

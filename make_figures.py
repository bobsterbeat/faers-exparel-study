"""
Manuscript figures for the Exparel disproportionality study.
Produces:
  - fig1_composite_forest.png: 2x2 era-stratified forest plots
      top: peroneal / prolonged-block composite (the main signal)
      bottom: LAST-spectrum composite (the Weber-effect contrast)
  - fig2_pt_forest.png: PT-level forest for individual signals
  - fig3_era_trend.png: Exparel cases per year with composite rates over time

Usage:
    python make_figures.py --out fig/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _forest(ax, df, label_col, title, xlim=(0.2, 10), sig_col="signal",
            err_cols=("ROR", "ROR_lower", "ROR_upper")):
    """Horizontal forest plot with log-x, reference at ROR=1."""
    df = df.reset_index(drop=True)
    ror, lo, hi = err_cols
    y = np.arange(len(df))
    # Error bars as asymmetric from point estimate
    err = np.array([df[ror] - df[lo], df[hi] - df[ror]])
    # Color by signal if column present, else uniform
    colors = ["#c0392b" if (sig_col in df.columns and df.loc[i, sig_col]) else "#7f8c8d"
              for i in range(len(df))]
    for i in range(len(df)):
        ax.errorbar(df.loc[i, ror], y[i],
                    xerr=[[err[0][i]], [err[1][i]]],
                    fmt="o", color=colors[i], capsize=4, markersize=7, lw=1.4)
    ax.axvline(1.0, ls="--", color="#888", lw=0.8, alpha=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(df[label_col])
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_xlabel("ROR (95% CI), log scale")
    ax.set_title(title, fontsize=11, loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    # Counts next to each row
    for i in range(len(df)):
        ax.text(xlim[1] * 0.92, y[i],
                f"{int(df.loc[i,'a'])}/{int(df.loc[i,'b'])}",
                va="center", ha="right", fontsize=8, color="#555")


def fig1_composite(out_dir: Path) -> None:
    peroneal = pd.read_csv("out/peroneal/composite_by_era.csv")
    last = pd.read_csv("out/last_composite.csv")

    # last_composite.csv has an extra "tier" column and different ordering;
    # filter to the era rows we want and sort by era order.
    era_order = ["all", "pre_2018", "2018_2020", "2021_plus"]
    era_labels = {
        "all": "All time",
        "pre_2018": "Pre-2018",
        "2018_2020": "2018-2020",
        "2021_plus": "2021+",
    }

    def _prep(df, comparator, add_signal=True):
        d = df[df["comparator"] == comparator].copy()
        d["era_sort"] = d["era"].map({e: i for i, e in enumerate(era_order)})
        d = d.sort_values("era_sort").reset_index(drop=True)
        d["label"] = d["era"].map(era_labels)
        if add_signal and "signal" not in d.columns:
            d["signal"] = (
                (d["a"] >= 3)
                & (d["ROR_lower"] > 1)
                & (d["PRR"] >= 2)
                & (d["IC025"] > 0)
            )
        return d

    per_bupi = _prep(peroneal, "bupivacaine_plain")
    per_ropi = _prep(peroneal, "ropivacaine")
    last_bupi = _prep(last, "bupivacaine_plain", add_signal=False)
    last_ropi = _prep(last, "ropivacaine", add_signal=False)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
    _forest(axes[0, 0], per_bupi, "label",
            "A. Prolonged-block composite — vs plain bupivacaine",
            xlim=(0.3, 12))
    _forest(axes[0, 1], per_ropi, "label",
            "B. Prolonged-block composite — vs ropivacaine",
            xlim=(0.3, 12))
    _forest(axes[1, 0], last_bupi, "label",
            "C. LAST-spectrum composite — vs plain bupivacaine",
            xlim=(0.3, 12))
    _forest(axes[1, 1], last_ropi, "label",
            "D. LAST-spectrum composite — vs ropivacaine",
            xlim=(0.3, 12))

    for ax in axes.flat:
        ax.tick_params(axis="y", length=0)

    fig.suptitle(
        "Era-stratified disproportionality: Exparel vs active comparators\n"
        "FAERS 2012 Q4–2025 Q4, case-level composite endpoints",
        fontsize=12, y=0.99, x=0.02, ha="left"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "fig1_composite_forest.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")


def fig2_pt_forest(out_dir: Path) -> None:
    tbl = pd.read_csv("out/peroneal/pt_level_vs_bupivacaine_plain.csv")
    # Keep only PTs that meet a >= 3 and have IC025 > -1 (loose filter to include marginals)
    tbl = tbl[(tbl["a"] >= 3) & (tbl["IC025"] > -1)].copy()
    tbl = tbl.sort_values("ROR", ascending=False).reset_index(drop=True)
    # Drop columns with insane CIs that will break the x-scale
    tbl = tbl[tbl["ROR_upper"] < 20].reset_index(drop=True)
    tbl["signal"] = (
        (tbl["a"] >= 3)
        & (tbl["ROR_lower"] > 1)
        & (tbl["PRR"] >= 2)
        & (tbl["IC025"] > 0)
    )
    tbl["label"] = tbl["pt"].str.title() + "  [" + tbl["group"] + "]"

    fig, ax = plt.subplots(figsize=(9, 0.4 * len(tbl) + 2))
    _forest(ax, tbl, "label",
            "PT-level signals — Exparel vs plain bupivacaine",
            xlim=(0.3, 12))
    fig.tight_layout()
    p = out_dir / "fig2_pt_forest.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")


def fig3_era_trend(out_dir: Path) -> None:
    """Real per-year data — no interpolation. Counts come from
    ds.groupby(event_year) via case_meta.parquet join, and composite
    rates are computed per year from the analysis dataset."""
    cases = pd.read_csv("out/exparel_cases.csv")
    ds = pd.read_parquet("out/analysis_dataset.parquet")
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    # Join case_meta event_year onto every ds row so we can groupby year
    # regardless of drug_group (cases.csv only has Exparel).
    case_meta = pd.read_parquet("out/case_meta.parquet")
    year_map = case_meta.set_index("primaryid")["event_year"]
    ds["event_year"] = ds["primaryid"].map(year_map)

    ex_year = ds.loc[ds["drug_group"] == "exparel", ["primaryid", "event_year"]].drop_duplicates()
    by_year = ex_year.groupby("event_year")["primaryid"].nunique()
    by_year = by_year[(by_year.index >= 2012) & (by_year.index <= 2025)]

    from duration_analysis import PATHOLOGY_PTS, EXPECTED_PTS, SYSTEMIC_ABSORPTION_PTS
    from last_composite import LAST_PTS
    PERONEAL_PTS = PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS

    def _rate_by_year(pts):
        rows = []
        for y in sorted(by_year.index):
            year_ids = set(
                ds.loc[(ds["drug_group"] == "exparel") & (ds["event_year"] == y), "primaryid"]
            )
            hit_ids = set(
                ds.loc[(ds["drug_group"] == "exparel") & (ds["event_year"] == y) & ds["pt_lower"].isin(pts), "primaryid"]
            )
            rows.append({"year": int(y), "total": len(year_ids), "hit": len(hit_ids),
                         "rate": len(hit_ids) / len(year_ids) if year_ids else np.nan})
        return pd.DataFrame(rows)

    per = _rate_by_year(PERONEAL_PTS)
    lst = _rate_by_year(LAST_PTS)

    # Save the per-year table for Supplementary S8
    combined = per[["year", "total", "hit"]].rename(
        columns={"total": "exparel_cases", "hit": "prolonged_block_cases"}
    ).merge(
        lst[["year", "hit"]].rename(columns={"hit": "last_cases"}),
        on="year",
    )
    combined["prolonged_block_pct"] = 100 * combined["prolonged_block_cases"] / combined["exparel_cases"]
    combined["last_pct"] = 100 * combined["last_cases"] / combined["exparel_cases"]
    combined.to_csv(Path("out") / "annual_composite_rates.csv", index=False)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1.3]})

    ax1.bar(by_year.index, by_year.values, color="#2c3e50", alpha=0.85)
    ax1.set_ylabel("Exparel cases\n(reports per year)")
    ax1.set_title("Exparel FAERS reporting volume and composite signal rate, 2012-2025",
                  fontsize=11, loc="left")
    ax1.spines[["top", "right"]].set_visible(False)

    ax2.plot(per["year"], per["rate"] * 100, "-o", color="#c0392b",
             label="Prolonged-block composite", lw=2, markersize=6)
    ax2.plot(lst["year"], lst["rate"] * 100, "-o", color="#2980b9",
             label="LAST-spectrum composite", lw=2, markersize=6)
    ax2.set_ylabel("% of Exparel cases with ≥1 PT\nin composite")
    ax2.set_xlabel("Year")
    ax2.legend(loc="upper right", frameon=False)
    ax2.spines[["top", "right"]].set_visible(False)

    # Key event annotations
    events = [
        (2018.25, "ISB label"),
        (2021.1, "Anesthesiology\npublication"),
        (2025.0, "NOPAIN"),
    ]
    for x, txt in events:
        for ax in (ax1, ax2):
            ax.axvline(x, ls="--", color="#888", lw=0.7, alpha=0.5)
        ax2.text(x, ax2.get_ylim()[1] * 0.95, txt, fontsize=8,
                 rotation=90, va="top", ha="right", color="#555")

    fig.tight_layout()
    p = out_dir / "fig3_era_trend.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")


def fig4_pathology_pharmacology(out_dir: Path) -> None:
    """Three-panel duration/permanence stratification figure.
      A: forest plot (pathology vs expected vs systemic, vs both comparators)
      B: serious-outcome distribution by stratum
      C: event-to-report lag histogram
    """
    # --- Data ---
    split = pd.read_csv("out/duration/pt_split_pathology_vs_expected.csv")
    outcomes = pd.read_csv("out/duration/outcomes_by_stratum.csv")
    lag_stats = pd.read_csv("out/duration/event_to_report_lag.csv").iloc[0]

    # Forest rows
    forest_rows = []
    for _, r in split.iterrows():
        stratum = r["stratum"].replace("_", " ").title()
        comp = "bupivacaine" if r["comparator"] == "bupivacaine_plain" else "ropivacaine"
        forest_rows.append({
            "label": f"{stratum}\nvs {comp}",
            "a": int(r["a"]), "b": int(r["b"]),
            "ROR": r["ROR"], "ROR_lower": r["ROR_lower"], "ROR_upper": r["ROR_upper"],
            "signal": bool((r["a"] >= 3) and (r["ROR_lower"] > 1)),
        })
    forest = pd.DataFrame(forest_rows)

    # Outcome bar data (pathology vs expected-pharm)
    p_row = outcomes[outcomes["stratum"] == "pathology_implying"].iloc[0]
    e_row = outcomes[outcomes["stratum"] == "expected_pharmacology"].iloc[0]
    outc_pcts = {
        "Disability":      (p_row["pct_DS"], e_row["pct_DS"]),
        "Death":           (p_row["pct_DE"], e_row["pct_DE"]),
        "Life-\nthreatening": (p_row["pct_LT"], e_row["pct_LT"]),
        "Hospitalized":    (p_row["pct_HO"], e_row["pct_HO"]),
    }

    # Lag histogram — use the same stratum definition as duration_analysis
    # (post-adjustment sets without 'anaesthesia', with 'paralysis'/'monoparesis')
    # so Panel C n matches §3.3a's "n=97".
    cases = pd.read_csv("out/exparel_cases.csv")
    ds = pd.read_parquet("out/analysis_dataset.parquet")
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()
    from duration_analysis import PATHOLOGY_PTS, EXPECTED_PTS, SYSTEMIC_ABSORPTION_PTS
    PB_PTS = PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS
    hit_pids = set(ds.loc[
        (ds["drug_group"] == "exparel") & ds["pt_lower"].isin(PB_PTS),
        "primaryid"
    ])
    sub = cases[cases["primaryid"].isin(hit_pids)].copy()
    sub["event_dt_p"] = pd.to_datetime(sub["event_dt"], format="%Y%m%d", errors="coerce")
    sub["init_fda_dt_p"] = pd.to_datetime(sub["init_fda_dt"], format="%Y%m%d", errors="coerce")
    lag_days = (sub["init_fda_dt_p"] - sub["event_dt_p"]).dt.days.dropna()

    # --- Figure layout: A on top, B+C on bottom ---
    fig = plt.figure(figsize=(11.5, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.2], hspace=0.45, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # Panel A — forest
    _forest(ax_a, forest, "label",
            "A. Stratum-specific disproportionality — Exparel vs comparator",
            xlim=(0.2, 12))

    # Panel B — grouped bars
    cats = list(outc_pcts.keys())
    path_vals = [outc_pcts[k][0] for k in cats]
    exp_vals  = [outc_pcts[k][1] for k in cats]
    x = np.arange(len(cats))
    w = 0.38
    b1 = ax_b.bar(x - w/2, path_vals, w, color="#c0392b", label=f"Pathology-implying (n={int(p_row['n_cases'])})")
    b2 = ax_b.bar(x + w/2, exp_vals,  w, color="#34495e", label=f"Expected-pharmacology (n={int(e_row['n_cases'])})")
    for bar, val in zip(b1, path_vals):
        ax_b.text(bar.get_x() + bar.get_width()/2, val + 0.2,
                  f"{val:.1f}%", ha="center", fontsize=8, color="#c0392b")
    for bar, val in zip(b2, exp_vals):
        ax_b.text(bar.get_x() + bar.get_width()/2, val + 0.2,
                  f"{val:.1f}%", ha="center", fontsize=8, color="#34495e")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(cats, fontsize=9)
    ax_b.set_ylabel("% of stratum cases")
    ax_b.set_title("B. Serious-outcome distribution by stratum", loc="left", fontsize=11)
    ax_b.legend(fontsize=8, loc="upper right", frameon=False)
    ax_b.spines[["top", "right"]].set_visible(False)

    # Panel C — lag histogram
    bins = [0, 30, 90, 180, 365, max(366, int(lag_days.max()) + 1)]
    labels = ["<30 d", "30–90 d", "90–180 d", "180–365 d", ">365 d"]
    counts, _ = np.histogram(lag_days, bins=bins)
    ax_c.bar(labels, counts, color=["#95a5a6", "#7f8c8d", "#34495e", "#c0392b", "#8e44ad"])
    for i, c in enumerate(counts):
        ax_c.text(i, c + 0.5, f"{int(c)}\n({100*c/len(lag_days):.0f}%)",
                  ha="center", fontsize=8)
    ax_c.set_ylabel(f"Cases (of {len(lag_days)} with event + report dates)")
    ax_c.set_title("C. Event-to-report lag", loc="left", fontsize=11)
    ax_c.set_ylim(0, max(counts) * 1.25)
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.text(0.98, 0.97,
              f"Median {int(lag_stats['median_lag_days'])} d "
              f"(IQR {int(lag_stats['q25_lag_days'])}–{int(lag_stats['q75_lag_days'])})\n"
              f"{lag_stats['pct_lag_gt_90d']:.0f}% reported > 3 months after event",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=8,
              bbox=dict(boxstyle="round,pad=0.3", fc="#fef9e7", ec="#c0392b"))

    fig.suptitle("Figure 4. Duration and permanence stratification of the prolonged-block composite",
                 fontsize=12, y=0.995, x=0.02, ha="left")
    p = out_dir / "fig4_pathology_pharmacology.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=Path("fig"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)

    fig1_composite(args.out)
    fig2_pt_forest(args.out)
    fig3_era_trend(args.out)
    fig4_pathology_pharmacology(args.out)
    print(f"\nAll figures → {args.out}/")


if __name__ == "__main__":
    main()

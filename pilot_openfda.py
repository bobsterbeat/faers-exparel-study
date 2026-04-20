"""
FAERS Exparel Disproportionality Analysis — Pilot (openFDA API version)
========================================================================

This is the PILOT script. It uses the openFDA API for rapid iteration.
For the primary analysis, switch to FAERS quarterly ASCII files (see
`faers_ascii_pipeline.py`).

Author: R. Aldwinckle
Status: Working skeleton — run end-to-end for a rapid signal scan.
License: MIT

Quick start:
    python -m pip install pandas numpy scipy requests matplotlib seaborn
    python pilot_openfda.py

Outputs:
    /out/exparel_signals_pilot.csv
    /out/exparel_top_pts.png
    /out/exparel_timetrend.png
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from scipy.stats import chi2_contingency, norm

# Minimal .env loader — no external dependency. Reads KEY=VALUE lines.
_env_path = Path(__file__).parent / ".env"
if _env_path.exists():
    for _line in _env_path.read_text().splitlines():
        _line = _line.strip()
        if not _line or _line.startswith("#") or "=" not in _line:
            continue
        _k, _v = _line.split("=", 1)
        os.environ.setdefault(_k.strip(), _v.strip())

OUT = Path(__file__).parent / "out"
OUT.mkdir(exist_ok=True)

OPENFDA_BASE = os.environ.get("OPENFDA_BASE_URL", "https://api.fda.gov/drug/event.json")
OPENFDA_API_KEY = os.environ.get("OPENFDA_API_KEY", "").strip() or None
HEADERS = {
    "User-Agent": os.environ.get(
        "HTTP_USER_AGENT", "faers-exparel-study/1.0 (academic research)"
    )
}

# ---------------------------------------------------------------------------
# 1. Drug definitions and search strings
# ---------------------------------------------------------------------------

# openFDA drug field search — returns counts per MedDRA PT
# Note: openFDA API has rate limits (240/min, 1000/day without API key;
# register a free key at open.fda.gov/apis/authentication for higher limits).

EXPAREL_QUERY = (
    '(patient.drug.medicinalproduct:"EXPAREL" '
    'OR patient.drug.openfda.brand_name:"EXPAREL" '
    'OR patient.drug.openfda.generic_name:"BUPIVACAINE LIPOSOME" '
    'OR patient.drug.medicinalproduct.exact:"BUPIVACAINE LIPOSOME")'
)

# Plain bupivacaine — PRIMARY comparator. Isolates the Exparel-specific
# differential (formulation + use pattern) from the molecule itself.
BUPIVACAINE_QUERY = (
    '(patient.drug.medicinalproduct:"BUPIVACAINE" '
    'OR patient.drug.openfda.generic_name:"BUPIVACAINE") '
    'AND NOT patient.drug.medicinalproduct:"LIPOSOME" '
    'AND NOT patient.drug.medicinalproduct:"EXPAREL"'
)

# Ropivacaine — SECONDARY comparator. Alternative long-acting amide LA;
# places findings in context of other perineural agents. Different
# pharmacology (lower cardiac tox, less motor block) so interpret with care.
ROPIVACAINE_QUERY = (
    '(patient.drug.medicinalproduct:"ROPIVACAINE" '
    'OR patient.drug.openfda.generic_name:"ROPIVACAINE")'
)

COMPARATORS = [
    ("bupivacaine", "primary", BUPIVACAINE_QUERY),
    ("ropivacaine", "secondary", ROPIVACAINE_QUERY),
]

# Pre-specified AE categories — MedDRA PTs for pilot scan
# (Full list lives in the protocol; this is a minimal working set)
AE_CATEGORIES: dict[str, list[str]] = {
    "prolonged_block": [
        "Hypoaesthesia", "Paraesthesia", "Muscular weakness",
        "Neurological symptom", "Hemiparesis", "Hypoaesthesia oral",
    ],
    "LAST_spectrum": [
        "Cardiac arrest", "Ventricular tachycardia", "Seizure",
        "Loss of consciousness", "Bradycardia", "Convulsion",
    ],
    "neurotoxicity": [
        "Peripheral nerve injury", "Neuralgia", "Neuropathy peripheral",
        "Nerve injury", "Brachial plexopathy",
    ],
    "wound": [
        "Wound dehiscence", "Wound", "Wound healing impaired",
        "Wound infection", "Incision site complication",
    ],
    "chondrolysis": [
        "Cartilage injury", "Joint destruction", "Chondromalacia",
    ],
    "admin_error": [
        "Product administered to patient of inappropriate age",
        "Incorrect route of product administration", "Medication error",
        "Off label use", "Wrong drug administered",
        "Product administered at inappropriate site",
    ],
    "hypersensitivity": [
        "Anaphylactic reaction", "Hypersensitivity", "Drug eruption",
        "Angioedema", "Rash",
    ],
    "cardiac": [
        "Hypotension", "Cardiac arrest", "Atrial fibrillation",
        "Ventricular tachycardia", "Cardiovascular insufficiency",
    ],
    "cns": [
        "Somnolence", "Confusional state", "Dizziness", "Tremor",
        "Visual impairment",
    ],
}

# ---------------------------------------------------------------------------
# 2. openFDA fetchers
# ---------------------------------------------------------------------------

def _with_key(params: dict) -> dict:
    if OPENFDA_API_KEY:
        params = {**params, "api_key": OPENFDA_API_KEY}
    return params


# openFDA caps limit=1000 on count queries and skip up to ~25,000.
# Pagination below retrieves the full PT distribution, avoiding the
# b=0 truncation artifact seen when limit is fixed at 500.
def fetch_pt_counts(
    drug_query: str,
    per_call: int = 1000,
    max_total: int = 25000,
) -> tuple[pd.DataFrame, bool]:
    """
    Return (df, truncated) where df has columns [pt, n] covering every PT
    returned by the count endpoint, and `truncated` is True iff pagination
    hit max_total (i.e., there may be more PTs we did not retrieve).
    """
    frames = []
    skip = 0
    truncated = False
    while skip < max_total:
        params = _with_key({
            "search": drug_query,
            "count": "patient.reaction.reactionmeddrapt.exact",
            "limit": per_call,
            "skip": skip,
        })
        r = requests.get(OPENFDA_BASE, params=params, headers=HEADERS, timeout=60)
        # openFDA returns 404 when no results and 400 when skip exceeds
        # the available result set — both mean "we're done paginating".
        if r.status_code in (400, 404):
            break
        r.raise_for_status()
        results = r.json().get("results", [])
        if not results:
            break
        frames.append(pd.DataFrame(results))
        if len(results) < per_call:
            break
        skip += per_call
        time.sleep(0.1)
    else:
        truncated = True

    if not frames:
        return pd.DataFrame(columns=["pt", "n"]), truncated
    df = pd.concat(frames, ignore_index=True).rename(
        columns={"term": "pt", "count": "n"}
    )
    return df, truncated


def fetch_total_reports(drug_query: str) -> int:
    """Return total number of reports matching drug query."""
    params = _with_key({"search": drug_query, "limit": 1})
    r = requests.get(OPENFDA_BASE, params=params, headers=HEADERS, timeout=60)
    r.raise_for_status()
    data = r.json()
    return data.get("meta", {}).get("results", {}).get("total", 0)


def fetch_time_trend(drug_query: str) -> pd.DataFrame:
    """Return counts by receivedate year/month for time-trend plot."""
    params = _with_key({
        "search": drug_query,
        "count": "receivedate",
        "limit": 1000,
    })
    r = requests.get(OPENFDA_BASE, params=params, headers=HEADERS, timeout=60)
    r.raise_for_status()
    data = r.json()
    if "results" not in data:
        return pd.DataFrame(columns=["date", "n"])
    df = pd.DataFrame(data["results"])
    df["date"] = pd.to_datetime(df["time"], format="%Y%m%d", errors="coerce")
    df["n"] = df["count"]
    df = df.dropna(subset=["date"]).sort_values("date")
    df["year_month"] = df["date"].dt.to_period("M")
    return df.groupby("year_month", as_index=False)["n"].sum().assign(
        year_month=lambda d: d["year_month"].dt.to_timestamp()
    )

# ---------------------------------------------------------------------------
# 3. Disproportionality calculations
# ---------------------------------------------------------------------------

@dataclass
class Signal:
    pt: str
    a: int          # Exparel + AE
    b: int          # Comparator + AE
    c: int          # Exparel + other AE
    d: int          # Comparator + other AE
    ROR: float
    ROR_lower: float
    ROR_upper: float
    PRR: float
    PRR_chi2: float
    IC: float
    IC025: float

    @property
    def signal(self) -> bool:
        return (
            self.a >= 3
            and self.ROR_lower > 1
            and self.PRR >= 2
            and self.PRR_chi2 >= 4
            and self.IC025 > 0
        )

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["signal"] = self.signal
        return d


def compute_signal(pt: str, a: int, b: int, c: int, d: int) -> Signal:
    """Compute ROR, PRR, and IC (BCPNN) for a single drug-event pair."""
    # Haldane-Anscombe correction for zero cells
    aa, bb, cc, dd = [x + 0.5 if x == 0 else x for x in (a, b, c, d)]

    # ROR
    ror = (aa / cc) / (bb / dd)
    se_ln_ror = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    ln_ror = np.log(ror)
    ror_lower = np.exp(ln_ror - 1.96 * se_ln_ror)
    ror_upper = np.exp(ln_ror + 1.96 * se_ln_ror)

    # PRR
    prr = (aa / (aa + cc)) / (bb / (bb + dd))

    # Chi-squared (Yates-corrected)
    try:
        chi2, _, _, _ = chi2_contingency(
            [[a, b], [c, d]], correction=True
        )
    except ValueError:
        chi2 = 0.0

    # Information Component (BCPNN, per Bate et al.)
    N = a + b + c + d
    if N == 0 or (a + c) == 0 or (a + b) == 0:
        ic = ic025 = -np.inf
    else:
        exp_a = (a + b) * (a + c) / N
        ic = np.log2((a + 0.5) / (exp_a + 0.5))
        # Approximate 95% CI for IC
        var_ic = (1 / np.log(2) ** 2) * (
            (N - a + 0.5) / ((a + 0.5) * (1 + N))
            + (N - (a + b) + 0.5) / (((a + b) + 0.5) * (1 + N))
            + (N - (a + c) + 0.5) / (((a + c) + 0.5) * (1 + N))
        )
        ic025 = ic - 1.96 * np.sqrt(var_ic)

    return Signal(
        pt=pt, a=a, b=b, c=c, d=d,
        ROR=ror, ROR_lower=ror_lower, ROR_upper=ror_upper,
        PRR=prr, PRR_chi2=chi2,
        IC=ic, IC025=ic025,
    )

# ---------------------------------------------------------------------------
# 4. Build contingency tables and iterate
# ---------------------------------------------------------------------------

def build_contingency(
    target_pts: pd.DataFrame,
    comparator_pts: pd.DataFrame,
    target_total: int,
    comparator_total: int,
    comparator_truncated: bool = False,
) -> pd.DataFrame:
    """
    For every PT appearing in either drug's reports, build the 2x2 and
    compute signal metrics. Adds a `b_is_truncation_candidate` flag for
    rows where b=0 AND the comparator fetch was paginated to its cap —
    i.e., we cannot rule out that the true b is just outside our window.
    """
    t_map = dict(zip(target_pts["pt"].str.upper(), target_pts["n"]))
    c_map = dict(zip(comparator_pts["pt"].str.upper(), comparator_pts["n"]))
    pts = set(t_map).union(c_map)

    rows = []
    for pt in pts:
        a = int(t_map.get(pt, 0))
        b = int(c_map.get(pt, 0))
        c = target_total - a
        d = comparator_total - b
        if a + b < 3:
            continue
        sig = compute_signal(pt, a, b, c, d)
        r = sig.as_dict()
        r["b_is_truncation_candidate"] = bool(b == 0 and comparator_truncated)
        rows.append(r)

    df = pd.DataFrame(rows).sort_values(
        ["signal", "ROR"], ascending=[False, False]
    )
    return df


def aggregate_by_category(
    target_pts: pd.DataFrame,
    comparator_pts: pd.DataFrame,
    target_total: int,
    comparator_total: int,
    ae_categories: dict[str, list[str]],
) -> pd.DataFrame:
    """
    Composite-endpoint analysis: sum a and b across all PTs within each
    pre-specified clinical category, then run the same disproportionality
    calc on the aggregated 2x2. Aggregation is pre-specified (categories
    fixed in the protocol, not data-driven) — disclose this in methods.
    """
    t_map = dict(zip(target_pts["pt"].str.upper(), target_pts["n"]))
    c_map = dict(zip(comparator_pts["pt"].str.upper(), comparator_pts["n"]))

    rows = []
    for cat, pts in ae_categories.items():
        pts_up = [p.upper() for p in pts]
        a = sum(int(t_map.get(p, 0)) for p in pts_up)
        b = sum(int(c_map.get(p, 0)) for p in pts_up)
        c = target_total - a
        d = comparator_total - b
        n_matched = sum(1 for p in pts_up if t_map.get(p, 0) or c_map.get(p, 0))
        if a + b < 3:
            continue
        sig = compute_signal(cat, a, b, c, d)
        r = sig.as_dict()
        r["category"] = cat
        r["n_prespec_pts"] = len(pts_up)
        r["n_matching_pts"] = n_matched
        rows.append(r)

    return pd.DataFrame(rows).sort_values(
        ["signal", "ROR"], ascending=[False, False]
    )

# ---------------------------------------------------------------------------
# 5. Main pilot pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"openFDA API key: {'loaded' if OPENFDA_API_KEY else 'NOT set (anonymous, 1k/day)'}")

    pt_to_category = {
        pt.upper(): cat for cat, pts in AE_CATEGORIES.items() for pt in pts
    }

    print("Fetching Exparel reports…")
    ex_total = fetch_total_reports(EXPAREL_QUERY)
    ex_pts, ex_trunc = fetch_pt_counts(EXPAREL_QUERY)
    print(f"  {ex_total:,} reports, {len(ex_pts):,} distinct PTs"
          f"{' [TRUNCATED]' if ex_trunc else ''}")

    all_pt_signals: list[pd.DataFrame] = []
    all_category_signals: list[pd.DataFrame] = []

    for comp_name, tier, comp_query in COMPARATORS:
        print(f"\nFetching {comp_name} ({tier} comparator)…")
        comp_total = fetch_total_reports(comp_query)
        time.sleep(0.1)
        comp_pts, comp_trunc = fetch_pt_counts(comp_query)
        print(f"  {comp_total:,} reports, {len(comp_pts):,} distinct PTs"
              f"{' [TRUNCATED]' if comp_trunc else ''}")

        print(f"  Computing PT-level signals vs {comp_name}…")
        pt_sigs = build_contingency(
            ex_pts, comp_pts, ex_total, comp_total,
            comparator_truncated=comp_trunc,
        )
        pt_sigs["comparator"] = comp_name
        pt_sigs["tier"] = tier
        pt_sigs["prespecified_category"] = (
            pt_sigs["pt"].str.upper().map(pt_to_category)
        )
        all_pt_signals.append(pt_sigs)

        print(f"  Aggregating by pre-specified category…")
        cat_sigs = aggregate_by_category(
            ex_pts, comp_pts, ex_total, comp_total, AE_CATEGORIES,
        )
        cat_sigs["comparator"] = comp_name
        cat_sigs["tier"] = tier
        all_category_signals.append(cat_sigs)

        n_pt_sig = pt_sigs["signal"].sum()
        n_pt_sig_b_zero = ((pt_sigs["signal"]) & (pt_sigs["b"] == 0)).sum()
        n_cat_sig = cat_sigs["signal"].sum() if not cat_sigs.empty else 0
        print(f"  → {n_pt_sig} PT signals  "
              f"({n_pt_sig - n_pt_sig_b_zero} with b≥1, "
              f"{n_pt_sig_b_zero} b=0)  "
              f"| {n_cat_sig} category signals")

    pt_all = pd.concat(all_pt_signals, ignore_index=True)
    cat_all = pd.concat(all_category_signals, ignore_index=True)

    # ----- three-tier outputs -----
    primary_pt = pt_all[(pt_all["signal"]) & (pt_all["b"] >= 1)].copy()
    primary_pt.to_csv(OUT / "primary_pt_signals.csv", index=False)

    cat_all[cat_all["signal"]].to_csv(OUT / "category_signals.csv", index=False)
    cat_all.to_csv(OUT / "category_signals_all.csv", index=False)

    supp_b_zero = pt_all[(pt_all["signal"]) & (pt_all["b"] == 0)].copy()
    supp_b_zero.insert(
        0, "caveat",
        "Comparator count is zero after full pagination; ROR is inflated "
        "by Haldane-Anscombe correction. Interpret with caution — hypothesis-"
        "generating only, do not quote effect sizes without sensitivity analysis."
    )
    supp_b_zero.to_csv(OUT / "supplementary_b_zero.csv", index=False)

    # full PT tables per comparator for reference / sensitivity work
    pt_all.to_csv(OUT / "pt_signals_full.csv", index=False)

    print("\n" + "=" * 70)
    print("TIER SUMMARY")
    print("=" * 70)
    for comp_name, tier, _ in COMPARATORS:
        slc = primary_pt[primary_pt["comparator"] == comp_name]
        print(f"\n{tier.upper()} — Exparel vs {comp_name}")
        print(f"  Primary PT signals (b≥1): {len(slc)}")
        if len(slc):
            print(slc.sort_values("ROR", ascending=False).head(15)[
                ["pt", "a", "b", "ROR", "ROR_lower", "PRR", "IC025",
                 "prespecified_category"]
            ].to_string(index=False))
        cats = cat_all[(cat_all["comparator"] == comp_name) & (cat_all["signal"])]
        print(f"\n  Category signals: {len(cats)} of {cat_all[cat_all['comparator']==comp_name].shape[0]} testable")
        if len(cats):
            print(cats[["category", "a", "b", "ROR", "ROR_lower", "PRR",
                         "IC025", "n_matching_pts"]].to_string(index=False))

    print("\nStep: time trend…")
    try:
        tt = fetch_time_trend(EXPAREL_QUERY)
        tt.to_csv(OUT / "exparel_time_trend.csv", index=False)
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(tt["year_month"], tt["n"], linewidth=1.2)
            ax.set_title("Exparel — FAERS reports per month")
            ax.set_xlabel("Month")
            ax.set_ylabel("Reports")
            events = [
                ("2018-04-01", "ISB label"),
                ("2021-02-01", "Anesthesiology pub"),
                ("2021-04-01", "Pacira v. ASA"),
                ("2025-01-01", "NOPAIN"),
            ]
            for dt, label in events:
                ax.axvline(pd.Timestamp(dt), color="#888", linestyle="--", linewidth=0.8)
                ax.text(pd.Timestamp(dt), ax.get_ylim()[1] * 0.95, label,
                        rotation=90, fontsize=8, va="top", color="#555")
            plt.tight_layout()
            plt.savefig(OUT / "exparel_timetrend.png", dpi=150)
            print(f"  saved plot to {OUT / 'exparel_timetrend.png'}")
        except ImportError:
            print("  matplotlib not installed — skipping plot")
    except Exception as e:
        print(f"  time trend fetch failed: {e}")

    print(f"\nPilot complete. Results in {OUT}")


if __name__ == "__main__":
    main()

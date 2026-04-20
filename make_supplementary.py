"""
Build the remaining supplementary artifacts that require computation:

  S6: Representative case narratives (privacy-respecting — ages redacted
      to 5-year bands, concomitant drugs and reporter country omitted).

  S7: Drug-name regex audit — random sample of matched and unmatched
      free-text drug names to show the regex catches intended variants
      and is not over-broad.

Usage: python make_supplementary.py
"""

from __future__ import annotations

import re
import random
from pathlib import Path

import pandas as pd

from ascii_pipeline import LIPOSOMAL_PATTERN, PLAIN_BUPI_PATTERN, ROPI_PATTERN


def age_band(age) -> str:
    if pd.isna(age):
        return "unknown"
    a = int(age)
    lo = (a // 5) * 5
    return f"{lo}–{lo+4}"


def s6_case_narratives(out_dir: Path) -> None:
    """Twenty representative cases (ages redacted to 5-year bands, no
    concomitant-drug lists, no reporter country, no free-text drugname).
    Covers era and sub-phenotype: peroneal-specific, prolonged-block,
    pathology, systemic-absorption."""
    cases = pd.read_csv("out/exparel_cases.csv")
    ds = pd.read_parquet("out/analysis_dataset.parquet")
    ds["pt_lower"] = ds["pt"].astype(str).str.lower()

    from duration_analysis import PATHOLOGY_PTS, EXPECTED_PTS, SYSTEMIC_ABSORPTION_PTS
    all_block = PATHOLOGY_PTS | EXPECTED_PTS | SYSTEMIC_ABSORPTION_PTS

    pids_any = set(ds.loc[
        (ds["drug_group"] == "exparel") & ds["pt_lower"].isin(all_block),
        "primaryid",
    ])
    cohort = cases[cases["primaryid"].isin(pids_any)].copy()

    # Pick 20 spanning era + stratum
    random.seed(42)
    selected = []

    def _pick_from(mask, n):
        sub = cohort[mask]
        if len(sub) == 0:
            return []
        return sub.sample(min(n, len(sub)), random_state=42).to_dict("records")

    # 3 pre-2018, 7 from 2018-2020, 10 from 2021+
    selected += _pick_from(cohort["event_year"].between(2012, 2017), 3)
    selected += _pick_from(cohort["event_year"].between(2018, 2020), 7)
    selected += _pick_from(cohort["event_year"] >= 2021, 10)

    rows = []
    for c in selected:
        rows.append({
            "era": ("pre-2018" if c["event_year"] <= 2017
                    else "2018–2020" if c["event_year"] <= 2020
                    else "2021+"),
            "age_band": age_band(c.get("age")),
            "sex": c.get("sex") if pd.notna(c.get("sex")) else c.get("gndr_cod", "unknown"),
            "reporter_type": c.get("occp_cod", "unknown"),
            "indication": c.get("indications", ""),
            "reactions": c.get("reactions", ""),
            "outcomes": c.get("outcomes", ""),
        })

    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "S6_case_narratives.csv", index=False)
    print(f"S6: {len(out)} cases written")


def s7_regex_audit(out_dir: Path) -> None:
    """Random sample of matched and unmatched drugnames to show the regex
    is neither too strict nor too permissive. Reviewers can inspect."""
    data = Path("data/raw")
    quarter_dirs = sorted(
        [p for p in data.iterdir() if p.is_dir() and "faers" in p.name.lower()]
    )
    # Sample a recent quarter only — representative and fast
    recent = [q for q in quarter_dirs if "2024" in q.name or "2025" in q.name]
    frames = []
    for q in recent[:4]:
        for sub in ("ASCII", "ascii"):
            if (q / sub).exists():
                ascii_dir = q / sub
                break
        else:
            ascii_dir = q
        drug_f = [p for p in (list(ascii_dir.glob("*.txt")) + list(ascii_dir.glob("*.TXT")))
                  if p.stem.upper().startswith("DRUG")]
        if drug_f:
            df = pd.read_csv(drug_f[0], sep="$", encoding="latin-1",
                             low_memory=False, on_bad_lines="skip",
                             usecols=lambda c: c.lower() in {"drugname", "prod_ai", "role_cod"})
            df.columns = [c.lstrip("\ufeff").lstrip("ï»¿").lower() for c in df.columns]
            frames.append(df)
    drug = pd.concat(frames, ignore_index=True)
    drug = drug[drug["role_cod"] == "PS"]

    def _classify(s):
        if not isinstance(s, str):
            return None
        if LIPOSOMAL_PATTERN.search(s):
            return "exparel"
        if ROPI_PATTERN.search(s):
            return "ropivacaine"
        if PLAIN_BUPI_PATTERN.search(s):
            return "bupivacaine_plain"
        return None

    drug["match"] = drug["drugname"].apply(_classify)

    # 1. Top 20 matched drugnames per group — show what the regex captures
    matched_summary = []
    for grp in ("exparel", "bupivacaine_plain", "ropivacaine"):
        top = (drug[drug["match"] == grp]["drugname"]
               .value_counts().head(20).reset_index())
        top.columns = ["drugname", "count"]
        top["drug_group"] = grp
        matched_summary.append(top)
    matched = pd.concat(matched_summary, ignore_index=True)
    matched.to_csv(out_dir / "S7_regex_matched_samples.csv", index=False)

    # 2. Random sample of 30 UNMATCHED drugnames that plausibly could be
    #    a local anesthetic — show the regex doesn't catch unrelated drugs.
    unmatched = drug[drug["match"].isna()].copy()
    # Filter to LA-adjacent keywords to make review tractable
    adj = unmatched[unmatched["drugname"].astype(str).str.contains(
        "CAINE|ANAES|ANEST|BLOCK|LIPOSOM|NERVE", case=False, na=False
    )]
    sample_n = min(30, len(adj))
    sample = (adj["drugname"].value_counts().head(sample_n).reset_index()
              if sample_n else pd.DataFrame(columns=["drugname", "count"]))
    sample.columns = ["drugname", "count"]
    sample.to_csv(out_dir / "S7_regex_unmatched_la_adjacent.csv", index=False)

    print(f"S7: {len(matched)} matched + {len(sample)} unmatched-adjacent rows written")


def main() -> None:
    out_dir = Path("supplementary")
    out_dir.mkdir(exist_ok=True)
    s6_case_narratives(out_dir)
    s7_regex_audit(out_dir)
    print(f"\nSupplementary files → {out_dir}/")


if __name__ == "__main__":
    main()

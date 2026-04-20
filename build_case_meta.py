"""
Build a per-case metadata table (primaryid, event_year) across ALL cases,
not just Exparel. Needed so that era-stratified disproportionality uses
era-matched comparator pools rather than the full-period background —
which is what a reviewer will (correctly) demand.

Fastest path: scan only the DEMO tables across all quarters (smaller than
DRUG/REAC). Uses the same dedup and event-year-fallback logic as the main
pipeline.

Output: out/case_meta.parquet with columns [primaryid, event_year].

Usage:  python build_case_meta.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


DATA_DIR = Path("data/raw")
OUT = Path("out/case_meta.parquet")


def _read_demo(quarter_dir: Path) -> pd.DataFrame | None:
    for sub in ("ASCII", "ascii"):
        if (quarter_dir / sub).exists():
            ascii_dir = quarter_dir / sub
            break
    else:
        ascii_dir = quarter_dir
    txt = list(ascii_dir.glob("*.txt")) + list(ascii_dir.glob("*.TXT"))
    demo_f = [p for p in txt if p.stem.upper().startswith("DEMO")]
    if not demo_f:
        return None
    df = pd.read_csv(demo_f[0], sep="$", encoding="latin-1",
                     low_memory=False, on_bad_lines="skip")
    df.columns = [c.lstrip("\ufeff").lstrip("ï»¿").lower() for c in df.columns]
    keep = [c for c in ("primaryid", "caseid", "fda_dt", "event_dt",
                        "init_fda_dt") if c in df.columns]
    return df[keep]


def main() -> None:
    quarter_dirs = sorted(
        [p for p in DATA_DIR.iterdir() if p.is_dir() and "faers" in p.name.lower()]
    )
    print(f"Loading DEMO from {len(quarter_dirs)} quarters…")
    frames = []
    for q in quarter_dirs:
        d = _read_demo(q)
        if d is not None:
            frames.append(d)
    demo = pd.concat(frames, ignore_index=True)
    print(f"  total rows: {len(demo):,}")

    # FDA-standard dedup: per caseid, keep latest fda_dt, tie-break highest primaryid.
    before = len(demo)
    demo = demo.sort_values(["caseid", "fda_dt", "primaryid"],
                            ascending=[True, False, False])
    demo = demo.drop_duplicates(subset=["caseid"], keep="first")
    print(f"  after dedup: {len(demo):,}  (-{before-len(demo):,})")

    # Event year with fallback: event_dt > init_fda_dt > fda_dt.
    def _year(row) -> int | pd._libs.missing.NAType:
        for col in ("event_dt", "init_fda_dt", "fda_dt"):
            v = row.get(col)
            if pd.notna(v):
                y = pd.to_datetime(v, format="%Y%m%d", errors="coerce")
                if pd.notna(y):
                    return y.year
        return pd.NA

    print("  deriving event_year with fallback…")
    demo["event_year"] = demo.apply(_year, axis=1)
    out = demo[["primaryid", "event_year"]].copy()
    out["primaryid"] = out["primaryid"].astype("int64")

    OUT.parent.mkdir(exist_ok=True, parents=True)
    out.to_parquet(OUT, index=False)
    covered = out["event_year"].notna().sum()
    print(f"  saved {len(out):,} rows → {OUT}")
    print(f"  event_year populated: {covered:,}/{len(out):,} ({100*covered/len(out):.1f}%)")


if __name__ == "__main__":
    main()

# FAERS Exparel Disproportionality Study

Real-world adverse event signal detection for liposomal bupivacaine (Exparel®) compared with plain bupivacaine and ropivacaine, using the FDA Adverse Event Reporting System (FAERS), 2012 Q4 – 2025 Q4.

Manuscript: *Liposomal Bupivacaine (Exparel®) and Prolonged Sensory/Motor Block: A Disproportionality Analysis of FAERS 2012–2025* (in preparation).

## Headline findings

- **Prolonged sensory/motor block composite** (n=2,127 Exparel cases): all-time ROR 2.91 (95% CI 2.40–3.52) vs plain bupivacaine; signal intensifies after the 2018 interscalene block indication (2018–2020 ROR 4.52) and persists into 2021+ (ROR 2.26).
- **Pathology-implying sub-stratum** (PTs whose MedDRA definitions cannot represent expected <72 h pharmacology): ROR 2.16 (95% CI 1.58–2.95), holds under serious-outcome-restricted sensitivity (ROR 2.20). 8.5% of flagged cases carry a Disability outcome code.
- **LAST-spectrum composite** previously reported in FAERS: confined to pre-2018 (ROR 3.69), inverted in the modern era (2018–2020 ROR 0.67), consistent with a Weber-effect reporting artifact.
- **prod_ai misclassification audit**: 2,145 of 2,209 Exparel primary-suspect records (97.1%) would be silently reclassified as plain bupivacaine by naive `prod_ai.fillna(drugname)` logic — a pre-registered methodological refinement of prior FAERS work.

## Reproduction

### Prerequisites
- Python 3.11+ (tested on 3.14)
- `pip install -r requirements.txt`
- ~60 GB free disk for raw FAERS + derived outputs

### One-command (if you've downloaded raw FAERS)
```bash
make all       # (if Makefile added) — OR run the steps below
```

### Step by step

```bash
# 1. Fetch FAERS quarterly ASCII releases (12–15 GB, ~45 min)
bash download_faers.sh
bash verify_downloads.sh

# 2. Unzip
cd data/raw && for f in faers_ascii_*.zip; do
  d="${f%.zip}"; [ -d "$d" ] || unzip -q "$f" -d "$d"
done && cd ..

# 3. Build per-case metadata (primaryid → event_year, with fallback)
python3 build_case_meta.py

# 4. Main ingest and case export
python3 ascii_pipeline.py --data-dir data/raw --out out/

# 5. Composite-endpoint and subgroup analyses
python3 last_composite.py   --dataset out/analysis_dataset.parquet --cases out/exparel_cases.csv --out out/last_composite.csv
python3 peroneal_palsy.py   --dataset out/analysis_dataset.parquet --cases out/exparel_cases.csv --out-dir out/peroneal/
python3 duration_analysis.py --dataset out/analysis_dataset.parquet --cases out/exparel_cases.csv --out-dir out/duration/
python3 sensitivity_analysis.py --dataset out/analysis_dataset.parquet --cases out/exparel_cases.csv --out out/sensitivity.csv
python3 audit_misclassification.py

# 6. Presentation artifacts
python3 make_figures.py --out fig/
python3 make_tables.py  --out out/tables/
python3 make_report.py             # produces report.html
```

Every step is idempotent. Rerunning after a FAERS update is safe.

## Repository layout

```
faers-exparel-study/
├── download_faers.sh            # fetch FAERS ASCII releases with retry/resume
├── verify_downloads.sh          # CRC + content check for every ZIP
├── ascii_pipeline.py            # ingest 53 quarters, build analysis dataset
├── build_case_meta.py           # per-case event_year with fallback
├── last_composite.py            # LAST spectrum composite + era stratification
├── peroneal_palsy.py            # prolonged-block subgroup deep dive
├── duration_analysis.py         # pathology/expected-pharm/systemic strata
├── sensitivity_analysis.py      # 4 pre-registered restrictions
├── audit_misclassification.py   # verify prod_ai misclassification count
├── make_figures.py              # Figures 1–4
├── make_tables.py               # Tables 1–3
├── make_report.py               # self-contained HTML report
├── pilot_openfda.py             # (pilot) openFDA API version for rapid iteration
├── manuscript.md                # manuscript source
├── report.html                  # rendered manuscript + figures + tables
├── report.pdf                   # PDF export of report.html
├── fig/                         # Figures 1–4 PNG
├── out/
│   ├── analysis_dataset.parquet # (primaryid, drug_group, pt) case-level
│   ├── case_meta.parquet        # primaryid → event_year (17M rows)
│   ├── exparel_cases.csv        # Exparel case-level extract (2,127 cases)
│   ├── primary_signals.csv      # PT-level signals vs plain bupivacaine
│   ├── secondary_ropivacaine_signals.csv
│   ├── last_composite.csv       # LAST composite by era × comparator
│   ├── sensitivity.csv          # 4 sensitivity analyses
│   ├── annual_composite_rates.csv # per-year rates, supports Figure 3
│   ├── peroneal/                # subgroup deep-dive outputs
│   ├── duration/                # pathology/expected/systemic strata
│   └── tables/                  # Tables 1–3 CSV + PNG
└── supplementary/               # READUS-PV supplementary appendices
```

Raw FAERS data (`data/raw/`) is gitignored — it is public domain and fetched by `download_faers.sh` directly from FDA.

## Methodological notes for re-users

1. **prod_ai precedence.** FAERS `prod_ai` (Product Active Ingredient) normalizes Exparel to `"BUPIVACAINE"`. A naive `prod_ai.fillna(drugname)` classifier will silently misclassify ~97% of Exparel records from 2014 Q3 onward. This pipeline uses drugname-first, most-specific-label-wins precedence. See `audit_misclassification.py` for verification.
2. **BOM handling.** FAERS 2012 Q4 files are UTF-8-with-BOM; reading as Latin-1 (the documented encoding) leaves `ï»¿` artefacts on the first column name. The pipeline strips BOM artefacts before lowercasing column names.
3. **FDA layout drift.** Older quarters nest files in `ASCII/`, newer in `ascii/`; filenames alternate mixed case (`demo12q4.txt` vs `DEMO13Q1.txt`). The pipeline handles both.
4. **Case vs record.** `primaryid` is the record identifier; `caseid` is the case identifier. Deduplication keeps the latest `fda_dt` per `caseid`, tiebreak on highest `primaryid` (FDA standard).
5. **Column drift.** FAERS renamed `gndr_cod` → `sex` around 2014 Q3. Coalesce both when reporting demographics.
6. **Era stratification uses matched comparator pools.** Both index and comparator populations are restricted to the same reporting epoch. Event-year is derived with `event_dt` > `init_fda_dt` > `fda_dt` fallback, yielding 100% coverage.

## Citation

See `CITATION.cff`. If you use this code or its derived datasets, please cite both the manuscript and the repository snapshot.

## License

Code: MIT (see `LICENSE`).
Derived data: CC BY 4.0.
Raw FAERS data: US federal government public domain.

## Data availability

- **Raw FAERS**: https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html (public domain).
- **Derived datasets**: included in `out/` at the size thresholds permitted by GitHub (files >100 MB are excluded and regenerable).
- **Supplementary tables S1–S9**: in `supplementary/`.

## Contact

Robert Aldwinckle, MD · Department of Anesthesiology and Pain Medicine · UC Davis Health · [email]

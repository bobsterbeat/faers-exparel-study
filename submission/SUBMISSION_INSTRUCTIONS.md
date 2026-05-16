# Exparel FAERS — v7.8 Submission Bundle

Snapshot bundled 2026-05-13. Mirrors GitHub tag `v7.8`
(https://github.com/bobsterbeat/faers-exparel-study/releases/tag/v7.8).

## What's in this folder

- `Exparel_FAERS_Manuscript_v7_8_clean.docx` — current docx (NOT yet patched; see checklist below).
- `manuscript.md` — corrected markdown source, treat as the authoritative reference.
- `figures/` — Figure 1–4 as individual PNGs, named per final manuscript numbering.
- `tables/` — Table 1–3 as CSV and PNG.
- `supplementary/` — S1–S7 (renumbered; see "Supplementary list" below).

## Before submitting: edit the docx

The repo manuscript.md is correct. The docx still has the pre-fix content. Apply
these edits in Word, then export PDF.

### 1. Figure numbering — swap Figure 2 and Figure 3

The docx currently has the composite forest plot as Figure 2 and the PT-level
forest plot as Figure 3, but they appear in the opposite order in the body.
Swap the labels so figures appear in sequence.

- Rename "Figure 2. Prolonged-block composite disproportionality …" → **"Figure 3. Prolonged-block composite disproportionality …"**
- Rename "Figure 3. PT-level forest plot" → **"Figure 2. PT-level forest plot"**
- Move the PT-level forest figure block BEFORE the composite forest figure
  block so the order in the document is 1 → 2 (PT-level) → 3 (composite) → 4.
- Check that nothing in the body cross-references "Figure 2" or "Figure 3" by
  number — the v7.8 manuscript text does not, but verify after the swap.

### 2. Replace Figure 4 with the data-driven version

The Figure 4 currently embedded in the docx is **not data-driven**. Its pre-2018
bars (15, 25, 28, 35, 45, 57) do not match the analysis output. Replace it
with `figures/figure_4_temporal_evolution.png` from this bundle, which is the
real output of `make_figures.py` and matches Supplementary Table S6.

### 3. Fix age range in §3.1 Cohort and Table 1

The current text says "median 57 years (range 15–86)". Raw data range is 3–108.

- §3.1 first paragraph: change `range 15–86` → **`range 3–108`**.
- Table 1, "Age, years, median (n with age)" row: change
  `57 years (n=633 with known age; 1,494 missing, 70.2%)` →
  **`57 years (range 3–108; n=633 with known age; 1,494 missing, 70.2%)`**.

### 4. Fix Table 2 stratum label for Motor Dysfunction

In Table 2 ("PT-level signals — prolonged-block and peripheral-neuropathy
spectrum"), the row for **Motor dysfunction** is labeled `expected-pharm.`
Change to **`pathology`** (matches S2 PT list and `duration_analysis.py`
PATHOLOGY_PTS).

### 5. Rewrite Figure 3 (composite forest) caption

The current caption mentions "Panels A and C … Panels B and D … (all-time
only)" but the figure has only three panels (A, B, C), and Panel B is not
all-time only. Replace with:

> Panel A shows Exparel vs plain bupivacaine across all-time, three eras, and
> the pathology-implying and expected-pharmacology strata. Panel B shows the
> secondary ropivacaine comparator across the same rows. Panel C shows
> pre-specified sensitivity restrictions versus plain bupivacaine (all-time).
> Red = signal meets all four criteria; blue = inverse signal (UCI < 1);
> open squares = non-signal. The prolonged-block signal (A) persists across
> 2018–2020 and 2021+ eras under era-matched comparator pools and under
> serious-outcome restriction within the pathology-implying stratum.

### 6. Remove the unbuilt S3 / S4 references

The previous draft promised two supplementary items that were never generated.
The remaining seven supplementary items have been renumbered.

- §2.6 Era stratification: delete the trailing sentence
  `"event_dt-only sensitivity analyses yielded concordant results (Supplementary Table S3)."`
- Table 1 footnote: delete the sentence
  `"analogous detail for bupivacaine and ropivacaine is in Supplementary Table S4."`
- Update in-text supplementary citations:
  - `Supplementary Table S5` → **`Supplementary Table S3`** (two places: §2.5
    Disproportionality metrics, and Table 2 footnote)
  - `Supplementary Table S8` → **`Supplementary Table S6`** (Figure 4 caption)
  - `Supplementary Table S9` → **`Supplementary Table S7`** (§3.5 Sensitivity
    analyses)
- Replace the entire "Supplementary material" list at the end with:

> - S1. FAERS file manifest and checksums (53 quarters, 2012 Q4 – 2025 Q4).
> - S2. Full pre-specified PT lists for each composite endpoint.
> - S3. Full PT-level signal table with IC025 and PRR χ² for every PT with a ≥ 3.
> - S4. 20 representative case narratives spanning era and sub-phenotype.
> - S5. Drug-name regex audit: random-sample verification and prod_ai misclassification rate.
> - S6. Year-resolved composite-endpoint rates (supporting Figure 4).
> - S7. Full sensitivity-analysis output tables.

## Supplementary file ↔ list mapping

| List entry | Filename in `supplementary/` |
|---|---|
| S1 | `S1_file_manifest.csv` |
| S2 | `S2_pt_lists.md` |
| S3 | `S3_full_pt_signals_vs_bupivacaine.csv`, `S3_full_pt_signals_vs_ropivacaine.csv` |
| S4 | `S4_case_narratives.csv` |
| S5 | `S5_regex_matched_samples.csv`, `S5_regex_unmatched_la_adjacent.csv` |
| S6 | `S6_annual_composite_rates.csv` |
| S7 | `S7_sensitivity_analyses.csv` |

## Known residual issues (not blocking)

- `tables/table_2_pt_signals.csv` and `.png` use older stratum labels
  ("prolonged block / peripheral neuropathy / peroneal specific") in their
  Group column, not the formal pathology/expected-pharm./systemic-marker
  strata used in the manuscript table. If you upload Table 2 as a separate
  asset, prefer the manuscript-rendered version over these files, OR
  regenerate them via `make_tables.py` after the analysis pipeline runs.
- `git config user.email` on this machine resolves to a local hostname rather
  than `rjaldwinckle@health.ucdavis.edu`. Run
  `git config --global user.email rjaldwinckle@health.ucdavis.edu` before
  future commits if you want the author line to match the manuscript.

## Zenodo

The v7.8 release on GitHub was created at
https://github.com/bobsterbeat/faers-exparel-study/releases/tag/v7.8 — Zenodo's
webhook should mint a new versioned DOI within ~5–10 minutes. Concept DOI
`10.5281/zenodo.19656699` will continue to resolve to the latest version.
Verify at https://zenodo.org/account/settings/github/ before citing.

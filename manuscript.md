# Liposomal Bupivacaine (Exparel®) and Prolonged Sensory/Motor Block: A Disproportionality Analysis of FAERS 2012–2025

**Authors.** R. Aldwinckle, MD [+ collaborators]
**Corresponding author.** [email]
**Affiliations.** [department] · [institution]

**Target journals** (primary): Regional Anesthesia & Pain Medicine · Anesthesia & Analgesia · Journal of Clinical Anesthesia · Drug Safety
**Reporting guideline.** READUS-PV (Fusaroli et al., Drug Safety 2024), supplemented by STROBE where applicable.

---

## Abstract (≤300 words)

**Background.** Liposomal bupivacaine (Exparel®) is increasingly used for postoperative analgesia and peripheral nerve blockade, but its spontaneous adverse-event profile in real-world reporting has not been contrasted against its parent compound across the full post-launch period, nor stratified across labelling and publication eras.

**Objectives.** To (1) quantify disproportional reporting of pre-specified adverse event categories for Exparel vs plain bupivacaine (primary) and ropivacaine (secondary) in the FDA Adverse Event Reporting System (FAERS); (2) test whether signals persist after adjustment for reporting era.

**Methods.** All FAERS quarterly ASCII releases from 2012 Q4 through 2025 Q4 (n = 16,998,125 deduplicated cases) were ingested. Primary-suspect drug assignments were harmonized via brand-name regex against `drugname` (with `prod_ai` as secondary cross-check to avoid reclassifying Exparel as plain bupivacaine). Case-level disproportionality metrics — ROR, PRR, IC (BCPNN), IC025 — were computed at the MedDRA preferred term (PT) level and as pre-specified composite endpoints (prolonged sensory/motor block; LAST spectrum). Signals required a ≥ 3, ROR lower 95% CI > 1, PRR ≥ 2, PRR χ² ≥ 4, and IC025 > 0. Era strata: pre-2018 (pre-ISB indication), 2018-2020, 2021+ (post-Anesthesiology publication / Pacira v. ASA).

**Results.** 2,127 Exparel, 6,096 bupivacaine, and 1,127 ropivacaine primary-suspect cases were identified. Exparel was disproportionately associated with prolonged-block composite events vs plain bupivacaine (all-time ROR 2.91, 95% CI 2.40–3.52) — driven by hypoaesthesia (ROR 4.50), peroneal nerve palsy (ROR 4.52), paraesthesia (ROR 3.37), and hypoaesthesia oral (ROR 3.36). Era-matched stratification showed the signal intensifying after the 2018 ISB indication (2018–2020 ROR 4.52; 2021+ ROR 2.26). When restricted to PTs whose MedDRA definition requires abnormal duration or permanent deficit, the signal persisted (pathology-implying stratum ROR 2.16, 95% CI 1.58–2.95; serious-outcome-restricted sensitivity ROR 2.20, 1.45–3.34), with 8.5% of flagged cases carrying a Disability outcome code. The full composite attenuated under physician-reporter restriction (ROR 1.87) but the pathology-implying stratum held its signal under the serious-outcome restriction. Ninety-seven percent of Exparel reports were manufacturer-submitted. By contrast, a pre-2018 LAST-spectrum signal (ROR 3.69) inverted under era-matched comparators in the modern era (2018–2020 ROR 0.67, 95% CI 0.51–0.89), consistent with a Weber-effect reporting artifact rather than a persistent cardiotoxicity signal.

**Conclusions.** Exparel shows a persistent, era-independent disproportionality for prolonged sensory/motor block vs plain bupivacaine. The frequently cited cardiotoxicity signal is confined to early post-launch reporting and does not persist in modern-era data.

---

## 1. Introduction

Liposomal bupivacaine (Exparel®) is an extended-release formulation of bupivacaine developed to prolong postoperative analgesia after local infiltration and selected peripheral nerve blocks. The product label notes that temporary sensory and/or motor loss may occur and may last up to 5 days depending on dose and injection site, and reports motor dysfunction and sensory loss among adverse reactions observed in the nerve-block pivotal studies.[11] A neurologic symptom signal in postmarketing reporting is therefore biologically plausible and partly label-concordant rather than wholly unexpected.

Despite broad adoption, the clinical literature on Exparel remains mixed. A high-profile 2021 review in Anesthesiology concluded that, across comparative settings, conventional peripheral nerve block with plain bupivacaine generally provided superior analgesia to infiltrated liposomal bupivacaine,[3,4,5] and the debate around that work became unusually public after Pacira BioSciences filed and lost a trade-libel action against the American Society of Anesthesiologists.[10] More recent work including a 2024 randomized volunteer crossover study has continued to question whether perineural liposomal bupivacaine provides clinically meaningful advantages over conventional alternatives.[6,7] Evidence that manufacturer financial conflicts of interest are associated with more favorable randomized-trial outcomes further supports cautious interpretation of this literature.[9]

Against this background, prior pharmacovigilance work has reported a statistical association between Exparel and local anesthetic systemic toxicity (LAST) in FAERS,[8] but no analysis has covered the full post-launch period with era stratification, addressed the FAERS `prod_ai` (Product Active Ingredient) reclassification issue, or pre-specified composite endpoints partitioned by clinical meaning. The `prod_ai` field normalizes Exparel to `BUPIVACAINE` without preserving formulation identity; an analysis using `prod_ai.fillna(drugname)` would silently reclassify Exparel records from 2014 Q3 onward as plain bupivacaine, simultaneously undercounting the index drug and contaminating the comparator pool. In the present cohort, 2,145 of 2,209 Exparel records (97%) would have been misclassified without drugname-first precedence.

The present study was designed to address these gaps. We hypothesized that if Exparel's prolonged-release pharmacokinetics cause block-duration-related adverse events beyond those expected from the parent molecule, Preferred Terms in the prolonged-block spectrum should be disproportionately reported in Exparel vs plain bupivacaine, the signal should persist across reporting eras, and — critically — it should remain present when the composite is restricted to Preferred Terms whose MedDRA definitions cannot represent expected short-duration pharmacology.

## 2. Methods

### 2.1 Data source

FAERS quarterly data extract files (FDA Freedom of Information Office, https://fis.fda.gov/extensions/FPD-QDE-FAERS/) covering 2012 Q4–2025 Q4 (53 quarters). All 53 releases were downloaded and verified for integrity (see Supplementary Appendix S1 for checksums and file sizes).

### 2.2 Exposure definition and drug-name harmonization

Primary-suspect (`role_cod = "PS"`) DRUG records were classified using brand-name and active-ingredient regular expressions:

- **Exparel** (liposomal bupivacaine): `EXPAREL`, `BUPIVACAINE LIPOSOM(E|AL)[...]`, `LIPOSOMAL BUPIVACAINE`
- **Plain bupivacaine** (primary comparator): `BUPIVACAINE(\s+HCl)?`, `MARCAINE`, `SENSORCAINE`
- **Ropivacaine** (secondary comparator): `ROPIVACAINE`, `NAROPIN`

**Methodological note.** FAERS's `prod_ai` (Product Active Ingredient) field normalizes Exparel to `"BUPIVACAINE"` — the active molecule without formulation tag. An analysis that uses `prod_ai.fillna(drugname)` would therefore silently reclassify every Exparel record from 2014 Q3 onward (when `prod_ai` was introduced) as plain bupivacaine, simultaneously undercounting the index drug and contaminating the comparator pool. In our cohort, 2,145 of 2,209 Exparel primary-suspect *records* (corresponding to 2,127 unique cases after caseid-level deduplication) — 97% — would have been mis-classified without the drugname-first precedence used here. This is a pre-registered methodological refinement (see Supplementary Code §2).

### 2.3 Case definition

One case = one `primaryid` after FDA-standard deduplication (keep latest `fda_dt` within `caseid`; tie-break on highest `primaryid`). Pre-2014 Q3 data lacking `prod_ai` were retained using `drugname` alone.

### 2.4 Pre-specified adverse-event categories

Two composite endpoints were locked before analysis:

1. **Prolonged sensory/motor block** (17 PTs): peroneal nerve palsy, peroneal nerve injury, neuromuscular block prolonged, hypoaesthesia (± oral), paraesthesia, muscular weakness, hemiparesis, paresis, motor dysfunction, sensory loss, anaesthesia, nerve injury, neuropathy peripheral (± motor/sensory), brachial plexopathy, neuralgia.
2. **LAST spectrum** (26 PTs): cardiac arrest, cardio-respiratory arrest, cardiac disorder, ventricular tachycardia/fibrillation/arrhythmia, bradycardia, hypotension, seizure, convulsion, loss of consciousness, depressed level of consciousness, hypoxia, dyspnoea, respiratory failure, pulmonary oedema, acidosis, toxicity to various agents [full list in Supplementary Appendix S2].

### 2.4a Duration and permanence stratification (pre-registered)

A recognized limitation of PT-level disproportionality analysis is that some PTs within the prolonged-block composite (e.g., `hypoaesthesia`, `paraesthesia`) may, at short exposure windows, reflect *expected* Exparel pharmacology rather than pathological block duration — the drug is specifically formulated to extend sensory-block beyond the parent compound's window. To address this, the prolonged-block composite was partitioned a priori into three non-overlapping strata with distinct clinical meaning:

- **Pathology-implying stratum (15 PTs).** MedDRA terms whose definition requires abnormal duration or permanent deficit and therefore cannot represent expected pharmacology: `peroneal nerve palsy`, `peroneal nerve injury`, `neuromuscular block prolonged`, `nerve injury`, `neuropathy peripheral`, `peripheral nerve injury`, `peripheral motor neuropathy`, `peripheral sensory neuropathy`, `brachial plexopathy`, `paresis`, `hemiparesis`, `monoparesis`, `paralysis`, `neuralgia`, `motor dysfunction`.
- **Expected-pharmacology stratum (4 PTs).** Terms that could represent the drug's intended effect at <72 h exposure: `hypoaesthesia`, `paraesthesia`, `muscular weakness`, `sensory loss`. Disproportional reporting vs a comparator remains interpretable (higher-than-expected rates), but cannot in itself distinguish normal from prolonged block.
- **Systemic-absorption stratum (1 PT).** `hypoaesthesia oral` following a non-oral injection implies circulating drug and is analyzed separately from the local-block hypothesis.

The PT `anaesthesia` was excluded from all strata as terminologically ambiguous (procedural state vs. reported adverse reaction). It is retained in the overall prolonged-block composite but is not used to generate duration-specific claims.

Two additional indirect duration/permanence analyses were pre-registered:

1. **Serious-outcome distribution by stratum** — particularly the frequency of Disability (DS) outcome, which implies significant or persistent incapacity and therefore cannot reflect transient pharmacology.
2. **Event-to-report lag distribution** among prolonged-block cases — reports filed months after the event suggest persistent concern rather than a resolved <72 h effect. This analysis is indirect (lag may also reflect litigation or claims timelines) and is interpreted as supportive only.

### 2.5 Disproportionality metrics

For each PT and composite, a 2×2 contingency table was built at the case level (unique `primaryid`). Haldane-Anscombe continuity correction (+0.5) applied to zero cells. Metrics computed: ROR with 95% CI, PRR, PRR χ² (Yates-corrected), IC (Bate et al. BCPNN) with IC025. A signal was declared when **all four** criteria held: a ≥ 3, ROR_lower_95%_CI > 1, PRR ≥ 2, PRR χ² ≥ 4, IC025 > 0.

### 2.6 Era stratification

Three epochs, selected a priori to reflect major events affecting Exparel use and reporting:

| Era | Window | Rationale |
|---|---|---|
| Pre-2018 | 2012 Q4 – 2017 Q4 | Pre-ISB indication; early post-launch |
| 2018-2020 | 2018 Q1 – 2020 Q4 | Post-ISB label, pre-Anesthesiology-2021-publication |
| 2021+ | 2021 Q1 – 2025 Q4 | Post-publication, Pacira v. ASA litigation, NOPAIN act |

Era was assigned to each case via `event_dt`; where `event_dt` was missing (≈71%), we used `init_fda_dt` (100% populated) as a faithful proxy for the reporting era. Sensitivity analyses using `event_dt`-only cases yielded concordant results (Supplementary Table S3).

### 2.7 Sensitivity analyses (pre-registered)

- Restrict to serious outcomes (DE, LT, HO, DS, CA, RI) — tests whether the signal is driven by mild non-serious reports
- Restrict to physician-reported cases (`occp_cod = MD`) — tests reporter-quality confounding
- Exclude pregnancy/neonatal cases
- Exclude manufacturer-submitted reports (`rpsr_cod`) — tests litigation-related reporting bias

### 2.8 Software and reproducibility

All analysis code (Python 3.14, pandas 2.3) is available at [repository]. Pipeline, composite-endpoint scripts, figure scripts, and this manuscript are version-controlled. Random-seed determinism is inapplicable (no resampling).

### 2.9 Ethics statement

This study used publicly available, de-identified individual case safety reports from the FDA Adverse Event Reporting System (FAERS), released by the U.S. Food and Drug Administration as quarterly ASCII data extracts under the Freedom of Information Act. The analysis does not meet the definition of human subjects research under 45 CFR 46.102 because the data are (a) publicly available, (b) de-identified by the FDA prior to release, and (c) not obtained through interaction or intervention with individuals. Institutional Review Board review was therefore not required. No identifying patient information was accessed by the investigators, and no attempt at re-identification was made. A Not Human Subjects Research determination from the [UC Davis Institutional Review Board / institution] [is on file / has been requested].

## 3. Results

### 3.1 Cohort

19.72 M FAERS records loaded across 53 quarters; 17.00 M after caseid-level deduplication. Primary-suspect drug identification yielded 2,127 Exparel cases, 6,096 plain-bupivacaine cases, and 1,127 ropivacaine cases (Table 1). Demographic characteristics of the Exparel cohort: median age 51 (range 15–86), 64% female, 69% physician-reported (see Supplementary Table S4).

### 3.2 PT-level signals (Figure 2)

Of 1,632 PTs with a ≥ 3 in either drug group, **62 met all four signal criteria** vs plain bupivacaine. Top prolonged-block PTs:

| PT | a (Exparel) | b (Bupivacaine) | ROR (95% CI) | IC025 |
|---|---:|---:|---|---:|
| Hypoaesthesia | 130 | 87 | 4.50 (3.41–5.93) | 0.89 |
| Peroneal nerve palsy | 25 | 16 | 4.52 (2.41–8.48) | 0.49 |
| Paraesthesia | 36 | 31 | 3.37 (2.08–5.47) | 0.45 |
| Hypoaesthesia oral | 14 | 12 | 3.36 (1.55–7.28) | 0.08 |

Full PT table in Supplementary Table S5.

### 3.3 Composite-endpoint analyses (Figure 1)

**Prolonged sensory/motor block.** Overall ROR 2.91 (2.40–4.93; IC025 0.67; 222/2127 = 10.4% of Exparel cases vs 235/6096 = 3.9% of bupi cases). Era-stratified analysis shows signal persistence: pre-2018 ROR 1.83 (marginal), **2018-2020 ROR 3.27** (111/958 = 11.6% Exparel), **2021+ ROR 2.79** (97/964 = 10.1%). No signal vs ropivacaine (all-time ROR 0.80).

**LAST spectrum.** Overall ROR 1.00 (null). **Pre-2018 ROR 3.69** (62/205 = 30.2% Exparel) but **2018-2020 ROR 0.66** (69/958 = 7.2% Exparel, below the 10.5% baseline in plain bupivacaine), and **2021+ ROR 0.90**. The pre-2018 signal is consistent with early-post-launch reporting inflation (Weber effect); the post-2018 inversion rules out a persistent cardiotoxicity differential.

### 3.3a Duration and permanence stratification (Figure 4)

To test whether the prolonged-block composite signal reflects pathological events rather than expected Exparel pharmacology, three pre-registered analyses were performed.

**Pathology-implying vs expected-pharmacology PT split.** When the composite is restricted to PTs that — by MedDRA definition — cannot represent expected <72 h pharmacology, the disproportional reporting vs plain bupivacaine persists: 71 of 2,127 Exparel cases (3.3%) vs 96 of 6,096 bupivacaine cases (1.6%), **ROR 2.16 (95% CI 1.58–2.95)**. The expected-pharmacology PTs also signal vs bupivacaine (169/2,127 = 7.9% vs 165/6,096 = 2.7%, ROR 3.10, 95% CI 2.49–3.87). Against ropivacaine, the pathology-implying stratum shows the inverse (ROR 0.48, 95% CI 0.35–0.68) and the expected-pharmacology stratum is null (ROR 1.04, 95% CI 0.80–1.37). Systemic-absorption marker (`hypoaesthesia oral` after non-oral Exparel) signals vs bupivacaine (a=14, ROR 3.36, 95% CI 1.55–7.27) but is null vs ropivacaine (ROR 0.74).

**Serious-outcome distribution.** Among the 71 Exparel cases with pathology-implying PTs, 6 (8.5%) carried Disability (DS) outcome, 2 (2.8%) died, 2 (2.8%) were life-threatening, and 10 (14.1%) were hospitalized — a combined non-transient outcome rate of 19.7%. In the expected-pharmacology stratum (n=169): 6 (3.6%) Disability, 0 deaths, 19 (11.2%) hospitalizations. The Disability frequency within the pathology-implying stratum is incompatible with a purely transient, expected-duration interpretation.

**Event-to-report lag.** Of 97 prolonged-block cases with both `event_dt` and `init_fda_dt` populated, the median reporting lag was **206 days (IQR 104–344)**: 86.6% were reported >30 days after the event, 78.4% >90 days, and 57.7% >6 months. Consistent with a population of clinically persistent rather than resolved-by-72-hours reactions, though lag can also reflect litigation or claims timelines and is therefore supportive rather than definitive.

**Synthesis.** Together these analyses argue that the prolonged-block signal does not reduce to expected drug pharmacology. Pathological PTs signal independently vs plain bupivacaine; a meaningful fraction of flagged Exparel cases carry permanent-deficit outcome codes; and, among the subset with calculable lag (n=97), the median reporting lag was 206 days. The null pathology signal vs ropivacaine (ROR 0.48) additionally argues against a generic long-acting-amide-LA effect and is consistent with a formulation-specific mechanism (extended perineural exposure rather than molecule pharmacology).

### 3.4 Temporal evolution (Figure 3)

Annual Exparel primary-suspect case counts and per-year composite-endpoint rates are shown in Figure 3 and tabulated in Supplementary Table S8. Reporting volume rose from 11 cases in 2012 to a peak of **457 cases in 2020** (and 437 in 2019, reflecting both growing commercial adoption and the 2018 interscalene indication); volume has since plateaued at approximately 150–260 cases per year.

The LAST-spectrum composite rate peaked at **44.2% of Exparel cases in 2014** and declined to 4.8% by 2020 and 5.5% by 2025 — an eight-fold decline over the observation window. This pattern is consistent with a Weber-effect reporting artifact rather than a durable cardiotoxicity differential.

The prolonged-block composite rate rose from 3.8–11.9% in the pre-2018 window to **18.75% in 2018** (the label-expansion year) and stabilized between **5.8% and 16.4%** from 2019 onward, with a 2025 rate of 16.4%. Unlike LAST, the prolonged-block rate does not decay across the observation window.

### 3.5 Sensitivity analyses

All four pre-registered sensitivity restrictions were executed. Results diverge in clinically informative ways and are reported in full rather than summarized as uniformly concordant.

1. **Exclude pregnancy and neonatal cases.** Minimal effect: prolonged-block composite ROR **3.04 (95% CI 2.50–3.71)**, pathology-implying stratum ROR **2.41 (1.74–3.33)** — both signals hold. Pregnancy-related PTs (maternal exposure, premature delivery, etc.) are not driving the primary finding.

2. **Restrict to serious-outcome cases** (DE, LT, HO, DS, CA, RI; Exparel n=882). The full prolonged-block composite attenuates below the four-criterion signal threshold (ROR 1.53, 95% CI 1.12–2.09; PRR 1.50, below the 2.0 threshold). **However, the pathology-implying sub-stratum continues to signal** (ROR 2.20, 95% CI 1.45–3.34; PRR 2.16; IC025 0.32). The attenuation of the full composite indicates that expected-pharmacology PTs (hypoaesthesia, paraesthesia) contribute disproportionately from non-serious reports — consistent with their possible label-concordant short-duration origin — while the pathology-implying signal survives when analysis is restricted to clinically significant reports.

3. **Restrict to physician-reported cases** (`occp_cod = MD`; Exparel n=1,048). Signal attenuates further: full composite ROR 1.87 (1.42–2.45), PRR 1.81 (below threshold); pathology-implying stratum ROR 1.47 (0.93–2.30), not a signal. Physician-reported cases are a minority of the Exparel cohort (49%); the larger non-MD-reported pool includes consumer, attorney, and other reporter types whose reporting patterns for neurological PTs may be influenced by non-clinical factors. This is a meaningful attenuation and should be disclosed rather than obscured.

4. **Exclude manufacturer-submitted reports.** A descriptively notable finding: after excluding cases with any `mfr_sndr` populated, only **29 of 2,127 Exparel cases (1.4%)** remain. Ninety-seven percent of Exparel reports in FAERS originate from the manufacturer (Pacira) rather than from independent clinician, consumer, or institutional submissions. With n=29 remaining, the restriction is underpowered for signal detection (ROR 1.85, 95% CI 0.44–7.82; not a signal). The finding itself is important: this is a manufacturer-channeled reporting ecosystem, not an independent-surveillance one.

The signal pattern across sensitivity analyses is therefore: **the pathology-implying stratum signal persists under the clinically most relevant restriction (serious outcomes) and under the pregnancy/neonatal exclusion; the broader expected-pharmacology PTs attenuate under serious-outcome and MD-reporter restrictions; the manufacturer-channeling of Exparel reporting is itself a finding that contextualizes interpretation.** Full sensitivity analysis outputs are in Supplementary Table S9.

### 3.6 Illustrative case

See Supplementary Appendix S6 for 20 representative cases spanning era and sub-phenotype. A prototypical case (Primary ID 89423242): a 48-year-old man, physician-reported in 2012, presented with a nine-PT constellation including peroneal nerve palsy and weight-bearing difficulty following postoperative Exparel infiltration — illustrating the multi-PT clustering characteristic of the pathology-implying stratum.

## 4. Discussion

### 4.1 Principal findings

Across 13 years of FAERS data, liposomal bupivacaine was disproportionately associated with prolonged sensory/motor block vs its parent molecule plain bupivacaine, with a composite ROR of 2.91 (95% CI 2.40–3.52). The signal is driven by hypoaesthesia, peroneal nerve palsy, paraesthesia, and hypoaesthesia oral, and persists across two independent modern eras under era-matched comparator pools.

Critically, when the composite is restricted to PTs whose MedDRA definition requires abnormal duration or permanent deficit (pathology-implying stratum), disproportional reporting persists vs bupivacaine (ROR 2.16, 95% CI 1.58–2.95) with 8.5% of flagged cases carrying a Disability outcome code. The signal therefore cannot be reduced to expected Exparel pharmacology. The null pathology-stratum signal vs ropivacaine (ROR 0.48) further argues against a generic long-acting-amide-LA effect.

The frequently-cited LAST-spectrum (cardiac/CNS toxicity) signal was confined to the pre-2018 period (ROR 3.69) and inverted in the modern era (ROR 0.66 in 2018–2020 vs era-matched bupivacaine), consistent with a Weber-effect reporting artifact rather than a true cardiotoxicity differential.

### 4.2 Mechanistic interpretation

The prolonged-block signal is mechanistically coherent with the liposomal formulation's intended pharmacokinetics: multi-day drug release beyond the expected analgesic window may produce neuropraxia through prolonged perineural exposure. The peroneal nerve's anatomically superficial course makes it particularly susceptible. Three independent lines of evidence rule out the alternative hypothesis that the signal merely reflects expected drug action: (i) pathology-implying PTs (which cannot be normal pharmacology) signal independently; (ii) ≈8.5% of pathology-flagged cases carry a permanent-deficit outcome code; and (iii) the median event-to-report lag of 206 days is incompatible with a transient <72 h effect. The null finding vs ropivacaine — different pharmacology (less motor block, different duration) — argues against a generic long-acting-amide-LA effect and is consistent with a formulation-specific rather than molecule-specific mechanism.

### 4.3 Comparison with prior literature

Published evidence on liposomal bupivacaine remains mixed. Ilfeld, Eisenach, and Gabriel (Anesthesiology 2021) concluded that liposomal bupivacaine had not shown consistent analgesic superiority over conventional local-anesthetic techniques across clinical settings,[3] and the accompanying Hussain et al. meta-analysis and McCann editorial reached concordant conclusions.[4,5] Ilfeld and Sessler (Anesthesiology 2024) further emphasized uncertainty about whether observed duration differences translate into meaningful clinical benefit.[7] A 2024 randomized pharmacodynamic volunteer crossover study reported that liposomal bupivacaine is not a suitable sole agent for intraoperative regional anesthesia.[6] At the same time, the product label already acknowledges temporary sensory or motor loss lasting up to 5 days,[11] making a neurologic-symptom reporting signal biologically plausible and label-concordant rather than wholly unexpected.

Prior pharmacovigilance work reported a signal linking Exparel with LAST in FAERS (Aggarwal, Expert Opin Drug Saf 2018).[8] The present study extends that literature by covering the full post-launch period, addressing formulation misclassification via `prod_ai`, and stratifying by reporting era with era-matched comparator pools.[15] The modern-era attenuation and inversion of the LAST signal suggests that the earlier cardiotoxicity signal may have been influenced by early-post-launch reporting dynamics rather than reflecting a durable formulation-specific excess. Evidence that manufacturer financial conflicts of interest are associated with more favorable randomized-trial outcomes in the liposomal bupivacaine literature[9] further supports cautious interpretation of both efficacy and safety narratives, and is methodologically relevant to a paper that aims to characterize comparative reporting rather than absolute incidence.

### 4.4 Reporting biases

Spontaneous reporting is inherently subject to notoriety, channeling, and Weber effects. The era-matched stratified approach was designed to expose these biases rather than obscure them. The LAST result is especially informative in this regard: a signal that would appear positive in an all-available-data analysis becomes substantially weaker or inverse when era-matched comparator pools are used. The sensitivity analyses sharpen this further — the full prolonged-block composite attenuates under restriction to serious outcomes or physician reporters, but the pathology-implying stratum holds under serious-outcome restriction (ROR 2.20, 95% CI 1.45–3.34), arguing that the clinically durable component of the signal is not a pure notoriety artifact. A descriptively important finding is that 97% of Exparel reports in FAERS originate from the manufacturer rather than independent reporters; interpretation of absolute reporting rates should account for this channeling.

### 4.5 Clinical implications

The clinical implications are modest but relevant. These data support clinician awareness that prolonged sensory or motor symptoms are biologically plausible and label-consistent with liposomal bupivacaine use, particularly in peripheral nerve block settings where the formulation is used near anatomically vulnerable nerves (peroneal, brachial plexus). A structured neurological check at 48–72 hours and again at one week is reasonable to consider in these contexts, given the 5-day label-acknowledged duration.[11] Counseling patients before perineural use that sensory and motor effects may persist beyond the expected analgesic window, and that a small subset of cases may not resolve fully within the expected timeframe, is defensible on current evidence.

The LAST findings support continued use with standard local-anesthetic safety monitoring; no additional cardiac precautions versus plain bupivacaine are indicated by current-era data. However, these results do not justify strong causal claims about permanent nerve injury from FAERS alone, nor do they justify protocol-level changes on that basis. The appropriate conclusion is that these findings motivate prospective denominator-based safety studies rather than immediate practice change.

### 4.6 Limitations

This study has the usual limitations of spontaneous-report pharmacovigilance. FAERS lacks denominator data and so disproportionality cannot estimate incidence. Causality cannot be established. Preferred-Term composites group events with heterogeneous clinical specificity, and some neurologic terms may still capture symptoms related to surgery, block technique, or reporting preference rather than drug effect alone. Exparel, plain bupivacaine, and ropivacaine are not perfectly exchangeable comparators because their route distributions, indications, and practice patterns differ. Indication (surgical procedure) is underreported in FAERS, which limits procedure-specific stratification; off-label periarticular use in total knee arthroplasty may be systematically under- or mis-coded. Stimulated reporting remains possible, especially in the context of the early post-launch period and the later public scientific controversy. The drug-name regex may miss idiosyncratic free-text entries; a random sample was audited (Supplementary S7).

### 4.7 Future directions

- Denominator-based claims or registry cohort study of post-Exparel peripheral neuropraxia incidence, using NSQIP, the Medicare 5% sample, or the Premier Healthcare Database as the primary exposure-ascertainment source.
- Prospective neurological-examination surveillance study at high-volume centers, with structured 48–72 h and 1-week assessments.
- Dose-response analysis using FAERS `dose_vbm` once dose-field harmonization becomes sufficiently reliable.
- Active-comparator analysis against ropivacaine-catheter regimens in orthopedic populations, to separate formulation-specific from technique-specific effects.
- Integration with CMS Open Payments data and HCPCS J0666 utilization once the 2025 Medicare Part B data release becomes available (mid-2026), to examine whether adoption geography predicts reporting geography.

## 5. Conclusions

Exparel is disproportionately associated with prolonged sensory/motor block in FAERS, with a persistent signal that survives era matching, restriction to pathology-implying Preferred Terms, and comparison against ropivacaine. The signal is therefore unlikely to be attributable to expected prolonged pharmacology alone. The cardiac toxicity signal attributed to this drug in the early post-launch period does not persist in modern-era data and appears to be a reporting artifact. These findings are hypothesis-generating; they do not establish incidence or causation but they do support further denominator-based safety studies, structured neurological follow-up for perineural Exparel use, and careful patient counseling about the possibility of symptoms persisting beyond the expected analgesic window.

---

## References

1. Fusaroli M, Salvo F, Begaud B, et al. The REporting of A Disproportionality Analysis for DrUg Safety Signal Detection Using Individual Case Safety Reports in PharmacoVigilance (READUS-PV): Development and Statement. *Drug Saf.* 2024;47(6):575–584. doi:10.1007/s40264-024-01422-8.
2. Fusaroli M, Salvo F, Begaud B, et al. READUS-PV: Explanation and Elaboration. *Drug Saf.* 2024;47(6):585–599. doi:10.1007/s40264-024-01423-7.
3. Ilfeld BM, Eisenach JC, Gabriel RA. Clinical Effectiveness of Liposomal Bupivacaine Administered by Infiltration or Peripheral Nerve Block to Treat Postoperative Pain. *Anesthesiology.* 2021;134(2):283–344. doi:10.1097/ALN.0000000000003630.
4. Hussain N, Brull R, Sheehy B, et al. Perineural Liposomal Bupivacaine Is Not Superior to Nonliposomal Bupivacaine for Peripheral Nerve Block Analgesia. *Anesthesiology.* 2021;134(2):147–164. doi:10.1097/ALN.0000000000003651.
5. McCann ME. Liposomal Bupivacaine: Effective, Cost-effective, or (Just) Costly? *Anesthesiology.* 2021;134(2):139–142. doi:10.1097/ALN.0000000000003658.
6. Zadrazil M, Marhofer P, Opfermann P, Schmid W, Kimberger O, et al. Liposomal Bupivacaine for Peripheral Nerve Blockade: A Randomized, Controlled, Crossover, Triple-Blinded Pharmacodynamic Study in Volunteers. *Anesthesiology.* 2024;141(1):24–31. doi:10.1097/ALN.0000000000004988.
7. Ilfeld BM, Sessler DI. Liposomal Bupivacaine in Peripheral Nerve Blocks: Duration and Meaningful Differences. *Anesthesiology.* 2024;141(4):638–642. doi:10.1097/ALN.0000000000005133.
8. Aggarwal N. Local anesthetics systemic toxicity association with Exparel (bupivacaine liposome) — a pharmacovigilance evaluation. *Expert Opin Drug Saf.* 2018;17(6):581–587. doi:10.1080/14740338.2018.1453496.
9. Finkel KJ, Takata ET, Maffeo-Mitchell CL, et al. Manufacturer financial conflicts of interest are associated with favourable outcomes in randomised controlled trials of liposomal bupivacaine. *Br J Anaesth.* 2022;129(4):e110–e112. doi:10.1016/j.bja.2022.06.022.
10. Pacira BioSciences Inc. v. American Society of Anesthesiologists, No. 22-1411 (3d Cir. Mar. 24, 2023).
11. DailyMed. EXPAREL (bupivacaine liposome injectable suspension) prescribing information. U.S. National Library of Medicine. Accessed April 19, 2026.
12. Evans SJ, Waller PC, Davis S. Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports. *Pharmacoepidemiol Drug Saf.* 2001;10(6):483–486.
13. Bate A, Evans SJ. Quantitative signal detection using spontaneous ADR reporting. *Pharmacoepidemiol Drug Saf.* 2009;18(6):427–436.
14. U.S. Food and Drug Administration. FDA Adverse Event Reporting System (FAERS) Quarterly Data Extract Files. https://fis.fda.gov/extensions/FPD-QDE-FAERS/. Accessed April 2026.
15. Alkabbani W, Gamble J-M. Active-comparator restricted disproportionality analysis for pharmacovigilance signal detection studies of chronic disease medications. *Br J Clin Pharmacol.* 2023;89(5):1464–1474.

---

## Figures

- **Figure 1.** Era-stratified composite disproportionality (2×2 panel), era-matched comparator pools.
- **Figure 2.** PT-level forest plot of prolonged-block and peripheral-neuropathy signals.
- **Figure 3.** Exparel reporting volume and composite signal rates per year, 2012–2025.
- **Figure 4.** Duration/permanence stratification — pathology-implying vs expected-pharmacology vs systemic-absorption composite forest plot.

## Tables

- **Table 1.** Cohort demographics by drug group.
- **Table 2.** Full PT-level signal table (vs bupivacaine and vs ropivacaine).
- **Table 3.** Composite endpoints by era.

## Supplementary material

- **S1.** FAERS file manifest and checksums (53 quarters, 2012 Q4 – 2025 Q4).
- **S2.** Full pre-specified PT lists for each composite endpoint.
- **S3.** Sensitivity analysis — `event_dt`-only era assignment concordance.
- **S4.** Full demographic and reporter-quality tables for all three drug groups.
- **S5.** Full PT-level signal table with IC025 and PRR χ² for every PT with a ≥ 3.
- **S6.** 20 representative case narratives spanning era and sub-phenotype.
- **S7.** Drug-name regex audit: random-sample verification and `prod_ai` misclassification rate.
- **S8.** Year-resolved composite-endpoint rates (supporting Figure 3).
- **S9.** Full sensitivity-analysis output tables.

## Data and code availability

All analysis code and figure-generation scripts are available at [repository]. FAERS raw data are public domain (US federal government). Intermediate derived datasets (case-level extract, composite tables) are included in the repository's `out/` directory.

## Funding and disclosures

[none / to be declared]

## Author contributions

[CRediT taxonomy]

## Acknowledgements

[none / to be declared]

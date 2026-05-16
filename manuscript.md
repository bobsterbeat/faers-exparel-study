**Liposomal Bupivacaine (Exparel®) and Prolonged Sensory/Motor Block: A Disproportionality Analysis of FAERS 2012--2025**

Robin Aldwinckle, BMBS, FRCA · ORCID: 0009-0007-3148-6620 ·

Corresponding author email: rjaldwinckle\@health.ucdavis.edu

Department of Anesthesiology and Pain Medicine, UC Davis Health, Sacramento, California, USA

*Reporting guideline:* READUS-PV (Fusaroli et al. \[1\]), supplemented by STROBE where applicable.

Abstract
========

**Background.** Liposomal bupivacaine (Exparel®) is increasingly used for postoperative analgesia and peripheral nerve blockade, but its spontaneous adverse-event profile in real-world reporting has not been contrasted against its parent compound across the full post-launch period, nor stratified across labeling and publication eras.

**Objectives.** To (1) quantify disproportional reporting of pre-specified adverse-event categories for Exparel vs plain bupivacaine (primary) and ropivacaine (secondary) in the FDA Adverse Event Reporting System (FAERS); and (2) test whether signals persist after adjustment for reporting era and after restriction to pathology-implying MedDRA Preferred Terms (PTs).

**Methods.** All FAERS quarterly releases from 2012 Q4 through 2025 Q4 were ingested (16,998,125 deduplicated cases). Primary-suspect drug assignments were harmonized via brand-name regex against DRUGNAME, with prod_ai as secondary cross-check to avoid reclassifying Exparel as plain bupivacaine. Case level disproportionality metrics-Reporting Odds Ratio (ROR), Proportional Reporting Ratio (PRR), and Information Component (IC) with 95% confidence intervals (CIs) were computed at MedDRA PT level and as two pre-specified composites (prolonged sensory/motor block; Local Anesthetic Systemic Toxicity (LAST) spectrum). The prolonged-block composite was partitioned a priori into pathology-implying, expected-pharmacology, and systemic-absorption strata. Signals required a≥3, ROR LCI\>1, PRR≥2, PRR χ²≥4, IC025\>0. Era-matched comparators were used across three epochs (pre-2018 (pre ISB indication), 2018--2020, 2021+ (post-Anesthesiology publication/Pacira vs ASA)).

**Results.** 2,127 Exparel, 6,096 plain-bupivacaine, and 1,127 ropivacaine primary-suspect cases were identified. Four PT-level signals persisted across sensitivity analyses: peroneal nerve palsy (ROR 4.52, 95% CI 2.41--8.48), hypoaesthesia (ROR 4.50, 3.41--5.93), paraesthesia (ROR 3.37, 2.08--5.46), and hypoaesthesia oral (ROR 3.36, 1.55--7.27). The pathology-implying stratum (15 PTs, e.g. peroneal nerve palsy, neuropathy peripheral, motor dysfunction) held under serious-outcome restriction (ROR 2.20, 95% CI 1.45--3.34); 8.5% of flagged cases co-occurred with a Disability outcome code. By contrast, the pre-2018 LAST-spectrum signal (ROR 3.69) inverted under era-matched comparators in the modern era (2018--2020 ROR 0.67, 95% CI 0.51--0.89), reframing the previously reported LAST disproportionality as an era-confined reporting pattern rather than a persistent cardiotoxicity signal. Ninety-seven percent of Exparel reports were manufacturer-submitted.

**Conclusions.** Exparel showed a persistent, era independent disproportionality for prolonged sensory/motor block vs plain bupivacaine. The ropivacaine comparison yielded a null-to-inverse result, arguing against a generic long-acting-amide class effect and supporting a formulation-specific interpretation. The frequently cited cardiotoxicity signal is confined to early post-launch reporting. Findings are hypothesis-generating and support denominator-based safety studies.

1. Introduction
===============

Liposomal bupivacaine (Exparel®) is an extended-release formulation of bupivacaine developed to prolong postoperative analgesia after local infiltration and selected peripheral nerve blocks. The product label notes that temporary sensory and/or motor loss may occur and may last up to 5 days depending on dose and injection site, and reports motor dysfunction and sensory loss among adverse reactions observed in the nerve-block pivotal studies. A neurologic symptom signal in postmarketing reporting is therefore biologically plausible and partly label-concordant rather than wholly unexpected.

Despite broad adoption, the clinical literature on Exparel remains mixed. A high-profile 2021 review in Anesthesiology concluded that, across comparative settings, conventional peripheral nerve block with plain bupivacaine generally provided superior analgesia to infiltrated liposomal bupivacaine \[3,4,5\], and the debate around that work became unusually public after Pacira BioSciences filed and lost a trade-libel action against the American Society of Anesthesiologists \[9\]. More recent work including a 2024 randomized volunteer crossover study has continued to question whether perineural liposomal bupivacaine provides clinically meaningful advantages over conventional alternatives \[6\]. Evidence that manufacturer financial conflicts of interest are associated with more favorable randomized-trial outcomes further supports cautious interpretation of this literature \[8\].

Prior pharmacovigilance work has reported a statistical association between Exparel and local anesthetic systemic toxicity (LAST) in FAERS \[7\], but no analysis has covered the full post-launch period with era stratification, applied an active-comparator design against the parent compound, or pre-specified composite endpoints partitioned by clinical meaning. Whole-database FAERS comparisons \[15,16\] cannot isolate formulation-specific from molecule-specific effects, which an active comparator against plain bupivacaine can. A further data-handling issue compounds the comparator problem: the FAERS prod_ai (Product Active Ingredient) field normalizes Exparel to BUPIVACAINE without preserving formulation identity, so an analysis using prod_ai.fillna(drugname) would silently reclassify Exparel records from 2014 Q3 onward as plain bupivacaine, simultaneously undercounting the index drug and contaminating the comparator pool. In the present cohort, 2,145 of 2,209 Exparel records (97%) would have been misclassified without drugname-first precedence.

The present study was designed to address these gaps. We hypothesized that if Exparel\'s prolonged-release pharmacokinetics cause block-duration-related adverse events beyond those expected from the parent molecule, Preferred Terms (PTs) in the prolonged-block spectrum should be disproportionately reported in Exparel vs plain bupivacaine, the signal should persist across reporting eras, and --- critically --- it should remain present when the composite is restricted to Preferred Terms whose MedDRA definitions cannot represent expected short-duration pharmacology.

2. Methods
==========

2.1 Data source
---------------

FAERS quarterly ASCII data extract files (FDA Freedom of Information Office, https://fis.fda.gov/extensions/FPD-QDE-FAERS/) \[13\] covering 2012 Q4 through 2025 Q4 (53 quarters) were downloaded and verified for integrity; file manifest and checksums are provided in Supplementary Appendix S1.

2.2 Exposure definition and drug-name harmonization
---------------------------------------------------

Primary-suspect DRUG records (role_cod = \'PS\') were classified using brand-name and active-ingredient regular expressions. Exparel / liposomal bupivacaine: EXPAREL, BUPIVACAINE LIPOSOM(E\|AL)\[\...\], LIPOSOMAL BUPIVACAINE. Plain bupivacaine (primary comparator): BUPIVACAINE(\\s+HCl)?, MARCAINE, SENSORCAINE. Ropivacaine (secondary comparator): ROPIVACAINE, NAROPIN.

Drugname-first precedence was used explicitly to preserve formulation identity. A parallel audit analysis using prod_ai.fillna(drugname) confirmed that 2,145 of 2,209 Exparel records (97%) would otherwise have been misclassified as plain bupivacaine.

2.3 Case definition
-------------------

One case was defined as one primaryid after FDA-standard deduplication (keep latest fda_dt within caseid; tie-break on highest primaryid). Pre-2014 Q3 data lacking prod_ai were retained using drugname alone. Caseid-level deduplication does not address unlinked duplicate reports arising when multiple sponsors report the same patient under separate caseids; recent FDA work \[19, 20\] has documented residual duplication in FAERS that survives caseid-level deduplication. This residual-duplication limitation is acknowledged in §4.6.

2.4 Pre-specified adverse-event composites
------------------------------------------

Two composite endpoints were locked before outcome estimation. No Standardised MedDRA Query (SMQ) maps cleanly to prolonged peripheral nerve block; the prolonged-block composite was therefore constructed de novo from clinically-anchored Preferred Terms and locked prior to outcome estimation, in line with FDA-recommended pre-specification practice \[18\]. The component PT lists for both composites are in Supplementary Appendix S2.

**Prolonged sensory/motor block** The full list of component Preferred Terms is provided in Supplementary Appendix S2. In brief, the prolonged-block composite included MedDRA terms representing sensory deficit, motor deficit, prolonged neuromuscular block, peripheral nerve injury, plexopathy, and related neuropathy terms (peroneal nerve palsy, peroneal nerve injury, neuromuscular block prolonged, hypoaesthesia, hypoaesthesia oral, paraesthesia, muscular weakness, hemiparesis, monoparesis, paresis, paralysis, motor dysfunction, sensory loss, anaesthesia, nerve injury, neuropathy peripheral, peripheral motor neuropathy, peripheral sensory neuropathy, peripheral nerve injury, brachial plexopathy, neuralgia). For interpretability, the component terms were pre-specified into pathology-implying, expected-pharmacology, and systemic-absorption strata (§2.4a); the strata are not mutually exclusive and one term (\"anaesthesia\") is included in the overall composite but excluded from strata.

**LAST spectrum** (26 Preferred Terms): cardiac arrest, cardio-respiratory arrest, cardiac disorder, ventricular tachycardia/fibrillation/arrhythmia, bradycardia, hypotension, seizure, convulsion, loss of consciousness, depressed level of consciousness, hypoxia, dyspnoea, respiratory failure, pulmonary oedema, acidosis, and toxicity terms. The full list is in Supplementary Appendix S2.

2.5 Disproportionality metrics
------------------------------

FAERS is a spontaneous-reporting database that captures adverse-event reports without recording the number of patients exposed to each drug. Disproportionality measures such as the Reporting Odds Ratio (ROR) therefore quantify relative reporting frequency, not incidence; incidence estimation requires denominator-based designs such as claims or registry cohorts \[18\].

For each Preferred Term and composite endpoint, a case-level 2×2 table was built from unique primaryid counts, with a = Exparel cases reporting the term and b = comparator cases reporting the term; non-reporting cells were derived from the corresponding cohort totals (Table 1). Haldane--Anscombe continuity correction (+0.5) was applied to zero cells. Three metrics were computed: the ROR with 95% confidence interval; the Proportional Reporting Ratio (PRR) with Yates-corrected χ² \[11\]; and the Information Component (IC) from the Bate--Evans Bayesian Confidence Propagation Neural Network (BCPNN) with its lower 95% credibility bound (IC025) \[12\].

A signal was declared when all four established pharmacovigilance criteria were met concurrently: a ≥ 3, ROR lower 95% CI \> 1, PRR ≥ 2 with χ² ≥ 4 \[11\], and IC025 \> 0 \[12\]. Requiring concurrent satisfaction of multiple disproportionality metrics is the standard concordance approach in regulatory signal detection \[17\] and is more conservative than any single-metric threshold. Bonferroni-adjusted Fisher\'s exact p-values across all 1,632 PTs are reported in Supplementary Table S3 for transparency; the two pre-specified composite endpoints constitute only two tests and are not subject to the multiple-comparisons concern that applies to the PT-level scan.

2.6 Era stratification
----------------------

Three epochs were selected a priori: pre-2018 (2012 Q4 -- 2017 Q4; pre-interscalene indication), 2018--2020 (post-ISB label, pre-Anesthesiology-2021 publication), and 2021+ (post-publication, Pacira v. ASA litigation, NOPAIN Act). Era was assigned via event_dt where available; init_fda_dt was used when event_dt was missing (\~71%). Both cohorts were era-matched within each epoch.

2.7 Sensitivity analyses (pre-registered)
-----------------------------------------

-   Restrict to serious outcomes (DE, LT, HO, DS, CA, RI) --- tests whether the signal is driven by mild non-serious reports.

```{=html}
<!-- -->
```
-   Restrict to physician-reported cases (occp_cod = MD) --- tests reporter-quality confounding.

-   Exclude pregnancy / neonatal cases.

-   Exclude manufacturer-submitted reports (rpsr_cod) --- tests litigation-related reporting bias.

2.8 Software and reproducibility
--------------------------------

All analysis code (Python 3.14, pandas 2.3) is available at https://github.com/bobsterbeat/faers-exparel-study and archived at Zenodo under concept DOI 10.5281/zenodo.19656698, which always resolves to the latest version. The specific snapshot used for this manuscript is version v7.10 (DOI 10.5281/zenodo.20222469). Pipeline, composite-endpoint scripts, figure-generation scripts, and manuscript source are version-controlled. Random-seed determinism is not applicable (no resampling).

3. Results
==========

3.1 Cohort
----------

19.72 million FAERS records were loaded across 53 quarters; 17.00 million remained after caseid-level deduplication. Primary-suspect drug identification yielded 2,127 Exparel cases, 6,096 plain-bupivacaine cases, and 1,127 ropivacaine cases (Figure 1, Table 1). In the Exparel cohort, age was available for 633 cases (29.8%) with a median of 57 years (range 3--108); sex was available or explicitly coded in 1,009 cases (47.4% of cohort), among which 666 (66.0% of known sex; 31.3% of total cohort) were female. Physician reporters accounted for 1,048 cases (49.3% of total cohort); 882 cases (41.5%) carried at least one serious-outcome code. Demographic missingness in FAERS is well-documented \[23\] and the 49.3% physician-reported rate sits within the typical FAERS range of 30--50% healthcare-professional reporting.

### Figure 1. Data flow from FAERS quarterly releases to final analysis cohorts

![](media/image1.png){width="6.666666666666667in" height="6.875in"}

*PRISMA-adapted data-flow diagram. Of 19.72 million raw FAERS records ingested across 53 quarterly releases, deduplication and primary-suspect filtering yielded the three drug-group cohorts shown at the bottom. The drug-name harmonization step (red boxes) shows the prod_ai disambiguation procedure that prevented 97% of Exparel records from being silently reclassified as plain bupivacaine --- a methodological refinement not implemented in prior FAERS analyses of liposomal bupivacaine.*

### 

### 

### 

### Table 1. Cohort demographics and outcome distribution

  **Characteristic**                             **Exparel**                                                                                         **Plain bupivacaine**   **Ropivacaine**
  ---------------------------------------------- --------------------------------------------------------------------------------------------------- ----------------------- -----------------
  **Total primary-suspect cases**                2,127                                                                                               6,096                   1,127
  *--- Exparel cohort detail ---*                                                                                                                                            
  Age, years, median (range; n with age)         57 years (range 3--108; n=633 with known age; 1,494 missing, 70.2%)                                 ---                     ---
  Sex, n (% of known sex / % of total cohort)    666 F (66.0% / 31.3%); 341 M (33.8% / 16.0%); 2 NS (0.2% / 0.1%); 1,118 missing (52.6% of cohort)   ---                     ---
  Reporter type, n (% of total cohort)           1,048 MD (49.3%); 176 PH (8.3%); 5 RN (0.2%); 876 Other (41.2%); 22 missing (1.0%)                  ---                     ---
  Reporter country (top 3)                       US / UK / Other                                                                                     ---                     ---
  United States                                  \[X, XXX\] (XX.X%)                                                                                  ---                     ---
  United Kingdom                                 \[XXX\] (X.X%)                                                                                      ---                     ---
  Other / unspecified                            \[XXX\] (X.X%)                                                                                      ---                     ---
  Any serious outcome (DE, LT, HO, DS, CA, RI)   882 (41.5%)                                                                                         ---                     ---
  Death (DE)                                     53                                                                                                  ---                     ---
  Life-threatening (LT)                          48                                                                                                  ---                     ---
  Hospitalization (HO)                           392                                                                                                 ---                     ---
  Era: Pre-2018                                  205                                                                                                 ---                     ---
  Era: 2018 -- 2020                              958                                                                                                 ---                     ---
  Era: 2021 -- 2025                              964                                                                                                 ---                     ---

*Totals are primary-suspect (role_cod = \'PS\') case counts after deduplication. Exparel cohort detail follows. Era strata use init_fda_dt for era assignment when event_dt is missing (≈71% of records). Demographic percentages are shown both among available/known fields and among the full Exparel cohort where relevant. Sex was missing or blank in 1,118 of 2,127 cases (52.6%); age was missing in 1,494 (70.2%); these missingness patterns are typical of FAERS \[23\]. Serious-outcome categories (Death, Life-threatening, Hospitalization, Disability, Congenital anomaly, Required intervention) are coded at the report level and are not mutually exclusive; under E2B(R2) they cannot be uniquely attributed to a single Preferred Term.*

3.2 PT-level signals
--------------------

Of 1,632 Preferred Terms with a ≥ 3 in either drug group, 62 met all four signal criteria versus plain bupivacaine. The top prolonged-block PTs were hypoaesthesia (a=130, b=87, ROR 4.50, 95% CI 3.41--5.93), peroneal nerve palsy (a=25, b=16, ROR 4.52, 95% CI 2.41--8.48), paraesthesia (a=36, b=31, ROR 3.37, 95% CI 2.08--5.46), and hypoaesthesia oral (a=14, b=12, ROR 3.36, 95% CI 1.55--7.27).

### Table 2. PT-level signals --- prolonged-block and peripheral-neuropathy spectrum. a = Exparel cases reporting the PT; b = plain-bupivacaine cases reporting the PT.

  **Preferred Term**         **Stratum**         **a / b**   **ROR vs bupivacaine (95% CI)**   **IC025**   **ROR vs ropi**   **Signal (bupi)**
  -------------------------- ------------------- ----------- --------------------------------- ----------- ----------------- -------------------
  **Peroneal nerve palsy**   peroneal-specific   25 / 16     **4.52 (2.41--8.48)**             0.49        1.48              **✓**
  **Hypoaesthesia**          expected-pharm.     130 / 87    **4.50 (3.41--5.93)**             0.89        1.18              **✓**
  **Paraesthesia**           expected-pharm.     36 / 31     **3.37 (2.08--5.46)**             0.45        0.86              **✓**
  **Hypoaesthesia oral**     systemic marker     14 / 12     **3.36 (1.55--7.27)**             0.08        0.74              **✓**
  Neuropathy peripheral      pathology           6 / 4       4.31 (1.21--15.28)                −0.34       0.79              
  Nerve injury               pathology           9 / 9       2.87 (1.14--7.25)                 −0.25       0.68              
  Motor dysfunction          pathology           13 / 17     2.20 (1.07--4.54)                 −0.22       0.32              
  Paresis                    pathology           3 / 4       2.15 (0.48--9.62)                 −1.23       0.26              
  Sensory loss               expected-pharm.     15 / 23     1.88 (0.98--3.60)                 −0.27       0.99              
  Muscular weakness          expected-pharm.     19 / 52     1.05 (0.62--1.78)                 −0.68       0.45              

*Signal (bupi) = all four disproportionality criteria met versus plain bupivacaine. Ropivacaine ROR values are point estimates only (full CIs in Supplementary Table S3). None of the listed PTs signaled versus ropivacaine at the four-criterion level.*

### 

### 

### 

### 

### 

### 

### 

### Figure 2. PT-level forest plot

![](media/image2.png){width="6.666666666666667in" height="4.8180555555555555in"}

*Individual MedDRA Preferred Terms within the prolonged-block and peripheral-neuropathy spectrum, Exparel versus plain bupivacaine (log-scale). Bold entries meet all four signal criteria. Peroneal nerve palsy (ROR 4.52) is the clinically most specific signal. Hypoaesthesia, paraesthesia, and hypoaesthesia oral are also four-criterion signals; the remaining PTs fall below signal threshold but are shown for completeness.*

### 

3.3 Composite-endpoint analyses
-------------------------------

**Prolonged sensory/motor block.** The all-time ROR versus plain bupivacaine was 2.91 (95% CI 2.40--3.52; IC025 0.67), with 222/2,127 Exparel cases (10.4%) versus 235/6,096 bupivacaine cases (3.9%). Era-matched stratification showed the signal intensifying after the 2018 interscalene indication: pre-2018 ROR 1.76 (95% CI 0.98--3.17; marginal, lower CI crossing 1), 2018--2020 ROR 4.52 (95% CI 3.24--6.30), and 2021+ ROR 2.26 (95% CI 1.69--3.03). No signal was observed versus ropivacaine (all-time ROR 0.80, 95% CI 0.64--1.00).

**LAST spectrum.** The all-time ROR was 1.00 (95% CI 0.85--1.17). Era-matched stratification revealed a pre-2018 signal of ROR 3.69 (95% CI 2.66--5.13) that inverted in subsequent eras: 2018--2020 ROR 0.67 (95% CI 0.51--0.89, statistically below the comparator rate) and 2021+ ROR 0.88 (95% CI 0.68--1.14, null). The pre-2018 elevation is consistent with early post-launch reporting inflation; the modern-era inversion argues against a persistent cardiotoxicity reporting differential under the chosen comparator design. Chen et al. \[24\] reported cardiogenic shock as a positive off-label IME-list signal in their pooled 2004--2024 whole-database analysis; that finding is methodologically consistent with the present era-stratified result, in that an aggregate signal is plausibly carried by the early-era reports and would not necessarily survive era-matched comparator restriction.

### Table 3. Composite-endpoint disproportionality by reporting era. a = Exparel cases reporting the composite endpoint; b = comparator cases reporting the composite endpoint.

  **Composite**                       **Era**      **Exparel cases**   **a**   **b**   **ROR vs bupi (95% CI)**   **Signal**   **ROR vs ropi (95% CI)**
  ----------------------------------- ------------ ------------------- ------- ------- -------------------------- ------------ --------------------------
  **Prolonged sensory/motor block**   All time     2,127               222     235     **2.91 (2.40--3.52)**      **✓**        0.80 (0.64--1.00)
  Prolonged sensory/motor block       Pre-2018     205                 14      85      1.76 (0.98--3.17)                       0.51 (0.27--0.94)
  Prolonged sensory/motor block       2018--2020   958                 111     55      **4.52 (3.24--6.30)**      **✓**        0.95 (0.63--1.45)
  Prolonged sensory/motor block       2021+        964                 97      95      **2.26 (1.69--3.03)**      **✓**        0.74 (0.53--1.04)
  **LAST spectrum**                   All time     2,127               223     641     1.00 (0.85--1.17)                       0.30 (0.25--0.37)
  LAST spectrum                       Pre-2018     205                 62      224     **3.69 (2.66--5.13)**      **✓**        0.95 (0.65--1.36)
  LAST spectrum                       2018--2020   958                 69      202     0.67 (0.51--0.89)                       0.27 (0.18--0.39)
  LAST spectrum                       2021+        964                 92      215     0.88 (0.68--1.14)                       0.28 (0.21--0.37)

*All era comparisons use era-matched comparator pools (both Exparel and plain-bupivacaine case sets restricted to the same epoch). Signal column denotes all four criteria met versus plain bupivacaine in that era. Era-specific denominators (from last_composite.csv): Exparel pre-2018 n=205, 2018--2020 n=958, 2021+ n=964 (Table 1; sum=2,127); era-matched plain-bupivacaine pools pre-2018 n=2,131, 2018--2020 n=1,950, 2021+ n=2,015 (sum=6,096); era-matched ropivacaine pools pre-2018 n=388, 2018--2020 n=265, 2021+ n=474 (sum=1,127). Era assignment uses event_dt where available; init_fda_dt is used as proxy for the \~71% of cases without event_dt.*

### 

### Figure 3. Prolonged-block composite disproportionality --- era-stratified, comparator-stratified, and sensitivity-restricted (unified forest plot)

![](media/image3.png){width="6.665972222222222in" height="7.070833333333334in"}

*Panel A shows Exparel vs plain bupivacaine across all-time, three eras, and the pathology-implying and expected-pharmacology strata. Panel B shows the secondary ropivacaine comparator across the same rows. Panel C shows pre-specified sensitivity restrictions versus plain bupivacaine (all-time). Red = signal meets all four criteria; blue = inverse signal (UCI < 1); open squares = non-signal. The prolonged-block signal (A) persists across 2018--2020 and 2021+ eras under era-matched comparator pools and under serious-outcome restriction within the pathology-implying stratum.*

3.4 Temporal evolution
----------------------

Figure 4 shows annual Exparel reporting volume and the time course of both composite endpoints. Reporting volume rose from single-digit counts immediately post-launch to a peak of \~450 cases in 2020, reflecting both growing commercial adoption and label expansion; reporting volume has since plateaued around 200 cases per year. The LAST composite rate peaked in the 2013--2014 reporting window and collapsed after 2017, consistent with a Weber-effect reporting artifact rather than a durable cardiotoxicity differential. The prolonged-block composite rate rose progressively and stabilized at approximately 10--12% of Exparel cases from 2018 onward.

### Figure 4. Temporal evolution of Exparel FAERS reporting, 2012--2025

![](media/image4.png){width="6.666666666666667in" height="5.434027777777778in"}

*Panel A: annual Exparel case counts. Panel B: annual rate of each composite endpoint among Exparel cases, with exact per-year values reported in Supplementary Table S6. Era shading denotes the three pre-specified epochs. Dashed vertical lines mark key events: ISB label (April 2018), Anesthesiology meta-analysis (Feb 2021), Pacira v. ASA (April 2021), and NOPAIN Act in effect (January 2025).*

3.5 Sensitivity analyses
------------------------

All four pre-specified sensitivity restrictions were executed (full outputs in Supplementary Table S7). Pregnancy/neonatal exclusion had minimal effect (composite ROR 3.04; pathology stratum ROR 2.41). Serious-outcome restriction attenuated the full composite below threshold (ROR 1.53, PRR 1.50) but the pathology-implying sub-stratum continued to signal (ROR 2.20, 95% CI 1.45--3.34). Physician-reporter restriction attenuated both the full composite (ROR 1.87) and the pathology stratum (ROR 1.47). Excluding manufacturer-submitted reports left only 29 of 2,127 Exparel cases (1.4%) --- underpowered for signal detection and a descriptively important finding: 97% of Exparel reports originate from the manufacturer channel.

4. Discussion
=============

4.1 Principal findings
----------------------

Across 13 years of FAERS data, liposomal bupivacaine was disproportionately associated with a pre-specified prolonged sensory/motor block composite compared with its parent compound plain bupivacaine (all-time ROR 2.91, 95% CI 2.40--3.52). The signal is driven by hypoaesthesia, peroneal nerve palsy, paraesthesia, and hypoaesthesia oral, and it persists across two independent modern eras under era-matched comparator pools.

Restricting the prolonged-block composite to pathology-implying PTs --- those whose MedDRA definition implies abnormal duration, anatomic injury, or persistent deficit --- preserved the signal versus plain bupivacaine (ROR 2.16, 95% CI 1.58--2.95); 8.5% of these cases also carried a Disability outcome code on the same report. The same stratum was inverse versus ropivacaine (ROR 0.48, 95% CI 0.35--0.68), arguing against a generic long-acting amide class effect and supporting a formulation-specific mechanism. Hypoaesthesia and peroneal palsy were independently top-ranked signals in Chen et al. \[24\] --- convergent evidence across two methodologically distinct designs.

By contrast, the frequently-cited LAST-spectrum signal was confined to the pre-2018 period (ROR 3.69) and inverted in the modern era (ROR 0.67 in 2018--2020 and 0.88 in 2021+ vs era-matched plain bupivacaine). This pattern is consistent with a Weber-effect reporting artifact rather than a durable cardiotoxicity differential.

4.2 Mechanistic interpretation
------------------------------

The prolonged-block signal is plausible given the liposomal formulation\'s intended pharmacokinetics: multi-day drug release beyond the expected analgesic window may produce neuropraxia through prolonged perineural exposure. The common peroneal nerve\'s anatomically superficial course around the fibular head makes it particularly susceptible. The most likely exposure pathway, given how widely it is practiced, is off-label periarticular infiltration in total knee arthroplasty --- a surgeon-administered technique in which posterior capsule injectate can reach the peroneal nerve. Popliteal sciatic block, in which the local anesthetic directly bathes the sciatic nerve at its peroneal-tibial bifurcation, is anatomically the most direct route but was off-label for most of the study period (approved on-label only in November 2023).

Three lines of evidence argue against a pure expected-pharmacology interpretation: (i) pathology-implying PTs signal independently of expected-pharmacology PTs versus plain bupivacaine; (ii) the median event-to-report lag of 206 days is more compatible with persistent than with transient \< 72 h reactions, although lag may also reflect litigation or claims timelines; and (iii) approximately 8.5% of pathology-flagged cases carry a Disability outcome code on the same report (event-level attribution caveat in §4.6). The inverse-direction finding versus ropivacaine is also informative --- ropivacaine differs from bupivacaine in producing less motor block and having a shorter elimination half-life, but it shares the long-acting amide class. The absence of a comparable pathology signal argues against a generic long-acting-amide class effect.

4.3 Comparison with prior literature
------------------------------------

Published evidence on liposomal bupivacaine remains mixed. Ilfeld et al. \[3\] concluded that liposomal bupivacaine had not shown consistent analgesic superiority over conventional local-anesthetic techniques across clinical settings, and Ilfeld and Sessler \[6\] emphasized uncertainty about whether observed duration differences translate into meaningful clinical benefit. At the same time, the product label already acknowledges temporary sensory or motor loss lasting up to 5 days, making a neurologic-symptom reporting signal biologically plausible and label-concordant rather than wholly unexpected.

Prior pharmacovigilance work reported a signal linking Exparel with LAST in FAERS \[7\]. The present study extends that literature by covering the full post-launch period, addressing formulation misclassification via prod_ai, and stratifying by reporting era with era-matched comparator pools. The modern-era attenuation and inversion of the LAST signal suggests that the earlier cardiotoxicity signal may have been influenced by early-post-launch reporting dynamics rather than reflecting a durable formulation-specific excess. A 2022 letter by Finkel et al. in the British Journal of Anaesthesia \[8\] reported an association between manufacturer financial conflicts of interest and more favorable outcomes in randomised controlled trials of liposomal bupivacaine --- supporting cautious interpretation of both efficacy and safety narratives, and reinforcing the methodological relevance of an analysis aimed at comparative reporting rather than absolute incidence.

Three design features distinguish the present analysis from typical single-drug FAERS scans \[15, 16\]. First, plain bupivacaine serves as the active comparator, separating formulation-specific effects from molecule-class effects \[14\]. Second, era-matched stratification reveals the LAST inversion, which is hidden in all-time pooled analyses. Third, a pre-specified pathology-vs-pharmacology PT partition addresses the objection --- unique to extended-release formulations --- that any apparent signal merely reflects expected drug action. Chen et al. \[24\], a 2025 whole-database FAERS scan of liposomal bupivacaine reporting 58 positive signals over 2004 Q1--2024 Q2, applied none of these features and included seven years of data preceding the drug\'s commercial availability.

4.4 Reporting biases
--------------------

Spontaneous reporting is subject to notoriety, channeling, and Weber effects; the era-matched stratified design was chosen to expose rather than obscure these biases. Two sensitivity restrictions tested signal robustness. The pathology-implying sub-stratum continued to signal under serious-outcome restriction (ROR 2.20, 95% CI 1.45--3.34) but fell below threshold under physician-reporter-only restriction (ROR 1.47). The full prolonged-block composite fell below the four-criterion threshold under both restrictions (serious-outcome ROR 1.53; physician-reporter ROR 1.87). Ninety-seven percent of Exparel reports originate from the manufacturer (Pacira) channel, so absolute reporting rates should be interpreted with this in mind \[21\]. The ropivacaine comparator, which is unaffected by Pacira-related reporting dynamics, did not show concurrent 2018--2020 intensification --- partial reassurance against pure stimulated reporting, but not full exclusion. Comparator choice also affects ROR magnitude: peroneal nerve palsy yielded ROR 4.52 here versus ROR 63.98 in Chen et al. \[24\] under a whole-database reference. This \~14-fold difference is not a substantive contradiction but a direct consequence of comparator-pool sensitivity: whole-database designs systematically inflate magnitudes for procedure-specific drugs whose reference pool contains millions of unrelated indications \[22\].

4.5 Clinical implications
-------------------------

The clinical implications are modest. These data support clinician awareness that prolonged sensory or motor symptoms are biologically plausible and label-consistent with liposomal bupivacaine use, particularly in peripheral nerve block settings where the formulation is used near anatomically vulnerable nerves (peroneal, brachial plexus). Whether structured neurological reassessment at 48--72 hours and again at one week would meaningfully improve detection of persistent symptoms is a question for prospective study; current data do not justify a protocol-level recommendation. Counseling patients before perineural use that sensory and motor effects may persist beyond the expected analgesic window, and that a small subset of cases may not resolve fully within the expected timeframe, is defensible on current evidence.

The LAST findings support continued use with standard local-anesthetic safety monitoring; no additional cardiac precautions versus plain bupivacaine are indicated by current-era data. However, these results do not justify strong causal claims about permanent nerve injury from FAERS alone, nor do they justify protocol-level changes on that basis. These findings call for prospective denominator-based safety studies, not immediate practice change.

4.6 Limitations
---------------

This study has the usual limitations of spontaneous-report pharmacovigilance \[18\]. FAERS lacks denominator data, so disproportionality cannot estimate incidence, prevalence, relative risk, or causality. The dabigatran experience \[25\] illustrates this concern: a large early-post-launch FAERS bleeding signal versus warfarin did not survive subsequent claims-based cohort analysis. Whether the present prolonged-block signal would survive analogous assessment is an open question.

Several specific issues also apply. Under the ICH E2B(R2) reporting standard in force during the study period, seriousness outcomes are coded at the report level rather than the event level; the transition to E2B(R3) by April 2026 should resolve this. Caseid-level deduplication does not address unlinked duplicates from multi-sponsor reporting \[19\]. Exparel, plain bupivacaine, and ropivacaine are not perfectly exchangeable comparators. Stimulated reporting remains possible: the Munoz and Dal Pan \[21\] framework on litigation-associated reporting bias is directly relevant to the 2018--2020 period. Finally, the plain-bupivacaine cohort (n=6,096) is \~3.8-fold larger per unit time than that in Chen et al. \[24\] (n=1,589), reflecting differences in drug-name capture.

4.7 Future directions
---------------------

-   Denominator-based claims or registry cohort study of post-Exparel peripheral neuropraxia incidence, using the National Surgical Quality Improvement Program (NSQIP), the Medicare 5% sample, or the Premier Healthcare Database as the primary exposure-ascertainment source.

-   Prospective neurological-examination surveillance study at high-volume centers, with structured 48--72 h and 1-week assessments.

-   Dose-response analysis using FAERS dose_vbm once dose-field harmonization becomes sufficiently reliable.

-   Active-comparator analysis against ropivacaine-catheter regimens in orthopedic populations, to separate formulation-specific from technique-specific effects.

-   Integration with Centers for Medicare & Medicaid Services (CMS) Open Payments data and J0666 utilization once the 2025 Medicare Part B data release becomes available (mid-2026), to examine whether adoption geography predicts reporting geography.

5. Conclusions
==============

Exparel is disproportionately associated with prolonged sensory/motor block in FAERS. The signal persists versus plain bupivacaine across both modern reporting eras (2018--2020 and 2021+) and survives restriction to pathology-implying Preferred Terms. Versus ropivacaine the result is null to inverse, which argues against a generic long-acting-amide class effect and is more consistent with a formulation-specific mechanism. The signal is therefore unlikely to reduce to expected prolonged pharmacology alone. The cardiac toxicity signal attributed to this drug in the early post-launch period does not persist in modern-era data and appears to be a reporting artifact. These findings are hypothesis-generating; they do not establish incidence or causation but they do support further denominator-based safety studies, prospective neurological follow-up research for perineural Exparel use, and careful patient counseling about the possibility of symptoms persisting beyond the expected analgesic window.

Statements and Declarations
===========================

Funding
-------

The author received no external funding for this work. No financial support of any kind was received from any organization for the submitted work, and no compensation was received by the author for conducting this study.

Competing interests
-------------------

The author declares no commercial or financial relationship with Pacira BioSciences Inc. (manufacturer of liposomal bupivacaine / Exparel), with any manufacturer of plain bupivacaine or ropivacaine, or with any party with a financial interest in the outcome of this analysis. The author has no patents, consulting arrangements, speaker fees, stock holdings, or other financial interests related to liposomal bupivacaine or any comparator agent within the past three years.

Ethics approval
---------------

This study used the U.S. FDA Adverse Event Reporting System (FAERS) Quarterly Data Extract Files, which are publicly available, de-identified post-marketing surveillance data. Under 45 CFR 46.104(d)(4)(ii), research involving the use of publicly available, de-identified data does not constitute human subjects research and is exempt from Institutional Review Board oversight. The UC Davis IRB has confirmed that analyses of FAERS data of this type do not require IRB review or approval.

Consent to participate
----------------------

Not applicable. This study used publicly available, de-identified post-marketing surveillance data; no individual participants were enrolled, contacted, or directly studied.

Consent to publish
------------------

Not applicable. No identifiable individual data are reported.

Data availability
-----------------

FAERS raw data are publicly available from the U.S. FDA Freedom of Information Office (https://fis.fda.gov/extensions/FPD-QDE-FAERS/). The analytic dataset used in this study, including all intermediate derived datasets, is archived at Zenodo under concept DOI 10.5281/zenodo.19656698 (always-latest); the v7.10 snapshot corresponding to this manuscript is permanently archived at DOI 10.5281/zenodo.20222469.

Code availability
-----------------

All analysis code, the FAERS-ingestion pipeline, composite-endpoint construction scripts, figure-generation scripts, and the manuscript source are version-controlled and publicly available at https://github.com/bobsterbeat/faers-exparel-study. The code was run under Python 3.14 with pandas 2.3; full dependency manifest is in the repository.

Author contributions
--------------------

R. Aldwinckle is the sole author and contributed to all aspects of this work: conceptualization, methodology, data curation, formal analysis, investigation, visualization, manuscript drafting, and manuscript revision.

AI assistance
-------------

Generative AI tools (Anthropic Claude) were used to assist with figure rendering (matplotlib code generation), manuscript copy-editing for grammar and style, and structured discussion of methodological design choices. No autonomous content generation was performed; the author conceived the study design, validated all analytic outputs against source data, and accepts full responsibility for the manuscript content, accuracy, and conclusions. AI tools did not contribute to authorship and do not meet ICMJE authorship criteria.

References
==========

> 1\. Fusaroli M, Salvo F, Begaud B, et al. The REporting of A Disproportionality Analysis for DrUg Safety Signal Detection Using Individual Case Safety Reports in PharmacoVigilance (READUS-PV): Development and Statement. Drug Saf. 2024;47(6):575--584. https://doi.org/10.1007/s40264-024-01422-8.
>
> 2\. Fusaroli M, Salvo F, Begaud B, et al. READUS-PV: Explanation and Elaboration. Drug Saf. 2024;47(6):585--599. https://doi.org/10.1007/s40264-024-01423-7.
>
> 3\. Ilfeld BM, Eisenach JC, Gabriel RA. Clinical Effectiveness of Liposomal Bupivacaine Administered by Infiltration or Peripheral Nerve Block to Treat Postoperative Pain. Anesthesiology. 2021;134(2):283--344. https://doi.org/10.1097/ALN.0000000000003630.
>
> 4\. Hussain N, Brull R, Sheehy B, et al. Perineural Liposomal Bupivacaine Is Not Superior to Nonliposomal Bupivacaine for Peripheral Nerve Block Analgesia. Anesthesiology. 2021;134(2):147--164. https://doi.org/10.1097/ALN.0000000000003651.
>
> 5\. McCann ME. Liposomal Bupivacaine: Effective, Cost-effective, or (Just) Costly? Anesthesiology. 2021;134(2):139--142. https://doi.org/10.1097/ALN.0000000000003658.
>
> 6\. Ilfeld BM, Sessler DI. Liposomal Bupivacaine in Peripheral Nerve Blocks: Duration and Meaningful Differences. Anesthesiology. 2024;141(4):638--642. https://doi.org/10.1097/ALN.0000000000005133.
>
> 7\. Aggarwal N. Local anesthetics systemic toxicity association with Exparel (bupivacaine liposome) --- a pharmacovigilance evaluation. Expert Opin Drug Saf. 2018;17(6):581--587. https://doi.org/10.1080/14740338.2018.1453496.
>
> 8\. Finkel KJ, Takata ET, Maffeo-Mitchell CL, et al. Manufacturer financial conflicts of interest are associated with favourable outcomes in randomised controlled trials of liposomal bupivacaine. Br J Anaesth. 2022;129(4):e90--e93. https://doi.org/10.1016/j.bja.2022.06.032.
>
> 9\. Pacira BioSciences Inc. v. American Society of Anesthesiologists, No. 22-1411 (3d Cir. Mar. 24, 2023).
>
> 10\. DailyMed. EXPAREL (bupivacaine liposome injectable suspension) prescribing information. U.S. National Library of Medicine. Accessed April 19, 2026.
>
> 11\. Evans SJ, Waller PC, Davis S. Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports. Pharmacoepidemiol Drug Saf. 2001;10(6):483--486.
>
> 12\. Bate A, Evans SJ. Quantitative signal detection using spontaneous ADR reporting. Pharmacoepidemiol Drug Saf. 2009;18(6):427--436.
>
> 13\. U.S. Food and Drug Administration. FDA Adverse Event Reporting System (FAERS) Quarterly Data Extract Files. https://fis.fda.gov/extensions/FPD-QDE-FAERS/. Accessed April 2026.
>
> 14\. Alkabbani W, Gamble J-M. Active-comparator restricted disproportionality analysis for pharmacovigilance signal detection studies of chronic disease medications. Br J Clin Pharmacol. 2023;89(5):1464--1474.
>
> 15\. Zou J, Zhan C, Feng Y, Liu Q, Peng G, Cai Z, Luo L, Zou J. Drug-induced movement disorder: A disproportionality analysis using the FDA adverse event reporting system (FAERS) from 2004 to 2024. PLoS One. 2025;20(10):e0335449. https://doi.org/10.1371/journal.pone.0335449.
>
> 16\. Yang C, Zhao W, Chen H, Yao Y, Zhang J. Cardiac adverse events associated with lacosamide: a disproportionality analysis of the FAERS database. Sci Rep. 2024;14(1):16202. https://doi.org/10.1038/s41598-024-67209-0.
>
> 17\. Hauben M, Bate A. Decision support methods for the detection of adverse events in post-marketing data. Drug Discov Today. 2009;14(7-8):343--357.
>
> 18\. Potter E, Reyes M, Naples J, Dal Pan G. FDA Adverse Event Reporting System (FAERS) Essentials: A Guide to Understanding, Applying, and Interpreting Adverse Event Data Reported to FAERS. Clin Pharmacol Ther. 2025;118(3):567--582. https://doi.org/10.1002/cpt.3701.
>
> 19\. Kreimeyer K, Dang O, Spiker J, et al. Increased confidence in deduplication of drug safety reports with natural language processing of narratives at the US Food and Drug Administration. Front Drug Saf Regul. 2022;2:918897. https://doi.org/10.3389/fdsfr.2022.918897.
>
> 20\. Hauben M, Zou C, Bright S, Hung E. More extreme duplication in the U.S. FDA FAERS database and a suggested check point for disproportionality analysis. Pharmacoepidemiol Drug Saf. 2023;32(4):387--391.
>
> 21\. Muñoz MA, Dal Pan GJ. The impact of litigation-associated reports on signal identification in the US FDA\'s adverse event reporting system. Drug Saf. 2019;42(10):1199--1201. https://doi.org/10.1007/s40264-019-00834-1.
>
> 22\. Mouffak A, Lepelley M, Revol B, Bernardeau C, Salvo F, Pariente A, Roustit M, Cracowski J-L, Khouri C. High prevalence of spin was found in pharmacovigilance studies using disproportionality analyses to detect safety signals: a meta-epidemiological study. J Clin Epidemiol. 2021;138:73--79.
>
> 23\. Pham P, Cheng C, Wu E, Kim I, Zhao Y, Dal Pan G. Leveraging case narratives to enhance patient age ascertainment from adverse event reports. Pharmaceut Med. 2021;35(5):307--316.
>
> 24\. Chen DX, Chen XM, Chen SM, Wang YD. Real-World Pharmacovigilance Analysis of Adverse Events Associated with Liposomal Bupivacaine and Bupivacaine. J Pain Res. 2025;18:1805--1816. https://doi.org/10.2147/JPR.S519504.
>
> 25\. Southworth MR, Reichman ME, Unger EF. Dabigatran and post-marketing reports of bleeding. N Engl J Med. 2013;368(14):1272--1274.

Supplementary material
======================

-   S1. FAERS file manifest and checksums (53 quarters, 2012 Q4 -- 2025 Q4).

-   S2. Full pre-specified PT lists for each composite endpoint.

-   S3. Full PT-level signal table with IC025 and PRR χ² for every PT with a ≥ 3.

-   S4. 20 representative case narratives spanning era and sub-phenotype.

-   S5. Drug-name regex audit: random-sample verification and prod_ai misclassification rate.

-   S6. Year-resolved composite-endpoint rates (supporting Figure 4).

-   S7. Full sensitivity-analysis output tables.

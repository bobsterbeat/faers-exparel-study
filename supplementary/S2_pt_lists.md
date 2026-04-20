# Supplementary Appendix S2 — Pre-specified MedDRA Preferred Term lists

All PT lists below were locked in the analysis code prior to outcome estimation. See `duration_analysis.py` and `last_composite.py` in the code repository for the authoritative definitions.

## Prolonged sensory/motor block composite (20 PTs, case-level)

Union of the three mutually-exclusive strata plus the ambiguous term `anaesthesia`.

### Pathology-implying stratum (15 PTs)

PTs whose MedDRA definition requires abnormal duration or permanent deficit and therefore cannot represent expected <72 h pharmacology.

- peroneal nerve palsy
- peroneal nerve injury
- neuromuscular block prolonged
- nerve injury
- neuropathy peripheral
- peripheral nerve injury
- peripheral motor neuropathy
- peripheral sensory neuropathy
- brachial plexopathy
- paresis
- hemiparesis
- monoparesis
- paralysis
- neuralgia
- motor dysfunction

### Expected-pharmacology stratum (4 PTs)

PTs that could reflect Exparel's intended effect at <72 h exposure. Disproportional reporting remains interpretable but cannot in itself distinguish normal from prolonged block.

- hypoaesthesia
- paraesthesia
- muscular weakness
- sensory loss

### Systemic-absorption marker (1 PT)

- hypoaesthesia oral

### Excluded as terminologically ambiguous (1 PT)

- anaesthesia — excluded from all three strata. Retained in the overall prolonged-block composite for completeness but not used for duration-specific claims.

## LAST-spectrum composite (26 PTs, case-level)

### Cardiac

- cardiac arrest
- cardio-respiratory arrest
- cardiac disorder
- ventricular tachycardia
- ventricular fibrillation
- ventricular arrhythmia
- arrhythmia
- bradycardia
- hypotension
- blood pressure decreased
- cardiovascular insufficiency

### CNS

- seizure
- convulsion
- generalised tonic-clonic seizure
- loss of consciousness
- depressed level of consciousness
- unresponsive to stimuli

### Respiratory (often CNS-depression-induced hypoventilation)

- hypoxia
- dyspnoea
- hypercapnia
- respiratory depression
- respiratory failure
- pulmonary oedema
- acute pulmonary oedema

### Metabolic / nonspecific toxicity marker

- acidosis
- toxicity to various agents

## Version pinning

MedDRA version at analysis: whichever version was current in each FAERS quarterly release (FDA does not normalize across releases). Reviewers wishing to reproduce should note that PT harmonization across MedDRA versions is a known limitation of multi-year FAERS studies and is discussed in §4.6 of the manuscript.

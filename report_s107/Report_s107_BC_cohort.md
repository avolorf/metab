# Bladder-cancer cohort of 107 urine samples: reprocessing with the s266 workflow, and re-evaluation of previously reported markers

*Supplement to the s266 analysis (BC, PC and controls). All numbers below were produced by the same code (`METAB_DATASET=s107 python run_all.py`, `python report_s107.py`); tables are in `report_s107/Tables_s107.xlsx` (version v2, also saved as `Tables_s107_v2.xlsx`), figures in `report_s107/`.*

*Version v2 (12 Nov 2026). Changes from v1: (1) sample No. 32 (bladder leukoplakia) excluded, as in the ABC article; (2) sex of samples Nos. 35 and 86 corrected to female from the primary records; (3) revised finding rule — a compound not detected in any control passes the high- vs low-index control criterion, which is undefined for it; (4) categories transferred from the annotations of the previous report and the ABC article are kept only when confirmed by the data. A v1 → v2 comparison is in `report_s107/v2_changes/`.*

## 1. Methods (differences from the s266 workflow)

**Samples.** The MS-DIAL export contained 107 urine samples (67 labelled BC, 40 controls) and 21 pooled-QC injections; no solvent blanks were exported. Three BC-labelled samples were excluded because their histology was outside the study design: a ureteral tumour (No. 54), bladder leukoplakia (No. 32) and Brunn's nest hyperplasia without carcinoma (No. 105). Urine was collected before any diagnostic or surgical manipulation; intact urine was aliquoted into 4-mL tubes and immediately stored at −20 °C, and each sample was mixed by pulse-vortexing after thawing. The analysed cohort thus comprised 64 BC patients (54 men, 10 women; age 65.8 ± 13.4 years, range 26–87) and 40 controls (16 men, 24 women; 57.3 ± 14.8 years, range 18–79). Age (p = 0.002, Mann–Whitney U test) and sex (p = 4.2 × 10⁻⁶, Fisher's exact test) differed between the groups. The sex distribution agrees with the primary records; Table 2 of the ABC article gives 12 men / 28 women for the controls, which is a typographical error (16 / 24). Eighty-two of the 107 samples are part of the 84-sample set analysed in the previous report and in the ABC article.

**Processing.** Feature filtering, PQN normalization (median of the pooled QC as reference), imputation, log₂ transformation, grouping of ions into compounds, the adjusted linear model (group + sex + age + exogenous-signal index), the finding rule, cross-validated classification (20 × 5-fold CV, top-5 panel, PLS-DA, 1,000 permutations) and all sensitivity analyses were identical to the s266 analysis. The following settings were adapted to this dataset:

- **Batches.** Samples were analysed in three batches (Nos. 1–36, 37–72 and 73–107). In each batch the pooled QC was injected at the start, after every six samples and at the end (seven QC injections per batch). Solvent blanks were injected but not exported; blank filtering was done in MS-DIAL. The groups were spread similarly over the batches (BC 23/20/21, controls 12/15/13). Before correction, 3,497 of 5,887 features differed between the QC injections of different batches (Kruskal–Wallis test, q < 0.05), and the median QC RSD was 15.7%. Each batch was therefore corrected feature by feature, by multiplying its values by the ratio of the median of all QC injections to the median of the QC injections of that batch. The correction was applied before the RSD filter. After correction, no feature differed between batches, and the median QC RSD fell to 10.4% (Figure S107-0).
- **Exogenous-signal index.** The index was built from the same homologous series. Five members were detected in this dataset, at retention times about 0.42 min shorter than in s266.
- **Splitting the index.** Unlike in s266, the index was higher in controls than in BC (mean z-score +0.69 vs −0.43). Splitting at the median of all samples would therefore have left only 8 low-index controls. Controls were instead split into low- and high-index halves (20/20) at the median of the controls.

**Annotation.** The 30 markers of the previous report (report Table 20, ABC) and the alternative annotations in `res_annotated_2.csv` were linked to features (±5 ppm, ±0.3 min). Each was then checked for agreement between formula, ion and m/z. The category (endogenous, food, drug or tobacco) proposed in the report was kept only when the data supported the identity, i.e. when the MS-DIAL annotation of the same feature named the same compound; this was the case for acetylcysteine, 4-hydroxy-3-(2-hydroxyethyl)acetophenone, isoferulic acid, venlafaxine, Thr-Phe and Asn-Phe. All other report categories were set to "ambiguous". For m/z 269.0880 (2.02 min) and m/z 362.1248 (8.79 min) the raw spectra were inspected (samples 4 and 56, QC17, blank 15): neither the hypoxanthine fragment m/z 137.0458 nor a consistent isotope pattern was found for the former, and no chlorine isotope or MS/MS spectrum for the latter, so both remain unidentified (MSI 4). As in s266, exogenous compounds could not qualify as findings.

**Finding rule (revision).** The criterion "no difference between controls with high and low exogenous-signal index" cannot be tested for a compound that is not detected in any control (all control values are the same imputed value). In v1 such compounds were rejected; in v2 the criterion is considered met for them, and they are flagged (`finding_via_ctrl_not_detected`). The same code was re-run for s266, where no compound is affected: all s266 results are unchanged.

## 2. Results

**Data quality.**
- **Features and compounds.** Of 12,930 aligned features, 5,887 passed the detection filter and, after the batch correction, 4,730 passed the QC RSD filter (median RSD in QC 8.1%). These were grouped into 3,603 compounds, 483 of them multi-ion.
- **Imputation.** 11.0% of the values were imputed.
- **Drift.** Within batches, the leave-one-QC-out test gave a median |Δlog₂| of 0.067 without correction and 0.065 with LOWESS, so no drift correction was applied.
- **Dilution.** Three samples had a PQN dilution factor < 0.2.

**Univariate analysis.**
- **Adjusted model.** 918 compounds had q < 0.05 (463 higher and 455 lower in BC). 582 of them met all finding criteria: 294 higher and 288 lower in BC. Twenty-five findings, all increased in BC, are not detected in any control and passed the index criterion under the revised rule.
- **Matched pairs.** The pair criterion used 23 pairs matched on age within sex.
- **Normalization.** With MSTUS instead of PQN, 503 of the 607 MSTUS findings were also PQN findings.
- **Most dilute samples.** Excluding the three most dilute samples (Nos. 24, 48, 106) kept 534 of the 582 findings.
- **Pathways.** Only leukotriene metabolism reached q < 0.05 (mummichog; q = 0.028). The mummichog result depends on random permutations (q = 0.046–0.068 in earlier runs), so it is not interpreted.

**Separation of the groups is almost complete.** The cross-validated top-5 panel gave AUC 0.993 ± 0.003 and PLS-DA 1.000 (both p = 0.001). The separation persisted in the subsets that control for the main confounders:

| Subset or model | Panel AUC |
|---|---|
| 23 age- and sex-matched pairs | 0.985 |
| 20 pairs matched on the exogenous-signal index | 1.000 |
| Compounds adjusted for the index within each training fold | 0.889 |
| Only the 697 candidate compounds detected in every sample | 0.981 (PLS-DA 0.991) |

The separation is therefore not produced by imputation or presence/absence patterns. Injection order did not explain it either: in the PCA of the 752 compounds detected in all samples, PC1 separated the groups (AUC 0.73) but did not follow injection order (ρ = −0.18, Figure S107-1).

The two baseline models were weaker but not uninformative:
- **Sex only:** AUC 0.68.
- **Exogenous-signal index only:** AUC 0.80. The index is higher in controls, the opposite of its behaviour in s266, where it was highest in PC. This again indicates that it reflects sample collection rather than disease.

**Composition of the BC signal.** Short peptides dominated the compounds increased in BC (Figure S107-4):
- **Enrichment.** Library matches to di- and tripeptides were enriched among compounds increased in BC (75 of 463 vs 202 of 3,140 other compounds; odds ratio 2.8, Fisher's exact test p = 3 × 10⁻¹¹).
- **Strongest signals.** Val-Leu (log₂FC 3.3, AUC 0.98) and m/z 244.1293 at 2.58 min (log₂FC 9.3, AUC 0.98, detected in 98% of BC vs 32% of control samples). The latter is annotated in report Table 20 as N-mevalonylcreatinine; its formula, C₁₀H₁₇N₃O₄, is also that of the dipeptide Gln-Pro, the alternative annotation given for the same feature in `res_annotated_2.csv`. Val-Leu formally fails the index criterion (controls with high and low index differ, p = 0.027) and is therefore not counted as a finding.
- **Other peptides.** Glu-Ile-Val at 7.95 min (report annotation; log₂FC 3.8, AUC 0.87), Thr-Phe and Asn-Phe. The Glu-Ile-Val of the ABC article (m/z 360.2132 at 2.57 min) is a separate peak: it is not detected in any control and is detected in 80% of BC samples (log₂FC 6.1, AUC 0.90); it is a finding under the revised rule.

Because urine was collected before any diagnostic or surgical manipulation, cystoscopy, catheterisation and perioperative treatment can be excluded as causes. An almost complete separation carried by peptides that are nearly absent from control urine is nevertheless more consistent with blood or protein in the urine of BC patients than with tumour metabolism: haematuria is the cardinal symptom of BC, and haemoglobin and plasma proteins are substrates for proteolysis. Storage at −20 °C rather than −80 °C does not exclude some proteolysis in frozen samples, and storage times of the two groups are not known. Haematuria and proteinuria were not recorded, so these explanations cannot be separated here. The AUC values of this cohort should therefore not be read as diagnostic accuracy.

**Replication of the s266 findings (Table 2 of the s266 analysis; T2 in Tables_s107.xlsx; Figure S107-3).** All three steroid glucuronides were increased in BC in the 107-sample cohort (q < 0.05 in the adjusted model). 11-Oxoetiocholanolone glucuronide and tetrahydrocortisone glucuronide isomer II met all finding criteria; isomer I missed the matched-pair criterion (p = 0.08 in 23 pairs) after the cohort revision:

| Compound | log₂FC | q | AUC | p (matched pairs) |
|---|---|---|---|---|
| 11-Oxoetiocholanolone glucuronide | 0.76 | 0.011 | 0.73 | 0.016 |
| Tetrahydrocortisone glucuronide, isomer I | 0.83 | 0.012 | 0.75 | 0.080 |
| Tetrahydrocortisone glucuronide, isomer II | 0.53 | 0.025 | 0.70 | 0.021 |

Each compound was found as the same group of ions (e.g., m/z 498.2697 with the aglycone ion 287.2007; m/z 541.2642 with its [M+NH₄]⁺ ion 558.2912), 0.30–0.33 min earlier than in s266. The other two s266 findings did not replicate:
- **Biliverdin** showed no difference (q = 0.45).
- **LysoPC(16:0)** did not pass the filters in this dataset.

The steroid glucuronides are thus the only annotated BC result observed independently in both cohorts (all three increased with q < 0.05; two of three meet all criteria in s107). Their effect size (AUC about 0.7) is much smaller than that of the peptide signal.

**Re-evaluation of the 30 markers of the previous report (T1 in Tables_s107.xlsx).** All 30 markers were found in the data, and the reported direction of change was confirmed for 28 of them in the adjusted model. However:

- **Exogenous compounds.** Only five of the 30 markers have an exogenous identity supported by the data (the MS-DIAL annotation of the same feature names the same compound): the drugs acetylcysteine and venlafaxine and the food-derived 4-hydroxy-3-(2-hydroxyethyl)acetophenone, isoferulic acid and fraxin (the last two are the same compound in this dataset). They reflect concomitant medication and diet, which differ systematically between patients and volunteers (samples were collected before any diagnostic or surgical manipulation, so perioperative treatment is not the cause), and cannot be regarded as BC markers. The other twelve markers that the report classed as drugs, food constituents or tobacco alkaloids (e.g. samuraciclib, the propofol–N-acetylcysteine adduct, 5-methylmyosmine, plant glycosides) are not supported by the data — MS-DIAL proposes a different compound or the mass does not fit — and are treated as unidentified ("ambiguous").
- **Findings.** Thirteen of the 30 markers met all finding criteria (T1): m/z 244.129 (Gln-Pro / N-mevalonylcreatinine), Thr-Phe, Asn-Phe, Glu-Ile-Val (7.95 min), m/z 269.088 (2.02 min; unidentified, see Annotation), the two oxylipin glucuronides at m/z 508.312 (one compound), the γ-CEHC glucuronide, and five compounds whose report identities are not supported (5-methylmyosmine, N-acetyl-S-(N-allylthiocarbamoyl)cysteine, vanillyl glycol glycoside, samuraciclib, the propofol–N-acetylcysteine adduct). The last group and m/z 269.088 are increased or decreased in BC, but their identities, and therefore their origin, remain open.
- **Identities that do not fit the m/z.**
  - m/z 445.243 is reported as grayanotoxin I, but grayanotoxin I does not fit this m/z as [M+H]⁺ (>77,000 ppm); in v2 it no longer meets the matched-pair criterion.
  - 5-methylmyosmine deviates by 63 ppm.
  - The reported formula of the peptide EFTPAVAH (C₄₀H₅₈N₁₀O₁₂) gives [M+2H]²⁺ m/z 436.218, not 436.254.
  - m/z 244.129 is more likely Gln-Pro than N-mevalonylcreatinine.
- **Not confirmed.** Three markers were not confirmed: 3-hydroxyoctanoylcarnitine (q = 0.37), N,N-dimethylsphingosine (opposite direction in the adjusted model), and the antipyrine glucuronide (q = 0.56, opposite direction).

**The 14 compounds of Table 5 of the ABC article** were checked separately (`report_s107/abc_table5/`). All 14 were detected and increased in BC. Twelve met all finding criteria. Val-Leu (4.14 min) failed the index criterion (controls with high and low index differ, p = 0.027). The N-acetyltryptophan feature (m/z 245.0923 at 6.34 min; 249.0923 in Table 5 is a typographical error) was removed by the QC RSD filter (31%) and does not fit the formula (MSI 4). Leu-Val (3.01 min) and m/z 360.2129 at 7.95 min are separate peaks, not retention-time shifts of Val-Leu and Glu-Ile-Val.

## 3. Prostate-cancer markers of the previous processing of s266 (report Tables 22, 28 and 29)

The 25 PC marker signals of report Table 28 were matched to the s266 features (±10 ppm, ±0.2 min). For each, we report its PC-vs-control statistics in our analysis (T3 in Tables_s107.xlsx; `results_s107_extra/PC_markers_old_report_vs_s266.csv`). All 25 signals were found.

- **Exogenous series.** Fourteen of the 25 signals (Nos. 11–24) are members of the homologous exogenous series used for the exogenous-signal index: alkyl diethanolamines and their oxygenated and ethoxylated forms, and the co-varying C₂₄H₃₀O₆ ions. In report Table 22, several of these signals are annotated as phytosphingosine, 2-amino-1,3,4,5-eicosanetetrol, 2-amino-12,16-dimethyleicosane-pentol, 15(S)-HETE, 2(R)-hydroxyeicosanoic acid, armillarin, nirantin and eplerenone.
  - Phytosphingosine (C₁₈H₃₉NO₃) is isobaric with C₁₄-DEA+O and elutes at its retention time within the series.
  - 15(S)-HETE (C₂₀H₃₂O₃, [M+H]⁺ 321.242) and 2(R)-hydroxyeicosanoic acid (C₂₀H₄₀O₃, [M+H]⁺ 329.305) do not match the observed m/z 302.306 and 346.332, which are C₁₄-DEA and C₁₆-DEA+O.
  - These signals correlate strongly with the index (ρ = 0.83–0.98). They are higher in PC because the index itself is higher in PC samples, and none of them differs after adjustment for the index (q ≥ 0.26).
- **Signal No. 25** (m/z 531.275, 18.1 min) also correlates with the index (ρ = 0.84) and is not significant after adjustment (q = 0.49).
- **Signals 1–10** reached nominal p < 0.05 in the unadjusted Mann–Whitney test, but none passed FDR control (q = 0.06–0.26). None was significant in the adjusted model (q ≥ 0.13) or in the comparison with male controls (q ≥ 0.61).
  - For kynurenic acid, Pro-Asp, Pro-Tyr and m/z 345.155, report Table 22 states an increase in PC, whereas the data show a decrease (log₂FC −1.0 to −3.8). The decrease agrees with the negative logFC values in report Table 29.
  - Several annotations in report Table 22 do not match the m/z. For Pro-Asp (C₉H₁₄N₂O₅, [M+H]⁺ 231.098 vs 211.072), Pro-Tyr (C₁₄H₁₈N₂O₄, [M+H]⁺ 279.134 vs 261.124; at most a water-loss ion) and adenosine ([M+H]⁺ 268.104; [M+Na]⁺ 290.086, 35 ppm from 290.096), the stated formula and ion do not give the observed mass. In several rows the formula column belongs to another compound.
- **Shared BC/PC signals.** The five signals of report Table 29 are Nos. 2, 3, 4, 8 and 25 of Table 28 and are covered above.

**Conclusion on the PC markers.** None of the 25 previously reported PC signals withstands adjustment for sex, age and the exogenous-signal index together with FDR control. More than half of them are the exogenous series itself. This agrees with our analysis, in which no compound met the finding criteria for PC vs controls. The apparent PC markers of the previous processing are explained by the higher exogenous-signal index of the PC samples (a pre-analytical, container- or collection-related signal), by annotation errors, and by the absence of FDR control.

## 4. Summary for integration with the s266 results

1. **Steroid glucuronides.** The increase of three urinary steroid glucuronides (11-oxoetiocholanolone glucuronide and two tetrahydrocortisone glucuronide isomers) in BC was found independently in two cohorts (s266: 29 BC / 45 controls; s107: 64 BC / 40 controls), with similar effect sizes (AUC 0.69–0.81). All three are increased with q < 0.05 in both cohorts; in s107, isomer I narrowly misses the matched-pair criterion. This is the most robust annotated BC result of the project.
2. **Peptide signal in s107.** The peptide-dominated signal of the 107-sample cohort gives almost perfect separation, but this is very likely influenced by haematuria/proteinuria or by pre-analytical differences between patients and volunteers. It should be reported with this caveat and not as diagnostic performance.
3. **Earlier marker lists.** Of the 30 markers of the previous report, five are exogenous compounds whose identity the data support; the exogenous identities of most others are not supported. Thirteen markers meet all finding criteria, among them the peptides and oxylipin glucuronides, but several identities need revision (Gln-Pro instead of N-mevalonylcreatinine; grayanotoxin I and EFTPAVAH inconsistent with m/z; m/z 269.088 unidentified).
4. **PC markers.** No valid PC markers were found in either processing once the exogenous series and multiple testing are accounted for.

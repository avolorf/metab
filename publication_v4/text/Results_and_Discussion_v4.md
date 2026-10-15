# 3. Results and Discussion

> Draft for *Analytica* (MDPI): "Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study". Version 4 (15 October 2026): streamlined. The homologous exogenous series is no longer a separate result; it is handled as a data-processing step (Section 2) and documented in the Supplementary Materials. [ABC] = Frolova et al., *Anal. Bioanal. Chem.* 2026, 418, 6169–6184, doi:10.1007/s00216-026-06678-w.
> Figures and tables for this version: folder `publication_v4/` (Figure 1 without panel d; Figures 2–5 renumbered; series → Figure S6; "exogenous-signal index" throughout; Figure 4 with the sex-only baseline only; index baseline in Table S2). Section 3.6 and Figure 5 / Table S3 are optional and can be removed together.
> Numbered references [1]–[7] are listed at the end; renumber them together with the Introduction when the manuscript is assembled.

## 3.1. Analytical performance and data processing

The study cohort comprised 94 urine samples (45 controls, 29 bladder cancer [BC] and 20 prostate cancer [PC] patients) analysed in one batch with six pooled QC injections. Of 14,693 aligned features, 8,867 were detected in ≥ 80% of QC injections and in ≥ 50% of the samples of at least one group. Of these, 7,403 had a QC relative standard deviation (RSD) ≤ 30% (median 8.6%; Figure 1b). QC injections clustered tightly in the PCA score plot, and their total ion intensity was stable over the sequence (Figure 1a, c), so no drift correction was applied. These values meet the commonly used acceptance criteria for pooled-QC-based untargeted assays [1].

Grouping co-eluting, correlated ions (ΔRT ≤ 0.02 min, Spearman ρ ≥ 0.8) reduced the 7,403 features to 5,792 compounds, 718 of them represented by two or more ions; all statistics and false-discovery-rate control were applied at the compound level. After PQN normalisation [2], 11.8% of values were imputed. Forty-one features behaved as presence/absence variables and were analysed separately (Table S5).

In line with the data-processing algorithm (Section 2), compounds of exogenous origin — drugs and their metabolites, food- and tobacco-related compounds, and background ions — were flagged and excluded from marker selection (Figure 2, crosses). This step included a homologous series of non-endogenous compounds (repeating C2H4O units; Figure S6, Table S8), whose level varied between samples; part of the series was also present in solvent blanks (Figure S5). Because such variation is unrelated to disease but can differ between groups, the summed level of the series (exogenous-signal index) was used as a covariate in all comparisons.

## 3.2. Cohort characteristics

The groups differed in age and sex (Table 1; Figure S2a). Controls were younger (53.3 ± 13.9 years) than BC (62.9 ± 9.7 years; p = 0.006) and PC patients (70.7 ± 6.6 years; p < 0.001). Controls were mostly female (18 male/25 female), BC patients mostly male (24/5; p < 0.001) and all PC patients were male. Age was missing for four controls and sex for two; these participants were excluded from the adjusted models only (complete-case analysis, n = 88).

Each comparison was therefore evaluated with an adjusted linear model (group + sex + age + exogenous-signal index) and confirmed in 1:1 age- and sex-matched pairs (18 BC–control pairs, 18 BC–PC pairs; Figure S2b). A compound was reported as a finding only if it passed both checks and the low-index control comparison (Section 2).

## 3.3. Bladder cancer versus controls

In the adjusted model, 295 of 5,792 compounds had q < 0.05, and 212 of them passed all robustness checks: 117 were higher and 95 lower in BC (Figure 2a, Figure S1; Table S1). Normalisation by MSTUS instead of PQN gave rank-identical fold changes (ρ = 1.00), and 120 of the 125 MSTUS findings were among the PQN findings (Figure S3; Table S7).

Seven findings were annotated (Table 2; Figure 3):
- **MSI level 2** (precursor within 5 ppm of the curated formula and MS/MS library match): dodecanoylcarnitine, biliverdin, LysoPC(16:0) and Val-Leu;
- **MSI level 3**: a C19H28O3 oxo-androstane steroid (putatively 11-oxoetiocholanolone; [M+H−H2O]+ and [M+H−2H2O]+ at the same retention time; log2FC 1.25, q = 0.006, AUC 0.80) whose MS/MS matched two isomeric library entries, and two isomeric tetrahydrocortisone glucuronides (C27H40O11; [M+H]+ 541.2642 at 11.06 min and [M+NH4]+ 558.2914 at 10.76 min; log2FC 1.14 and 0.87), each co-eluting with aglycone ions whose MS/MS matched tetrahydrocortisone.

The remaining 205 findings are at MSI level 4. All seven annotated findings were higher in BC (AUC 0.69–0.80). Isomers not distinguished by MS/MS are listed in Table 2 (e.g., Leu-Val, Val-Ile and Ile-Val for Val-Leu; sn-1/sn-2 regioisomers for LysoPC).

The annotated findings fall into three groups.

1. **Glucocorticoid and androgen metabolites** (tetrahydrocortisone glucuronides, C19H28O3 steroid). Altered urinary androgen and corticosteroid metabolites have recently been reported in BC [3]. The increase remained significant after adjustment for sex and age and in matched pairs (p = 0.011 and 0.006 for the two glucuronides; p < 0.001 for the steroid). Because glucocorticoid output also rises with physiological stress, a non-specific illness response cannot be excluded.
2. **Acylcarnitines** (dodecanoylcarnitine). Increased urinary carnitine and acylcarnitines in BC have been reported in several cohorts and linked to altered fatty-acid β-oxidation [4,5].
3. **Blood- and protein-derived compounds** (LysoPC(16:0), biliverdin, Val-Leu). Lysophosphatidylcholines are abundant plasma lipids, biliverdin is a heme degradation product, and short peptides reflect proteolysis or proteinuria. Macroscopic haematuria was absent or mild in this cohort; microscopic haematuria and proteinuria were not quantified and may contribute to these differences.

No pathway reached FDR < 0.05 in mummichog analysis (Figure S4), so pathway-level conclusions are not drawn.

## 3.4. Classification performance

Classifiers were evaluated with 20 × 5-fold cross-validation, with feature selection nested inside the training folds and significance tested with 1,000 label permutations. Exogenous compounds and compounds correlated with the exogenous-signal index (|ρ| ≥ 0.5) were removed from the candidate pool beforehand.

For BC versus controls, the five-compound panel reached AUC 0.81 ± 0.03 (p < 0.001) and PLS-DA 0.83 ± 0.02 (Figure 4a; Table S2). With features adjusted for the exogenous-signal index inside each training fold, the panel kept AUC 0.78 ± 0.05 (p = 0.002), above the sex-only baseline evaluated with the same cross-validation (0.65 ± 0.02). In 18 age- and sex-matched pairs the panel reached 0.71 ± 0.05 (p = 0.059). The estimate on all samples (0.81) should therefore be regarded as an upper bound, and confirmation in a larger matched cohort is needed.

## 3.5. Bladder versus prostate cancer, and prostate cancer versus controls

**BC vs PC.** Three unannotated compounds (MSI 4), all lower in BC, met all criteria. The panel separated BC from PC with AUC 0.86 ± 0.04 (p = 0.002), but after in-fold adjustment for the exogenous-signal index the separation was not significant (0.63 ± 0.07, p = 0.11; Figure 4b). PLS-DA on all retained compounds gave AUC 0.72 (p = 0.011; Table S2), but this model is not adjusted for the exogenous-signal index. A robust metabolic difference between the two cancers was therefore not established.

**PC vs controls.** No compound reached q < 0.05 after adjustment, and the panel did not discriminate PC from controls (AUC 0.54 ± 0.06, p = 0.35; Figure 4c). The higher AUC after in-fold adjustment for the exogenous-signal index (0.74; Figure 4c) is not interpreted as a PC-specific signal, because the index itself nearly separates PC from controls (Table S2), so that the adjusted estimate extrapolates beyond the overlap of the two groups; in the sex-homogeneous comparison with male controls, the adjustment did not increase the AUC. PC patients differed from controls in age and sex (all PC patients were male and the oldest group), and in the comparison of PC with male controls the panel AUC was 0.61–0.62 (Table S2). PLS-DA on all retained compounds reached higher values (0.78 vs all controls, 0.77 vs male controls; Table S2); because this model uses all compounds without adjustment for age, sex or the exogenous-signal index, and no compound passed the adjusted analysis, these values were not regarded as evidence of PC-specific metabolites. Endogenous PC-specific urinary compounds were thus not found under the conditions of this study; a sex- and age-balanced design will be required to assess them.

## 3.6. Comparison with previously reported bladder cancer markers (optional)

The 14 metabolites reported in our previous study ([ABC], Table 5) were searched within ± 5 ppm of the reported or theoretical m/z and ± 0.5 min of the reported retention time (Figure 5; Table S3); the retention-time window allowed for shifts between analytical campaigns performed at different times with the same LC–HRMS method (observed shifts −0.43 to +0.19 min).
- **Not detected:** seven markers.
- **Below the detection filter:** two (Leu-Leu-Arg, 2-iminobiotin).
- **Quantified:** five. Val-Leu (MS/MS-annotated feature at 3.72 min; MSI 2) was higher in BC (log2FC 1.99, AUC 0.76, p = 1.4 × 10⁻⁴), remained significant in the adjusted model (q = 0.048) and in age-matched pairs (p = 0.004), and was one of the seven annotated findings. Glutaminylproline and LysoPE(20:4) changed in the same direction but not significantly, Leu-Thr-Ser in the opposite direction (p = 0.035), and 5-hydroxyomeprazole, a drug metabolite, did not differ between groups (p = 0.84).

The limited overlap is consistent with the small size of the earlier age-matched subset, the presence/absence behaviour of several earlier markers (log2FC 6–11), which is sensitive to detection limits, and the putative level of annotation. The replication of Val-Leu, supported by MS/MS, points to peptide-related changes in BC urine, which may, however, reflect proteinuria or microscopic haematuria rather than tumour proteolysis.

## 3.7. Sensitivity analyses

The main conclusions were robust to the analytical choices (Table S7):
1. **Samples with a PQN factor < 0.2** (3 BC, 3 PC, 1 control): excluding them increased the number of findings to 371 (190 of the 212 retained) and left the panel AUC unchanged (0.80).
2. **Alternative index definition** (series members absent from the blanks only): 163 findings (159 shared); five of the seven annotated findings retained, biliverdin (q = 0.057) and Val-Leu (q = 0.073) just above the FDR threshold; panel AUC 0.79 (p < 0.001).
3. **MSTUS normalisation:** see Section 3.3.

## 3.8. Limitations

1. **Single batch, unbalanced groups.** Age and sex distributions differed between groups, so matched analyses were limited to 18 pairs for BC and were not possible for PC.
2. **Blanks.** Only solvent blanks were included; process blanks would further support the assignment of signal origin.
3. **Clinical covariates.** Microscopic haematuria, proteinuria, renal function, medication, smoking and diet were not recorded.
4. **Annotation.** No reference standards were analysed; most findings remain at MSI level 4, and level 2–3 annotations require confirmation with authentic standards [6].

## 4. Conclusions (draft)

- The DoE-optimised sample preparation and LC–HRMS workflow, applied to a new cohort, gave a median QC RSD of 8.6%.
- The data-processing algorithm — compound-level grouping, exclusion of exogenous compounds, adjustment for age and sex, matched-pair checks and nested cross-validation with permutation testing — provided conservative estimates of group differences.
- Urinary differences between BC patients and controls remained after adjustment and involved glucocorticoid and androgen metabolites, an acylcarnitine and blood- or protein-derived compounds; Val-Leu was confirmed from our previous study. The classification performance (AUC 0.81; ≈ 0.71 in matched pairs) warrants validation in a larger, demographically balanced cohort.
- No endogenous PC-specific compounds were found; PC-related comparisons require an age- and sex-balanced design.

## References

[1] Broadhurst, D.; Goodacre, R.; Reinke, S.N.; et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* **2018**, 14, 72. doi:10.1007/s11306-018-1367-3

[2] Dieterle, F.; Ross, A.; Schlotterbeck, G.; Senn, H. Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. Application in 1H NMR metabonomics. *Anal. Chem.* **2006**, 78, 4281–4290. doi:10.1021/ac051632c

[3] Saygın, H.; Bolat, S.; Kablan, D.; et al. Metabolic profiling reveals insights into bladder cancer pathogenesis and recurrence. *Adv. Urol.* **2026**, 2026, 6524384. doi:10.1155/aiu/6524384

[4] Kim, W.T.; Yun, S.J.; Yan, C.; et al. Metabolic pathway signatures associated with urinary metabolite biomarkers differentiate bladder cancer patients from healthy controls. *Yonsei Med. J.* **2016**, 57, 865–871. doi:10.3349/ymj.2016.57.4.865

[5] Wu, C.; Wei, X.; Huang, Z.; et al. Urinary microbiome dysbiosis is associated with an inflammatory environment and perturbed fatty acids metabolism in the pathogenesis of bladder cancer. *J. Transl. Med.* **2024**, 22, 628. doi:10.1186/s12967-024-05446-7

[6] Sumner, L.W.; Amberg, A.; Barrett, D.; et al. Proposed minimum reporting standards for chemical analysis. *Metabolomics* **2007**, 3, 211–221. doi:10.1007/s11306-007-0082-2

[7] Burton, C.; Ma, Y. Current trends in cancer biomarker discovery using urinary metabolomics: achievements and new challenges. *Curr. Med. Chem.* **2019**, 26, 5–28. doi:10.2174/0929867324666170914102236 [optional, for the Introduction]

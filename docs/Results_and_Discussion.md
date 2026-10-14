# 3. Results and Discussion

> Draft for *Analytica* (MDPI). [ABC] = Frolova et al., *Anal. Bioanal. Chem.* 2026, 418, 6169–6184, doi:10.1007/s00216-026-06678-w (stage 1).
> Numbered references [1]–[8] are listed at the end; insert them in the journal's reference style.
> Items in square brackets marked "authors" must be confirmed before submission.

## 3.1. Analytical performance and data quality

The independent cohort comprised 94 urine samples (45 controls, 29 bladder cancer [BC] and 20 prostate cancer [PC] patients) analysed in one batch with six pooled QC injections. Of 14,693 aligned features, 8,867 were detected in ≥ 80% of QC injections and in ≥ 50% of the samples of at least one group. Of these, 7,403 had a QC relative standard deviation (RSD) ≤ 30% (median QC RSD 10.4% before and 8.6% after the RSD filter; Figure 1b). The QC injections clustered tightly in the PCA score plot, and the total ion intensity of QC injections was stable over the sequence (Figure 1a, c). A leave-one-QC-out test showed that LOWESS drift correction would reduce the prediction error only marginally (0.081 vs 0.074 log2 units). With six QC injections such a correction cannot be validated independently, so it was not applied. These figures meet the commonly used acceptance criteria for pooled-QC-based untargeted assays [1].

Grouping co-eluting, correlated ions (ΔRT ≤ 0.02 min, Spearman ρ ≥ 0.8) reduced the 7,403 features to 5,792 compounds, 718 of them represented by two or more ions. All statistics and false-discovery-rate control were therefore applied at the compound level. After PQN normalisation [2], 11.8% of values were imputed. Forty-one features behaved as presence/absence variables and were analysed separately (Table S5).

## 3.2. Cohort characteristics: a confounded design

The groups differed in age and sex (Table 1; Figure S2a). Controls were younger (53.3 ± 13.9 years) than BC (62.9 ± 9.7 years; p = 0.006) and PC patients (70.7 ± 6.6 years; p < 0.001). Controls were mostly female (18 male/25 female), whereas BC patients were mostly male (24/5; p < 0.001) and all PC patients were male. Age was missing for four controls and sex for two. These participants were excluded from the adjusted models only (complete-case analysis, n = 88; Table 1 note).

Because of this imbalance, every comparison was evaluated in three ways:
- with an adjusted linear model (group + sex + age + contamination index);
- in 1:1 age- and sex-matched pairs: 18 BC–control pairs (13 male, 5 female), 18 BC–PC pairs and 9 PC–male-control pairs (Figure S2b);
- against controls with low contamination.

A compound was reported as a "finding" only if it passed all three checks (Section 2.x).

## 3.3. Plastic-derived surfactant contamination

### 3.3.1. Identification and association with groups

A homologous series of alkyl diethanolamines and their oxygenated and ethoxylated analogues was detected:
- mass increments of C2H4 (28.031), O (15.995) and C2H4O (44.026);
- mass errors < 1.5 ppm;
- systematic increase of retention time with chain length (Figure 2a).

Such surfactants are typical additives of laboratory and medical plastics, which can leach into samples and interfere with LC–MS analyses [3]. The contamination index (mean z-score of ten series members) was:
- highest in PC (mean +0.87);
- intermediate in BC (−0.02);
- lowest in controls (−0.37; Figure 2b).

On its own, in the same cross-validation as the metabolite models, the index separated PC from controls with AUC 0.86 ± 0.00 and from male controls with AUC 0.93 ± 0.01. It separated BC from PC with AUC 0.78 ± 0.02, but BC from controls only with AUC 0.59 ± 0.02 (p = 0.087; Figure 5). The index was also inversely related to the PQN dilution factor (ρ = −0.45).

### 3.3.2. Origin: solvent blanks versus samples

To establish where the series comes from, extracted-ion chromatograms of all series members were obtained from the raw data of batch 1. The comparison included seven methanol solvent blanks (two of them injected before the first sample), the six pooled QC injections and five study samples: the two controls with the lowest and the three PC patients with the highest index (Figure S5; Table S8). For C12-DEA, the peak areas from these chromatograms agreed with the MS-DIAL areas used in the statistical analysis (QC median 8.36 × 10⁶ vs 8.52 × 10⁶).

Two components were distinguished.

1. **A constant system/solvent background.** C12-DEA, C14-DEA, C18-DEA+O and C20-DEA+O were present in every solvent blank, with retention times matching those in the samples. Their blank areas were:
   - C12-DEA: 17% of the QC median;
   - C14-DEA: 46%;
   - C18-DEA+O: 36%;
   - C20-DEA+O: 61%.

   The two blanks injected before any sample contained the same amounts as the later blanks. Endogenous reference ions (a C19H28O3 steroid and tetrahydrocortisone) were absent from all blanks. Both observations exclude carry-over and point to a background from the solvent or the LC–MS system.

2. **A sample-borne component.** The oxygenated and ethoxylated members (C12-DEA+O, C14-DEA+O, C16-DEA+O, C14-DEA+C2H4O2 and C14-DEA+(C2H4O)2+O) were absent from, or < 3% of QC in, the solvent blanks. In raw, non-normalised areas, the five study samples differed by orders of magnitude:
   - C14-DEA+O: 0.07 × 10⁶ in the least contaminated control vs 16–19 × 10⁶ in two PC samples;
   - C12-DEA: 1.8 × 10⁶ in the least contaminated control (close to the blank level) vs 38–40 × 10⁶ in the same two PC samples.

   One of the most contaminated PC samples was strongly diluted (PQN factor 0.06). Its contamination is therefore not an artefact of normalisation.

A contamination index rebuilt only from the five members absent from the blanks was almost identical to the original index (ρ = 0.96). It reproduced the group pattern (controls −0.38, BC −0.01, PC +0.86) and the separation of PC from controls (AUC 0.84).

The series is therefore introduced with the samples, before or during sample preparation, and its level depends on the group. The blank design of this batch (solvent blanks only, without process or container blanks) cannot tell whether the source is the urine collection containers, storage or the laboratory consumables. Its strong association with PC suggests a pre-analytical difference between recruitment streams, such as containers, collection site or storage time [authors: compare collection procedures of the PC, BC and control groups].

## 3.4. Bladder cancer versus controls

In the adjusted model, 295 of 5,792 compounds had q < 0.05. Of these, 212 also passed the age-matched comparison and both contamination checks: 117 were higher and 95 lower in BC (Figure 3a, Figure S1; Table S1). Normalisation by MSTUS instead of PQN changed only a per-sample scaling: fold changes were rank-identical (ρ = 1.00), and 120 of the 125 MSTUS findings were among the PQN findings (Figure S3; Table S7).

Annotation confidence was limited (Table 2; Figure 4).
- **MSI level 2** (precursor within 5 ppm of the curated formula and MS/MS library match) — four findings: dodecanoylcarnitine, biliverdin, LysoPC(16:0) and Val-Leu.
- **MSI level 3** — three findings:
  - a C19H28O3 oxo-androstane steroid (putatively 11-oxoetiocholanolone; [M+H−H2O]+ and [M+H−2H2O]+ at the same retention time; log2FC 1.25, q = 0.006, AUC 0.80), whose MS/MS matched two isomeric library entries;
  - two isomeric tetrahydrocortisone glucuronides (C27H40O11; [M+H]+ 541.2642 at 11.06 min and [M+NH4]+ 558.2914 at 10.76 min; log2FC 1.14 and 0.87). Each co-eluted with aglycone ions whose MS/MS matched tetrahydrocortisone.
- **MSI level 4** — the remaining 205 findings, including one library match judged implausible in human urine (a synthetic thiophene building block matched to an in-source fragment).

All seven annotated findings were higher in BC, with AUC 0.69–0.80 (Table 2). The isomers that the MS/MS spectra cannot distinguish are listed with each compound in Table 2, for example Leu-Val, Val-Ile and Ile-Val for Val-Leu, and sn-1/sn-2 regioisomers for LysoPC.

The annotated findings fall into three groups.

1. **Glucocorticoid and androgen metabolites** (tetrahydrocortisone glucuronides, C19H28O3 steroid). Altered urinary steroid hormone metabolites, including androgen and corticosteroid metabolites, have recently been reported in BC [4]. Glucocorticoid output also rises with physiological stress and differs by sex and age. The increase survived adjustment for sex and age and remained significant in matched pairs (p = 0.011 and 0.006 for the two glucuronides; p < 0.001 for the C19H28O3 steroid). Even so, a non-specific stress or illness response cannot be excluded.
2. **Acylcarnitines** (dodecanoylcarnitine). Increased urinary carnitine and acylcarnitines in BC have been reported in independent cohorts and linked to altered fatty-acid β-oxidation [5,6].
3. **Blood- and protein-derived compounds** (LysoPC(16:0), biliverdin, Val-Leu). Lysophosphatidylcholines are abundant plasma lipids, biliverdin is a heme degradation product, and short peptides indicate proteolysis or proteinuria. Haematuria is the cardinal sign of BC. The most parsimonious explanation for these three compounds is therefore blood and protein in the urine of BC patients rather than tumour metabolism [authors: haematuria/proteinuria status, if recorded, should be added as a covariate].

No pathway reached FDR < 0.05 in mummichog analysis (smallest q = 0.57 for BC vs controls; Figure S4). Pathway-level conclusions are therefore not drawn.

## 3.5. Classification performance

Classifiers were evaluated with repeated cross-validation. Feature selection was nested inside the training folds, and significance was tested with 1,000 label permutations. Exposure-related compounds and compounds correlated with the contamination index (|ρ| ≥ 0.5) were removed from the candidate pool beforehand.

For BC versus controls:
- the five-compound panel reached AUC 0.81 ± 0.03 (p < 0.001), and PLS-DA reached 0.83 ± 0.02 (Figure 5a; Table S2);
- with the features residualised on the contamination index inside each training fold, the panel kept AUC 0.78 ± 0.05 (p = 0.002);
- single-covariate baselines, evaluated with exactly the same cross-validation, were weaker: sex 0.65 ± 0.02 and contamination index 0.59 ± 0.02.

The confounder-balanced designs gave lower and less certain estimates:
- in 18 age- and sex-matched pairs the panel reached 0.71 ± 0.05 (p = 0.059);
- in 15 pairs matched on the contamination index within sex, in which the index itself had no discriminating power (AUC 0.32), it reached 0.71 ± 0.07 (p = 0.073).

The two matched estimates agree (≈ 0.71) but, with 15–18 pairs, do not reach significance. A diagnostic value of the panel that is independent of age, sex and contamination is therefore not established. The 0.81 obtained on all samples should be regarded as an upper bound.

## 3.6. Bladder versus prostate cancer, and prostate cancer versus controls

**BC vs PC.** Three compounds met all criteria. All were lower in BC and unannotated (MSI 4). The panel separated BC from PC with AUC 0.86 ± 0.04 (p = 0.002). The separation collapsed when contamination was taken into account (Figure 5b):
- after in-fold adjustment for the index: 0.63 ± 0.07 (p = 0.11);
- in 13 pairs matched on the index within sex: 0.57 ± 0.06 (p = 0.28).

The index alone gave 0.78. The apparent BC–PC difference is therefore attributed mainly to the different contamination of the two groups rather than to tumour-specific metabolism.

**PC vs controls.** No compound reached q < 0.05 after adjustment, and the panel did not discriminate PC from controls (AUC 0.54 ± 0.06, p = 0.35; Figure 5c). The index alone gave 0.86 and sex alone 0.75. Matched analyses were not possible: nine age-matched and five contamination-matched male pairs. After residualisation on the contamination index, the panel AUC rose to 0.74 ± 0.07 (p = 0.012). We examined this increase in detail:
- residualisation coefficients were estimated in the training folds only, and the procedure gave AUC 0.45 on random features and 0.58 with permuted index values, so it does not create signal by itself;
- the increase was absent in the sex-homogeneous comparison of PC with male controls (0.62 → 0.61);
- it became larger, up to 0.85 (the level of the index alone), when the adjustment slope was estimated within groups.

When a covariate nearly separates two groups, the "adjusted" difference of every feature that co-varies with it is an extrapolation across the gap between the groups. Part of the covariate's group difference then re-enters the adjusted features. The adjusted AUC is therefore not interpreted as evidence of a PC-specific signal. More generally, statistical adjustment cannot rescue a comparison in which the groups barely overlap in a strong confounder. Such confounders must be balanced by design.

## 3.7. Replication of the stage-1 markers

The 14 metabolites reported in stage 1 ([ABC], Table 5) were searched in the full feature table within ± 5 ppm of the reported or theoretical m/z and ± 0.5 min of the reported retention time (Figure 6; Table S3).
- **Not detected:** seven markers.
- **Below the detection filter:** two (Leu-Leu-Arg, 2-iminobiotin).
- **Quantified:** five:
  - Val-Leu was replicated. Of two isomeric candidates, the feature annotated as Val-Leu by MS/MS (3.72 min; MSI 2) was higher in BC (log2FC 1.99, AUC 0.76, p = 1.4 × 10⁻⁴), remained significant in the adjusted model (q = 0.048) and in age-matched pairs (p = 0.004), and was one of the seven annotated findings;
  - glutaminylproline and LysoPE(20:4) changed in the same direction but not significantly;
  - Leu-Thr-Ser changed in the opposite direction (log2FC −1.03, p = 0.035);
  - 5-hydroxyomeprazole was found only at its theoretical m/z (362.1171, +0.5 ppm, ΔRT −0.06 min), indicating that the m/z listed in Table 5 of [ABC] is in error, and did not differ between groups (p = 0.84).

For Leu-Leu-Tyr and Leu-Thr-Ser, the formulas of the named peptides were used, because those printed in [ABC] do not correspond to these peptides.

One of 14 stage-1 markers was thus replicated. Possible reasons for the low replication rate:
- stage-1 markers were selected by VIP in a small age-matched subset;
- several had very large stage-1 fold changes (log2FC 6–11), suggesting presence/absence behaviour that is sensitive to detection limits;
- annotation was putative;
- urine collection and handling may have differed between the cohorts.

The retention-time window of ± 0.5 min still needs to be justified by the retention shift between the two stages for compounds identified by MS/MS in both datasets [authors: provide the list of MS/MS-identified compounds of stage 1; the pipeline computes the median and range of the shift]. The replication of Val-Leu, together with its consistent MS/MS annotation, supports peptide-related changes in BC urine. Section 3.4, however, suggests that these changes may reflect proteinuria or haematuria rather than tumour proteolysis.

## 3.8. Sensitivity analyses

The main conclusions were robust to the analytical choices (Table S7).

1. **Samples with a PQN factor < 0.2** (3 BC, 3 PC, 1 control) were retained in the main analysis. Excluding them increased the number of significant BC–control compounds to 490 (277 shared) and of findings to 371 (190 of the 212 retained), and left the panel AUC unchanged (0.80).
2. **Index built from the five members absent from the blanks.** With this index as covariate, 206 compounds had q < 0.05 (205 shared with the main analysis) and 163 were findings (159 shared). Five of the seven annotated findings were retained; biliverdin (q = 0.057) and Val-Leu (q = 0.073) fell just above the FDR threshold. The contamination-adjusted BC panel reached AUC 0.79 (p < 0.001).
3. **MSTUS normalisation:** see Section 3.4.

## 3.9. Limitations

1. **Single batch, unbalanced groups.** The analysis is based on a single batch with unbalanced age and sex distributions. Matched analyses were therefore limited to 15–18 pairs for BC and were not possible for PC.
2. **Blank design.** Only solvent blanks were available. Process and container blanks would be needed to locate the source of the sample-borne contamination. In addition, three background-affected series members (QC-to-blank ratio 1.6–2.8) are present in the exported table [authors: confirm whether the MS-DIAL sample/blank filter was applied to this export, and with which blanks].
3. **Unrecorded clinical covariates.** Haematuria, proteinuria, renal function, medication, smoking and diet were not recorded and may contribute to the BC findings.
4. **Annotation.** No reference standards were analysed. Most findings remain at MSI level 4, and the levels 2–3 annotations require confirmation with authentic standards [7]. MS-DIAL match scores were not exported; the MS/MS level relies on the MS-DIAL identification flag.

## Conclusions (draft)

- **Bladder cancer.** In an independent cohort, urinary metabolic differences between BC patients and controls were detectable after adjustment for age, sex and plastic-derived contamination. They involved glucocorticoid and androgen metabolites, an acylcarnitine and blood- or protein-derived compounds, with Val-Leu replicated from stage 1. Their discriminating ability in confounder-balanced pairs (AUC ≈ 0.71) was modest and not significant.
- **Prostate cancer.** The apparent metabolic separation of PC from BC and from controls was explained by contamination introduced with the samples.
- **Practical recommendations** for urinary metabolomics in multi-group clinical cohorts:
  - standardised collection containers;
  - process and container blanks;
  - monitoring of contaminant series;
  - design-based balancing of strong confounders.

## References

[1] Broadhurst, D.; Goodacre, R.; Reinke, S.N.; et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* **2018**, 14, 72. doi:10.1007/s11306-018-1367-3

[2] Dieterle, F.; Ross, A.; Schlotterbeck, G.; Senn, H. Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. Application in 1H NMR metabonomics. *Anal. Chem.* **2006**, 78, 4281–4290. doi:10.1021/ac051632c

[3] Schauer, K.L.; Broccardo, C.J.; Webb, K.M.; Covey, P.A.; Prenni, J.E. Mass spectrometry contamination from Tinuvin 770, a common additive in laboratory plastics. *J. Biomol. Tech.* **2013**, 24, 57–61. doi:10.7171/jbt.13-2402-004

[4] Saygın, H.; Bolat, S.; Kablan, D.; et al. Metabolic profiling reveals insights into bladder cancer pathogenesis and recurrence. *Adv. Urol.* **2026**, 2026, 6524384. doi:10.1155/aiu/6524384

[5] Kim, W.T.; Yun, S.J.; Yan, C.; et al. Metabolic pathway signatures associated with urinary metabolite biomarkers differentiate bladder cancer patients from healthy controls. *Yonsei Med. J.* **2016**, 57, 865–871. doi:10.3349/ymj.2016.57.4.865

[6] Wu, C.; Wei, X.; Huang, Z.; et al. Urinary microbiome dysbiosis is associated with an inflammatory environment and perturbed fatty acids metabolism in the pathogenesis of bladder cancer. *J. Transl. Med.* **2024**, 22, 628. doi:10.1186/s12967-024-05446-7

[7] Sumner, L.W.; Amberg, A.; Barrett, D.; et al. Proposed minimum reporting standards for chemical analysis. *Metabolomics* **2007**, 3, 211–221. doi:10.1007/s11306-007-0082-2

[8] Burton, C.; Ma, Y. Current trends in cancer biomarker discovery using urinary metabolomics: achievements and new challenges. *Curr. Med. Chem.* **2019**, 26, 5–28. doi:10.2174/0929867324666170914102236 [optional, for the Introduction]

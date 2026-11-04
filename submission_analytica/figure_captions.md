# Figure and table captions — submission version (v5) ("Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study")

## Main text

**Table 1.** Clinical and demographic characteristics of the study cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). Footnote: see table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the 30% acceptance threshold. (c) PCA score plot of QC and study samples (log₂, Pareto scaling; one representative feature per compound; the 1,366 compounds detected in all study samples, i.e., without imputed values).

**Figure 2.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + exogenous-signal index). Filled circles, findings, i.e., compounds that met all finding criteria (Section 2.5): q < 0.05 in the adjusted model and p < 0.05 with the same direction of change in age- and sex-matched pairs (Wilcoxon signed-rank test); in comparisons with controls (a, c), also p < 0.05 with the same direction against controls with a below-median exogenous-signal index and no difference between controls with high and low index (p ≥ 0.05); crosses, exogenous compounds excluded from marker selection (drugs, food- and tobacco-related compounds, the homologous exogenous series); dashed line, FDR 5%. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (BC vs. controls; MSI levels 2–3). log₂FC and q (Benjamini–Hochberg) are from the model adjusted for sex, age and the exogenous-signal index; AUC, all samples; p (matched), Wilcoxon signed-rank test in age- and sex-matched pairs. Finding criteria and MSI levels are defined in Section 2.5; full data are in Tables S1 and S3.

**Figure 3.** Annotated bladder-cancer findings (the compounds of Table 2) across groups (log₂ PQN-normalized intensities). Open symbols, values below the detection limit, imputed as half of the feature minimum.

**Figure 4.** Cross-validated ROC curves (20×5-fold CV; top-5 compounds selected by Mann–Whitney U inside each training fold, logistic regression). Exogenous compounds and compounds correlated with the exogenous-signal index (|ρ| ≥ 0.5) were excluded from the candidate pool before CV. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot (a, b), compounds adjusted for the exogenous-signal index inside each training fold (for PC vs Control, see Table S2 and its note); dashed, 1:1 age- and sex-matched pairs (shown only with ≥ 10 pairs); dotted, sex-only baseline evaluated with the same CV. AUC, mean ± SD over repeats; p, permutation test (1000 permutations).

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log₂ intensities). Column bars: group, sex, age and exogenous-signal index. Rows are labeled with the curated annotation where available (putative, MSI level 3) and otherwise by *m/z* and retention time.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) Optimal 1:1 matching on age within sex (caliper ±5 years, without replacement).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labeled samples have PQN factor < 0.2. (b) Number of compounds with q < 0.05 under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, 5 ppm; exogenous features excluded). No pathway reached Benjamini–Hochberg q < 0.05.

**Figure S5.** Members of the homologous exogenous series in methanol solvent blanks, pooled QC and urine samples (raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas in solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest index; † members with solvent-blank area ≥ 10% of the QC median; n.d., not detected. (b, c) Extracted-ion chromatograms of a member present in the blanks (C12) and of one absent from the blanks (C14+O).

**Figure S6.** Homologous exogenous series used for the exogenous-signal index. (a) Members ([M+H]⁺; mass increments C₂H₄, O and C₂H₄O) in the *m/z*–retention-time plane; open diamonds, co-varying ions not included in the index. (b) Exogenous-signal index (mean z-score of the ten members) by group (Mann–Whitney U test).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, 5th–95th percentiles, Q², permutation p), including the single-covariate baselines (sex only; exogenous-signal index only).

**Table S2b.** Composition of the cross-validated top-5 panels: number of the 100 training folds (20 repeats × 5 folds) in which each compound was selected, and the five compounds selected by the same rule on all samples of each comparison. Compounds selected in ≥ 10% of folds or included in the all-sample panel are listed.

**Table S3.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (same rule as in Table 2), supporting ions, detection by group, and correlation with the exogenous-signal index and the dilution factor; annotations with mass error > 5 ppm are listed separately (S3b). No reference standards were analyzed (no MSI level 1).

**Table S4.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S5.** Compound-level statistics for all comparisons (all models and sensitivity checks).

**Table S6.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < 0.2, with the exogenous-signal index built only from series members absent from solvent blanks, restricted to compounds without imputed values, and with imputed values treated as missing; overlap with the main analysis and cross-validated panel AUC.

**Table S7.** Solvent-blank check of the homologous exogenous series: median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.

**Table S8.** Keyword list used to flag drug-, food- and tobacco-related compounds: regular expressions (case-insensitive) matched against MS-DIAL library names (NIST HR-MS/MS), compiled before the statistical analysis, and the number of compounds flagged by each expression (library matches that met the MS-DIAL identification criteria: dot-product scores ≥ 600, ≥ 25% matched spectrum, ≥ 3 matched peaks; compounds whose drug- or tobacco-related match was labeled "low score" were excluded only from the candidate pool of the classification models).

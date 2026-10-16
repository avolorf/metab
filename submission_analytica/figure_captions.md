# Figure and table captions — submission version (v5) ("Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study")

## Main text

**Table 1.** Clinical and demographic characteristics of the study cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). Footnote: see table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the 30% acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound; the 1,366 compounds detected in all study samples, i.e., without imputed values).

**Figure 2.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + exogenous-signal index). Filled circles, findings that passed all robustness criteria; crosses, exogenous compounds excluded from marker selection (drugs, food- and tobacco-related compounds, the homologous exogenous series); dashed line, FDR 5%. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3), steroid glucuronides first, then other compounds, each ordered by q. log2FC and q, adjusted model (sex, age, exogenous-signal index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-index controls), BC vs controls with an exogenous-signal index below the median. Findings met all of the following: q < 0.05 in the adjusted model; p < 0.05 in age-matched pairs and against low-index controls, with the same direction of change; and no difference between controls with high and low index (p ≥ 0.05), which excludes compounds whose level itself depends on the index. MSI 2, precursor within 5 ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the feature itself; MSI 3, formula-consistent annotation without MS/MS confirmation, MS/MS match to several isomers, or MS/MS match without a curated formula. "Isomers not excluded" lists isomers that the MS/MS spectrum does not distinguish. All findings, including MSI level 4, are listed in Table S1.

**Figure 3.** Annotated bladder-cancer findings (the compounds of Table 2) across groups (log2 PQN-normalized intensities). Open symbols, values below the detection limit, imputed as half of the feature minimum. q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 4.** Cross-validated ROC curves (20×5-fold CV; top-5 compounds selected by Mann–Whitney U inside each training fold, logistic regression). Exogenous compounds and compounds correlated with the exogenous-signal index (|ρ| ≥ 0.5) were excluded from the candidate pool before CV. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot (a, b), compounds adjusted for the exogenous-signal index inside each training fold (for PC vs Control, see Table S2 and its note); dashed, 1:1 age- and sex-matched pairs (shown only with ≥ 10 pairs); dotted, sex-only baseline evaluated with the same CV. AUC, mean ± SD over repeats; p, permutation test (1000 permutations).

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and exogenous-signal index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±5 years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < 0.2. (b) Number of compounds with q < 0.05 under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, 5 ppm; exogenous features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Figure S5.** Members of the homologous exogenous series in methanol solvent blanks, pooled QC and urine samples (raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas in solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest index; † members with solvent-blank area ≥ 10% of the QC median; n.d., not detected. (b, c) Extracted-ion chromatograms of a member present in the blanks (C12) and of one absent from the blanks (C14+O).

**Figure S6.** Homologous exogenous series used for the exogenous-signal index. (a) Members ([M+H]+; mass increments C2H4, O and C2H4O) in the m/z–retention-time plane; open diamonds, co-varying ions not included in the index. (b) Exogenous-signal index (mean z-score of the ten members) by group (Mann–Whitney U test).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, 5th–95th percentiles, Q², permutation p), including the single-covariate baselines (sex only; exogenous-signal index only).

**Table S3.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (same rule as in Table 2), supporting ions, detection by group, and correlation with the exogenous-signal index and the dilution factor; annotations with mass error > 5 ppm are listed separately (S3b). No reference standards were analysed (no MSI level 1).

**Table S4.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S5.** Compound-level statistics for all comparisons (all models and robustness checks).

**Table S6.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < 0.2, with the exogenous-signal index built only from series members absent from solvent blanks, restricted to compounds without imputed values, and with imputed values treated as missing; overlap with the main analysis and cross-validated panel AUC.

**Table S7.** Solvent-blank check of the homologous exogenous series: median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.

**Table S8.** Keyword list used to flag drug-, food- and tobacco-related compounds: regular expressions (case-insensitive) matched against MS-DIAL library names (NIST 23 HRMS, MS-DIAL MSP), compiled before the statistical analysis, and the number of compounds flagged by each expression (high-scoring library matches; compounds with low-scoring drug- or tobacco-related matches were excluded only from the candidate pool of the classification models).

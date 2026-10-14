# Figure and table captions (draft)

## Main text

**Table 1.** Clinical and demographic characteristics of the independent cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). The footnote on missing age and sex is given in table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the 30% acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound). (d) The same scores coloured by the contamination index.

**Figure 2.** Plastic-derived surfactant contamination. (a) Homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+) used to build the contamination index; open diamonds, co-varying ions not included in the index. (b) Contamination index by group (Mann–Whitney U test).

**Figure 3.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + contamination index). Filled circles, findings that passed all robustness criteria; crosses, exposure-related compounds (contaminants, food, drugs); dashed line, FDR 5%. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3). log2FC and q, adjusted model (sex, age, contamination index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-contamination controls), BC vs controls with a contamination index below the median. MSI 2, precursor within 5 ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the feature itself; MSI 3, formula-consistent annotation without MS/MS confirmation, MS/MS match to several isomers, or MS/MS match without a curated formula. "Isomers not excluded" lists isomers that the MS/MS spectrum does not distinguish (e.g., Leu-Val, Val-Ile and Ile-Val for Val-Leu; sn-1/sn-2 regioisomers for LysoPC). For compounds without a curated annotation, the MS-DIAL library match is given for transparency. Library matches judged implausible in human urine (ethyl 2-amino-5-methyl-3-thiophenecarboxylate, a synthetic building block matched to an in-source fragment) were assigned MSI level 4. All findings, including MSI level 4, are listed in Table S1.

**Figure 4.** Annotated bladder-cancer findings (MSI levels 2–3; the compounds of Table 2) across groups (log2 PQN-normalized intensities). q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 5.** Cross-validated ROC curves (20×5-fold CV; top-5 features selected by Mann–Whitney U inside each training fold, logistic regression). Candidate pool defined before CV: exposure-related compounds (contaminants, food, drugs, tobacco) and compounds correlated with the contamination index (|ρ| ≥ 0.5) excluded. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot, features residualized on the contamination index inside each training fold; dashed, 1:1 age-matched pairs; dash-dot-dot, 1:1 contamination-matched pairs (within sex, caliper 0.25 SD of the index; sensitivity analysis); matched sets shown only with ≥ 10 pairs. p, permutation test (1000 permutations). Dotted lines, single-covariate baselines (sex only; contamination index only) evaluated with the same repeated CV (logistic regression on one covariate) and permutation test. AUC, mean ± SD over repeats. In (c) the contamination index alone nearly separates PC from controls; the contamination-adjusted curve is then an extrapolation beyond the overlap of the two groups and is not interpretable as a contamination-independent signal (see Section 3.6).

**Figure 6.** Stage-1 bladder-cancer markers (Frolova et al., Anal Bioanal Chem 2026, Table 5) in the independent cohort. Features were matched by m/z (±5 ppm) and retention time (±0.5 min); all isomeric candidates are shown. Points, unadjusted log2 fold change (BC vs Control) of quantified features; coloured, same direction and Mann–Whitney p < 0.05; text gives AUC, unadjusted p and FDR-adjusted q of the adjusted model. Markers not detected or removed by the detection/QC filters are labelled. Exact values are given in Table S3.

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and contamination index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±5 years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < 0.2. (b) Number of compounds with q < 0.05 under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, 5 ppm; exposure-related features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Figure S5.** Alkyl diethanolamine series in solvent blanks, pooled QC and urine samples (batch 1; raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas of each series member in methanol solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest contamination index; † members with solvent-blank area ≥ 10% of the QC median (system/solvent background); n.d., not detected. (b, c) Extracted-ion chromatograms of C12-DEA (present in blanks) and its oxygenated analogue C14-DEA+O (absent from blanks).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, Q², permutation p).

**Table S3.** Stage-1 markers: all candidate features within ±5 ppm of the reported or the theoretical m/z (formula and ion) and ±0.5 min of the reported retention time, detection/filter status and statistics in the independent cohort. For Leu-Leu-Tyr and Leu-Thr-Ser the formula of the named peptide was used to compute the theoretical m/z.

**Table S4.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (assigned by the same rule as in Table 2), supporting ions, detection by group, and correlation with the contamination index and the dilution factor; annotations with mass error > 5 ppm are listed separately (S4b). No reference standards were analysed (no MSI level 1). MS-DIAL match scores were not part of the exported table; the MS/MS match is the MS-DIAL identification flag.

**Table S5.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S6.** Compound-level statistics for all comparisons (all models and robustness checks).

**Table S7.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < 0.2, and with a contamination index built only from the series members that were absent from solvent blanks; overlap with the main analysis and cross-validated panel AUC.

**Table S8.** Solvent-blank check of the alkyl diethanolamine series (batch 1): median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.

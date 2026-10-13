# Figure captions (draft)

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the 30% acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound). (d) The same scores coloured by the contamination index.

**Figure 2.** Plastic-derived surfactant contamination. (a) Homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+) used to build the contamination index; open diamonds, co-varying ions not included in the index. (b) Contamination index by group (Mann–Whitney U test).

**Figure 3.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + contamination index). Filled circles, findings that passed all robustness criteria; crosses, exposure-related compounds (contaminants, food, drugs). Dashed line, FDR 5%.

**Figure 4.** Heatmap of bladder-cancer findings (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and contamination index.

**Figure 5.** Cross-validated ROC curves (20×5-fold CV; top-5 features selected by Mann–Whitney U inside each training fold, logistic regression). Shaded band, 10th–90th percentile across repeats; p, permutation test (1000 permutations). Dotted lines, single-covariate baselines (sex; contamination index).

**Figure 6.** Top bladder-cancer findings across groups (log2 PQN-normalized intensities). q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure S1.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±5 years).

**Figure S2.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < 0.2. (b) Number of compounds with q < 0.05 under PQN and MSTUS normalization and their overlap.

**Figure S3.** Pathway enrichment (mummichog, positive mode, 5 ppm; exposure-related features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

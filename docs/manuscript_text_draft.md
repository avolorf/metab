# Draft text for the validation-stage manuscript (Analytica, MDPI)

> **Status:** numbers below come from `run_all.py` with `metadata/metadata.csv` (ages available for 23/45 controls and 15/20 PC patients).
> They **will change** after re-running with the corrected ages (`metadata_stage2.csv`). Everything in [square brackets] must be completed by the authors.
> Stage 1 = Frolova et al., *Anal Bioanal Chem* 2026, 418:6169–6184, doi:10.1007/s00216-026-06678-w (cited below as [ABC]).

---

## 2. Materials and Methods (data processing and statistics)

### 2.x Study design of the independent cohort
First-morning mid-stream urine samples were collected from [n] patients with bladder cancer (BC), [n] patients with prostate cancer (PC) and [n] controls [collection site/period, ethics approval as in ABC or new approval number]. One patient with a benign bladder lesion and PC was assigned to the PC group; one patient with multiple primary tumours (PC, BC and skin cancer) and patients with kidney (n = 3) or pancreatic (n = 2) cancer were excluded. The final analysis set comprised 45 controls, 29 BC and 20 PC patients (Table 1). Group membership, sex, age and exclusions were read from a single metadata file; no sample subsets were defined manually.

### 2.x LC-HRMS analysis and feature extraction
Samples were prepared and analysed with the DoE-optimised protocol and the LC-HRMS method described in [ABC] (Acquity UPLC BEH C18, 23-min gradient, Orbitrap Fusion Lumos, ESI+, ddMS2). Pooled QC samples (n = 6) were injected at the start of the sequence and after every 20 study samples. Peak picking, alignment and library annotation were performed in [MS-DIAL v. …, parameters]. [Procedural blanks: if available, state number and that the blank filter (QC/blank ≥ 5) was applied; currently no blanks were included in the exported table.]

### 2.x Data pre-processing
Zero intensities were treated as missing. Features were retained when detected in ≥ 80% of QC injections and in ≥ 50% of samples of at least one group (14,693 → 8,867 features) and when the relative standard deviation in QC samples was ≤ 30% (→ 7,403 features). Drift correction was not applied: leave-one-QC-out prediction error was 0.081 log2 units without correction and 0.074 with LOWESS, and with six QC injections the correction cannot be validated independently; QC RSD is therefore reported without correction (median 10.4% before and 8.6% after the RSD filter). Urine dilution was corrected by probabilistic quotient normalisation (PQN; reference = median QC spectrum); MSTUS normalisation (1,900 features detected in all samples) was computed in parallel as a sensitivity analysis. Missing values (11.8%) were replaced by half of the feature minimum and data were log2-transformed. The fraction of imputed values per group is reported for every feature; features detected in < 50% of one group and ≥ 80% of another (n = 41) were analysed separately as presence/absence features (Fisher's exact test). Seven samples with a PQN factor < 0.2 were flagged for chromatogram inspection and analysed in a sensitivity analysis without them.

### 2.x Grouping of features into compounds
Features eluting within 0.02 min and correlated across samples (Spearman ρ ≥ 0.8) were grouped into compounds (connected components); adduct, isotope and in-source-fragment mass differences were annotated. The 7,403 features formed 5,792 compounds (718 with ≥ 2 ions). Statistics and false-discovery-rate (FDR) control were applied at the compound level using the most intense ion as representative.

### 2.x Contamination index
A homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+, mass increments C2H4 28.031, O 15.995 and C2H4O 44.026; ten features, mass error < 1.5 ppm) was identified as a plastic-derived surfactant contamination. A contamination index was defined as the mean z-score of these ten features and was used as a continuous covariate and dichotomised at the median (high/low). A co-varying unknown C24H30O6 ([M+H]+ and [M+Na]+) is reported separately and is not part of the index.

### 2.x Exposure-related compounds
Compounds annotated as plastic contaminants, drugs, tobacco or food constituents (curated list and keyword rules) were labelled as exposure-related. Low-score library matches were reported as "unknown" in figures (putative names only in the Supplementary Tables). Exposure-related compounds could not qualify as findings and were excluded from classification models and pathway analysis.

### 2.x Statistical analysis
For each comparison (BC vs Control, primary; BC vs PC; PC vs Control; PC vs male controls) the main model was a linear model of log2 intensity with group, sex, age and contamination index as covariates (complete cases), with Benjamini–Hochberg FDR. Robustness was assessed by (i) 1:1 nearest-age matching within sex (caliper ± 5 years, optimal assignment; Mann–Whitney U test on matched sets), (ii) comparison with controls of low contamination index, and (iii) comparison of controls with high vs low contamination index. A compound was considered a **finding** when q < 0.05 in the main model, p < 0.05 in the age-matched comparison in the same direction, p < 0.05 versus low-contamination controls in the same direction, and p ≥ 0.05 between high- and low-contamination controls; exposure-related compounds were excluded.

Classification was evaluated by 20 × 5-fold stratified cross-validation. Before cross-validation, exposure-related compounds and compounds correlated with the contamination index (|ρ| ≥ 0.5) were removed from the candidate pool (5,351 compounds remained). Within each training fold the five compounds with the smallest Mann–Whitney p were selected and a logistic regression was fitted; PLS-DA (two components, Pareto scaling within folds) on all pool compounds was used for comparison. A third model residualised all features on the contamination index within each training fold. Significance was assessed by 1,000 label permutations. Single-covariate models (sex; contamination index) served as baselines. Age-matched cross-validation was performed only with ≥ 10 pairs.

Pathway enrichment (mummichog, positive mode, 5 ppm; top 10% of features ranked by the main-model p; exposure-related features excluded; 500 permutations; BH correction across pathways) is reported in the Supplementary Material only.

### 2.x Comparison with stage-1 markers
The 14 metabolites reported as altered in BC in stage 1 ([ABC], Table 5) were searched in the full, unfiltered feature table by m/z (± 5 ppm) and retention time (± 0.5 min; a systematic retention shift of −0.4 to +0.2 min between batches was observed for matched compounds). All candidate features within the window were reported (isomers flagged). For N-acetyltryptophan the theoretical [M+H]+ (m/z 247.1077) was used because the m/z listed in [ABC] Table 5 (249.0923) does not correspond to the stated formula; [authors: check also 5-hydroxyomeprazole, m/z 362.1248 vs. theoretical 362.1169, and the formulas given for Leu-Leu-Tyr and Leu-Thr-Ser]. Replication criteria, defined before the analysis, were: same direction and Mann–Whitney p < 0.05 (nominal) and same direction with q < 0.05 in the adjusted model (replicated).

### 2.x Annotation confidence
Annotation levels follow the Metabolomics Standards Initiative: level 3 for compounds with a molecular formula consistent within 5 ppm and supporting ions, but without isomer-specific MS/MS or a reference standard; level 4 for features with low-score or no library match. Low-score library matches that are implausible in human urine (e.g., aplaviroc, apicidin, flunixin) are listed in Table 2 for transparency only.

### 2.x Software and code availability
Python 3 (pandas, numpy, scipy, statsmodels, scikit-learn, matplotlib/seaborn) and mummichog 2.7. The complete pipeline, parameters and random seed are available at [repository URL/DOI]; a single command reproduces all tables and figures from the feature table and the metadata file.

---

## 3. Results

### 3.1 Data quality
QC samples clustered tightly in PCA (Figure 1c) and total ion intensity of QC injections was stable over the sequence (Figure 1a). After filtering, 7,403 features (5,792 compounds) with median QC RSD 8.6% were retained (Figure 1b).

### 3.2 Cohort characteristics
The groups differed in age (Kruskal–Wallis p = 3·10⁻⁵) and sex (BC vs Control, Fisher p = 6·10⁻⁴; Table 1). [Update after age correction: number of controls with known age, matched pairs.]

### 3.3 Plastic-derived contamination
The contamination index was highest in PC (mean +0.87 vs −0.37 in controls and −0.02 in BC; Figure 2). As a single variable it discriminated PC from controls with AUC 0.87 and from male controls with AUC 0.94, and BC from PC with AUC 0.80.

### 3.4 Bladder cancer vs controls
In the adjusted model, 70 compounds had q < 0.05 and 34 met all robustness criteria (Figures 3a, 4, 6; Table 2). Only one finding had MSI level 3: a C19H28O3 oxo-androstane steroid (putatively 11-oxoetiocholanolone) detected as [M+H−H2O]+ m/z 287.2007 and [M+H−2H2O]+ m/z 269.1901 at the same retention time (ρ = 0.99; 0.4–0.5 ppm); library MS/MS matching returned two isomeric steroid names, so the isomer remains unassigned. The remaining 33 findings are MSI level 4 and are reported as candidates for confirmation, not as a diagnostic panel. In cross-validation, the five-compound panel reached AUC 0.81 ± 0.03 (permutation p < 0.001) and 0.78 ± 0.05 (p = 0.002) after in-fold adjustment for contamination; baseline AUC was 0.70 for sex and 0.63 for the contamination index (Figure 5a). [Age-matched estimate after age correction; with the current ages only 8 pairs were available and the model was not estimated.]

### 3.5 Bladder vs prostate cancer
No compound reached q < 0.05 in the adjusted model. The unadjusted panel AUC of 0.86 ± 0.04 decreased to 0.63 ± 0.07 (permutation p = 0.11) after in-fold adjustment for the contamination index (age-matched, 14 pairs: 0.84 → 0.67, p = 0.11), while the contamination index alone gave AUC 0.80 (Figure 5b). The apparent separation of BC and PC is therefore attributed mainly to contamination rather than to tumour-specific metabolism.

### 3.6 Prostate cancer vs controls
No compound reached q < 0.05 after adjustment for sex, age and contamination. After removal of contamination-correlated compounds the panel did not discriminate PC from controls (AUC 0.54 ± 0.06, p = 0.35; Figure 5c). [Re-evaluate after age correction, including PC vs age-matched male controls.]

### 3.7 Comparison with stage-1 markers
Of the 14 stage-1 markers, 8 were not detected within ± 5 ppm and ± 0.5 min, 2 were detected below the detection filter (Leu-Leu-Arg, 2-iminobiotin) and 4 were quantified (Table 3, Figure 7). Val-Leu had two isomeric candidates; the candidate annotated as Val-Leu (3.72 min) was higher in BC (log2FC 1.99, AUC 0.76, p = 1.4·10⁻⁴) but not significant after adjustment (q = 0.34). Glutaminylproline and LysoPE(20:4) changed in the same direction without significance, and Leu-Thr-Ser changed in the opposite direction (log2FC −1.03, p = 0.035). Thus, one of 14 stage-1 markers was nominally reproduced and none met the adjusted replication criterion. [Discuss: different cohort composition and matching (25 vs 25 in stage 1), univariate selection by VIP/p without FDR in stage 1, large stage-1 fold changes suggesting presence/absence behaviour, and annotation level.]

### 3.8 Sensitivity analyses
PQN and MSTUS differed only by a per-sample scaling: compound fold changes were perfectly rank-correlated (ρ = 1.00), but MSTUS yielded fewer significant compounds (19 vs 70, all contained in the PQN set; 14 vs 34 findings, 12 shared) and similar panel AUC (0.80 vs 0.81). Excluding seven samples with PQN factor < 0.2 increased the number of significant compounds (162; 68 shared with the main analysis). No pathway reached FDR < 0.05 (Supplementary Figure S3).

---

## Supplementary Materials (proposed list)
- **Figure S1** Age by group and sex; age-matched pairs. **Figure S2** PQN vs MSTUS. **Figure S3** Pathway enrichment.
- **Table S1** Classification results (all models, permutation p). **Table S2** Full compound-level statistics for all comparisons (`results.xlsx`, sheets `stats_*`). **Table S3** Curated candidates incl. mass-mismatch annotations. **Table S4** Presence/absence features. **Table S5** Stage-1 marker matching, all candidates.

## Data availability (draft)
The feature table, metadata (pseudonymised), processing code and parameters are available at [repository]. Raw LC-HRMS files are available from the corresponding author on reasonable request [or: deposited in MetaboLights/Metabolomics Workbench, accession …].

# 2. Materials and Methods (v4)

> Draft for *Analytica* (MDPI): "Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study". Terminology matches Results and Discussion v4 and the files in `publication_v4/figures_tables/`.
> [ABC] = Frolova et al., *Anal. Bioanal. Chem.* 2026, 418, 6169–6184, doi:10.1007/s00216-026-06678-w.
> Items in square brackets marked "authors" must be completed or confirmed. Section 2.9 is optional and is removed together with Section 3.6, Figure 5 and Table S3.

## 2.1. Study cohort

First-morning mid-stream urine samples were collected from patients with bladder cancer (BC), patients with prostate cancer (PC) and controls [authors: collection site and period; ethics committee and approval number; informed consent]. One patient with a benign bladder lesion and PC was assigned to the PC group. Excluded were:
- one patient with multiple primary tumours (PC, BC and skin cancer);
- patients with kidney (n = 3) or pancreatic (n = 2) cancer.

The final analysis set comprised 45 controls, 29 BC and 20 PC patients (Table 1). Age was missing for four controls and sex for two controls. Group membership, sex, age and exclusions were read from a single metadata file, and no sample subsets were defined manually.

## 2.2. Sample preparation and LC–HRMS analysis

Samples were prepared and analysed with the DoE-optimised protocol and the LC–HRMS method described previously [ABC]:
- urine–methanol 1:2 (v/v), vortexing and centrifugation;
- evaporation under nitrogen and reconstitution in methanol–water 50:50 (v/v);
- Acquity UPLC BEH C18 column (2.1 × 100 mm, 1.7 µm), 23-min gradient;
- Orbitrap Fusion Lumos, positive electrospray ionisation, data-dependent MS/MS.

Pooled QC samples (n = 6) were injected at the start of the sequence and after every 20 study samples. Methanol solvent blanks (n = 7) were injected at the start of the sequence and before each QC injection.

## 2.3. Feature extraction and annotation

Raw data were processed in MS-DIAL 5.5 with the following settings:
- **feature detection:** minimum peak height 10,000 counts, mass slice width 0.1 Da, Savitzky–Golay smoothing (5 scans), minimum peak width 5 scans;
- **deconvolution:** sigma window 0.5, amplitude cut-off 100 counts;
- **alignment:** a representative QC sample as reference, retention-time and m/z tolerances 0.2 min and 0.015 Da (retention-time and m/z factors 0.5/0.5);
- **adduct annotation:** default ESI+ adduct and dimer library, supplemented with sodium formate clusters.

[Authors: the blank-filter sentence needs to be checked. In the exported table, several background ions have QC-to-blank ratios of 1.6–2.8 (Table S8). Use either (A) "No blank subtraction was applied in MS-DIAL; background ions were handled downstream (Section 2.6)", or (B) the actual MS-DIAL setting, if a sample/blank filter was applied with other blanks or another threshold.]

[Authors: keep one of the two sentences: (A) "A QC-based LOWESS signal correction was applied within MS-DIAL." or (B) "No QC-based signal correction was applied in MS-DIAL."]

MS/MS spectra were matched against NIST 23 HRMS and the MS-DIAL metabolomics MSP spectral kit, with additional annotation in MS-FINDER 3.6. The aligned peak-area table was exported and processed further in Python (Section 2.10).

## 2.4. Data pre-processing

**Filtering.**
- Zero intensities were treated as missing.
- Detection filter: features were retained when detected in ≥ 80% of QC injections and in ≥ 50% of the samples of at least one group (14,693 → 8,867 features).
- QC repeatability filter: features with a relative standard deviation (RSD) in QC samples ≤ 30% were retained (→ 7,403 features).

**Drift.** No drift correction was applied in the processing pipeline. The leave-one-QC-out prediction error was 0.081 log2 units without correction and 0.074 with LOWESS. With six QC injections, such a correction cannot be validated independently. QC RSD values refer to the exported table [if (A) in Section 2.3: "after the MS-DIAL QC correction"]. The median QC RSD was 10.4% for features passing the detection filter and 8.6% after the RSD filter.

**Normalisation and imputation.**
- Urine dilution was corrected by probabilistic quotient normalisation (PQN; reference, median QC spectrum).
- MSTUS normalisation (1,900 features detected in all samples) was computed in parallel as a sensitivity analysis.
- Missing values (11.8%) were replaced by half of the feature minimum, and the data were log2-transformed. The fraction of imputed values per group is reported for every feature.
- The data-quality PCA (Figure 1c) was computed on compounds detected in all samples (no imputed values; Pareto scaling).

**Special cases.**
- Features detected in < 50% of one group and ≥ 80% of another (n = 41) were analysed separately as presence/absence features (Fisher's exact test).
- Seven samples with a PQN factor < 0.2 were retained in the main analysis and excluded in a sensitivity analysis.
- The influence of imputation was assessed by repeating the analysis (i) on the compounds detected in all samples and (ii) with imputed values treated as missing, each compound being analysed on its detected values only (Table S7).

## 2.5. Grouping of features into compounds

Features were grouped into compounds (connected components) when they eluted within 0.02 min and were correlated across samples with Spearman ρ ≥ 0.8. Adduct, isotope and in-source-fragment mass differences were annotated. The 7,403 features formed 5,792 compounds, 718 of them with two or more ions. Statistics and false-discovery-rate (FDR) control were applied at the compound level. Each compound was represented by its curated intact-molecule ion ([M+H]+, [M+Na]+, [M+NH4]+ or [M+K]+) when such an ion was annotated in the group, and otherwise by its most intense ion, so that in-source fragments did not represent a compound.

## 2.6. Exogenous compounds and the exogenous-signal index

Compounds of exogenous origin were flagged in three ways:

1. **Curated annotations.** A manually curated annotation (Table S4) classified the compound as a drug or drug metabolite, or as a food- or tobacco-related constituent.
2. **Keyword list.** The compound's MS-DIAL library name matched a pre-specified keyword list of 21 regular expressions covering:
   - drugs commonly prescribed to urological and elderly patients, grouped by therapeutic class;
   - dietary methylxanthines, polyphenols, alkaloids and other food constituents;
   - nicotine and its metabolites.

   The list was compiled before the statistical analysis and applied to the library names returned by MS-DIAL; no external database was queried automatically.
3. **Background ions.** The compound belonged to a homologous series of non-endogenous compounds that differ by repeating C2H4 and C2H4O units and include their oxygenated analogues:
   - ten [M+H]+ ions, mass error < 1.5 ppm, retention times increasing with chain length (Figure S6a);
   - members of the series were also detected in methanol solvent blanks (Section 2.7).

Exogenous compounds could not qualify as findings and were excluded from the classification models and the pathway analysis. Low-score library matches are shown as "unknown" in the figures; their putative names are given only in the Supplementary Tables.

Because the level of the homologous series varied between samples and could differ between groups, it was entered as a covariate. The exogenous-signal index was defined as the mean z-score of the ten series members (Figure S6b). It was used as a continuous covariate and, for one robustness check, dichotomised at its median (high/low).

## 2.7. Solvent-blank check

Extracted-ion chromatograms (± 5 ppm) of all series members were obtained from the raw data (converted to mzML with ProteoWizard msconvert, vendor centroiding, MS1 only) for the following injections:
- the seven solvent blanks;
- the six QC injections;
- five study samples: the two controls with the lowest and the three PC samples with the highest index.

Peak areas were integrated around the expected retention time after per-injection retention-time alignment (Figure S5, Table S8). A series member was considered present in the blanks when its median blank area was ≥ 10% of the median QC area.

As a sensitivity analysis, the index was recalculated from the five members absent from the blanks, and all analyses were repeated (Table S7).

## 2.8. Statistical analysis

**Comparisons.** BC vs controls (primary), BC vs PC, PC vs controls and PC vs male controls.

**Main model.** For each comparison, a linear model of log2 intensity was fitted with group, sex, age and the exogenous-signal index as covariates (complete cases), with Benjamini–Hochberg FDR control.

**Robustness checks.**
1. 1:1 nearest-age matching within sex: caliper ± 5 years, optimal assignment (Hungarian algorithm), Mann–Whitney U test on the matched sets. Participants without age or sex were not matched.
2. Comparison with controls whose index was below the median.
3. Comparison of controls with high vs low index.

**Definition of a finding.** A compound was considered a finding when all of the following held:
- q < 0.05 in the main model;
- p < 0.05 in the age-matched comparison, in the same direction;
- p < 0.05 versus low-index controls, in the same direction;
- p ≥ 0.05 between high- and low-index controls.

Exogenous compounds could not qualify.

**Classification.**
- Validation scheme: 20 × 5-fold stratified cross-validation.
- Candidate pool, defined before cross-validation: exogenous compounds and compounds correlated with the index (|ρ| ≥ 0.5) were removed, leaving 5,351 compounds.
- Main model (panel): within each training fold, the five compounds with the smallest Mann–Whitney p were selected and a logistic regression was fitted.
- Index-adjusted panel: the same procedure after residualising all compounds on the exogenous-signal index within each training fold.
- PLS-DA (two components, Pareto scaling within folds) on all pool compounds, for comparison.
- Single-covariate baselines (sex only; index only) were evaluated with the same cross-validation.
- Significance of all models: 1,000 label permutations.
- Age-matched cross-validation was performed only with ≥ 10 pairs.
- Supplementary sensitivity analysis (Table S2): cases and references matched 1:1 within sex on the index (caliper 0.25 SD).

**Pathways.** Pathway enrichment used mummichog (positive mode, 5 ppm; top 10% of features ranked by the main-model p; exogenous features excluded; 500 permutations; Benjamini–Hochberg correction across pathways). It is reported in the Supplementary Material only.

## 2.9. Comparison with previously reported bladder cancer markers (optional)

The 14 metabolites reported as altered in BC in our previous study ([ABC], Table 5) were searched in the full, unfiltered feature table:
- m/z: ± 5 ppm of the reported m/z or of the theoretical m/z calculated from the formula and ion;
- retention time: ± 0.5 min of the reported value, to allow for retention-time shifts between analytical campaigns performed at different times with the same LC–HRMS method.

All candidate features within these windows were reported, with isomers flagged. For Leu-Leu-Tyr and Leu-Thr-Ser, the formula of the named peptide was used. Two replication criteria were defined before the analysis:
- **nominal:** same direction and Mann–Whitney p < 0.05;
- **replicated:** same direction and q < 0.05 in the adjusted model.

## 2.10. Annotation confidence

Annotation levels follow the Metabolomics Standards Initiative [6], consistent with [ABC]:
- **Level 2:** precursor within 5 ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the same feature. Isomers that the MS/MS spectrum does not distinguish (e.g., Leu-Val, Val-Ile and Ile-Val for Val-Leu; sn-1/sn-2 regioisomers of lysophosphatidylcholines) are listed with each compound.
- **Level 3:** a consistent formula without MS/MS confirmation, MS/MS matches to several isomers, or MS/MS matches not verified against a curated formula.
- **Level 4:** all other features, including annotations with a mass error > 5 ppm and library matches judged implausible in human urine.

No reference standards were analysed, so no level-1 annotations are reported. The same rule and the same source of MS/MS flags were used for Table 2, Figures 2–3 and Table S4.

## 2.11. Software and data availability

**Software.**
- MS-DIAL 5.5 and MS-FINDER 3.6.
- Python 3.13 with numpy 2.5, pandas 3.0, scipy 1.18, scikit-learn 1.9, statsmodels 0.15, matplotlib 3.11, seaborn 0.13, joblib 1.6, adjustText 1.4, pyteomics 4.7 and mummichog 2.7 (exact versions in requirements.txt).

**Reproducibility.** All random steps (cross-validation splits, permutations, bootstrap, mummichog permutations) were seeded from a single fixed seed (20261013). All parameters are set in one configuration file. A single command regenerates every table and figure from the MS-DIAL feature table and the sample metadata file.

**Data availability (draft).** The MS-DIAL feature table, sample metadata (age, sex and group; no personal identifiers), curation files and the processing code are available at [authors: repository and DOI]. Raw LC–HRMS files are available from [authors: MetaboLights/MassIVE accession, or the corresponding author on reasonable request].

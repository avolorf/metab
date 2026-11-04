# Response to Reviewer

We thank the Reviewer for these comments. Both points led to changes in the analysis and in the manuscript. All results affected by the changes were recalculated, and the revised values are given below and in the revised manuscript.

---

**Comment 1.** *The authors performed age- and sex-matching but used the Mann–Whitney U test for the matched groups. This test does not account for the paired design.*

**Response.** We agree. The matched-pair criterion is now evaluated with the Wilcoxon signed-rank test on within-pair differences (log2 intensity, case minus matched participant), and the direction of change is taken from the median paired difference. Pairs with a zero difference (both members below the detection limit and therefore assigned the same imputed value) are omitted, as in the standard Wilcoxon procedure. The exact null distribution is used when the non-zero differences contain no ties, and the normal approximation is used otherwise. The pairs themselves are unchanged: optimal 1:1 matching on age within sex (caliper ±5 years, without replacement), giving 20 BC–control, 18 BC–PC and 11 PC–control pairs. The two comparisons that are not paired (BC vs controls with a below-median exogenous-signal index, and controls with high vs low index) still use the Mann–Whitney U test.

The whole finding procedure and all sensitivity analyses were recalculated with the paired test. The adjusted linear model is not affected by this change (340 compounds with q < 0.05 in BC vs Control, as before). The effect on the findings is as follows:

- **BC vs Control:** 188 compounds met all finding criteria (104 higher and 84 lower in BC), compared with 210 previously; 176 of the previous 210 findings were retained, 34 no longer met the matched-pair criterion, and 12 were added.
- **Annotated findings (Table 2):** all five were retained. Their paired p-values are 0.006 (11-oxoetiocholanolone glucuronide), 0.021 and 0.033 (tetrahydrocortisone glucuronide isomers I and II), 0.033 (biliverdin) and 0.007 (LysoPC(16:0)). All 20 pairs had a non-zero difference for each of these compounds.
- **BC vs PC:** all three findings were retained (paired p = 0.00002–0.0009).
- **PC vs Control:** no findings, as before.

The cross-validated classification results do not depend on this test and are unchanged (BC vs Control: top-5 panel AUC 0.81 ± 0.03, PLS-DA 0.83 ± 0.02). The sensitivity analyses were updated accordingly. With MSTUS normalisation, 109 of the 122 MSTUS findings were also PQN findings. After excluding the most dilute samples, 176 of the 188 findings were kept. An index built only from series members absent from the blanks gave 165 findings, 157 of them shared, and retained all five annotated compounds. All 76 findings among the compounds detected in all samples were retained. When imputed values were treated as missing, 170 findings were obtained, 113 of them shared with the main analysis.

*Changes in the manuscript:* Section 2.5 (description of the matched-pair test); Section 3 (numbers of findings and sensitivity analyses); captions of Table 2 and Figure 2; Tables S1, S5 and S6. Table S5 additionally reports the number of pairs with a non-zero difference and, for reference, the previous Mann–Whitney p-value, labelled as superseded.

---

**Comment 2.** *The five compounds were selected within each training fold. It is therefore unclear whether the reported panel contains the same five compounds across all models.*

**Response.** The Reviewer is correct that the panel is not a fixed set of compounds. To avoid selection bias, the five compounds are re-selected by the Mann–Whitney test within each training fold, so the cross-validated AUC estimates the performance of the selection-and-classification procedure rather than of one fixed set of five compounds. We have now stated this explicitly in Section 2.6. To show how stable the selection is, we added Table S2b. For each comparison and for both panels (unadjusted and index-adjusted), it reports:

- the number of the 100 training folds (20 repeats × 5 folds) in which each compound was selected;
- the five compounds selected by the same rule on all samples, with their rank;
- whether each compound is also a univariate finding.

For BC vs Control, the selection has a stable core. The compounds at m/z 304.1182 and 301.1186 were selected in 88 and 72 of the 100 training folds, and the all-sample panel consists of these two and three further compounds (m/z 819.4420, 306.1738 and 438.2351; selected in 35–48 folds). All five compounds of the all-sample panel are univariate findings, as are 10 of the 11 compounds selected in at least 10% of the folds. The index-adjusted panel shares three of its five all-sample compounds with the unadjusted panel. For BC vs PC, two compounds were selected in 98 and 90 folds, and both belong to the three univariate BC–PC findings. For PC vs Control, the selection was unstable (no compound in more than 53% of folds), and none of the selected compounds is a univariate finding. This is consistent with our conclusion that no robust PC-specific signal was found. The compounds of the BC panels are not annotated (MSI level 4), so the panel is reported as evidence of a separable BC profile, not as a candidate diagnostic test.

*Changes in the manuscript:* Section 2.6 (new sentence on per-fold re-selection and reference to Table S2b); new Table S2b in the Supplementary Material.

#!/usr/bin/env python3
"""Анализ чувствительности: индекс загрязнения только из членов ряда, отсутствующих в бланках растворителя.

    python sensitivity_index.py        # после run_all.py и blank_check.py

Правило отбора (то же, что в blank_check.py): член ряда входит в «чистый» индекс, если медиана площади в бланках
растворителя < BLANK_FRAC медианы QC. Затем для всех сравнений повторяются основная модель (с «чистым» индексом
вместо исходного), правило находок, базовая модель «только индекс» и панель с поправкой на индекс (та же CV и
перестановки). Результат: results/tables/sensitivity_contam_index.csv.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from blank_check import BLANK_FRAC
from pipeline import config as C
from pipeline import features, multivariate
from pipeline import stats as st

warnings.filterwarnings("ignore")


def main():
    t = Path(C.OUT) / "tables"
    bsum = pd.read_csv(Path(C.OUT) / "blank_check" / "blank_check_summary.csv")
    meta = pd.read_csv(t / "samples.csv", index_col=0)
    meta.index = meta.index.astype(str)
    L = pd.read_csv(t / "processed_log2_matrix.csv", index_col=0)
    info = L[["rt_min", "mz"]]
    L = L[[c for c in L.columns if c in meta.index]].T
    L.index = L.index.astype(str)
    feat = pd.read_csv(t / "features.csv", index_col=0)
    series = pd.read_csv(t / "contaminant_series.csv", index_col=0)
    study = meta.index[meta.study]

    s = bsum[bsum.in_index == 1]
    clean_lab = s.label[s.solvent_blank_frac_of_QC < BLANK_FRAC].tolist()
    lab_of = {}
    for f, r in series[series.in_index == 1].iterrows():
        lab_of[f] = r.label if "label" in r else f
    clean = [f for f, l in lab_of.items() if l in clean_lab and f in L.columns]
    print(f"Члены ряда без фона в бланках (< {BLANK_FRAC:.0%} QC): {clean_lab} -> {len(clean)} признаков")
    z = L.loc[study, clean]
    ci = ((z - z.mean()) / z.std()).mean(axis=1)
    m2 = meta.copy()
    m2.loc[study, "contam_index"] = ci
    m2.loc[study, "contam_level"] = np.where(ci > ci.median(), "high", "low")
    rho = stats.spearmanr(ci, meta.loc[study, "contam_index"]).statistic

    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    comp_cat = feat.loc[reps, "compound_category"].fillna("")
    cat = comp_cat
    pool = list(reps[~features.panel_excluded(cat) & ~(feat.loc[reps, "rho_contam_index"].abs() >= C.PANEL_MAX_ABS_RHO_CONTAM)])
    rows = []
    for k, (name, spec) in enumerate(C.COMPARISONS.items()):
        r0 = pd.read_csv(t / f"stats_{name}.csv", index_col=0)
        r1, _, _ = st.run_comparison(L.loc[study, reps], m2, name, spec, comp_cat)
        f0, f1 = set(r0.index[r0.finding]), set(r1.index[r1.finding])
        q0, q1 = set(r0.index[r0.q_model < C.FDR]), set(r1.index[r1.q_model < C.FDR])
        row = dict(comparison=name, members_clean_index=len(clean), rho_with_main_index=rho,
                   q05_main=len(q0), q05_clean=len(q1), q05_overlap=len(q0 & q1),
                   findings_main=len(f0), findings_clean=len(f1), findings_overlap=len(f0 & f1),
                   rho_log2FC=r0.beta_model.corr(r1.beta_model.reindex(r0.index), method="spearman"))
        if spec["primary"] or name == "PC_vs_Control":
            case, ref = st.comparison_sets(meta, spec)
            y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
            zz = m2.loc[case + ref, "contam_index"].to_numpy(float)
            a, *_ = multivariate.repeated_cv(zz[:, None], y, "covariate", C.SEED + k)
            p, _ = multivariate.perm_test(zz[:, None], y, "covariate", a.mean(), C.SEED + k + 1, C.N_PERM)
            X = L.loc[case + ref, pool].to_numpy()
            b, *_ = multivariate.repeated_cv(X, y, "panel_adj", C.SEED + k, zz)
            pb, _ = multivariate.perm_test(X, y, "panel_adj", b.mean(), C.SEED + k + 1, C.N_PERM, zz)
            row.update(AUC_index_only=a.mean(), p_index_only=p, AUC_panel_adj=b.mean(), AUC_panel_adj_sd=b.std(), p_panel_adj=pb)
        rows.append(row)
        print(row)
    out = pd.DataFrame(rows)
    grp = meta.loc[study, "group"]
    out.attrs["means"] = {g: ci[grp == g].mean() for g in C.STUDY_GROUPS}
    out.to_csv(t / "sensitivity_contam_index.csv", index=False)
    print("Средний «чистый» индекс:", {g: round(ci[grp == g].mean(), 2) for g in C.STUDY_GROUPS},
          f"; ρ с PQN-фактором {stats.spearmanr(ci, meta.loc[study, 'dilution_PQN']).statistic:+.2f}")


if __name__ == "__main__":
    main()

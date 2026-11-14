#!/usr/bin/env python3
"""Анализ чувствительности к порогам объединения ионов в соединения (ΔRT и ρ Спирмена).

    python sensitivity_grouping.py                     # s266 (results/)
    METAB_DATASET=s107 python sensitivity_grouping.py  # s107 (results_s107/)

Работает по готовым результатам run_all.py (обработанная матрица, метаданные проб, категории признаков), исходная
таблица MS-DIAL не нужна. Для каждой пары порогов повторяются группировка, выбор представителя (курированный ион
целой молекулы, иначе самый интенсивный), категория соединения (наиболее «экзогенная» среди его ионов) и
одномерная статистика с правилом находок. Совпадение с основным результатом считается по признакам-представителям
находок основной обработки: находка «сохранилась», если её представитель входит в соединение, которое является
находкой при новых порогах.
Результат: <OUT>/tables/sensitivity_grouping.csv.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import config as C
from pipeline import features
from pipeline import stats as st

warnings.filterwarnings("ignore")
T = Path(C.OUT) / "tables"
PRIO = ["contaminant", "drug", "tobacco", "food", "possible drug (low score)", "possible tobacco (low score)",
        "possible food (low score)", "endogenous", "ambiguous", "annotated (high score)", "putative (low score)", "unknown"]


def main():
    L = pd.read_csv(T / "processed_log2_matrix.csv", index_col=0)
    L = L.drop(columns=["rt_min", "mz", "annotation", "annotation_conf"]).T.astype(float)
    L.index = L.index.astype(str)
    meta = pd.read_csv(T / "samples.csv", index_col=0)
    meta.index = meta.index.astype(str)
    meta["study"] = meta.study.astype(bool)
    feat = pd.read_csv(T / "features.csv", index_col=0).loc[L.columns]
    info = feat[["rt_min", "mz"]]
    study = meta.index[meta.study]
    cur = pd.read_csv(C.CURATION)
    main_res = {n: pd.read_csv(T / f"stats_{n}.csv", index_col=0) for n in C.COMPARISONS}
    rows = []
    for tol in [0.01, 0.02, 0.03, 0.05]:
        for rho in [0.7, 0.8, 0.9]:
            C.GROUP_RT_TOL, C.GROUP_MIN_CORR = tol, rho
            g = features.group_features(L.loc[study], info, feat.median_int)
            for _, r in cur[cur.ion.isin(features.MOLECULAR_IONS)].iterrows():
                f = next(iter(info[((info.mz - r.mz).abs() < 0.002) & ((info.rt_min - r.rt_min).abs() < 0.02)].index), None)
                if f is not None and f in g.index:
                    mem = g.index[g.compound == g.loc[f, "compound"]]
                    g.loc[mem, "representative"] = f
                    g.loc[mem, "is_representative"] = mem == f
            reps = g.index[g.is_representative]
            rank = feat.category.map({c: i for i, c in enumerate(PRIO)}).fillna(len(PRIO))
            best = pd.DataFrame({"compound": g.compound, "rank": rank, "cat": feat.category}) \
                .sort_values("rank").drop_duplicates("compound").set_index("compound")["cat"]
            comp_cat = g.compound.map(best).loc[reps]
            for name, spec in C.COMPARISONS.items():
                r, info_, _ = st.run_comparison(L.loc[study, reps], meta, name, spec, comp_cat)
                fnd = set(r.index[r.finding])
                fmem = set(g.index[g.compound.isin(g.loc[list(fnd), "compound"])]) if fnd else set()
                m0 = main_res[name]
                main_f = list(m0.index[m0.finding.astype(bool)])
                kept = sum(f in fmem for f in main_f)
                rows.append(dict(comparison=name, rt_tol_min=tol, min_rho=rho, n_compounds=len(reps),
                                 n_multi_ion=int((g.groupby("compound").size() > 1).sum()),
                                 n_q05=int((r.q_model < C.FDR).sum()), n_findings=len(fnd),
                                 main_findings=len(main_f), main_findings_retained=kept,
                                 retained_pct=round(100 * kept / max(1, len(main_f)), 1)))
                print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(T / "sensitivity_grouping.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

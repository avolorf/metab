"""Сверка с маркерами первого этапа и статьи ABC (апробация на независимой выборке).

Входной файл curation/prior_markers.csv (одна строка = один маркер из предыдущей работы):
    source       stage1 | ABC
    marker       название
    comparison   BC_vs_Control | PC_vs_Control | BC_vs_PC | PC_vs_maleControl
    direction    up | down (в case относительно ref, как в исходной работе)
    mz           наблюдавшийся m/z (обязательно)
    rt_min       RT, мин (если тот же метод; иначе пусто -> сопоставление только по массе)
    ion, formula, hmdb, comment  (необязательно)

Для каждого маркера: признаки в пределах PPM (и PRIOR_RT_TOL), соединение-представитель, его статистика
в текущем сравнении. Сводка: найдено / то же направление / то же направление и p<0.05 (все контроли и пары
по возрасту) и биномиальный тест согласия направлений (H0: 50 %).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .features import ppm


def run(path, info, groups, res, log):
    if not Path(path).exists():
        log(f"\nСверка с маркерами: файл {path} не найден — раздел пропущен (шаблон: curation/prior_markers_template.csv)")
        return None, None
    pm = pd.read_csv(path)
    rows = []
    for _, m in pm.iterrows():
        cand = info[ppm(info.mz, float(m.mz)).abs() <= C.PPM]
        rt_used = pd.notna(m.get("rt_min")) and str(m.get("rt_min")).strip() != ""
        if rt_used:
            cand = cand[(cand.rt_min - float(m.rt_min)).abs() <= C.PRIOR_RT_TOL]
        base = m.to_dict()
        base["match_mode"] = "m/z + RT" if rt_used else "только m/z"
        if not len(cand) or m.comparison not in res:
            rows.append({**base, "matched": False, "n_candidates": len(cand)})
            continue
        r = res[m.comparison]
        reps = groups.loc[cand.index, "representative"].unique()
        reps = [x for x in reps if x in r.index]
        if not reps:
            rows.append({**base, "matched": False, "n_candidates": len(cand)})
            continue
        best = r.loc[reps].sort_values("p_model").index[0]       # при нескольких кандидатах — наименьший p
        x = r.loc[best]
        exp = 1 if str(m.direction).lower().startswith("up") else -1
        rows.append({**base, "matched": True, "n_candidates": len(cand), "feature": best,
                     "mz_obs": info.loc[best, "mz"], "rt_obs": info.loc[best, "rt_min"],
                     "ppm": round(ppm(info.loc[best, "mz"], float(m.mz)), 2),
                     "log2FC_model": x.beta_model, "p_model": x.p_model, "q_model": x.q_model,
                     "log2FC_matched": x.get("log2FC_matched"), "p_matched": x.get("p_matched"),
                     "same_direction": bool(np.sign(x.beta_model) == exp),
                     "replicated_p05": bool(np.sign(x.beta_model) == exp and x.p_model < 0.05),
                     "replicated_matched_p05": bool(np.sign(x.get("log2FC_matched", np.nan)) == exp and x.get("p_matched", 1) < 0.05),
                     "category": x.compound_category, "finding_this_study": bool(x.finding)})
    tab = pd.DataFrame(rows)
    summ = []
    for (src, comp), d in tab.groupby(["source", "comparison"]):
        mm = d[d.matched == True]
        k, n = int(mm.same_direction.sum()), len(mm)
        summ.append(dict(source=src, comparison=comp, markers=len(d), matched=n, same_direction=k,
                         same_direction_p05=int(mm.replicated_p05.sum()),
                         same_direction_p05_age_matched=int(mm.replicated_matched_p05.sum()),
                         sign_test_p=stats.binomtest(k, n, 0.5, alternative="greater").pvalue if n else np.nan))
    summ = pd.DataFrame(summ)
    log("\nСверка с маркерами первого этапа / ABC:\n" + summ.to_string(index=False))
    return tab, summ

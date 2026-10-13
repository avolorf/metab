"""Сверка с маркерами первого этапа (статья ABC: Frolova et al., Anal Bioanal Chem 2026, Table 5).

Входной файл curation/prior_markers.csv (одна строка = один маркер):
    source, marker, comparison (BC_vs_Control ...), direction (up/down), mz, rt_min, ion, formula, hmdb, comment

Для каждого маркера ищутся ВСЕ признаки в ПОЛНОЙ таблице (до фильтров): |Δm/z| ≤ PPM и |ΔRT| ≤ PRIOR_RT_TOL
(изомеры в окне показываются все и помечаются как неоднозначные). Статус: not detected / below detection filter / failed RSD filter / quantified.
Для измеренных — статистика соединения-представителя в текущем сравнении и согласие направления.
Критерии (заданы до просмотра результатов):
    nominal    — то же направление и p (Манн–Уитни, все пробы) < 0.05;
    replicated — то же направление и q (модель с поправками на пол, возраст, индекс загрязнения) < FDR.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .features import ADDUCTS, formula_mass, ppm


def run(path, info_all, feat_table, res, log):
    if not Path(path).exists():
        log(f"\nСверка с маркерами: файл {path} не найден — раздел пропущен (шаблон: curation/prior_markers_template.csv)")
        return None, None
    pm = pd.read_csv(path)
    rows = []
    for _, m in pm.iterrows():
        d = m.to_dict()
        theo = formula_mass(m.formula) + ADDUCTS[m.ion] if isinstance(m.get("formula"), str) and m.get("ion") in ADDUCTS else np.nan
        d["theo_mz"] = round(theo, 4) if pd.notna(theo) else np.nan
        d["reported_ppm_vs_theo"] = round(ppm(m.mz, theo), 1) if pd.notna(theo) else np.nan
        cand = info_all[ppm(info_all.mz, float(m.mz)).abs() <= C.PPM].copy()
        if pd.notna(m.get("rt_min")):
            cand["drt"] = cand.rt_min - float(m.rt_min)
            cand = cand[cand.drt.abs() <= C.PRIOR_RT_TOL]
        if not len(cand):
            rows.append({**d, "status": "not detected", "candidate": 0, "n_candidates": 0})
            continue
        if "drt" in cand:
            cand = cand.reindex(cand.drt.abs().sort_values().index)
        # все кандидаты в окне (изомеры), а не только ближайший по RT
        for k, f in enumerate(cand.index, 1):
            row = dict(d, candidate=k, n_candidates=len(cand), feature=f, mz_obs=info_all.loc[f, "mz"],
                       rt_obs=info_all.loc[f, "rt_min"],
                       delta_rt=round(info_all.loc[f, "rt_min"] - float(m.rt_min), 2) if pd.notna(m.get("rt_min")) else np.nan,
                       ppm_obs=round(ppm(info_all.loc[f, "mz"], float(m.mz)), 1),
                       annotation_this_study=info_all.loc[f, "annotation"])
            if f not in feat_table.index or not feat_table.loc[f, "pass_detection"]:
                row["status"] = "below detection filter"
            elif not feat_table.loc[f, "kept"]:
                row["status"] = "failed QC RSD filter"
            else:
                row["status"] = "quantified"
                x = res[m.comparison].loc[feat_table.loc[f, "representative"]]
                exp = 1 if str(m.direction).lower().startswith("up") else -1
                row.update(log2FC=x.log2FC, AUC=x.AUC, p_MW=x.p_MW, beta_model=x.beta_model, q_model=x.q_model,
                           p_matched=x.get("p_matched"), p_lowindex=x.get("p_lowindex"), category=x.compound_category,
                           same_direction=bool(np.sign(x.log2FC) == exp),
                           nominal=bool(np.sign(x.log2FC) == exp and x.p_MW < 0.05),
                           replicated=bool(np.sign(x.beta_model) == exp and x.q_model < C.FDR),
                           finding_this_study=bool(x.finding))
            for g in C.STUDY_GROUPS:
                if f in feat_table.index:
                    row[f"det_{g}"] = feat_table.loc[f, f"det_{g}"]
            rows.append(row)
    tab = pd.DataFrame(rows)
    summ = []
    rank = {"quantified": 0, "failed QC RSD filter": 1, "below detection filter": 2, "not detected": 3}
    for (src, comp), t in tab.groupby(["source", "comparison"]):
        per = []
        for mk, u in t.groupby("marker", sort=False):
            q = u[u.status == "quantified"]
            per.append(dict(status=min(u.status, key=rank.get), n_q=len(q),
                            any_same=bool(q.same_direction.any()) if len(q) else False,
                            all_same=bool(q.same_direction.all()) if len(q) else False,
                            any_nominal=bool(q.nominal.any()) if len(q) else False,
                            any_rep=bool(q.replicated.any()) if len(q) else False))
        per = pd.DataFrame(per)
        qm = per[per.status == "quantified"]
        summ.append(dict(source=src, comparison=comp, markers=len(per),
                         not_detected=int((per.status == "not detected").sum()),
                         below_detection_filter=int((per.status == "below detection filter").sum()),
                         failed_rsd=int((per.status == "failed QC RSD filter").sum()), quantified=len(qm),
                         ambiguous_isomers=int((qm.n_q > 1).sum()),
                         same_direction_all_candidates=int(qm.all_same.sum()),
                         same_direction_any_candidate=int(qm.any_same.sum()),
                         nominal_p05_any_candidate=int(qm.any_nominal.sum()),
                         replicated_q05_any_candidate=int(qm.any_rep.sum())))
    summ = pd.DataFrame(summ)
    log(f"\nСверка с маркерами (m/z ± {C.PPM} ppm, RT ± {C.PRIOR_RT_TOL} мин):\n" + summ.to_string(index=False))
    return tab, summ

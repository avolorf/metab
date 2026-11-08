"""Сводка «было → стало» для перезапуска s107 (v1 -> v2) и проверка неизменности s266.

    python report_s107/compare_versions.py OLD_DIR
OLD_DIR — копия результатов до перезапуска: OLD_DIR/results_s107, OLD_DIR/results.
Выход: report_s107/v2_changes/Сводка_v1_v2.xlsx (листы) и печать основных чисел.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OLD = Path(sys.argv[1])
OUT = Path("report_s107/v2_changes")
OUT.mkdir(exist_ok=True)
EXO = ["contaminant", "drug", "tobacco", "food"]


def stats(d):
    return pd.read_csv(d / "tables/stats_BC_vs_Control.csv", index_col=0)


o, n = stats(OLD / "results_s107"), stats(Path("results_s107"))
fo, fn = o.finding.astype(bool), n.finding.astype(bool)
summary = []


def add(sect, item, a, b):
    summary.append(dict(section=sect, item=item, v1=a, v2=b))


add("находки", "всего", int(fo.sum()), int(fn.sum()))
add("находки", "↑ при РМП", int((fo & (o.beta_model > 0)).sum()), int((fn & (n.beta_model > 0)).sum()))
add("находки", "↓ при РМП", int((fo & (o.beta_model < 0)).sum()), int((fn & (n.beta_model < 0)).sum()))
add("находки", "по правилу «не обнаружено в контролях»", 0, int(n.get("finding_via_ctrl_not_detected", pd.Series(False)).sum()))
add("находки", "q < 0.05 (модель)", int((o.q_model < 0.05).sum()), int((n.q_model < 0.05).sum()))
cats = sorted(set(o.compound_category.dropna()) | set(n.compound_category.dropna()))
for c in cats:
    add("находки по категориям", c, int((fo & (o.compound_category == c)).sum()), int((fn & (n.compound_category == c)).sum()))
for c in cats:
    add("все соединения по категориям", c, int((o.compound_category == c).sum()), int((n.compound_category == c).sum()))

# какие соединения поменяли статус
cols = ["name", "compound_category", "beta_model", "q_model", "p_matched", "p_lowindex", "p_ctrl_hi_vs_lo", "AUC",
        "det_Control", "det_BC"]
new = n.loc[fn & ~fo.reindex(n.index, fill_value=False), [c for c in cols if c in n] +
            (["finding_via_ctrl_not_detected"] if "finding_via_ctrl_not_detected" in n else [])].copy()
new.insert(1, "category_v1", o.compound_category.reindex(new.index))
new["why_new"] = np.where(new.get("finding_via_ctrl_not_detected", False), "правило «не обнаружено в контролях»",
                          np.where(new.category_v1.isin(EXO), "категория: " + new.category_v1.astype(str) + " -> " +
                                   new.compound_category.astype(str), "другое"))
lost = o.loc[fo & ~fn.reindex(o.index, fill_value=False), [c for c in cols if c in o]]
cat_change = pd.DataFrame({"category_v1": o.compound_category, "category_v2": n.compound_category.reindex(o.index)})
cat_change = cat_change[cat_change.category_v1 != cat_change.category_v2].join(o[["name", "beta_model", "q_model"]])
cat_change["finding_v1"], cat_change["finding_v2"] = fo.reindex(cat_change.index), fn.reindex(cat_change.index)

# классификация
mo, mn = (pd.read_csv(d / "tables/multivariate_CV.csv") for d in (OLD / "results_s107", Path("results_s107")))
key = ["comparison", "set", "model"]
mv = mo.merge(mn, on=key, how="outer", suffixes=("_v1", "_v2"))
mv = mv[key + [c for c in ["n_features_v1", "n_features_v2", "AUC_mean_v1", "AUC_mean_v2", "Q2_v1", "Q2_v2",
                           "perm_p_v1", "perm_p_v2", "most_selected_v1", "most_selected_v2"] if c in mv]]
po, pn = (pd.read_csv(d / "tables/panel_selection.csv") for d in (OLD / "results_s107", Path("results_s107")))
ps = po.merge(pn, on=["comparison", "set", "method", "feature"], how="outer", suffixes=("_v1", "_v2"))
ps = ps[["comparison", "set", "method", "feature", "compound_v1", "compound_v2", "selection_frequency_v1",
         "selection_frequency_v2", "rank_all_sample_panel_v1", "rank_all_sample_panel_v2"]]
ps = ps.sort_values(["set", "method", "selection_frequency_v2"], ascending=[True, True, False])

# 30 маркеров отчёта
ro, rn = (pd.read_csv(d / "tables/replication_markers.csv") for d in (OLD / "results_s107", Path("results_s107")))
rcols = ["marker", "feature", "category", "beta_model", "q_model", "replicated", "finding_this_study"]
rep = ro[rcols].merge(rn[rcols], on=["marker", "feature"], how="outer", suffixes=("_v1", "_v2"))
rep["changed"] = (rep.category_v1.astype(str) != rep.category_v2.astype(str)) | \
                 (rep.finding_this_study_v1 != rep.finding_this_study_v2)
so, sn = (pd.read_csv(d / "tables/replication_summary.csv") for d in (OLD / "results_s107", Path("results_s107")))
rsum = pd.concat([so.assign(version="v1"), sn.assign(version="v2")])
best = lambda r: r.sort_values("q_model").drop_duplicates("marker")
for v, r in (("v1", ro), ("v2", rn)):
    pass
add("30 маркеров отчёта", "находки (≥1 признак маркера)", int(ro.groupby("marker").finding_this_study.any().sum()),
    int(rn.groupby("marker").finding_this_study.any().sum()))
add("30 маркеров отчёта", "экзогенная категория (≥1 признак)", int(ro.groupby("marker").category.apply(lambda c: c.isin(EXO).any()).sum()),
    int(rn.groupby("marker").category.apply(lambda c: c.isin(EXO).any()).sum()))

# s266: таблицы до/после
s266 = []
for f in sorted((OLD / "results/tables").glob("*.csv")):
    g = Path("results/tables") / f.name
    if not g.exists():
        s266.append(dict(table=f.name, status="нет в новой версии"))
        continue
    a, b = pd.read_csv(f), pd.read_csv(g)
    newcols = [c for c in b.columns if c not in a.columns]
    b2 = b[[c for c in b.columns if c in a.columns]]
    same_shape = a.shape == b2.shape
    if same_shape:
        num = a.select_dtypes("number").columns
        dmax = float((a[num] - b2[num]).abs().max().max()) if len(num) else 0.0
        txt = int((a.drop(columns=num).astype(str) != b2.drop(columns=num).astype(str)).to_numpy().sum())
    else:
        dmax, txt = np.nan, np.nan
    s266.append(dict(table=f.name, same_shape=same_shape, max_abs_numeric_diff=dmax, n_text_cells_diff=txt,
                     new_columns=", ".join(newcols)))
s266 = pd.DataFrame(s266)
for c in ["BC_vs_Control", "BC_vs_PC", "PC_vs_Control", "PC_vs_maleControl"]:
    fa, fb = OLD / f"results/tables/stats_{c}.csv", Path(f"results/tables/stats_{c}.csv")
    if fa.exists() and fb.exists():
        add("s266", f"находки {c}", int(pd.read_csv(fa).finding.sum()), int(pd.read_csv(fb).finding.sum()))

S = pd.DataFrame(summary)
print(S.to_string(index=False))
print(mv[["set", "model", "AUC_mean_v1", "AUC_mean_v2", "Q2_v1", "Q2_v2"]].to_string(index=False))
print(s266.to_string(index=False))
with pd.ExcelWriter(OUT / "Сводка_v1_v2.xlsx") as w:
    S.to_excel(w, sheet_name="summary", index=False)
    new.to_excel(w, sheet_name="new_findings")
    lost.to_excel(w, sheet_name="lost_findings")
    cat_change.to_excel(w, sheet_name="category_changes")
    mv.to_excel(w, sheet_name="classification", index=False)
    ps.to_excel(w, sheet_name="panel_selection", index=False)
    rep.to_excel(w, sheet_name="report_30_markers", index=False)
    rsum.to_excel(w, sheet_name="report_30_summary", index=False)
    s266.to_excel(w, sheet_name="s266_tables_check", index=False)

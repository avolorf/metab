#!/usr/bin/env python3
"""Полная обработка данных ВЭЖХ-МС одним запуском.

    python run_all.py                       # параметры из pipeline/config.py, метаданные из metadata/metadata.csv
    python run_all.py --metadata metadata/metadata_v2.csv
    python run_all.py --quick               # быстрый прогон для проверки (мало перестановок, без путей)

Каталог результатов (results/) пересоздаётся при каждом запуске.
"""
import argparse
import random
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_curve

from pipeline import config as C
from pipeline import annotate, features, io, multivariate, plots, preprocess
from pipeline import stats as st

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=C.DATA)
    ap.add_argument("--metadata", default=C.METADATA)
    ap.add_argument("--out", default=C.OUT)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    if a.quick:
        C.N_PERM, C.CV_REPEATS, C.N_BOOT, C.RUN_PATHWAYS = 20, 3, 200, False
    np.random.seed(C.SEED)
    random.seed(C.SEED)

    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    for sub in ["tables", "figures", "external", "pathways"]:
        (out / sub).mkdir(parents=True)
    fig = out / "figures"
    log_lines, sheets = [], {}

    def log(s=""):
        print(s)
        log_lines.append(str(s))

    # ------------------------------------------------------------ 1. данные и метаданные
    X_all, kind, info_all, ndup = io.load_features(a.data)
    meta = io.load_metadata(a.metadata, kind)
    meta["inj_order"] = io.injection_order(meta)
    study = meta.index[meta.study]
    qc = meta.index[meta.kind == "QC"]
    log(f"Таблица: {len(X_all)} инъекций ({(kind == 'QC').sum()} QC, {(kind == 'blank').sum()} бланков), "
        f"{X_all.shape[1]} признаков; дублированных строк удалено: {ndup}")
    log(f"Метаданные: {a.metadata}; включено {len(study)}: " +
        ", ".join(f"{g} {(meta.loc[study, 'group'] == g).sum()}" for g in C.STUDY_GROUPS))
    exc = meta[(meta.kind == "sample") & (meta.include == 0)]
    log("Исключены: " + "; ".join(f"{i} ({r.get('exclude_reason', '')})" for i, r in exc.iterrows()))
    demo, demo_tests = st.demographics(meta)
    sheets["demographics"] = demo
    sheets["demographics_tests"] = demo_tests
    log("\nВозраст и пол по группам:\n" + demo.round(1).to_string(index=False))
    log(demo_tests.assign(p=demo_tests.p.map(lambda x: f"{x:.3g}")).to_string(index=False))

    # ------------------------------------------------------------ 2. предобработка
    use = meta.index[meta.study | (meta.kind.isin(["QC", "blank"]))]
    X1, feat_table = preprocess.filter_features(X_all.loc[use], meta.loc[use], log)
    loo = preprocess.drift_loo(X1, meta.inj_order, qc)
    log(f"Дрейф: leave-one-QC-out median |Δlog2| без коррекции {loo['none']:.3f}, LOWESS {loo['lowess']:.3f}. "
        f"Коррекция: {C.DRIFT}")
    if C.DRIFT == "lowess":
        X1 = preprocess.drift_correct(X1, meta.inj_order, qc)
    ana = meta.index[meta.study | (meta.kind == "QC")]
    norm = {}
    for m in C.COMPARE_NORMALIZATIONS:
        Xn, f, txt = preprocess.normalize(X1.loc[ana], meta.loc[ana], m)
        norm[m] = (preprocess.impute_log(Xn, ana), f)
        meta[f"dilution_{m}"] = f
        log(f"Нормализация {txt}: факторы образцов {f[study].min():.2f}–{f[study].max():.2f}")
    L = norm[C.NORMALIZATION][0]
    log(f"Основная нормализация: {C.NORMALIZATION}. Импутировано {X1.loc[ana].isna().to_numpy().mean():.1%} значений "
        "(½ минимума признака)")
    low = meta.loc[study][meta.loc[study, "dilution_PQN"] < C.LOW_DILUTION_FACTOR]
    nfeat = X1.loc[study].notna().sum(axis=1)
    low_tab = pd.DataFrame({"group": low.group, "PQN_factor": low.dilution_PQN, "MSTUS_factor": low.dilution_MSTUS,
                            "n_detected_features": nfeat[low.index], "median_n_detected": nfeat.median(),
                            "TIC_raw": X_all.loc[low.index].sum(axis=1)})
    sheets["low_dilution_samples"] = low_tab
    log(f"Пробы с PQN-фактором < {C.LOW_DILUTION_FACTOR}: {', '.join(low.index)} — требуют просмотра хроматограмм; "
        "ниже — анализ чувствительности без них")

    # детекция / импутация по группам, «есть/нет»
    lo, hi = C.PRESENCE_ABSENCE
    dcols = [f"det_{g}" for g in C.STUDY_GROUPS]
    for g in C.STUDY_GROUPS:
        feat_table[f"imputed_{g}"] = 1 - feat_table[f"det_{g}"]
    feat_table["presence_absence"] = (feat_table[dcols].min(axis=1) < lo) & (feat_table[dcols].max(axis=1) >= hi)

    # ------------------------------------------------------------ 3. соединения, индекс загрязнения, категории
    info = info_all.loc[L.columns]
    groups = features.group_features(L.loc[study], info, X1.loc[study].median())
    reps = groups.index[groups.is_representative]
    log(f"\nГруппировка признаков (ΔRT ≤ {C.GROUP_RT_TOL} мин, ρ ≥ {C.GROUP_MIN_CORR}): {L.shape[1]} признаков → "
        f"{len(reps)} соединений; многоионных соединений {(groups.groupby('compound').size() > 1).sum()}")
    series = features.match_series(info)
    ci, lvl, thr, idx_feats = features.contamination_index(L, series, study)
    meta["contam_index"], meta["contam_level"] = ci, lvl
    log(f"Индекс загрязнения: {len(idx_feats)} признаков ряда диэтаноламинов; граница high/low = {thr:.2f} (медиана)")
    log("  средний индекс: " + ", ".join(f"{g} {meta.contam_index[study][meta.group[study] == g].mean():+.2f}" for g in C.STUDY_GROUPS))
    cur_raw = pd.read_csv(C.CURATION)
    cur_raw["feature"] = [next(iter(info[((info.mz - r.mz).abs() < 0.002) & ((info.rt_min - r.rt_min).abs() < 0.02)].index), None)
                          for _, r in cur_raw.iterrows()]
    kw = features.exposure_category(info)
    cats = features.categorize(info, groups, series, cur_raw, kw)
    names = features.display_name(info, cats.category)
    comp_cat = cats.loc[reps, "compound_category"]
    feat_table = feat_table.join(info).join(groups.drop(columns=["rt_min", "mz"]), how="left").join(cats[["category", "compound_category"]])
    feat_table["display_name"] = names
    log("Категории соединений: " + ", ".join(f"{k} {v}" for k, v in comp_cat.value_counts().items()))

    # ------------------------------------------------------------ 4. одномерная статистика
    Lc = L.loc[study, reps]
    res, infos, pairs_all = {}, [], []
    for name, spec in C.COMPARISONS.items():
        r, inf, pairs = st.run_comparison(Lc, meta, name, spec, comp_cat)
        r.insert(0, "name", names[r.index])
        r.insert(1, "n_ions", groups.loc[r.index, "n_ions"])
        res[name] = r
        infos.append(inf)
        pairs_all.append(pairs)
    summ = pd.DataFrame(infos)
    sheets["summary_univariate"] = summ
    sheets["age_matched_pairs"] = pd.concat(pairs_all)
    log("\nОдномерная статистика (уровень соединений):\n" + summ.to_string(index=False))

    # признаки «есть/нет»: точный тест Фишера по детекции
    pa = feat_table[feat_table.kept & feat_table.presence_absence].copy()
    for name, spec in C.COMPARISONS.items():
        case, ref = st.comparison_sets(meta, spec)
        dc = X1.loc[case, pa.index].notna().sum()
        dr = X1.loc[ref, pa.index].notna().sum()
        pa[f"p_fisher_{name}"] = [stats.fisher_exact([[a_, len(case) - a_], [b_, len(ref) - b_]]).pvalue for a_, b_ in zip(dc, dr)]
    sheets["presence_absence"] = pa
    log(f"Признаки «есть/нет» (детекция <{lo:.0%} в одной группе и ≥{hi:.0%} в другой): {len(pa)}")

    # чувствительность к нормализации
    nc_rows, res_norm = [], {}
    for m, (Lm, _) in norm.items():
        res_norm[m] = {}
        for name, spec in C.COMPARISONS.items():
            r, _, _ = st.run_comparison(Lm.loc[study, reps], meta, name, spec, comp_cat)
            res_norm[m][name] = r
    a_, b_ = C.COMPARE_NORMALIZATIONS[:2]
    for name in C.COMPARISONS:
        s1 = set(res_norm[a_][name].index[res_norm[a_][name].q_model < C.FDR])
        s2 = set(res_norm[b_][name].index[res_norm[b_][name].q_model < C.FDR])
        f1 = set(res_norm[a_][name].index[res_norm[a_][name].finding])
        f2 = set(res_norm[b_][name].index[res_norm[b_][name].finding])
        nc_rows.append(dict(comparison=name, **{f"q<FDR_{a_}": len(s1), f"q<FDR_{b_}": len(s2)}, overlap=len(s1 & s2),
                            jaccard=len(s1 & s2) / max(1, len(s1 | s2)),
                            **{f"findings_{a_}": len(f1), f"findings_{b_}": len(f2)}, findings_overlap=len(f1 & f2),
                            rho_log2FC=res_norm[a_][name].beta_model.corr(res_norm[b_][name].beta_model, method="spearman"),
                            log2FC_shift=float((res_norm[b_][name].beta_model - res_norm[a_][name].beta_model).median())))
    # чувствительность: без проб с низким фактором разведения
    meta_s = meta.copy()
    meta_s.loc[low.index, "study"] = False
    for row in nc_rows:
        name = row["comparison"]
        r, _, _ = st.run_comparison(L.loc[meta_s.index[meta_s.study], reps], meta_s, name, C.COMPARISONS[name], comp_cat)
        s0 = set(res[name].index[res[name].q_model < C.FDR])
        s3 = set(r.index[r.q_model < C.FDR])
        row.update(q_without_lowdilution=len(s3), overlap_with_main=len(s0 & s3),
                   findings_without_lowdilution=int(r.finding.sum()),
                   findings_overlap_main=len(set(res[name].index[res[name].finding]) & set(r.index[r.finding])))
    ncs = pd.DataFrame(nc_rows)
    log("\nЧувствительность (PQN vs MSTUS; без проб с низким фактором разведения):\n" + ncs.round(2).to_string(index=False))

    # ------------------------------------------------------------ 5. многомерный анализ
    rho_ci = Lc.corrwith(meta.loc[study, "contam_index"], method="spearman")
    feat_table["rho_contam_index"] = rho_ci.reindex(feat_table.index)
    ex_cat = features.panel_excluded(comp_cat)
    ex_rho = rho_ci.abs() >= C.PANEL_MAX_ABS_RHO_CONTAM
    pool = list(Lc.columns[~ex_cat & ~ex_rho])
    log(f"\nПул для моделей: {len(pool)} соединений (исключены до CV: экспозиции {int(ex_cat.sum())}, "
        f"|ρ с индексом загрязнения| ≥ {C.PANEL_MAX_ABS_RHO_CONTAM} ещё {int((ex_rho & ~ex_cat).sum())}: "
        + ", ".join(names[f] for f in Lc.columns[ex_rho & ~ex_cat]) + ")")
    mv_rows, curves = [], {}
    for k, (name, spec) in enumerate(C.COMPARISONS.items()):
        case, ref = st.comparison_sets(meta, spec)
        pr = pairs_all[k]
        for label, (cs, rf) in {"all": (case, ref), "matched": (list(pr.case), list(pr.ref))}.items():
            need = C.MIN_PAIRS_CV if label == "matched" else C.CV_FOLDS
            if min(len(cs), len(rf)) < need:
                mv_rows.append(dict(comparison=name, set=label, n_case=len(cs), n_ref=len(rf),
                                    model=f"не оценивается: {min(len(cs), len(rf))} < {need} на группу"))
                continue
            rows, cv = multivariate.evaluate(Lc, meta, cs, rf, pool, label, C.SEED + k, names)
            for r_ in rows:
                r_["comparison"] = name
            mv_rows += rows
            curves.setdefault(name, {}).update(cv)
    mv = pd.DataFrame(mv_rows)
    cols = ["comparison", "set", "model", "n_case", "n_ref", "n_features", "AUC_mean", "AUC_sd", "AUC_p5", "AUC_p95",
            "Q2", "perm_p", "n_perm", "most_selected"]
    mv = mv[[c for c in cols if c in mv]]
    sheets["multivariate_CV"] = mv
    log("\nМногомерные модели (повторная CV, отбор признаков внутри фолдов):\n" +
        mv.drop(columns=["most_selected"], errors="ignore").round(3).to_string(index=False))
    for m in C.COMPARE_NORMALIZATIONS:
        for name in ["BC_vs_Control", "BC_vs_PC"]:
            case, ref = st.comparison_sets(meta, C.COMPARISONS[name])
            X = norm[m][0].loc[case + ref, pool].to_numpy()
            y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
            aucs, *_ = multivariate.repeated_cv(X, y, "panel", C.SEED)
            for row in nc_rows:
                if row["comparison"] == name:
                    row[f"panel_AUC_{m}"] = aucs.mean()
    sheets["normalization_sensitivity"] = pd.DataFrame(nc_rows)

    # ------------------------------------------------------------ 6. кандидаты и находки
    Lraw_all = np.log2(X_all.loc[study])
    cur_ok, cur_bad = annotate.curate(C.CURATION, info_all, info, Lraw_all, L, meta, groups, feat_table)
    for name in C.COMPARISONS:
        r = res[name]
        cmap = groups.compound
        for i, row in cur_ok.iterrows():
            f = row.feature
            if isinstance(f, str) and f in cmap.index and cmap[f] in set(groups.loc[r.index, "compound"]):
                rep = groups.index[(groups.compound == cmap[f]) & groups.is_representative][0]
                for c in ["beta_model", "q_model", "p_matched", "p_lowindex", "p_ctrl_hi_vs_lo", "finding"]:
                    if c in r:
                        cur_ok.loc[i, f"{name}:{c}"] = r.loc[rep, c]
    sheets["candidates_curated"] = cur_ok
    sheets["candidates_mass_mismatch"] = cur_bad
    log(f"\nКурированные кандидаты: {len(cur_ok)} с массой ≤{C.PPM} ppm; {len(cur_bad)} вынесены на лист "
        "candidates_mass_mismatch")
    find = []
    for name, r in res.items():
        f = r[r.finding].copy()
        f.insert(0, "comparison", name)
        find.append(f)
        sheets[f"stats_{name}"] = r.join(groups[["compound", "ion_relations"]]).join(
            feat_table[["rt_min", "mz", "annotation", "annotation_conf"] + [f"det_{g}" for g in C.STUDY_GROUPS] +
                       [f"imputed_{g}" for g in C.STUDY_GROUPS]]).sort_values("p_model")
    findings = pd.concat(find)
    sheets["findings"] = findings
    log("Находки (правило в README): " + ", ".join(f"{n} {int(r.finding.sum())}" for n, r in res.items()))

    # ------------------------------------------------------------ 6б. сверка с маркерами первого этапа / статьи ABC
    from pipeline import replication
    rep_tab, rep_sum = replication.run(C.PRIOR_MARKERS, info_all, feat_table, res, log)
    if rep_tab is not None:
        sheets["replication_markers"] = rep_tab
        sheets["replication_summary"] = rep_sum

    # ------------------------------------------------------------ 7. пути (приложение)
    if C.RUN_PATHWAYS:
        from pipeline import pathways
        excl = list(feat_table.index[feat_table.kept & feat_table.compound_category.isin(features.EXPOSURE)])
        pws, sets = [], []
        for k, name in enumerate(["BC_vs_Control", "BC_vs_PC", "PC_vs_Control"]):
            rf, _, _ = st.run_comparison(L.loc[study], meta, name, C.COMPARISONS[name], cats.compound_category)
            pw, s_ = pathways.run(rf, info, name, out / "pathways" / name, excl, C.SEED + k)
            pws.append(pw)
            sets.append(s_)
            plots.pathways(pw, fig / f"A1_pathways_{name}", f"mummichog: {name}")
        sheets["pathways_appendix"] = pd.concat(pws)
        sheets["pathways_settings"] = pd.DataFrame(sets)
        log("\nПути (mummichog, без экспозиций, BH по путям):\n" + pd.DataFrame(sets).to_string(index=False))

    # ------------------------------------------------------------ 8. рисунки
    plots.qc_overview(X_all, meta, meta.inj_order, feat_table.RSD_QC_raw[feat_table.pass_detection],
                      pd.DataFrame({"PQN": meta.dilution_PQN, "MSTUS": meta.dilution_MSTUS}), fig / "01_qc_overview")
    S, var = multivariate.pca_scores(L.loc[:, reps])
    plots.pca_panels(S, var, meta, meta.inj_order, fig / "02_pca")
    plots.contamination(meta, series, fig / "03_contamination_index")
    for k, (name, spec) in enumerate(C.COMPARISONS.items()):
        r = res[name]
        plots.volcano(r, name, names, fig / f"04_volcano_{name}")
        plots.heatmap(Lc, meta, r, names, comp_cat, spec["case"], spec["ref"], fig / f"05_heatmap_{name}", name)
        fnd = list(r[r.finding].sort_values("p_model").index[:16])
        plots.boxplots(L.loc[study], meta, fnd, names, fig / f"06_findings_{name}", f"{name}: находки")
        singles = []
        case, ref = st.comparison_sets(meta, spec)
        y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
        endo = r[r.finding & r.compound_category.isin(["endogenous", "annotated (high score)"])].sort_values("p_model")
        for f in endo.index[:6]:
            s = Lc.loc[case + ref, f].to_numpy()
            auc, lo_, hi_ = multivariate.boot_auc(y, s, C.SEED)
            fpr, tpr, _ = roc_curve(y, s if np.median(s[y == 1]) >= np.median(s[y == 0]) else -s)
            singles.append((names[f], fpr, tpr, auc, lo_, hi_))
        if name in curves:
            plots.roc(curves[name], singles, name, fig / f"07_roc_{name}")
    plots.norm_concordance(res_norm, fig / "08_normalization_PQN_vs_MSTUS")

    # ------------------------------------------------------------ 9. экспорт
    export_external(L, Lc, meta, res, feat_table, cur_ok, out / "external", cats)
    sheets["samples"] = meta.drop(columns=["exclude_reason", "comment"], errors="ignore")
    sheets["features"] = feat_table
    sheets["contaminant_series"] = series
    sheets["log"] = pd.DataFrame({"log": log_lines})
    order = ["log", "demographics", "demographics_tests", "summary_univariate", "findings", "multivariate_CV",
             "replication_summary", "replication_markers",
             "normalization_sensitivity", "candidates_curated", "candidates_mass_mismatch", "low_dilution_samples",
             "presence_absence", "age_matched_pairs", "samples", "contaminant_series"] + \
            [f"stats_{n}" for n in C.COMPARISONS] + ["features", "pathways_appendix", "pathways_settings"]
    with pd.ExcelWriter(out / "results.xlsx") as w:
        for k in order:
            if k in sheets:
                sheets[k].to_excel(w, sheet_name=k[:31], index=not isinstance(sheets[k].index, pd.RangeIndex))
                sheets[k].to_csv(out / "tables" / f"{k}.csv", index=not isinstance(sheets[k].index, pd.RangeIndex))
    L.T.join(info).to_csv(out / "tables" / "processed_log2_matrix.csv")
    (out / "log.txt").write_text("\n".join(log_lines), encoding="utf-8")
    import subprocess
    import sys
    subprocess.run([sys.executable, "make_figures.py", "--res", str(out), "--data", a.data], check=True)
    print(f"\nГотово → {out}")


def export_external(L, Lc, meta, res, feat_table, cur, out, cats):
    study = meta.index[meta.study]
    nm = feat_table.mz.round(4).astype(str) + "__" + feat_table.rt_min.round(3).astype(str)
    keep = Lc.columns[cats.loc[Lc.columns, "compound_category"] != "contaminant"]
    iqr = Lc.loc[study, keep].quantile(0.75) - Lc.loc[study, keep].quantile(0.25)
    keep = iqr.sort_values(ascending=False).index[:5000]   # лимит MetaboAnalyst; отбор без учёта групп
    X = (2 ** Lc[keep]).rename(columns=nm)
    for n, groups in {"BC_vs_Control": ["BC", "Control"], "BC_vs_PC": ["BC", "PC"], "PC_vs_Control": ["PC", "Control"],
                      "BC_PC_Control": ["BC", "PC", "Control"]}.items():
        s = study[meta.loc[study, "group"].isin(groups)]
        t = X.loc[s].copy()
        t.insert(0, "Label", meta.loc[s, "group"])
        t.index.name = "Sample"
        t.to_csv(out / f"1_stats_{n}.csv", float_format="%.5g")
    s = study[meta.loc[study, "sex"].isin(["M", "F"])]
    X.loc[s].rename_axis("Sample").to_csv(out / "2_metadata_module_data.csv", float_format="%.5g")
    md = pd.DataFrame({"Class": meta.loc[s, "group"], "Sex": meta.loc[s, "sex"], "Age": meta.loc[s, "age"],
                       "ContamIndex": meta.loc[s, "contam_index"].round(3), "ContamLevel": meta.loc[s, "contam_level"]})
    md.rename_axis("Sample").to_csv(out / "2_metadata_module_metadata.csv")
    pd.DataFrame({"compound_column": nm[keep], "category": cats.loc[keep, "compound_category"]}).to_csv(
        out / "0_compound_categories.csv")
    expo = features.panel_excluded(feat_table.compound_category.fillna(""))
    for name in ["BC_vs_Control", "BC_vs_PC", "PC_vs_Control"]:
        r, _, _ = st.run_comparison(L.loc[study, feat_table.index[feat_table.kept & ~expo]], meta, name,
                                    C.COMPARISONS[name], feat_table.compound_category)
        pd.DataFrame({"m.z": feat_table.mz[r.index], "p.value": r.p_model,
                      "t.score": np.sign(r.beta_model) * -np.log10(r.p_model),
                      "r.t": feat_table.rt_min[r.index] * 60}).dropna().to_csv(out / f"3_peaks_{name}.txt", sep="\t", index=False)
    endo = cur[cur.category.str.startswith("эндоген")]
    for focus in ["BC", "PC"]:
        g = endo[endo.focus == focus]
        (out / f"4_compounds_{focus}_HMDB.txt").write_text("\n".join(g.hmdb[g.hmdb != ""].drop_duplicates()))
        (out / f"4_compounds_{focus}_names.txt").write_text("\n".join(g.name_en[g.name_en != ""].drop_duplicates()))
        k = g[g.kegg != ""]
        col = f"{focus}_vs_Control:beta_model"
        if col in k:
            (out / f"6_kegg_mapper_{focus}.txt").write_text(
                "\n".join(f"{a} {'red' if b > 0 else 'blue'}" for a, b in zip(k.kegg, k[col])))
    sig = pd.concat([r[r.q_model < C.FDR] for r in res.values()])
    sig = feat_table.loc[sig.index.unique()]
    sig = sig[~sig.compound_category.isin(["contaminant"])]
    (out / "7_metlin_hmdb_mz.txt").write_text("\n".join(f"{m:.5f}" for m in sig.mz.sort_values()))
    sig[["mz", "rt_min", "annotation", "compound_category"]].sort_values("mz").to_csv(out / "7_metlin_hmdb_mz_rt.csv")


if __name__ == "__main__":
    main()

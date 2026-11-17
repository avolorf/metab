"""Матрица соединений s107 — та же, что в основном анализе (run_all.py, METAB_DATASET=s107).

Шаги повторяют run_all.py дословно: загрузка, фильтр детекции, межсерийная коррекция по QC, фильтр RSD,
PQN, импутация ½ минимума, log2, группировка ионов в соединения (представитель — курированный ион молекулы),
индекс экзогенного сигнала, категории. Код конвейера берётся из ветки основного анализа
(claude/metabolomics-hplc-ms-processing-ktvooo, коммит 97dffcb): путь к её копии — переменная S107_PIPELINE.
Проверка: числа на каждом шаге и результат основного сравнения РМП vs контроль совпадают с results_s107/log.txt.
"""
import os
import sys
from pathlib import Path

PIPE = Path(os.environ.get("S107_PIPELINE", "")).resolve()
if not (PIPE / "pipeline" / "config.py").exists():
    sys.exit("Укажите S107_PIPELINE — папку с копией ветки основного анализа s107 (там должны быть pipeline/ и data/s107/)")
os.environ["METAB_DATASET"] = "s107"
os.chdir(PIPE)
sys.path.insert(0, str(PIPE))

import numpy as np                                   # noqa: E402
import pandas as pd                                  # noqa: E402

from pipeline import config as C, features, io, preprocess, stats as st   # noqa: E402

# числа основного анализа (results_s107/log.txt) — должны совпасть
EXPECTED = {"injections": 128, "features": 12930, "after_detection": 5887, "after_rsd": 4730, "compounds": 3603,
            "study": 104, "q_model": 918, "findings": 582}


def build(log=print):
    np.random.seed(C.SEED)
    X_all, kind, info_all, _ = io.load_features(C.DATA)
    meta = io.load_metadata(C.METADATA, kind)
    meta["inj_order"] = io.injection_order(meta)
    meta["batch"] = io.batch_of(meta)
    study = meta.index[meta.study]
    use = meta.index[meta.study | (meta.kind.isin(["QC", "blank"]))]
    steps = []
    X1, feat_table = preprocess.filter_features(X_all.loc[use], meta.loc[use], steps.append)
    ana = meta.index[meta.study | (meta.kind == "QC")]
    Xn, f, _ = preprocess.normalize(X1.loc[ana], meta.loc[ana], C.NORMALIZATION)
    L = preprocess.impute_log(Xn, ana)
    info = info_all.loc[L.columns]
    groups = features.group_features(L.loc[study], info, X1.loc[study].median())
    cur_rep = pd.read_csv(C.CURATION)
    for _, r in cur_rep[cur_rep.ion.isin(features.MOLECULAR_IONS)].iterrows():
        ft = next(iter(info[((info.mz - r.mz).abs() < 0.002) & ((info.rt_min - r.rt_min).abs() < 0.02)].index), None)
        if ft is not None and ft in groups.index:
            members = groups.index[groups.compound == groups.loc[ft, "compound"]]
            groups.loc[members, "representative"] = ft
            groups.loc[members, "is_representative"] = members == ft
    reps = groups.index[groups.is_representative]
    series = features.match_series(info)
    ci, lvl, thr, _ = features.contamination_index(
        L, series, study, study[meta.loc[study, "group"] == "Control"] if C.CONTAM_SPLIT == "control_median" else None)
    meta["contam_index"], meta["contam_level"] = ci, lvl
    cur_raw = pd.read_csv(C.CURATION)
    cur_raw["feature"] = [next(iter(info[((info.mz - r.mz).abs() < 0.002) & ((info.rt_min - r.rt_min).abs() < 0.02)].index), None)
                          for _, r in cur_raw.iterrows()]
    cats = features.categorize(info, groups, series, cur_raw, features.exposure_category(info))
    names = features.display_name(info, cats.category)
    comp_cat = cats.loc[reps, "compound_category"]
    # проверка: основное сравнение РМП vs контроль
    spec = C.COMPARISONS["BC_vs_Control"]
    main, inf, _ = st.run_comparison(L.loc[study, reps], meta, "BC_vs_Control", spec, comp_cat)
    got = {"injections": len(X_all), "features": X_all.shape[1], "after_detection": int(feat_table.pass_detection.sum()),
           "after_rsd": L.shape[1], "compounds": len(reps), "study": len(study),
           "q_model": inf["n_q_model"], "findings": inf["n_findings"]}
    check = pd.DataFrame({"expected": EXPECTED, "reproduced": got})
    check["match"] = check.expected == check.reproduced
    for s in steps:
        log("  " + s)
    log(check.to_string())
    if not check.match.all():
        raise SystemExit("Матрица не совпадает с основным анализом s107 — расчёт остановлен")
    return dict(L=L, Lrep=L.loc[study, reps], X1=X1, meta=meta, info=info, groups=groups, reps=reps, names=names,
                cats=cats, comp_cat=comp_cat, main=main, cur=cur_raw, series=series, check=check,
                pipeline_commit=os.popen("git -C . rev-parse --short HEAD").read().strip())


if __name__ == "__main__":
    build()

"""Эталонный расчёт тем же алгоритмом на Python (модули pipeline/ из корня репозитория) для сверки с UroMetabR.
    python3 tests/crosscheck_python.py <каталог контрольного набора> <каталог вывода>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline import config as C, features, io, preprocess  # noqa: E402
from pipeline import stats as st                               # noqa: E402

d, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
X_all, kind, info_all, _ = io.load_features(d / "feature_table.csv")
meta = io.load_metadata(d / "metadata.csv", kind)
study = meta.index[meta.study]
use = meta.index[meta.study | meta.kind.isin(["QC", "blank"])]
X1, ft = preprocess.filter_features(X_all.loc[use], meta.loc[use], lambda s: None)
ana = meta.index[meta.study | (meta.kind == "QC")]
Xn, f, _ = preprocess.normalize(X1.loc[ana], meta.loc[ana], "PQN")
L = preprocess.impute_log(Xn, ana)
_, fm, _ = preprocess.normalize(X1.loc[ana], meta.loc[ana], "MSTUS")
pd.DataFrame({"pqn": f, "mstus": fm}).to_csv(out / "py_factors.csv")
ft[["RSD_QC_raw", "kept"]].to_csv(out / "py_filter.csv")
g = features.group_features(L.loc[study], info_all.loc[L.columns], X1.loc[study].median())
g[["compound", "is_representative"]].to_csv(out / "py_groups.csv")
C.CONTAMINANT_SERIES = str(Path(__file__).resolve().parents[1] / "example" / "contaminant_series.csv")
ser = features.match_series(info_all.loc[L.columns], C.CONTAMINANT_SERIES)
ci, lvl, _, _ = features.contamination_index(L, ser, study)
meta["contam_index"], meta["contam_level"] = ci, lvl
meta["contam_index"].to_csv(out / "py_contam.csv")
case, ref = st.comparison_sets(meta, C.COMPARISONS["BC_vs_Control"])
Ls = L.loc[study]
p, auc, fc = st.mw(Ls, case, ref)
b, pm, nc, nr, used = st.ols_group(Ls, meta, case, ref, C.COVARIATES)
pairs = st.age_match(meta, case, ref)
pd.DataFrame({"p_MW": p, "AUC": auc, "beta": b, "p_model": pm}, index=L.columns).to_csv(out / "py_stats.csv")
pd.DataFrame({"n_pairs": [len(pairs)], "sum_diff": [pairs["diff"].sum()]}).to_csv(out / "py_pairs.csv", index=False)
print("python reference done:", L.shape)

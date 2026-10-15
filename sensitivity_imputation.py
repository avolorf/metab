#!/usr/bin/env python3
"""Анализ чувствительности к импутации (после run_all.py).

    python sensitivity_imputation.py

(A) только соединения, обнаруженные во всех пробах (без импутированных значений): та же модель и правило находок,
    FDR внутри этого набора;
(B) все соединения, но импутированные значения заменены на пропуски (анализ по фактически обнаруженным значениям):
    модель, подбор по возрасту и сравнения с контролями считаются по доступным значениям каждого соединения.
Результат: results/tables/sensitivity_imputation.csv (сводка) и sensitivity_imputation_compounds.csv (по соединениям).
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pipeline import config as C
from pipeline import io
from pipeline import stats as st

warnings.filterwarnings("ignore")


def ols_one(y, D):
    m = ~np.isnan(y)
    if m.sum() < D.shape[1] + 3 or len(set(D[m, 1])) < 2:
        return np.nan, np.nan
    X, yy = D[m], y[m]
    XtX = np.linalg.pinv(X.T @ X)
    b = XtX @ X.T @ yy
    res = yy - X @ b
    dof = len(yy) - X.shape[1]
    se = np.sqrt((res ** 2).sum() / dof * XtX[1, 1])
    return b[1], 2 * stats.t.sf(abs(b[1] / se), dof)


def mw_one(a, b):
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    return stats.mannwhitneyu(a, b).pvalue, np.mean(a) - np.mean(b)


def complete_case(Lna, meta, spec, pairs, comp_cat):
    case, ref = st.comparison_sets(meta, spec)
    idx = case + ref
    D = pd.DataFrame({"g": [1.0 if s in set(case) else 0.0 for s in idx]}, index=idx)
    for c in C.COVARIATES:
        D[c] = meta.loc[idx, "sex"].map({"M": 1.0, "F": 0.0}) if c == "sex" else meta.loc[idx, c].astype(float)
    D = D.dropna()
    used = [c for c in C.COVARIATES if D[c].std() > 0]
    Dm = np.column_stack([np.ones(len(D)), D[["g"] + used].to_numpy()])
    lo = [s for s in ref if meta.contam_level[s] == "low"]
    hi = [s for s in ref if meta.contam_level[s] == "high"]
    rows = {}
    for f in Lna.columns:
        y = Lna.loc[D.index, f].to_numpy(float)
        b, p = ols_one(y, Dm)
        pm, fm = mw_one(Lna.loc[list(pairs.case), f].to_numpy(float), Lna.loc[list(pairs.ref), f].to_numpy(float))
        r = dict(beta_model=b, p_model=p, p_matched=pm, log2FC_matched=fm)
        if spec["ref"] == "Control":
            r["p_lowindex"], r["log2FC_lowindex"] = mw_one(Lna.loc[case, f].to_numpy(float), Lna.loc[lo, f].to_numpy(float))
            r["p_ctrl_hi_vs_lo"], _ = mw_one(Lna.loc[hi, f].to_numpy(float), Lna.loc[lo, f].to_numpy(float))
        rows[f] = r
    out = pd.DataFrame(rows).T
    out["q_model"] = st.bh(out.p_model.to_numpy(float))
    out["compound_category"] = comp_cat.reindex(out.index)
    out["finding"] = st.finding_rule(out, spec)
    return out


def main():
    t = Path(C.OUT) / "tables"
    meta = pd.read_csv(t / "samples.csv", index_col=0)
    meta.index = meta.index.astype(str)
    L = pd.read_csv(t / "processed_log2_matrix.csv", index_col=0)
    L = L[[c for c in L.columns if c in meta.index]].T
    L.index = L.index.astype(str)
    feat = pd.read_csv(t / "features.csv", index_col=0)
    pairs_all = pd.read_csv(t / "age_matched_pairs.csv", dtype={"case": str, "ref": str})
    X_all, *_ = io.load_features(C.DATA)
    study = meta.index[meta.study]
    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    comp_cat = feat.loc[reps, "compound_category"].fillna("")
    raw = X_all.reindex(index=study, columns=reps)
    detected = raw.notna() & (raw > 0)
    full = list(reps[detected.all().to_numpy()])
    Lna = L.loc[study, reps].where(detected)
    print(f"Соединений без пропусков: {len(full)} из {len(reps)}")
    rows, per = [], []
    for name, spec in C.COMPARISONS.items():
        r0 = pd.read_csv(t / f"stats_{name}.csv", index_col=0)
        f0 = set(r0.index[r0.finding])
        rA, _, _ = st.run_comparison(L.loc[study, full], meta, name, spec, comp_cat)
        rB = complete_case(Lna, meta, spec, pairs_all[pairs_all.comparison == name], comp_cat)
        fA, fB = set(rA.index[rA.finding]), set(rB.index[rB.finding])
        rows.append(dict(comparison=name, compounds_without_imputation=len(full), main_findings=len(f0),
                         main_findings_without_imputation=len(f0 & set(full)),
                         A_q05=int((rA.q_model < C.FDR).sum()), A_findings=len(fA), A_shared=len(fA & f0),
                         B_q05=int((rB.q_model < C.FDR).sum()), B_findings=len(fB), B_shared=len(fB & f0)))
        if name == "BC_vs_Control":
            x = rB[["beta_model", "p_model", "q_model", "p_matched", "p_lowindex", "finding"]].add_prefix("B_")
            x["imputed_fraction"] = 1 - detected.mean()
            x["in_full_set"] = x.index.isin(full)
            x = x.join(rA[["q_model", "finding"]].add_prefix("A_"))
            x["main_finding"] = x.index.isin(f0)
            per.append(x[x.main_finding])
        print(rows[-1])
    pd.DataFrame(rows).to_csv(t / "sensitivity_imputation.csv", index=False)
    pd.concat(per).to_csv(t / "sensitivity_imputation_compounds.csv")


if __name__ == "__main__":
    main()

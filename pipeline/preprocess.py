"""Фильтры, коррекция дрейфа (опционально), нормализация (PQN / MSTUS), импутация."""
import numpy as np
from scipy import stats
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from . import config as C


def rsd(df):
    return df.std(ddof=1) / df.mean()


def detection(X, meta):
    """Доля детекции (не-NaN до импутации) в QC и в каждой исследуемой группе."""
    out = pd.DataFrame(index=X.columns)
    out["det_QC"] = X.loc[meta.kind == "QC"].notna().mean()
    for g in C.STUDY_GROUPS:
        idx = meta.index[meta.study & (meta.group == g)]
        out[f"det_{g}"] = X.loc[idx].notna().mean()
    return out


def filter_features(X, meta, log):
    det = detection(X, meta)
    qc = meta.index[meta.kind == "QC"]
    keep = (det.det_QC >= C.MIN_QC_DETECT) & (det[[f"det_{g}" for g in C.STUDY_GROUPS]].max(axis=1) >= C.MIN_GROUP_DETECT)
    log(f"Фильтр детекции (QC ≥ {C.MIN_QC_DETECT:.0%}, ≥ {C.MIN_GROUP_DETECT:.0%} хотя бы в одной группе): "
        f"{X.shape[1]} → {keep.sum()}")
    blanks = meta.index[meta.kind == "blank"]
    det["blank_ratio"] = np.nan
    if len(blanks):
        ratio = X.loc[qc].median() / X.loc[blanks].fillna(0).median().replace(0, np.nan)
        ratio = ratio.fillna(np.inf)
        det["blank_ratio"] = ratio
        kb = ratio >= C.BLANK_RATIO
        log(f"Фильтр по бланкам ({len(blanks)} проб, QC/бланк ≥ {C.BLANK_RATIO}): удалено {(keep & ~kb).sum()}")
        keep &= kb
    else:
        log("Бланков в экспортированной таблице нет: по методике фильтр по бланкам (проба/бланк ≥ 5) выполнен в MS-DIAL; "
            "повторно не применяется.")
    if C.BATCHES and C.BATCH_CORRECTION and "batch" in meta:
        X = batch_correct(X, meta, keep, log)
    det["RSD_QC_raw"] = rsd(X.loc[qc])
    k2 = keep & (det.RSD_QC_raw <= C.MAX_QC_RSD)
    log(f"Фильтр RSD_QC ≤ {C.MAX_QC_RSD:.0%} (без коррекции дрейфа): {keep.sum()} → {k2.sum()}. Медианный RSD в QC: "
        f"{det.RSD_QC_raw[keep].median():.1%} до RSD-фильтра, {det.RSD_QC_raw[k2].median():.1%} после")
    det["pass_detection"] = keep
    det["kept"] = k2
    return X.loc[:, k2], det


def _trend(xq, yq, xall):
    fit = lowess(yq, xq, frac=1.0, it=0, delta=0, xvals=np.asarray(xall, float))
    lo, hi = xq.min(), xq.max()
    fit = np.where(xall < lo, fit[np.argmin(np.abs(xall - lo))], fit)
    return np.where(xall > hi, fit[np.argmin(np.abs(xall - hi))], fit)


def batch_correct(X, meta, keep, log):
    """Межсерийный сдвиг: для каждого признака значения серии k умножаются на (медиана всех QC) / (медиана QC серии k).
    Оценка до/после: число признаков с различием QC между сериями (Краскел–Уоллис, BH q < 0.05) и RSD QC."""
    from statsmodels.stats.multitest import multipletests
    qc = meta.index[meta.kind == "QC"]
    b = meta.batch
    feats = X.columns[keep]
    def n_diff(Y):
        Q = np.log2(Y.loc[qc, feats])
        p = []
        for f in feats:
            g = [Q.loc[qc[b[qc] == k], f].dropna() for k in sorted(b[qc].dropna().unique())]
            p.append(stats.kruskal(*g).pvalue if min(len(x) for x in g) >= 3 else np.nan)
        p = np.nan_to_num(np.array(p), nan=1.0)
        return int((multipletests(p, method="fdr_bh")[1] < 0.05).sum())
    before, rsd_b = n_diff(X), rsd(X.loc[qc, feats]).median()
    Y = X.copy()
    ref = X.loc[qc].median()
    for k in sorted(b[qc].dropna().unique()):
        fac = ref / X.loc[qc[b[qc] == k]].median()
        rows = meta.index[(b == k) & meta.index.isin(X.index)]
        Y.loc[rows] = X.loc[rows] * fac.replace([np.inf, -np.inf], np.nan).fillna(1.0)
    after, rsd_a = n_diff(Y), rsd(Y.loc[qc, feats]).median()
    log(f"Коррекция межсерийного сдвига по QC ({len(C.BATCHES)} серии): признаков с различием QC между сериями "
        f"(q < 0.05) {before} → {after} из {len(feats)}; медианный RSD QC {rsd_b:.1%} → {rsd_a:.1%}")
    return Y


def drift_loo(X, order, qc, batch=None):
    """Если заданы серии, тест выполняется внутри каждой серии (QC той же серии), ошибки объединяются."""
    if batch is not None and C.BATCHES:
        errs = {"none": [], "lowess": []}
        for k in sorted(batch[qc].dropna().unique()):
            r = _drift_loo_errors(X, order, qc[batch[qc] == k])
            for m in errs:
                errs[m] += r[m]
        return {k: float(np.median(v)) for k, v in errs.items()}
    return {k: float(np.median(v)) for k, v in _drift_loo_errors(X, order, qc).items()}


def _drift_loo_errors(X, order, qc):
    """Leave-one-QC-out: RSD внутренних QC, спрогнозированных без них (none vs lowess)."""
    L = np.log2(X.loc[qc].to_numpy())
    xq = order.loc[qc].to_numpy(float)
    inner = np.argsort(xq)[1:-1]
    err = {"none": [], "lowess": []}
    for j in range(L.shape[1]):
        y = L[:, j]
        for k in inner:
            m = ~np.isnan(y)
            m[k] = False
            if np.isnan(y[k]) or m.sum() < 4:
                continue
            err["none"].append(abs(y[k] - np.median(y[m])))
            err["lowess"].append(abs(y[k] - _trend(xq[m], y[m], np.array([xq[k]]))[0]))
    return err


def drift_correct(X, order, qc):
    xq = order.loc[qc].to_numpy(float)
    xall = order.loc[X.index].to_numpy(float)
    L = np.log2(X.to_numpy())
    Lq = np.log2(X.loc[qc].to_numpy())
    for j in range(L.shape[1]):
        m = ~np.isnan(Lq[:, j])
        if m.sum() >= 4:
            L[:, j] = L[:, j] - _trend(xq[m], Lq[m, j], xall) + np.median(Lq[m, j])
    return pd.DataFrame(2 ** L, index=X.index, columns=X.columns)


def pqn_factors(X, ref_idx):
    ref = X.loc[ref_idx].median()
    return X.div(ref, axis=1).median(axis=1, skipna=True)


def mstus_factors(X, idx, qc):
    """MSTUS: сумма интенсивностей признаков, обнаруженных во ВСЕХ пробах (исследуемых и QC),
    отнесённая к медиане QC."""
    full = X.loc[idx].notna().all()
    tot = X.loc[:, full].sum(axis=1)
    return tot / tot.loc[qc].median(), int(full.sum())


def normalize(X, meta, method):
    qc = meta.index[meta.kind == "QC"]
    idx = meta.index[meta.study | (meta.kind == "QC")]
    if method == "PQN":
        f = pqn_factors(X, qc)
        info = "PQN (референс — медианный спектр QC)"
    elif method == "MSTUS":
        f, n = mstus_factors(X, idx, qc)
        info = f"MSTUS ({n} признаков, обнаруженных во всех пробах)"
    else:
        raise ValueError(method)
    return X.div(f, axis=0), f, info


def impute_log(Xn, idx):
    """½ минимума признака по включённым пробам и QC, затем log2."""
    Xi = Xn.loc[idx]
    filled = Xi.fillna(Xi.min() / 2)
    return np.log2(filled)

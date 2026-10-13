"""Многомерный анализ без утечки: повторная CV, отбор признаков внутри фолда, перестановочные тесты."""
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold

from . import config as C

FPR = np.linspace(0, 1, 101)


def pareto_fit(X):
    mu, sd = X.mean(0), np.sqrt(X.std(0, ddof=1))
    sd[sd == 0] = 1
    return mu, sd


def residualize(Xtr, Xte, ztr, zte):
    """Убирает линейную связь признаков с ковариатой z; коэффициенты — только по обучающему фолду."""
    zc = ztr - ztr.mean()
    beta = (zc @ (Xtr - Xtr.mean(0))) / (zc @ zc)
    return Xtr - np.outer(zc, beta), Xte - np.outer(zte - ztr.mean(), beta)


def cv_predict(X, y, method, seed, z=None):
    """Одна 5-fold CV. method: 'panel' (top-k Манн–Уитни в фолде + логистическая регрессия) | 'plsda' |
    'panel_adj' (как panel, но признаки предварительно скорректированы на ковариату z внутри фолда)."""
    skf = StratifiedKFold(C.CV_FOLDS, shuffle=True, random_state=seed)
    pred = np.zeros(len(y))
    chosen = []
    for tr, te in skf.split(X, y):
        Xtr, Xte = X[tr], X[te]
        if method == "panel_adj":
            Xtr, Xte = residualize(Xtr, Xte, z[tr], z[te])
        if method in ("panel", "panel_adj"):
            p = stats.mannwhitneyu(Xtr[y[tr] == 1], Xtr[y[tr] == 0], axis=0).pvalue
            sel = np.argsort(p)[:C.PANEL_K]
            chosen.extend(sel)
            mu, sd = Xtr[:, sel].mean(0), Xtr[:, sel].std(0) + 1e-9
            m = LogisticRegression(C=1.0, max_iter=2000).fit((Xtr[:, sel] - mu) / sd, y[tr])
            pred[te] = m.predict_proba((Xte[:, sel] - mu) / sd)[:, 1]
        else:
            mu, sd = pareto_fit(Xtr)
            m = PLSRegression(2, scale=False).fit((Xtr - mu) / sd, y[tr])
            pred[te] = m.predict((Xte - mu) / sd).ravel()
    return pred, chosen


def repeated_cv(X, y, method, seed, z=None):
    rng = np.random.RandomState(seed)
    aucs, tprs, chosen, q2 = [], [], [], []
    for _ in range(C.CV_REPEATS):
        pred, ch = cv_predict(X, y, method, rng.randint(1e9), z)
        aucs.append(roc_auc_score(y, pred))
        fpr, tpr, _ = roc_curve(y, pred)
        tprs.append(np.interp(FPR, fpr, tpr))
        chosen.extend(ch)
        q2.append(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return np.array(aucs), np.array(tprs), np.array(chosen, int), float(np.mean(q2))


def _perm_one(X, y, method, seed, z=None):
    rng = np.random.RandomState(seed)
    yp = rng.permutation(y)          # z остаётся привязан к пробам
    pred, _ = cv_predict(X, yp, method, rng.randint(1e9), z)
    return roc_auc_score(yp, pred)


def perm_test(X, y, method, observed, seed, n_perm=C.N_PERM, z=None):
    seeds = np.random.RandomState(seed).randint(1e9, size=n_perm)
    null = Parallel(n_jobs=C.N_JOBS)(delayed(_perm_one)(X, y, method, s, z) for s in seeds)
    null = np.array(null)
    return (np.sum(null >= observed) + 1) / (n_perm + 1), null


def boot_auc(y, s, seed, n=C.N_BOOT):
    rng = np.random.RandomState(seed)
    a = roc_auc_score(y, s)
    sign = 1 if a >= 0.5 else -1
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    b = [_boot1(y, sign * s, i1, i0, rng) for _ in range(n)]
    return max(a, 1 - a), np.percentile(b, 2.5), np.percentile(b, 97.5)


def _boot1(y, s, i1, i0, rng):
    i = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
    return roc_auc_score(y[i], s[i])


def evaluate(L, meta, case, ref, pool, label, seed, names=None, n_perm=None):
    """Панель и PLS-DA на наборе case/ref по признакам pool (без экспозиций)."""
    idx = case + ref
    y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
    X = L.loc[idx, pool].to_numpy()
    z = meta.loc[idx, "contam_index"].to_numpy(float)
    rows, curves = [], {}
    for method, nm in [("panel", f"панель top-{C.PANEL_K} + логистическая регрессия (отбор внутри CV)"),
                       ("panel_adj", f"панель top-{C.PANEL_K}, признаки скорректированы на индекс загрязнения внутри фолда"),
                       ("plsda", "PLS-DA, 2 компоненты (все признаки без экспозиций)")]:
        aucs, tprs, chosen, q2 = repeated_cv(X, y, method, seed, z)
        p, _ = perm_test(X, y, method, aucs.mean(), seed + 1, n_perm or C.N_PERM, z) if (n_perm or C.N_PERM) else (np.nan, None)
        sel = ""
        if len(chosen):
            vc = pd.Series(chosen).value_counts().head(10)
            sel = "; ".join(f"{names[pool[i]] if names is not None else pool[i]} ({n / (C.CV_REPEATS * C.CV_FOLDS):.0%})"
                            for i, n in vc.items())
        rows.append(dict(set=label, model=nm, n_case=len(case), n_ref=len(ref), n_features=len(pool),
                         AUC_mean=aucs.mean(), AUC_sd=aucs.std(), AUC_p5=np.percentile(aucs, 5),
                         AUC_p95=np.percentile(aucs, 95), Q2=q2, perm_p=p, n_perm=C.N_PERM, most_selected=sel))
        curves[(label, method)] = (tprs.mean(0), aucs.mean())
    # базовые модели: только пол / только индекс загрязнения
    for cov in ["sex", "contam_index"]:
        v = meta.loc[idx, cov]
        v = v.map({"M": 1.0, "F": 0.0}) if cov == "sex" else v.astype(float)
        m = v.notna().to_numpy()
        if m.sum() > 5 and v[m].std() > 0 and len(set(y[m])) == 2:
            a = roc_auc_score(y[m], v[m])
            rows.append(dict(set=label, model=f"только {cov} (базовая линия; {'выше' if a >= 0.5 else 'ниже'} у case)",
                             n_case=int(y[m].sum()), n_ref=int((1 - y[m]).sum()), AUC_mean=max(a, 1 - a)))
            fpr, tpr, _ = roc_curve(y[m], v[m] if a >= 0.5 else -v[m])
            curves[(label, cov)] = (np.interp(FPR, fpr, tpr), max(a, 1 - a))
    return rows, curves


def pca_scores(L, n=3):
    X = L.to_numpy()
    mu, sd = pareto_fit(X)
    p = PCA(n).fit((X - mu) / sd)
    return pd.DataFrame(p.transform((X - mu) / sd), index=L.index), p.explained_variance_ratio_

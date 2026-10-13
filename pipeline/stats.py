"""Одномерная статистика: подбор по возрасту, основная модель с ковариатами, проверки устойчивости."""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
from statsmodels.stats.multitest import multipletests

from . import config as C


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    if m.any():
        q[m] = multipletests(p[m], method="fdr_bh")[1]
    return q


# ------------------------------------------------------------------ группы сравнения
def comparison_sets(meta, spec):
    s = meta[meta.study]
    if spec["sex"]:
        s = s[s.sex == spec["sex"]]
    return list(s.index[s.group == spec["case"]]), list(s.index[s.group == spec["ref"]])


def age_match(meta, case, ref):
    """1:1 по ближайшему возрасту (оптимальное назначение), допуск ±caliper, внутри пола."""
    pairs = []
    strata = ["M", "F"] if C.MATCH_WITHIN_SEX else [None]
    for sx in strata:
        a = [s for s in case if pd.notna(meta.age[s]) and (sx is None or meta.sex[s] == sx)]
        b = [s for s in ref if pd.notna(meta.age[s]) and (sx is None or meta.sex[s] == sx)]
        if not a or not b:
            continue
        D = np.abs(meta.age[a].to_numpy()[:, None] - meta.age[b].to_numpy()[None, :])
        cost = np.where(D <= C.MATCH_CALIPER_YEARS, D, 1e6)
        ri, ci = linear_sum_assignment(cost)
        for i, j in zip(ri, ci):
            if cost[i, j] < 1e6:
                pairs.append(dict(case=a[i], ref=b[j], sex=sx or "", age_case=meta.age[a[i]],
                                  age_ref=meta.age[b[j]], diff=D[i, j]))
    return pd.DataFrame(pairs, columns=["case", "ref", "sex", "age_case", "age_ref", "diff"])


def demographics(meta):
    s = meta[meta.study]
    rows = []
    for g in C.STUDY_GROUPS:
        d = s[s.group == g]
        rows.append(dict(group=g, n=len(d), male=int((d.sex == "M").sum()), female=int((d.sex == "F").sum()),
                         sex_unknown=int((~d.sex.isin(["M", "F"])).sum()), age_n=int(d.age.notna().sum()),
                         age_missing=int(d.age.isna().sum()), age_mean=d.age.mean(), age_sd=d.age.std(),
                         age_median=d.age.median(), age_min=d.age.min(), age_max=d.age.max(),
                         male_age_ge60=int(((d.sex == "M") & (d.age >= 60)).sum())))
    tab = pd.DataFrame(rows)
    ages = [s.age[(s.group == g)].dropna() for g in C.STUDY_GROUPS]
    tests = [dict(test="возраст, Краскел–Уоллис (3 группы)", p=stats.kruskal(*ages).pvalue)]
    for a, b in [("BC", "Control"), ("PC", "Control"), ("BC", "PC")]:
        x, y = s.age[s.group == a].dropna(), s.age[s.group == b].dropna()
        tests.append(dict(test=f"возраст {a} vs {b}, Манн–Уитни", p=stats.mannwhitneyu(x, y).pvalue))
        ct = pd.crosstab(s.group[s.group.isin([a, b]) & s.sex.isin(["M", "F"])], s.sex)
        tests.append(dict(test=f"пол {a} vs {b}, точный тест Фишера", p=stats.fisher_exact(ct.to_numpy()).pvalue))
    xm = s.age[(s.group == "PC") & (s.sex == "M")].dropna()
    ym = s.age[(s.group == "Control") & (s.sex == "M")].dropna()
    if len(ym) > 1:
        tests.append(dict(test="возраст PC vs мужчины-контроли, Манн–Уитни", p=stats.mannwhitneyu(xm, ym).pvalue))
    return tab, pd.DataFrame(tests)


# ------------------------------------------------------------------ модели
def mw(L, a, b):
    if len(a) < 2 or len(b) < 2:
        n = L.shape[1]
        return np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    A, B = L.loc[a].to_numpy(), L.loc[b].to_numpy()
    U = stats.mannwhitneyu(A, B, axis=0, alternative="two-sided")
    return U.pvalue, U.statistic / (len(a) * len(b)), A.mean(0) - B.mean(0)


def ols_group(L, meta, case, ref, covariates):
    """log2 ~ группа + ковариаты (complete case). Возвращает коэффициент группы, p, n и использованные ковариаты."""
    idx = case + ref
    D = pd.DataFrame({"group": [1.0 if s in set(case) else 0.0 for s in idx]}, index=idx)
    used = []
    for c in covariates:
        if c == "sex":
            v = meta.loc[idx, "sex"].map({"M": 1.0, "F": 0.0})
        else:
            v = meta.loc[idx, c].astype(float)
        D[c] = v
    D = D.dropna()
    for c in covariates:
        if D[c].std() > 0:
            used.append(c)
    D = D[["group"] + used]
    X = np.column_stack([np.ones(len(D)), D.to_numpy()])
    y = L.loc[D.index].to_numpy()
    XtX = np.linalg.pinv(X.T @ X)
    beta = XtX @ X.T @ y
    res = y - X @ beta
    dof = len(D) - X.shape[1]
    se = np.sqrt((res ** 2).sum(0) / dof * XtX[1, 1])
    p = 2 * stats.t.sf(np.abs(beta[1] / se), dof)
    n_case = int(D.group.sum())
    return beta[1], p, n_case, len(D) - n_case, used


def run_comparison(L, meta, name, spec, comp_cat):
    """Все варианты для одного сравнения на уровне соединений (L — представители соединений)."""
    case, ref = comparison_sets(meta, spec)
    out = pd.DataFrame(index=L.columns)
    info = {"comparison": name, "n_case": len(case), "n_ref": len(ref)}
    # 1. все контроли: Манн–Уитни и основная модель
    p, auc, fc = mw(L, case, ref)
    out["log2FC"], out["AUC"], out["p_MW"], out["q_MW"] = fc, auc, p, bh(p)
    b, p, nc, nr, used = ols_group(L, meta, case, ref, C.COVARIATES)
    out["beta_model"], out["p_model"], out["q_model"] = b, p, bh(p)
    info.update(model_covariates="+".join(used), model_n_case=nc, model_n_ref=nr)
    # 2. подбор по возрасту
    pairs = age_match(meta, case, ref)
    pairs.insert(0, "comparison", name)
    a, r = list(pairs.case), list(pairs.ref)
    p, auc, fc = mw(L, a, r)
    out["log2FC_matched"], out["AUC_matched"], out["p_matched"], out["q_matched"] = fc, auc, p, bh(p)
    info["matched_pairs"] = len(pairs)
    # 3. контроли с низким индексом загрязнения (только если референс — контроль)
    if spec["ref"] == "Control":
        lo = [s for s in ref if meta.contam_level[s] == "low"]
        hi = [s for s in ref if meta.contam_level[s] == "high"]
        p, _, fc = mw(L, case, lo)
        out["log2FC_lowindex"], out["p_lowindex"] = fc, p
        p, _, fc = mw(L, hi, lo)
        out["log2FC_ctrl_hi_vs_lo"], out["p_ctrl_hi_vs_lo"] = fc, p
        info.update(n_ref_lowindex=len(lo), n_ref_highindex=len(hi))
    out["compound_category"] = comp_cat.reindex(out.index)
    out["finding"] = finding_rule(out, spec)
    info["n_q_model"] = int((out.q_model < C.FDR).sum())
    info["n_findings"] = int(out.finding.sum())
    return out, info, pairs


def finding_rule(r, spec):
    """Находка: q(модель с ковариатами) < FDR И p(подбор по возрасту) < CONFIRM_P в том же направлении;
    для сравнений с контролем дополнительно: p(против контролей с низким индексом) < CONFIRM_P в том же
    направлении И нет различия между контролями с высоким и низким индексом (p ≥ CTRL_INDEX_P).
    Экспозиции (контаминанты, пища, лекарства) находками не считаются."""
    sgn = np.sign(r.beta_model)
    ok = (r.q_model < C.FDR) & (r.p_matched < C.CONFIRM_P) & (np.sign(r.log2FC_matched) == sgn)
    if "p_lowindex" in r:
        ok &= (r.p_lowindex < C.CONFIRM_P) & (np.sign(r.log2FC_lowindex) == sgn) & (r.p_ctrl_hi_vs_lo >= C.CTRL_INDEX_P)
    from .features import EXPOSURE
    ok &= ~r.compound_category.isin(EXPOSURE)
    return ok.fillna(False)

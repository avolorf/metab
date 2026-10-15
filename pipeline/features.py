"""Группировка признаков в соединения, индекс загрязнения, категории экспозиций, поиск подтверждающих ионов."""
import re

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C

PROTON = 1.007276
MOLECULAR_IONS = ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+"]
ADDUCTS = {"[M+H]+": PROTON, "[M+Na]+": 22.989218, "[M+K]+": 38.963158, "[M+NH4]+": 18.033823,
           "[M+H-H2O]+": PROTON - 18.010565, "[M+H-2H2O]+": PROTON - 2 * 18.010565}
# разности масс между ионами одного соединения
RELATIONS = {"13C": 1.003355, "Na-H": 21.981942, "K-H": 37.955882, "NH4-H": 17.026547, "H2O": 18.010565,
             "2H2O": 36.021130, "CO2": 43.989829, "HCOOH": 46.005479, "CO": 27.994915}
ELEM = {"C": 12.0, "H": 1.00782503, "N": 14.0030740, "O": 15.9949146, "P": 30.9737620,
        "S": 31.9720707, "Cl": 34.9688527, "Br": 78.9183371, "Na": 22.9897693, "K": 38.9637064}


def formula_mass(f):
    return sum(ELEM[e] * (int(n) if n else 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", f))


def ppm(obs, theo):
    return (obs - theo) / theo * 1e6


# ------------------------------------------------------------------ группировка
def group_features(L, info, raw_median):
    """Соединения = связные компоненты графа: |ΔRT| ≤ tol и корреляция Спирмена ≥ порога.
    Представитель — признак с наибольшей медианной интенсивностью."""
    feats = L.columns.to_numpy()
    rt = info.loc[feats, "rt_min"].to_numpy()
    o = np.argsort(rt)
    R = L.rank().to_numpy()
    R = (R - R.mean(0)) / (R.std(0) + 1e-12)
    n = len(feats)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edges = []
    for a in range(n):
        i = o[a]
        b = a + 1
        while b < n and rt[o[b]] - rt[i] <= C.GROUP_RT_TOL:
            j = o[b]
            r = float(R[:, i] @ R[:, j]) / len(R)
            if r >= C.GROUP_MIN_CORR:
                edges.append((i, j, r))
                pi, pj = find(i), find(j)
                if pi != pj:
                    parent[pi] = pj
            b += 1
    comp = np.array([find(i) for i in range(n)])
    df = pd.DataFrame({"feature": feats, "root": comp, "rt_min": rt, "mz": info.loc[feats, "mz"].to_numpy(),
                       "median_int": raw_median.loc[feats].to_numpy()})
    rep = df.sort_values("median_int", ascending=False).drop_duplicates("root").set_index("root")["feature"]
    df["representative"] = df.root.map(rep)
    ids = {r: f"C{k:05d}" for k, r in enumerate(df.sort_values("rt_min").root.unique(), 1)}
    df["compound"] = df.root.map(ids)
    df["n_ions"] = df.groupby("compound").feature.transform("size")
    # отношения масс внутри соединения
    rel = {f: [] for f in feats}
    for i, j, r in edges:
        d = abs(df.mz[i] - df.mz[j])
        for nm, v in RELATIONS.items():
            if abs(d - v) < 0.003:
                lo, hi = (i, j) if df.mz[i] < df.mz[j] else (j, i)
                rel[feats[hi]].append(f"{nm} от {df.mz[lo]:.4f}")
                rel[feats[lo]].append(f"{df.mz[hi]:.4f} (+{nm})")
    df["ion_relations"] = [";".join(rel[f]) for f in feats]
    df["is_representative"] = df.feature == df.representative
    return df.drop(columns="root").set_index("feature")


# ------------------------------------------------------------------ индекс загрязнения
def match_series(info, path=C.CONTAMINANT_SERIES):
    s = pd.read_csv(path)
    rows = []
    for _, r in s.iterrows():
        theo = formula_mass(r.formula) + ADDUCTS[r.ion]
        d = info[(ppm(info.mz, theo).abs() <= C.PPM) & ((info.rt_min - r.rt_min).abs() <= 0.05)]
        for f, x in d.iterrows():
            rows.append(dict(feature=f, label=r.label, formula=r.formula, ion=r.ion, rt_min=x.rt_min, mz=x.mz,
                             ppm=round(ppm(x.mz, theo), 2), in_index=int(r.in_index)))
    return pd.DataFrame(rows).drop_duplicates("feature").set_index("feature")


def contamination_index(L, series, study_idx):
    f = series.index[(series.in_index == 1) & series.index.isin(L.columns)]
    Z = L.loc[:, f]
    Z = (Z - Z.loc[study_idx].mean()) / Z.loc[study_idx].std()
    idx = Z.mean(axis=1)
    thr = idx.loc[study_idx].median()
    return idx, (idx > thr).map({True: "high", False: "low"}), thr, list(f)


# ------------------------------------------------------------------ категории экспозиций
def exposure_category(info, path=C.EXPOSURE_KEYWORDS):
    kw = pd.read_csv(path)
    cat = pd.Series("", index=info.index)
    for _, r in kw.iterrows():
        m = info.annotation.str.contains(r.pattern, flags=re.I, regex=True) & (cat == "")
        cat[m] = r.category
    return cat


def categorize(info, groups, series, curated, kw_cat):
    """Категория признака и соединения: contaminant > drug > tobacco > food > endogenous/unknown."""
    cat = pd.Series("unknown", index=info.index)
    cat[info.annotation_conf == "high"] = "annotated (high score)"
    cat[info.annotation_conf == "low"] = "putative (low score)"
    for f, c in kw_cat.items():
        if c:
            cat[f] = c if info.loc[f, "annotation_conf"] == "high" else f"possible {c} (low score)"
    cmap = {"food": "food", "drug": "drug", "endogenous": "endogenous", "ambiguous": "ambiguous"}
    for _, r in curated.iterrows():
        if isinstance(r.feature, str) and r.feature in cat.index:
            for k, v in cmap.items():
                if str(r.category).startswith(k):
                    cat[r.feature] = v
    cat[series.index.intersection(cat.index)] = "contaminant"
    out = groups[["compound"]].copy()
    out["category"] = cat.reindex(out.index)
    prio = ["contaminant", "drug", "tobacco", "food", "possible drug (low score)", "possible tobacco (low score)",
            "possible food (low score)", "endogenous", "ambiguous", "annotated (high score)", "putative (low score)",
            "unknown"]
    rank = out.category.map({c: i for i, c in enumerate(prio)})
    best = out.assign(rank=rank).sort_values("rank").drop_duplicates("compound").set_index("compound")["category"]
    out["compound_category"] = out.compound.map(best)
    return out


EXPOSURE = {"contaminant", "drug", "tobacco", "food"}


def panel_excluded(cat):
    ex = cat.isin(EXPOSURE)
    if C.EXCLUDE_DRUG_LIKE_UNKNOWN_FROM_PANEL:
        ex |= cat.str.startswith("possible drug") | cat.str.startswith("possible tobacco")
    return ex


def display_name(info, cat):
    """Имя для рисунков: только уверенные аннотации; low score (в т. ч. «лекарственные») -> unknown.
    Предполагаемое название low-score признаков остаётся в таблицах (столбец annotation)."""
    out = []
    for f, r in info.iterrows():
        tag = f"{r.mz:.4f}@{r.rt_min:.2f}"
        if r.annotation_conf == "high":
            out.append(f"{r.annotation[:32]} {tag}")
        else:
            out.append(f"unknown {tag}" + ("*" if r.annotation_conf == "low" else ""))
    return pd.Series(out, index=info.index)


# ------------------------------------------------------------------ подтверждающие ионы
def confirming_ions(feature, formula, info_all, Lraw):
    """Ищет в ПОЛНОЙ таблице (до фильтров) другие ионы того же M при том же RT (±0.03 мин, ≤5 ppm)."""
    rt = info_all.loc[feature, "rt_min"]
    M = formula_mass(formula)
    found = []
    near = info_all[(info_all.rt_min - rt).abs() <= 0.03]
    cands = {ion: M + sh for ion, sh in ADDUCTS.items()}
    obs = info_all.loc[feature, "mz"]
    cands["13C-изотоп"] = obs + 1.003355
    for ion, mz in cands.items():
        d = near[ppm(near.mz, mz).abs() <= C.PPM]
        for f2, r in d.iterrows():
            if f2 == feature:
                continue
            a, b = Lraw[feature], Lraw[f2]
            m = a.notna() & b.notna()
            rho = stats.spearmanr(a[m], b[m]).statistic if m.sum() > 5 else np.nan
            found.append(f"{ion} {r.mz:.4f} ({ppm(r.mz, mz):+.1f} ppm, ρ={rho:.2f})")
    return "; ".join(found)

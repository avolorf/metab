"""Этап 1 (s107): переобработка по клиническим подгруппам внутри РМП.

Матрица соединений — та же, что в основном анализе s107 (prepare_matrix.py проверяет совпадение чисел).
Выборка — только РМП (group = BC). Сравнения:
  S1 — стадия: Ta против T1 (T2plus, ambiguous, unknown, PUNLMP_uncertain — только описательно);
  S2 — степень злокачественности: LG против HG (G2_1973, not_stated, PUNLMP_uncertain исключены).
Модель для каждого соединения: log2 ~ подгруппа + пол + возраст + индекс экзогенного сигнала; для S2 — также
вариант с поправкой на стадию (только пробы Ta/T1); анализ чувствительности — только мужчины (log2 ~ подгруппа +
возраст + индекс).

Блок А (results_s107_clinical/steroids/) — три глюкуронида стероидов и ионы их агликонов.
Блок Б (results_s107_clinical/clinical/) — без стероидов: «кровяная» оценка (после согласования списка),
ненаправленное сравнение S1 и S2, классификация S1 (вложенная CV 5 × 20, перестановочный тест, базовая модель
«пол + возраст»).

Запуск (из этой папки):
  S107_PIPELINE=<копия ветки основного анализа> python run_clinical.py --clinical <s107_clinical_groups.xlsx> \
      [--blood_list blood_score_list.csv] [--n_perm 1000]
"""
import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
CWD = Path.cwd()                         # prepare_matrix переходит в папку конвейера: пути аргументов — от исходной папки
sys.path.insert(0, str(HERE))
import prepare_matrix  # noqa: E402  (переходит в папку конвейера и подключает pipeline)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.cross_decomposition import PLSRegression  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score, roc_curve  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

from pipeline import config as C, features  # noqa: E402

# ---------------------------------------------------------------- параметры расчёта
P = dict(
    seed=C.SEED,                         # то же начальное значение, что в основном анализе s107
    fdr=0.05,
    covariates=["sex", "age", "contam_index"],
    covariates_men=["age", "contam_index"],
    auc_bootstrap=2000,                  # стратифицированный бутстрэп для 95 % ДИ AUC
    jt_permutations=10000,               # перестановочный p для критерия Йонкхира
    cv_folds=5, cv_repeats=20,           # внешняя CV классификации S1
    inner_folds=5,                       # внутренняя CV: выбор размера панели / числа компонент PLS-DA
    panel_k_grid=[3, 5, 10, 20],
    plsda_comp_grid=[1, 2, 3],
    n_perm=1000,                         # перестановки классификации (вся процедура CV на каждую перестановку)
    pool_max_abs_rho_index=C.PANEL_MAX_ABS_RHO_CONTAM,
    top_n_if_no_q=50,
    blood_z_reference="все 104 пробы основного анализа (среднее и SD каждого соединения)",
)
S1 = dict(name="S1_stage", col="stage_group", case="T1", ref="Ta", label="T1 vs Ta",
          descriptive_only=["T2plus", "ambiguous", "unknown", "PUNLMP_uncertain"])
S2 = dict(name="S2_grade", col="grade_group", case="HG", ref="LG", label="HG vs LG",
          excluded=["G2_1973", "not_stated", "PUNLMP_uncertain"])
# три глюкуронида стероидов основного анализа s107 (report_s107.py, S266_T2): m/z и RT представителя
STEROIDS = [("11-Oxoetiocholanolone glucuronide", 498.2697, 10.73),
            ("Tetrahydrocortisone glucuronide, isomer I", 541.2642, 10.76),
            ("Tetrahydrocortisone glucuronide, isomer II", 558.2914, 10.43)]
# аннотации стероидов и стероидных конъюгатов (исключаются из блока Б вместе со всем соединением)
STEROID_RE = re.compile(
    r"andro?st|pregn|estr[aoi]|oestr|cortis|cortic|cortol|corton|testost|dehydroepiandro|\bDHEA|etiocholan|aldoster|"
    r"progest|steroid|sterone|\bchol(?:an|en|est)|cholic|sterol|-17-one|17-ol|-11,17-dione", re.I)

BILE_RE = re.compile(r"cholic|cholan|cholen|\bsterol\b|cholest", re.I)

LAB = {"Control": "Контроль", "Ta": "Ta", "T1": "T1", "LG": "LG", "HG": "HG", "T2plus": "T2+"}
COL = {"Control": "#8c8c8c", "Ta": "#6baed6", "T1": "#08519c", "LG": "#fdae6b", "HG": "#d94801", "T2plus": "#000000"}


# ---------------------------------------------------------------- статистика
def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    if m.sum():
        pv = p[m]
        o = np.argsort(pv)
        q = pv[o] * m.sum() / np.arange(1, m.sum() + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        r = np.empty_like(q)
        r[o] = np.minimum(q, 1)
        out[m] = r
    return out


def design(meta, case, ref, covariates, extra=None):
    idx = list(case) + list(ref)
    D = pd.DataFrame({"group": [1.0] * len(case) + [0.0] * len(ref)}, index=idx)
    for c in covariates:
        D[c] = meta.loc[idx, "sex"].map({"M": 1.0, "F": 0.0}) if c == "sex" else meta.loc[idx, c].astype(float)
    if extra is not None:
        D = D.join(extra)
    D = D.dropna()
    used = [c for c in D.columns[1:] if D[c].std() > 0]
    return D[["group"] + used], used


def ols(L, meta, case, ref, covariates, extra=None):
    """log2 ~ группа + ковариаты (полные случаи); коэффициент группы (log2FC) с 95 % ДИ, p."""
    D, used = design(meta, case, ref, covariates, extra)
    X = np.column_stack([np.ones(len(D)), D.to_numpy()])
    y = L.loc[D.index].to_numpy()
    XtX = np.linalg.pinv(X.T @ X)
    beta = XtX @ X.T @ y
    res = y - X @ beta
    dof = len(D) - X.shape[1]
    se = np.sqrt((res ** 2).sum(0) / dof * XtX[1, 1])
    t = stats.t.ppf(0.975, dof)
    p = 2 * stats.t.sf(np.abs(beta[1] / se), dof)
    n1 = int(D.group.sum())
    return pd.DataFrame({"log2FC": beta[1], "CI_low": beta[1] - t * se, "CI_high": beta[1] + t * se, "p": p,
                         "n_case": n1, "n_ref": len(D) - n1}, index=L.columns), used


def auc_ci(x_case, x_ref, rng, n_boot):
    """AUC = P(значение в группе case > значения в группе ref) по log2 без поправок; 95 % ДИ — стратифицированный бутстрэп."""
    y = np.r_[np.ones(len(x_case)), np.zeros(len(x_ref))]
    s = np.r_[x_case, x_ref]
    a = roc_auc_score(y, s)
    i1, i0 = np.arange(len(x_case)), np.arange(len(x_case), len(s))
    b = []
    for _ in range(n_boot):
        i = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
        b.append(roc_auc_score(y[i], s[i]))
    return a, np.percentile(b, 2.5), np.percentile(b, 97.5)


def auc_fast(L, case, ref):
    """AUC для всех столбцов сразу (через ранги, = Манн–Уитни U / n1n0)."""
    A = L.loc[list(case) + list(ref)].rank().to_numpy()
    n1, n0 = len(case), len(ref)
    return pd.Series((A[:n1].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0), index=L.columns)


def jonckheere(groups, rng, n_perm):
    """Критерий Йонкхира–Терпстры для упорядоченных групп (возрастающий тренд); p двусторонний —
    асимптотический (с поправкой на связи) и перестановочный."""
    def jt_stat(gs):
        s = 0.0
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                d = gs[j][None, :] - gs[i][:, None]
                s += (d > 0).sum() + 0.5 * (d == 0).sum()
        return s
    gs = [np.asarray(g, float) for g in groups]
    J = jt_stat(gs)
    n = np.array([len(g) for g in gs])
    N = n.sum()
    allv = np.concatenate(gs)
    _, t = np.unique(allv, return_counts=True)
    mean = (N ** 2 - (n ** 2).sum()) / 4
    var = (N * (N - 1) * (2 * N + 5) - (n * (n - 1) * (2 * n + 5)).sum() - (t * (t - 1) * (2 * t + 5)).sum()) / 72 \
        + ((n * (n - 1) * (n - 2)).sum() * (t * (t - 1) * (t - 2)).sum()) / (36 * N * (N - 1) * (N - 2)) \
        + ((n * (n - 1)).sum() * (t * (t - 1)).sum()) / (8 * N * (N - 1))
    z = (J - mean) / np.sqrt(var)
    cuts = np.cumsum(n)[:-1]
    null = np.array([jt_stat(np.split(rng.permutation(allv), cuts)) for _ in range(n_perm)])
    p_perm = (np.sum(np.abs(null - mean) >= abs(J - mean)) + 1) / (n_perm + 1)
    return J, z, 2 * stats.norm.sf(abs(z)), p_perm


# ---------------------------------------------------------------- классификация S1 (вложенная CV)
def _fit_predict(Xtr, ytr, Xte, method, k):
    if method == "panel":
        p = stats.mannwhitneyu(Xtr[ytr == 1], Xtr[ytr == 0], axis=0).pvalue
        sel = np.argsort(p)[:k]
        mu, sd = Xtr[:, sel].mean(0), Xtr[:, sel].std(0) + 1e-9
        m = LogisticRegression(C=1.0, max_iter=2000).fit((Xtr[:, sel] - mu) / sd, ytr)
        return m.predict_proba((Xte[:, sel] - mu) / sd)[:, 1], sel
    if method == "plsda":
        mu, sd = Xtr.mean(0), np.sqrt(Xtr.std(0, ddof=1))
        sd[sd == 0] = 1
        m = PLSRegression(k, scale=False).fit((Xtr - mu) / sd, ytr)
        return m.predict((Xte - mu) / sd).ravel(), []
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9           # covariate: логистическая регрессия на всех столбцах
    m = LogisticRegression(C=1.0, max_iter=2000).fit((Xtr - mu) / sd, ytr)
    return m.predict_proba((Xte - mu) / sd)[:, 1], []


def _inner_choice(Xtr, ytr, method, grid, seed):
    if len(grid) == 1:
        return grid[0]
    skf = StratifiedKFold(P["inner_folds"], shuffle=True, random_state=seed)
    best, best_auc = grid[0], -1
    for k in grid:
        pred = np.zeros(len(ytr))
        for a, b in skf.split(Xtr, ytr):
            pred[b], _ = _fit_predict(Xtr[a], ytr[a], Xtr[b], method, k)
        auc = roc_auc_score(ytr, pred)
        if auc > best_auc + 1e-12:
            best, best_auc = k, auc
    return best


def nested_cv(X, y, method, seed):
    """Внешняя CV 5 × 20; внутри каждого обучающего блока — внутренняя 5-fold CV для выбора размера панели
    (method = panel) или числа компонент (plsda); отбор признаков — только на обучающих пробах."""
    grid = P["panel_k_grid"] if method == "panel" else P["plsda_comp_grid"] if method == "plsda" else [0]
    rng = np.random.RandomState(seed)
    aucs, tprs, chosen, ks = [], [], [], []
    fpr_grid = np.linspace(0, 1, 101)
    for _ in range(P["cv_repeats"]):
        skf = StratifiedKFold(P["cv_folds"], shuffle=True, random_state=rng.randint(1e9))
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            k = _inner_choice(X[tr], y[tr], method, grid, rng.randint(1e9))
            pred[te], sel = _fit_predict(X[tr], y[tr], X[te], method, k)
            chosen.extend(sel)
            ks.append(k)
        aucs.append(roc_auc_score(y, pred))
        f, t, _ = roc_curve(y, pred)
        tprs.append(np.interp(fpr_grid, f, t))
    return np.array(aucs), np.array(tprs), chosen, ks


def _perm_one(X, y, method, seed):
    rng = np.random.RandomState(seed)
    yp = rng.permutation(y)
    aucs, *_ = nested_cv(X, yp, method, rng.randint(1e9))
    return aucs.mean()


def perm_test(X, y, method, observed, seed, n_perm):
    seeds = np.random.RandomState(seed).randint(1e9, size=n_perm)
    null = np.array(Parallel(n_jobs=-1)(delayed(_perm_one)(X, y, method, s) for s in seeds))
    return (np.sum(null >= observed) + 1) / (n_perm + 1), null


# ---------------------------------------------------------------- вспомогательные
def find_feature(info, mz, rt, tol_mz=0.003, tol_rt=0.05):
    f = info.index[((info.mz - mz).abs() < tol_mz) & ((info.rt_min - rt).abs() < tol_rt)]
    if len(f) != 1:
        raise SystemExit(f"m/z {mz} RT {rt}: найдено признаков {len(f)} (ожидался 1)")
    return f[0]


def save(fig, path):
    fig.savefig(str(path) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(path) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def fmt_p(p):
    return "" if pd.isna(p) else (f"{p:.2g}" if p < 0.001 else f"{p:.3f}")


def subgroup_sets(meta, study):
    bc = [s for s in study if meta.loc[s, "group"] == "BC"]
    ctrl = [s for s in study if meta.loc[s, "group"] == "Control"]
    def pick(col, val, men=False):
        return [s for s in bc if meta.loc[s, col] == val and (not men or meta.loc[s, "sex"] == "M")]
    return bc, ctrl, pick


def analyses(meta, pick):
    """Варианты модели: (имя, S, case, ref, ковариаты, дополнительный столбец, описание)."""
    stage_ind = meta["stage_group"].map({"T1": 1.0, "Ta": 0.0}).rename("stage_T1")
    return [
        ("S1", S1, pick("stage_group", "T1"), pick("stage_group", "Ta"), P["covariates"], None,
         "T1 vs Ta: подгруппа + пол + возраст + индекс"),
        ("S1_men", S1, pick("stage_group", "T1", True), pick("stage_group", "Ta", True), P["covariates_men"], None,
         "T1 vs Ta, только мужчины: подгруппа + возраст + индекс"),
        ("S2", S2, pick("grade_group", "HG"), pick("grade_group", "LG"), P["covariates"], None,
         "HG vs LG: подгруппа + пол + возраст + индекс"),
        ("S2_stage_adj", S2, [s for s in pick("grade_group", "HG") if meta.loc[s, "stage_group"] in ("Ta", "T1")],
         [s for s in pick("grade_group", "LG") if meta.loc[s, "stage_group"] in ("Ta", "T1")], P["covariates"],
         stage_ind.to_frame(), "HG vs LG, пробы Ta/T1: подгруппа + стадия (T1 vs Ta) + пол + возраст + индекс"),
        ("S2_men", S2, pick("grade_group", "HG", True), pick("grade_group", "LG", True), P["covariates_men"], None,
         "HG vs LG, только мужчины: подгруппа + возраст + индекс"),
    ]


# ---------------------------------------------------------------- описательная таблица
def descriptive(meta, study):
    rows = []
    order = [("Контроль", meta.group == "Control"), ("РМП, все", meta.group == "BC")]
    for v in ["Ta", "T1", "T2plus", "ambiguous", "unknown", "PUNLMP_uncertain"]:
        order.append((f"стадия {v}", (meta.group == "BC") & (meta.stage_group == v)))
    for v in ["LG", "HG", "G2_1973", "not_stated", "PUNLMP_uncertain"]:
        order.append((f"степень {v}", (meta.group == "BC") & (meta.grade_group == v)))
    m = meta.loc[study]
    for lab, mask in order:
        s = m[mask.reindex(m.index).fillna(False)]
        rows.append({"подгруппа": lab, "n": len(s), "мужчины": int((s.sex == "M").sum()), "женщины": int((s.sex == "F").sum()),
                     "мужчины, %": round(100 * (s.sex == "M").mean(), 1) if len(s) else np.nan,
                     "возраст, M ± SD": f"{s.age.mean():.1f} ± {s.age.std():.1f}" if len(s) > 1 else (f"{s.age.mean():.0f}" if len(s) else ""),
                     "возраст, Me [min–max]": f"{s.age.median():.0f} [{s.age.min():.0f}–{s.age.max():.0f}]" if len(s) else "",
                     "индекс, M ± SD": f"{s.contam_index.mean():+.2f} ± {s.contam_index.std():.2f}" if len(s) > 1 else (f"{s.contam_index.mean():+.2f}" if len(s) else ""),
                     "индекс, Me [IQR]": f"{s.contam_index.median():+.2f} [{s.contam_index.quantile(.25):+.2f}; {s.contam_index.quantile(.75):+.2f}]" if len(s) else ""})
    tests = []
    for nm, col, a, b in [("S1", "stage_group", "Ta", "T1"), ("S2", "grade_group", "LG", "HG")]:
        A = m[(m.group == "BC") & (m[col] == a)]
        B = m[(m.group == "BC") & (m[col] == b)]
        tests += [{"сравнение": f"{nm}: {a} vs {b}", "переменная": "возраст", "критерий": "Манн–Уитни",
                   "p": stats.mannwhitneyu(A.age, B.age).pvalue},
                  {"сравнение": f"{nm}: {a} vs {b}", "переменная": "пол", "критерий": "точный тест Фишера",
                   "p": stats.fisher_exact([[(A.sex == "M").sum(), (A.sex == "F").sum()], [(B.sex == "M").sum(), (B.sex == "F").sum()]]).pvalue},
                  {"сравнение": f"{nm}: {a} vs {b}", "переменная": "индекс экзогенного сигнала", "критерий": "Манн–Уитни",
                   "p": stats.mannwhitneyu(A.contam_index, B.contam_index).pvalue}]
    bc = m[m.group == "BC"]
    ct = pd.crosstab(bc.stage_group, bc.grade_group, margins=True, margins_name="всего")
    return pd.DataFrame(rows), pd.DataFrame(tests), ct


# ---------------------------------------------------------------- рисунок: значения по подгруппам
def fig_groups(meta, L, items, path, title, ylab="log2-интенсивность (PQN)"):
    order = ["Control", "Ta", "T1", "LG", "HG"]
    n = len(items)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.9 * nrow), squeeze=False)
    rng = np.random.RandomState(0)
    for ax, (lab, col, note) in zip(axes.ravel(), items):
        v = L[col]
        sets = {"Control": meta.index[meta.group == "Control"],
                "Ta": meta.index[(meta.group == "BC") & (meta.stage_group == "Ta")],
                "T1": meta.index[(meta.group == "BC") & (meta.stage_group == "T1")],
                "LG": meta.index[(meta.group == "BC") & (meta.grade_group == "LG")],
                "HG": meta.index[(meta.group == "BC") & (meta.grade_group == "HG")]}
        data = [v.reindex(sets[k]).dropna().to_numpy() for k in order]
        bp = ax.boxplot(data, positions=range(5), widths=0.55, showfliers=False, patch_artist=True)
        for patch, k in zip(bp["boxes"], order):
            patch.set_facecolor(COL[k])
            patch.set_alpha(0.35)
        for i, (k, d) in enumerate(zip(order, data)):
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(d)), d, s=9, color=COL[k], lw=0, zorder=3)
        t2 = v.reindex(meta.index[(meta.group == "BC") & (meta.stage_group == "T2plus")]).dropna()
        ax.scatter(np.full(len(t2), 5.0), t2, s=36, marker="D", facecolor="white", edgecolor="black", zorder=4, label="T2+")
        ax.axvline(2.5, color="#bbbbbb", lw=0.8, ls=":")
        ax.set_xticks(range(6))
        ax.set_xticklabels([f"{LAB[k]}\n(n={len(d)})" for k, d in zip(order, data)] + [f"T2+\n(n={len(t2)})"], fontsize=8)
        ax.set_title(lab, fontsize=9)
        if note:
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + 0.30 * (hi - lo))       # место для подписи над данными
            ax.text(0.02, 0.98, note, transform=ax.transAxes, va="top", fontsize=7,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        ax.set_ylabel(ylab, fontsize=8)
        ax.tick_params(axis="y", labelsize=8)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    save(fig, path)


def fig_volcano(res, path, title):
    fig, ax = plt.subplots(figsize=(6, 4.6))
    sig = res.q < P["fdr"]
    ax.scatter(res.log2FC[~sig], -np.log10(res.p[~sig]), s=6, c="#bdbdbd", lw=0)
    ax.scatter(res.log2FC[sig], -np.log10(res.p[sig]), s=10, c="#d94801", lw=0, label=f"q < {P['fdr']} (n = {int(sig.sum())})")
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("log2FC (модель с поправками)")
    ax.set_ylabel("−log10 p")
    ax.set_title(title, fontsize=10)
    if sig.any():
        ax.legend(fontsize=8, frameon=False)
    save(fig, path)


def classification(meta, LB, keepB, comp_cat, names, study, pick, dB, log):
    """Классификация S1 (T1 vs Ta): вложенная CV, перестановочный тест, базовая модель «пол + возраст»."""
    case, ref = pick("stage_group", "T1"), pick("stage_group", "Ta")
    idx = case + ref
    y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
    rho = LB.corrwith(meta.loc[study, "contam_index"], method="spearman")
    ex_cat = features.panel_excluded(comp_cat.reindex(keepB).fillna(""))
    pool = [c for c in keepB if not ex_cat[c] and abs(rho[c]) < P["pool_max_abs_rho_index"]]
    log(f"\nКлассификация S1: {len(case)} T1 / {len(ref)} Ta; пул {len(pool)} соединений (без стероидов, экспозиций "
        f"{int(ex_cat.sum())} и |ρ с индексом| ≥ {P['pool_max_abs_rho_index']})")
    X = LB.loc[idx, pool].to_numpy()
    Z = np.column_stack([meta.loc[idx, "sex"].map({"M": 1.0, "F": 0.0}).to_numpy(float), meta.loc[idx, "age"].to_numpy(float)])
    crow, curves, sel_rows = [], {}, None
    for k, (method, Xm, nm) in enumerate([
            ("panel", X, f"панель: top-k по Манн–Уитни внутри обучающего блока + логистическая регрессия; k ∈ {P['panel_k_grid']} — внутренняя CV"),
            ("plsda", X, f"PLS-DA (все соединения пула); число компонент ∈ {P['plsda_comp_grid']} — внутренняя CV"),
            ("covariate", Z, "базовая модель: пол + возраст (логистическая регрессия)")]):
        aucs, tprs, chosen, ks = nested_cv(Xm, y, method, P["seed"] + 10 * k)
        p_perm, null = perm_test(Xm, y, method, aucs.mean(), P["seed"] + 10 * k + 1, P["n_perm"]) if P["n_perm"] else (np.nan, None)
        crow.append({"model": nm, "n_T1": len(case), "n_Ta": len(ref), "n_features": Xm.shape[1],
                     "AUC_mean": aucs.mean(), "AUC_sd": aucs.std(), "AUC_p2.5": np.percentile(aucs, 2.5),
                     "AUC_p97.5": np.percentile(aucs, 97.5), "perm_p": p_perm, "n_perm": P["n_perm"],
                     "inner_choice": "; ".join(f"{v}: {c}" for v, c in pd.Series(ks).value_counts().sort_index().items())})
        curves[method] = (tprs.mean(0), aucs.mean())
        if method == "panel":
            cnt = pd.Series(chosen, dtype=int).value_counts()
            n_f = P["cv_folds"] * P["cv_repeats"]
            sel_rows = pd.DataFrame({"compound": [pool[i] for i in cnt.index], "name": [names.get(pool[i], pool[i]) for i in cnt.index],
                                     "n_selected": cnt.to_numpy(), "selection_frequency": cnt.to_numpy() / n_f})
        log(f"  {nm}: AUC {aucs.mean():.3f} (2,5–97,5 %: {np.percentile(aucs, 2.5):.3f}–{np.percentile(aucs, 97.5):.3f}), "
            f"перестановочный p = {fmt_p(p_perm)}")
    pd.DataFrame(crow).to_csv(dB / "classification_S1_nested_cv.csv", index=False, encoding="utf-8-sig")
    sel_rows.to_csv(dB / "classification_S1_panel_selection_frequency.csv", index=False, encoding="utf-8-sig")
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    fpr = np.linspace(0, 1, 101)
    for (m, (tpr, auc)), lab, c in zip(curves.items(), ["панель (top-k)", "PLS-DA", "пол + возраст"], ["#08519c", "#d94801", "#636363"]):
        ax.plot(fpr, tpr, color=c, label=f"{lab}: AUC {auc:.2f}")
    ax.plot([0, 1], [0, 1], ls=":", c="k", lw=0.7)
    ax.set_xlabel("1 − специфичность")
    ax.set_ylabel("чувствительность")
    ax.set_title(f"S1: T1 vs Ta, вложенная CV {P['cv_folds']} × {P['cv_repeats']}", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    save(fig, dB / "Figure_classification_S1_ROC")


# ---------------------------------------------------------------- основная программа
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical", required=True, help="s107_clinical_groups.xlsx")
    ap.add_argument("--out", default=str(HERE.parent / "results_s107_clinical"))
    ap.add_argument("--blood_list", default=None, help="согласованный список соединений «кровяной» оценки (csv, столбец compound)")
    ap.add_argument("--n_perm", type=int, default=P["n_perm"])
    ap.add_argument("--skip_classification", action="store_true", help="не пересчитывать классификацию S1 (файлы остаются)")
    a = ap.parse_args()
    a.clinical = str((CWD / a.clinical).resolve())
    a.out = str((CWD / a.out).resolve())
    a.blood_list = str((CWD / a.blood_list).resolve()) if a.blood_list else None
    P["n_perm"] = a.n_perm
    out = Path(a.out)
    dA, dB = out / "steroids", out / "clinical"
    for d in (dA, dB):
        d.mkdir(parents=True, exist_ok=True)
    log_lines = []

    def log(s=""):
        print(s)
        log_lines.append(str(s))

    # ------------------------------------------------------------ матрица основного анализа и клинические группы
    M = prepare_matrix.build(log)
    meta, L, info, groups, names, comp_cat = M["meta"], M["Lrep"], M["info"], M["groups"], M["names"], M["comp_cat"]
    Lall = M["L"]
    study = list(L.index)
    clin = pd.read_excel(a.clinical)
    clin["sample"] = clin.sample_no.astype(str)
    clin = clin.set_index("sample")
    meta["code"] = pd.read_csv(C.METADATA, dtype={"sample": str}).set_index("sample")["code"].reindex(meta.index)
    rec = meta.loc[study, ["group", "sex", "age", "code"]].join(clin[["code", "group", "stage_group", "grade_group", "CIS", "sex", "age"]],
                                                                   rsuffix="_clinical")
    rec["sex_clinical"] = rec.sex_clinical.map({"муж": "M", "жен": "F"})
    rec["code_match"] = rec.code.astype(str).str.strip() == rec.code_clinical.astype(str).str.strip()
    rec["group_match"] = rec.group == rec.group_clinical
    rec["age_match"] = rec.age == rec.age_clinical
    rec["sex_match"] = rec.sex == rec.sex_clinical
    if not (rec.code_match.all() and rec.group_match.all() and rec.age_match.all()):
        raise SystemExit("Клинический файл не совпадает с метаданными s107 по шифру, группе или возрасту")
    log(f"\nСопоставление по sample_no: {len(rec)} проб основного анализа; шифр (code), группа и возраст совпадают у всех.")
    mism = rec[~rec.sex_match]
    if len(mism):
        log("Пол расходится у проб " + ", ".join(f"{i} (файл: {r.sex_clinical}, основной анализ s107 v2: {r.sex})" for i, r in mism.iterrows())
            + " — используется пол основного анализа (исправлен в версии v2 после сверки)")
    meta["stage_group"] = clin.stage_group.reindex(meta.index)
    meta["grade_group"] = clin.grade_group.reindex(meta.index)
    meta["CIS"] = clin.CIS.reindex(meta.index)
    bc, ctrl, pick = subgroup_sets(meta, study)
    variants = analyses(meta, pick)
    log("Варианты модели: " + "; ".join(f"{v[0]} — {len(v[2])}/{len(v[3])}" for v in variants))

    # ------------------------------------------------------------ описательная таблица (в обе папки)
    desc, desc_tests, ct = descriptive(meta, study)
    for d in (dA, dB):
        desc.to_csv(d / "descriptive_subgroups.csv", index=False, encoding="utf-8-sig")
        desc_tests.to_csv(d / "descriptive_tests.csv", index=False, encoding="utf-8-sig")
        ct.to_csv(d / "stage_by_grade.csv", encoding="utf-8-sig")
    rec.to_csv(out / "sample_reconciliation.csv", encoding="utf-8-sig")
    log("\nОписательная таблица:\n" + desc.to_string(index=False))
    log(desc_tests.assign(p=desc_tests.p.map(fmt_p)).to_string(index=False))

    # ------------------------------------------------------------ блок А: стероиды
    rng = np.random.RandomState(P["seed"])
    ster = []
    for nm, mz, rt in STEROIDS:
        f = find_feature(info, mz, rt)
        rep = groups.loc[f, "representative"]
        members = list(groups.index[groups.compound == groups.loc[f, "compound"]])
        ster.append(dict(name=nm, feature=f, rep=rep, members=members))
    ster_reps = [s["rep"] for s in ster]
    ions = pd.DataFrame([dict(glucuronide=s["name"], feature=m, mz=info.loc[m, "mz"], rt_min=info.loc[m, "rt_min"],
                              is_representative=m == s["rep"], ion_relation=groups.loc[m, "ion_relations"],
                              msdial_annotation=info.loc[m, "annotation"])
                         for s in ster for m in s["members"]])
    ions.to_csv(dA / "steroid_ions.csv", index=False, encoding="utf-8-sig")
    rowsA, rows_ion = [], []
    for vname, S, case, ref, cov, extra, desc_v in variants:
        r, used = ols(L[ster_reps], meta, case, ref, cov, extra)
        r["q_BH_3_glucuronides"] = bh(r.p)
        for s in ster:
            x = r.loc[s["rep"]]
            auc, lo, hi = auc_ci(L.loc[case, s["rep"]].to_numpy(), L.loc[ref, s["rep"]].to_numpy(), rng, P["auc_bootstrap"])
            rowsA.append({"compound": s["name"], "analysis": vname, "comparison": S["label"], "model": desc_v,
                          "covariates": "+".join(used), "n_case": int(x.n_case), "n_ref": int(x.n_ref),
                          "log2FC": x.log2FC, "CI95_low": x.CI_low, "CI95_high": x.CI_high, "p": x.p,
                          "q_BH_3_glucuronides": x.q_BH_3_glucuronides, "AUC": auc, "AUC_CI95_low": lo, "AUC_CI95_high": hi})
        # ионы (агликон, аддукты) — проверка: та же модель для каждого иона отдельно
        ion_cols = [m for s in ster for m in s["members"]]
        ri, _ = ols(Lall.loc[study, ion_cols], meta, case, ref, cov, extra)
        for m in ion_cols:
            g = next(s["name"] for s in ster if m in s["members"])
            rows_ion.append({"glucuronide": g, "analysis": vname, "feature": m, "mz": info.loc[m, "mz"],
                             "rt_min": info.loc[m, "rt_min"], "ion_relation": groups.loc[m, "ion_relations"],
                             "representative": m == next(s["rep"] for s in ster if m in s["members"]),
                             "log2FC": ri.loc[m, "log2FC"], "CI95_low": ri.loc[m, "CI_low"], "CI95_high": ri.loc[m, "CI_high"],
                             "p": ri.loc[m, "p"]})
    tabA = pd.DataFrame(rowsA)
    tabA.to_csv(dA / "steroids_S1_S2_models.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(rows_ion).to_csv(dA / "steroids_ions_check.csv", index=False, encoding="utf-8-sig")
    log("\nБлок А — глюкурониды стероидов:\n" + tabA[["compound", "analysis", "n_case", "n_ref", "log2FC", "CI95_low",
                                                       "CI95_high", "p", "q_BH_3_glucuronides", "AUC"]].round(3).to_string(index=False))
    # тренд контроль → Ta → T1 (Йонкхир), все и только мужчины
    jt = []
    for men in (False, True):
        sx = (lambda s: meta.loc[s, "sex"] == "M") if men else (lambda s: True)
        grp = [[s for s in ctrl if sx(s)], [s for s in pick("stage_group", "Ta") if sx(s)], [s for s in pick("stage_group", "T1") if sx(s)]]
        for s in ster:
            J, z, p_as, p_perm = jonckheere([L.loc[g, s["rep"]].to_numpy() for g in grp], rng, P["jt_permutations"])
            jt.append({"compound": s["name"], "subset": "только мужчины" if men else "все", "order": "Контроль < Ta < T1",
                       "n": "/".join(str(len(g)) for g in grp), "JT": J, "z": z, "p_asymptotic": p_as, "p_permutation": p_perm,
                       "direction": "возрастающий" if z > 0 else "убывающий"})
    jt = pd.DataFrame(jt)
    jt.to_csv(dA / "steroids_trend_jonckheere.csv", index=False, encoding="utf-8-sig")
    log("\nТренд Контроль → Ta → T1 (Йонкхир):\n" + jt.drop(columns="order").round(4).to_string(index=False))
    # описательно: подгруппы вне S1/S2
    drows = []
    for col, vals in [("stage_group", ["Ta", "T1", "T2plus", "ambiguous", "unknown", "PUNLMP_uncertain"]),
                      ("grade_group", ["LG", "HG", "G2_1973", "not_stated", "PUNLMP_uncertain"])]:
        for v in vals:
            ids = pick(col, v)
            for s in ster:
                x = L.loc[ids, s["rep"]]
                drows.append({"compound": s["name"], "variable": col, "subgroup": v, "n": len(ids),
                              "median_log2": x.median(), "q25": x.quantile(.25), "q75": x.quantile(.75), "min": x.min(), "max": x.max()})
    for s in ster:
        x = L.loc[ctrl, s["rep"]]
        drows.append({"compound": s["name"], "variable": "group", "subgroup": "Control", "n": len(ctrl), "median_log2": x.median(),
                      "q25": x.quantile(.25), "q75": x.quantile(.75), "min": x.min(), "max": x.max()})
    pd.DataFrame(drows).to_csv(dA / "steroids_descriptive_by_subgroup.csv", index=False, encoding="utf-8-sig")
    items = []
    for s in ster:
        r1 = tabA[(tabA.compound == s["name"]) & (tabA.analysis == "S1")].iloc[0]
        r2 = tabA[(tabA.compound == s["name"]) & (tabA.analysis == "S2")].iloc[0]
        j = jt[(jt.compound == s["name"]) & (jt.subset == "все")].iloc[0]
        items.append((s["name"], s["rep"], f"S1: log2FC {r1.log2FC:+.2f}, p = {fmt_p(r1.p)}\nS2: log2FC {r2.log2FC:+.2f}, p = {fmt_p(r2.p)}\n"
                                          f"тренд К→Ta→T1: p = {fmt_p(j.p_permutation)}"))
    fig_groups(meta.loc[study], L, items, dA / "Figure_steroid_glucuronides_by_subgroup",
               "Глюкурониды стероидов: контроль, стадия (Ta, T1) и степень злокачественности (LG, HG); T2+ — отдельные точки")

    # ------------------------------------------------------------ блок Б: исключение стероидов
    member_ann = info.loc[groups.index, "annotation"].astype(str)
    hit = member_ann.str.contains(STEROID_RE)
    excl_comp = set(groups.loc[hit[hit].index, "compound"]) | {groups.loc[s["feature"], "compound"] for s in ster}
    ex_rows = []
    for c in sorted(excl_comp, key=str):
        mem = groups.index[groups.compound == c]
        rep = groups.loc[mem[0], "representative"]
        why = "глюкуронид стероида (блок А) и его ионы" if any(c == groups.loc[s["feature"], "compound"] for s in ster) \
            else "аннотация стероида / стероидного конъюгата у признака соединения"
        for m in mem:
            an = str(info.loc[m, "annotation"])
            cls = "" if not hit[m] else ("желчная кислота / стерол (стероидный скелет)" if BILE_RE.search(an)
                                         else "экдистероид" if re.search(r"cyasterone|turkesterone|ecdys", an, re.I)
                                         else "стероидный гормон / его метаболит или конъюгат")
            ex_rows.append({"compound_representative": rep, "feature": m, "mz": info.loc[m, "mz"], "rt_min": info.loc[m, "rt_min"],
                            "msdial_annotation": an, "steroid_annotation_match": bool(hit[m]), "steroid_class": cls,
                            "is_representative": m == rep, "reason": why})
    excl = pd.DataFrame(ex_rows)
    excl.to_csv(out / "excluded_steroid_features_for_block_B.csv", index=False, encoding="utf-8-sig")
    keepB = [c for c in L.columns if groups.loc[c, "compound"] not in excl_comp]
    log(f"\nБлок Б: исключено соединений {len(excl_comp)} (признаков {len(excl)}); осталось {len(keepB)} из {L.shape[1]}")

    # ------------------------------------------------------------ блок Б: ненаправленное сравнение S1 и S2
    LB = L[keepB]
    ann = pd.DataFrame({"name": names.reindex(keepB), "mz": info.loc[keepB, "mz"], "rt_min": info.loc[keepB, "rt_min"],
                        "msdial_annotation": info.loc[keepB, "annotation"], "annotation_confidence": info.loc[keepB, "annotation_conf"],
                        "compound_category": comp_cat.reindex(keepB), "n_ions": groups.loc[keepB, "n_ions"]})
    summ = []
    full = {}
    for vname, S, case, ref, cov, extra, desc_v in variants:
        r, used = ols(LB, meta, case, ref, cov, extra)
        r["q"] = bh(r.p)
        r["AUC"] = auc_fast(LB, case, ref)
        full[vname] = r
        summ.append({"analysis": vname, "model": desc_v, "covariates": "+".join(used), "n_case": int(r.n_case.iloc[0]),
                     "n_ref": int(r.n_ref.iloc[0]), "compounds": len(r), "q<0.05": int((r.q < P["fdr"]).sum()),
                     "p<0.05": int((r.p < 0.05).sum()), "min_q": r.q.min()})
    summ = pd.DataFrame(summ)
    summ.to_csv(dB / "untargeted_summary.csv", index=False, encoding="utf-8-sig")
    log("\nНенаправленное сравнение (BH по всем соединениям блока Б):\n" + summ.drop(columns="model").to_string(index=False))
    for Sname, prim, others in [("S1", "S1", ["S1_men"]), ("S2", "S2", ["S2_stage_adj", "S2_men"])]:
        r = full[prim]
        tab = ann.join(r[["n_case", "n_ref", "log2FC", "CI_low", "CI_high", "p", "q", "AUC"]])
        for o in others:
            tab = tab.join(full[o][["log2FC", "p", "q"]].add_suffix(f"_{o}"))
        tab = tab.sort_values("p")
        tab.to_csv(dB / f"untargeted_{Sname}_all_compounds.csv", encoding="utf-8-sig", index_label="compound")
        sig = tab[tab.q < P["fdr"]]
        sel = sig if len(sig) else tab.head(P["top_n_if_no_q"])
        sel.assign(selection=f"q < {P['fdr']}" if len(sig) else f"нет q < {P['fdr']}: первые {P['top_n_if_no_q']} по p") \
            .to_csv(dB / f"untargeted_{Sname}_top.csv", encoding="utf-8-sig", index_label="compound")
        fig_volcano(r, dB / f"Figure_volcano_{Sname}", f"{S1['label'] if Sname == 'S1' else S2['label']} (РМП; без стероидов)")
        log(f"{Sname}: q < {P['fdr']} — {len(sig)}; в таблицу: {'соединения с q < 0.05' if len(sig) else 'первые 50 по p'}")

    # ------------------------------------------------------------ блок Б: классификация S1 (вложенная CV)
    if not a.skip_classification:
        classification(meta, LB, keepB, comp_cat, names, study, pick, dB, log)


    # ------------------------------------------------------------ блок Б: «кровяная» оценка
    if a.blood_list:
        bl = pd.read_csv(a.blood_list)
        comps = [c for c in bl.compound if c in L.columns]
        miss = set(bl.compound) - set(comps)
        if miss:
            raise SystemExit(f"Соединения списка не найдены в матрице: {miss}")
        if set(comps) & {r for s in ster for r in s["members"]} or any(groups.loc[c, "compound"] in excl_comp for c in comps):
            raise SystemExit("В списке «кровяной» оценки есть стероиды")
        Zs = (L[comps] - L[comps].mean()) / L[comps].std()
        score = Zs.mean(axis=1).rename("blood_score")
        Ls = score.to_frame()
        rowsS = []
        for vname, S, case, ref, cov, extra, desc_v in variants:
            r, used = ols(Ls, meta, case, ref, cov, extra)
            auc, lo, hi = auc_ci(score[case].to_numpy(), score[ref].to_numpy(), rng, P["auc_bootstrap"])
            x = r.iloc[0]
            rowsS.append({"analysis": vname, "comparison": S["label"], "model": desc_v, "covariates": "+".join(used),
                          "n_case": int(x.n_case), "n_ref": int(x.n_ref), "difference_SD_units": x.log2FC,
                          "CI95_low": x.CI_low, "CI95_high": x.CI_high, "p": x.p, "AUC": auc, "AUC_CI95_low": lo, "AUC_CI95_high": hi})
        tabS = pd.DataFrame(rowsS)
        tabS.to_csv(dB / "blood_score_S1_S2_models.csv", index=False, encoding="utf-8-sig")
        jts = []
        for men in (False, True):
            sx = (lambda s: meta.loc[s, "sex"] == "M") if men else (lambda s: True)
            grp = [[s for s in ctrl if sx(s)], [s for s in pick("stage_group", "Ta") if sx(s)], [s for s in pick("stage_group", "T1") if sx(s)]]
            J, z, p_as, p_perm = jonckheere([score[g].to_numpy() for g in grp], rng, P["jt_permutations"])
            jts.append({"subset": "только мужчины" if men else "все", "order": "Контроль < Ta < T1", "n": "/".join(str(len(g)) for g in grp),
                        "JT": J, "z": z, "p_asymptotic": p_as, "p_permutation": p_perm})
        pd.DataFrame(jts).to_csv(dB / "blood_score_trend_jonckheere.csv", index=False, encoding="utf-8-sig")
        bl.assign(name=[names.get(c, c) for c in bl.compound]).to_csv(dB / "blood_score_compounds.csv", index=False, encoding="utf-8-sig")
        r1, r2 = tabS[tabS.analysis == "S1"].iloc[0], tabS[tabS.analysis == "S2"].iloc[0]
        fig_groups(meta.loc[study], Ls, [("«Кровяная» оценка", "blood_score",
                                          f"S1: {r1.difference_SD_units:+.2f} SD, p = {fmt_p(r1.p)}\nS2: {r2.difference_SD_units:+.2f} SD, p = {fmt_p(r2.p)}\n"
                                          f"тренд К→Ta→T1: p = {fmt_p(jts[0]['p_permutation'])}")],
                   dB / "Figure_blood_score_by_subgroup", f"«Кровяная» оценка ({len(comps)} соединений): среднее z-оценок",
                   ylab="средняя z-оценка")
        log("\n«Кровяная» оценка:\n" + tabS[["analysis", "n_case", "n_ref", "difference_SD_units", "CI95_low", "CI95_high", "p", "AUC"]].round(3).to_string(index=False))
    else:
        log("\n«Кровяная» оценка не рассчитана: список соединений ещё не согласован (--blood_list)")

    # ------------------------------------------------------------ параметры и журнал
    params = dict(P, run_at=datetime.now().isoformat(timespec="seconds"), pipeline_commit=M["pipeline_commit"],
                  matrix_check=M["check"].to_dict(), clinical_file=Path(a.clinical).name,
                  variants={v[0]: dict(description=v[6], n_case=len(v[2]), n_ref=len(v[3])) for v in variants},
                  steroid_regex=STEROID_RE.pattern, blood_list=Path(a.blood_list).name if a.blood_list else None,
                  sex_source="метаданные основного анализа s107 (v2)")
    for d in (out, dA, dB):
        fn = "parameters.json"
        pp = dict(params)
        if d == dB:
            pp.pop("steroid_regex", None)
        json.dump(pp, open(d / fn, "w", encoding="utf-8"), ensure_ascii=False, indent=2, default=str)
    (out / "log.txt").write_text("\n".join(log_lines), encoding="utf-8")


if __name__ == "__main__":
    main()

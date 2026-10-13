#!/usr/bin/env python3
"""Обработка данных ненаправленного метаболомного профилирования (ВЭЖХ-МС).

Шаги:
  1. Загрузка таблицы признаков и последовательности проб, сопоставление метаданных.
  2. Порядок инъекций (образцы + QC).
  3. Фильтрация признаков по детекции в QC и в образцах.
  4. Коррекция дрейфа сигнала по QC (QC-RLSC; метод выбирается по leave-one-QC-out).
  5. Фильтрация по RSD в QC.
  6. Нормализация на разведение мочи (PQN).
  7. Импутация пропусков, log2-преобразование.
  8. Многомерный анализ (PCA, PLS-DA с кросс-валидацией и перестановочным тестом).
  9. Одномерная статистика (Манн-Уитни + BH FDR, AUC, линейная модель с поправкой на пол).

Запуск:  python process.py --data data/s266_table_for_calc.csv --seq data/sequence.xlsx --out results
"""
import argparse
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# ----------------------------------------------------------------------------- параметры
QC_EVERY = 20               # QC вводился в начале и после каждых N образцов (QC1 ... QC6)
MIN_QC_DETECT = 0.8         # признак должен детектироваться в >= 80 % QC
MIN_GROUP_DETECT = 0.5      # ... и в >= 50 % образцов хотя бы одной группы (Cancer / Control)
MAX_QC_RSD = 0.30           # RSD в QC после коррекции дрейфа
DRIFT_METHODS = ("none", "linear", "lowess")
N_PERM = 200                # перестановок для PLS-DA
SEED = 1

GROUP_COLORS = {"QC": "#7f7f7f", "Control": "#2b7bba", "BC": "#d7301f", "PC": "#fc8d59",
                "KC": "#7b3294", "PanC": "#1a9850", "BC+PC": "#8c510a"}


# ----------------------------------------------------------------------------- 1. загрузка
def diag_from_text(t):
    t = str(t).strip().lower()
    if t.startswith("кг"):
        return "control"
    has_bc = "рмп" in t
    has_pc = "рпж" in t or "рак предст" in t
    if has_bc and has_pc:
        return "BC+PC"
    if has_bc:
        return "BC"
    if has_pc:
        return "PC"
    if "почк" in t:
        return "KC"
    if "поджелуд" in t:
        return "PanC"
    return "?"


def load(data_path, seq_path):
    raw = pd.read_csv(data_path, encoding="utf-8-sig", dtype={"Sample": str})
    meta_cols = ["Sample", "Group", "Diag"]
    feat_cols = [c for c in raw.columns if c not in meta_cols]

    # Пробы с двойным диагнозом (42, 71) продублированы в таблице с идентичными
    # интенсивностями -> оставляем одну инъекцию.
    dup = raw[raw.duplicated("Sample", keep=False)]
    merged_diag = dup.groupby("Sample")["Diag"].apply(lambda s: "+".join(sorted(s)))
    raw = raw.drop_duplicates("Sample").copy()
    raw.loc[raw.Sample.isin(merged_diag.index), "Diag"] = raw.Sample.map(merged_diag)

    X = raw.set_index("Sample")[feat_cols].astype(float)
    X = X.replace(0, np.nan)

    seq = pd.read_excel(seq_path, header=0)
    seq = seq.rename(columns={seq.columns[0]: "order", "Пол": "sex", "Возраст": "age",
                              "Диагноз / группа": "diag_text"})
    seq = seq.dropna(subset=["order"])
    seq["Sample"] = seq["order"].astype(int).astype(str)
    seq["sex"] = seq["sex"].astype(str).str.strip().replace({"nan": np.nan})
    seq["diag_seq"] = seq["diag_text"].map(diag_from_text)
    seq = seq[["Sample", "sex", "age", "diag_text", "diag_seq"]]

    meta = raw[meta_cols].set_index("Sample").join(seq.set_index("Sample"))
    meta["class"] = meta["Diag"].replace({"control": "Control"})
    meta.loc[meta["class"] == "BC+PC", "class"] = "BC+PC"
    meta["is_qc"] = meta["class"] == "QC"

    # проверка соответствия диагнозов в двух таблицах
    s = meta[~meta.is_qc]
    mism = s[(s["diag_seq"] != s["Diag"]) & ~((s["Diag"] == "BC+PC") & (s["diag_seq"].isin(["BC+PC", "PC"])))]
    feat = parse_features(feat_cols)
    return X, meta, feat, mism


def parse_features(cols):
    rows = []
    for c in cols:
        rt, mz, name = c.split("_", 2)
        low = name.startswith("low score:")
        clean = re.sub(r"^low score:\s*", "", name)
        rows.append(dict(feature=c, rt_min=float(rt), mz=float(mz), annotation=clean,
                         annotation_conf=("unknown" if clean == "Unknown" else ("low" if low else "high"))))
    return pd.DataFrame(rows).set_index("feature")


# ----------------------------------------------------------------------------- 2. порядок инъекций
def injection_order(meta):
    """QC1 – перед образцом 1, QC(k+1) – после образца k*QC_EVERY."""
    order = {}
    pos = 0
    qcs = sorted([s for s in meta.index if s.startswith("QC")], key=lambda s: int(s[2:]))
    samples = sorted([s for s in meta.index if not s.startswith("QC")], key=int)
    qi = iter(qcs)
    order[next(qi)] = pos
    for i, s in enumerate(samples, 1):
        pos += 1
        order[s] = pos
        if i % QC_EVERY == 0:
            q = next(qi, None)
            if q:
                pos += 1
                order[q] = pos
    return pd.Series(order, name="inj_order")


# ----------------------------------------------------------------------------- 4. дрейф
def _fit_predict(xq, yq, xall, method):
    """yq – log-интенсивности QC (без NaN). Возвращает прогноз log-тренда на xall."""
    if method == "none" or len(yq) < 3:
        return np.full(len(xall), np.nanmedian(yq) if len(yq) else np.nan)
    if method == "linear":
        b, a = np.polyfit(xq, yq, 1)
        return a + b * xall
    if method == "lowess":
        fit = lowess(yq, xq, frac=1.0, it=0, delta=0, xvals=np.asarray(xall, float))
        # за пределами диапазона QC – константа по крайним точкам
        lo, hi = xq.min(), xq.max()
        fit = np.where(xall < lo, fit[np.argmin(np.abs(xall - lo))], fit)
        fit = np.where(xall > hi, fit[np.argmin(np.abs(xall - hi))], fit)
        return fit
    raise ValueError(method)


def loo_drift_error(X, order, qc_idx, method):
    """Leave-one-QC-out: ошибка прогноза внутренних QC (median |log2 err| по признакам)."""
    L = np.log2(X.loc[qc_idx].to_numpy())
    xq_all = order.loc[qc_idx].to_numpy(float)
    inner = np.argsort(xq_all)[1:-1]  # крайние QC не прогнозируем (экстраполяция)
    errs = []
    for j in range(L.shape[1]):
        y = L[:, j]
        e = []
        for k in inner:
            m = np.ones(len(y), bool)
            m[k] = False
            m &= ~np.isnan(y)
            if np.isnan(y[k]) or m.sum() < 4:
                continue
            p = _fit_predict(xq_all[m], y[m], np.array([xq_all[k]]), method)[0]
            e.append(abs(y[k] - p))
        if e:
            errs.append(np.mean(e))
    return float(np.median(errs))


def drift_correct(X, order, qc_idx, method):
    """QC-RLSC: деление на тренд, построенный по QC, с возвратом к медиане QC."""
    Xc = X.copy()
    xq = order.loc[qc_idx].to_numpy(float)
    xall = order.loc[X.index].to_numpy(float)
    L = np.log2(X.to_numpy())
    Lq = np.log2(X.loc[qc_idx].to_numpy())
    out = L.copy()
    for j in range(L.shape[1]):
        y = Lq[:, j]
        m = ~np.isnan(y)
        if m.sum() < 4 or method == "none":
            continue
        trend = _fit_predict(xq[m], y[m], xall, method)
        out[:, j] = L[:, j] - trend + np.median(y[m])
    Xc[:] = 2 ** out
    return Xc


def rsd(df):
    return df.std(ddof=1) / df.mean()


# ----------------------------------------------------------------------------- 6. PQN
def pqn(X, ref_idx):
    ref = X.loc[ref_idx].median()
    q = X.div(ref, axis=1)
    factors = q.median(axis=1, skipna=True)
    return X.div(factors, axis=0), factors


# ----------------------------------------------------------------------------- 8–9. статистика
def pareto(L):
    c = L - L.mean()
    return c / np.sqrt(L.std(ddof=1))


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    if m.any():
        q[m] = multipletests(p[m], method="fdr_bh")[1]
    return q


def compare(L, meta, case, ctrl, feat, label):
    """L – log2-данные (образцы x признаки)."""
    a = L.loc[case]
    b = L.loc[ctrl]
    U = stats.mannwhitneyu(a, b, axis=0, alternative="two-sided")
    p_mw = U.pvalue
    auc = U.statistic / (len(a) * len(b))
    lfc = (a.mean() - b.mean()).to_numpy()
    lfc_med = (a.median() - b.median()).to_numpy()

    # линейная модель log2 ~ группа + пол (пол – главный конфаундер)
    idx = [s for s in list(case) + list(ctrl) if pd.notna(meta.loc[s, "sex"])]
    y = L.loc[idx].to_numpy()
    grp = np.array([1.0 if s in set(case) else 0.0 for s in idx])
    sex = (meta.loc[idx, "sex"] == "м").astype(float).to_numpy()
    p_adj = np.full(L.shape[1], np.nan)
    b_adj = np.full(L.shape[1], np.nan)
    if sex.std() > 0 and len(set(grp)) == 2:
        D = sm.add_constant(np.column_stack([grp, sex]))
        XtX_inv = np.linalg.inv(D.T @ D)
        beta = XtX_inv @ D.T @ y
        resid = y - D @ beta
        dof = len(idx) - D.shape[1]
        s2 = (resid ** 2).sum(0) / dof
        se = np.sqrt(s2 * XtX_inv[1, 1])
        t = beta[1] / se
        p_adj = 2 * stats.t.sf(np.abs(t), dof)
        b_adj = beta[1]

    res = pd.DataFrame({
        "log2FC_mean": lfc, "log2FC_median": lfc_med, "AUC": auc,
        "p_MW": p_mw, "q_MW": bh(p_mw),
        "log2FC_sex_adj": b_adj, "p_sex_adj": p_adj, "q_sex_adj": bh(p_adj),
    }, index=L.columns)
    res = feat.join(res, how="right")
    res.insert(0, "comparison", label)
    return res.sort_values("p_MW")


def vip_scores(pls):
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_
    p, h = w.shape
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total = s.sum()
    wn = (w / np.linalg.norm(w, axis=0)) ** 2
    return np.sqrt(p * (wn @ s).ravel() / total)


def plsda_cv(Z, y, ncomp=2, n_rep=5, seed=SEED):
    """Повторная стратифицированная 5-fold CV, возвращает средний AUC и Q2."""
    rng = np.random.RandomState(seed)
    aucs, q2s = [], []
    for r in range(n_rep):
        skf = StratifiedKFold(5, shuffle=True, random_state=rng.randint(1e6))
        pred = np.zeros(len(y))
        for tr, te in skf.split(Z, y):
            Ztr = Z[tr]
            mu, sd = Ztr.mean(0), np.sqrt(Ztr.std(0, ddof=1))
            sd[sd == 0] = 1
            m = PLSRegression(ncomp, scale=False).fit((Ztr - mu) / sd, y[tr])
            pred[te] = m.predict((Z[te] - mu) / sd).ravel()
        aucs.append(roc_auc_score(y, pred))
        q2s.append(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return float(np.mean(aucs)), float(np.mean(q2s))


# ----------------------------------------------------------------------------- графики
def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_pca(scores, var, meta, order, path, title):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for cls, col in GROUP_COLORS.items():
        m = meta.loc[scores.index, "class"] == cls
        if m.any():
            ax[0].scatter(scores.loc[m, 0], scores.loc[m, 1], c=col, label=f"{cls} (n={m.sum()})",
                          s=60 if cls == "QC" else 30, edgecolor="k", lw=0.4,
                          marker="s" if cls == "QC" else "o")
    ax[0].set_xlabel(f"PC1 ({var[0]:.1%})")
    ax[0].set_ylabel(f"PC2 ({var[1]:.1%})")
    ax[0].legend(fontsize=8)
    ax[0].set_title(title)
    sc = ax[1].scatter(scores[0], scores[1], c=order.loc[scores.index], cmap="viridis", s=30,
                       edgecolor="k", lw=0.4)
    plt.colorbar(sc, ax=ax[1], label="порядок инъекции")
    ax[1].set_xlabel(f"PC1 ({var[0]:.1%})")
    ax[1].set_ylabel(f"PC2 ({var[1]:.1%})")
    ax[1].set_title("Цвет = порядок инъекции")
    savefig(fig, path)


def plot_volcano(res, path, title):
    fig, ax = plt.subplots(figsize=(7, 6))
    x = res["log2FC_mean"]
    y = -np.log10(res["p_MW"])
    sig = (res["q_MW"] < 0.05) & (x.abs() >= 1)
    ax.scatter(x[~sig], y[~sig], s=5, c="#bbbbbb", lw=0)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=8, c="#d7301f", lw=0, label="↑ q<0.05, |log2FC|≥1")
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=8, c="#2b7bba", lw=0, label="↓ q<0.05, |log2FC|≥1")
    top = res[sig & (res.annotation_conf == "high")].head(10)
    for f, r in top.iterrows():
        ax.annotate(r.annotation[:30], (r.log2FC_mean, -np.log10(r.p_MW)), fontsize=6)
    ax.axvline(-1, ls=":", c="k", lw=0.6)
    ax.axvline(1, ls=":", c="k", lw=0.6)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p (Манн-Уитни)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    savefig(fig, path)


def plot_box(L, meta, feats, path, groups):
    n = len(feats)
    if n == 0:
        return
    cols = 4
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.3 * rows), squeeze=False)
    for ax, (f, title) in zip(axes.ravel(), feats):
        data = [L.loc[meta.index[meta["class"] == g].intersection(L.index), f] for g in groups]
        ax.boxplot(data, tick_labels=groups, showfliers=False)
        for i, d in enumerate(data, 1):
            ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=8, c=GROUP_COLORS.get(groups[i - 1]), alpha=0.7)
        ax.set_title(title, fontsize=7)
        ax.set_ylabel("log2 int.", fontsize=7)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    savefig(fig, path)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/s266_table_for_calc.csv")
    ap.add_argument("--seq", default="data/sequence.xlsx")
    ap.add_argument("--out", default="results")
    ap.add_argument("--drift", choices=("auto",) + DRIFT_METHODS, default="auto")
    args = ap.parse_args()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    X, meta, feat, mism = load(args.data, args.seq)
    say(f"Инъекций: {len(X)} (QC: {meta.is_qc.sum()}, образцов: {(~meta.is_qc).sum()}), признаков: {X.shape[1]}")
    say(f"Классы: {meta['class'].value_counts().to_dict()}")
    say(f"Несоответствий диагнозов CSV vs последовательность: {len(mism)}")
    if len(mism):
        say(mism[["Diag", "diag_seq", "diag_text"]].to_string())

    order = injection_order(meta)
    meta["inj_order"] = order
    qc = meta.index[meta.is_qc]
    smp = meta.index[~meta.is_qc]
    say("Позиции QC в последовательности: " + ", ".join(f"{q}={order[q]}" for q in qc))

    # --- QC/TIC обзор
    tic = X.sum(axis=1)
    fig, ax = plt.subplots(figsize=(11, 4))
    for cls, col in GROUP_COLORS.items():
        m = meta["class"] == cls
        ax.scatter(order[m.index[m]], tic[m] / 1e9, c=col, label=cls, s=50 if cls == "QC" else 20,
                   marker="s" if cls == "QC" else "o", edgecolor="k", lw=0.3)
    ax.set_yscale("log")
    ax.set_xlabel("порядок инъекции")
    ax.set_ylabel("суммарная интенсивность, ×1e9")
    ax.legend(fontsize=7, ncol=7)
    ax.set_title("Суммарный сигнал по порядку инъекции (до обработки)")
    savefig(fig, out / "figures/01_tic_by_order.png")

    # --- 3. фильтр по детекции
    n0 = X.shape[1]
    det_qc = X.loc[qc].notna().mean()
    det_cancer = X.loc[smp[meta.loc[smp, "Group"] == "Cancer"]].notna().mean()
    det_ctrl = X.loc[smp[meta.loc[smp, "Group"] == "Control"]].notna().mean()
    keep = (det_qc >= MIN_QC_DETECT) & ((det_cancer >= MIN_GROUP_DETECT) | (det_ctrl >= MIN_GROUP_DETECT))
    feat["det_QC"] = det_qc
    feat["det_Cancer"] = det_cancer
    feat["det_Control"] = det_ctrl
    feat["RSD_QC_raw"] = rsd(X.loc[qc])
    X1 = X.loc[:, keep]
    say(f"\nФильтр детекции (QC ≥ {MIN_QC_DETECT:.0%} и ≥ {MIN_GROUP_DETECT:.0%} в одной из групп): {n0} → {X1.shape[1]}")

    # --- 4. дрейф
    errs = {m: loo_drift_error(X1, order, qc, m) for m in DRIFT_METHODS}
    say("Leave-one-QC-out ошибка прогноза (median |Δlog2|): " + ", ".join(f"{k}={v:.4f}" for k, v in errs.items()))
    method = min(errs, key=errs.get) if args.drift == "auto" else args.drift
    say(f"Метод коррекции дрейфа: {method}")
    X2 = drift_correct(X1, order, qc, method)
    r_raw = rsd(X1.loc[qc])
    r_cor = rsd(X2.loc[qc])
    feat.loc[X2.columns, "RSD_QC_corr"] = r_cor
    say(f"Медианный RSD QC: до {r_raw.median():.1%}, после {r_cor.median():.1%}")

    # --- 5. RSD-фильтр
    keep2 = r_cor <= MAX_QC_RSD
    X3 = X2.loc[:, keep2]
    say(f"Фильтр RSD_QC ≤ {MAX_QC_RSD:.0%}: {X2.shape[1]} → {X3.shape[1]}")
    feat["kept"] = feat.index.isin(X3.columns)

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 51)
    ax.hist(r_raw.clip(upper=1), bins, alpha=0.5, label=f"до коррекции (медиана {r_raw.median():.1%})")
    ax.hist(r_cor.clip(upper=1), bins, alpha=0.5, label=f"после коррекции (медиана {r_cor.median():.1%})")
    ax.axvline(MAX_QC_RSD, c="k", ls="--")
    ax.set_xlabel("RSD в QC")
    ax.set_ylabel("число признаков")
    ax.legend()
    savefig(fig, out / "figures/02_qc_rsd.png")

    # --- 6. PQN (референс – медианный QC-спектр)
    X4, pqn_f = pqn(X3, qc)
    meta["PQN_factor"] = pqn_f
    say(f"PQN-факторы образцов: медиана {pqn_f[smp].median():.2f}, диапазон {pqn_f[smp].min():.2f}–{pqn_f[smp].max():.2f}")

    # --- 7. импутация (½ минимума признака) и log2
    minv = X4.min()
    n_missing = int(X4.isna().sum().sum())
    X5 = X4.fillna(minv / 2)
    L = np.log2(X5)
    say(f"Импутировано значений: {n_missing} ({n_missing / X4.size:.1%})")

    # --- 8. PCA
    Z = pareto(L)
    pca = PCA(5).fit(Z)
    S = pd.DataFrame(pca.transform(Z), index=Z.index)
    plot_pca(S, pca.explained_variance_ratio_, meta, order, out / "figures/03_pca_with_qc.png",
             "PCA (QC + образцы, log2, Pareto)")
    qc_spread = S.loc[qc, [0, 1]].std().mean() / S.loc[smp, [0, 1]].std().mean()
    say(f"Разброс QC относительно образцов на PC1/PC2: {qc_spread:.2f}")

    Ls = L.loc[smp]
    Zs = pareto(Ls)
    pca_s = PCA(5).fit(Zs)
    Ss = pd.DataFrame(pca_s.transform(Zs), index=Zs.index)
    plot_pca(Ss, pca_s.explained_variance_ratio_, meta, order, out / "figures/04_pca_samples.png",
             "PCA (только образцы)")
    rho = [stats.spearmanr(Ss[k], order[Ss.index]) for k in range(3)]
    say("Связь PC1–3 с порядком инъекции (Spearman): " + ", ".join(f"PC{k+1} ρ={r.statistic:.2f} p={r.pvalue:.3f}" for k, r in enumerate(rho)))
    sx = meta.loc[Ss.index, "sex"]
    for k in range(3):
        a, b = Ss.loc[sx == "м", k], Ss.loc[sx == "ж", k]
        say(f"PC{k+1}: м vs ж p={stats.mannwhitneyu(a, b).pvalue:.3g}; "
            f"Cancer vs Control p={stats.mannwhitneyu(Ss.loc[meta.loc[Ss.index,'Group']=='Cancer',k], Ss.loc[meta.loc[Ss.index,'Group']=='Control',k]).pvalue:.3g}")

    # --- 9. одномерная статистика
    cls = meta.loc[smp, "class"]
    sex_s = meta.loc[smp, "sex"]
    ctrl = smp[cls == "Control"]
    comps = {
        "Cancer_vs_Control": (smp[meta.loc[smp, "Group"] == "Cancer"], ctrl),
        "BC_vs_Control": (smp[cls == "BC"], ctrl),
        "PC_vs_Control": (smp[cls == "PC"], ctrl),
        "PC_vs_maleControl": (smp[cls == "PC"], ctrl[sex_s[ctrl] == "м"]),
        "Cancer_vs_Control_males": (smp[(meta.loc[smp, "Group"] == "Cancer") & (sex_s == "м")],
                                    ctrl[sex_s[ctrl] == "м"]),
    }
    results = {}
    summary = []
    for name, (a, b) in comps.items():
        r = compare(L, meta, list(a), list(b), feat[["rt_min", "mz", "annotation", "annotation_conf", "RSD_QC_corr"]], name)
        results[name] = r
        n_q = int((r.q_MW < 0.05).sum())
        n_q_fc = int(((r.q_MW < 0.05) & (r.log2FC_mean.abs() >= 1)).sum())
        has_adj = r.p_sex_adj.notna().any()  # в сравнениях только мужчин поправка на пол не нужна
        n_adj = int((r.q_sex_adj < 0.05).sum()) if has_adj else "n/a"
        n_both = int(((r.q_MW < 0.05) & (r.q_sex_adj < 0.05)).sum()) if has_adj else "n/a"
        summary.append(dict(comparison=name, n_case=len(a), n_control=len(b), q_MW_lt_0_05=n_q,
                            q_MW_lt_0_05_and_absFC_ge_2=n_q_fc, q_sex_adj_lt_0_05=n_adj,
                            both_q_lt_0_05=n_both))
        plot_volcano(r, out / f"figures/05_volcano_{name}.png", f"{name} (n={len(a)} vs {len(b)})")
    summary = pd.DataFrame(summary)
    say("\nОдномерная статистика (признаков):")
    say(summary.to_string(index=False))

    # --- PLS-DA: Cancer vs Control (и BC, PC)
    pls_rows = []
    rng = np.random.RandomState(SEED)
    for name in ["Cancer_vs_Control", "BC_vs_Control", "PC_vs_Control"]:
        a, b = comps[name]
        idx = list(a) + list(b)
        y = np.r_[np.ones(len(a)), np.zeros(len(b))]
        Zc = L.loc[idx].to_numpy()
        auc, q2 = plsda_cv(Zc, y)
        perm = []
        for _ in range(N_PERM):
            perm.append(plsda_cv(Zc, rng.permutation(y), n_rep=1, seed=rng.randint(1e6))[0])
        p_perm = (np.sum(np.array(perm) >= auc) + 1) / (N_PERM + 1)
        Zp = pareto(L.loc[idx]).to_numpy()
        pls = PLSRegression(2, scale=False).fit(Zp, y)
        vip = pd.Series(vip_scores(pls), index=L.columns)
        results[name]["VIP_PLSDA"] = vip
        pls_rows.append(dict(comparison=name, CV_AUC=auc, Q2=q2, perm_p=p_perm, n_perm=N_PERM))
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        sc = pd.DataFrame(pls.x_scores_, index=idx)
        for c, col in [(1, "#d7301f"), (0, "#2b7bba")]:
            ax[0].scatter(sc[y == c][0], sc[y == c][1], c=col, s=25, edgecolor="k", lw=0.3,
                          label=name.split("_vs_")[0] if c else "Control")
        ax[0].set_xlabel("t1")
        ax[0].set_ylabel("t2")
        ax[0].legend()
        ax[0].set_title(f"PLS-DA {name} (подгонка, без CV)")
        ax[1].hist(perm, 30, color="#bbbbbb")
        ax[1].axvline(auc, c="r")
        ax[1].set_title(f"CV AUC={auc:.2f}, Q2={q2:.2f}, perm p={p_perm:.3f}")
        ax[1].set_xlabel("CV AUC при перестановке меток")
        savefig(fig, out / f"figures/06_plsda_{name}.png")
    pls_df = pd.DataFrame(pls_rows)
    say("\nPLS-DA (2 компоненты, 5×5-fold CV, Pareto внутри фолда):")
    say(pls_df.to_string(index=False))

    # --- боксплоты топ-признаков (аннотированные, устойчивые к поправке на пол)
    r = results["Cancer_vs_Control"]
    top = r[(r.q_MW < 0.05) & (r.q_sex_adj < 0.05) & (r.annotation_conf != "unknown")].head(16)
    plot_box(L, meta, [(f, f"{t.annotation[:35]}\nRT {t.rt_min} m/z {t.mz}  q={t.q_MW:.1e}")
                       for f, t in top.iterrows()],
             out / "figures/07_top_Cancer_vs_Control.png", ["Control", "BC", "PC", "KC", "PanC", "BC+PC"])

    # --- сохранение
    meta_out = meta.drop(columns=["diag_text"]).copy()
    meta_out["TIC_raw"] = tic
    with pd.ExcelWriter(out / "metabolomics_results.xlsx") as w:
        pd.DataFrame({"log": log}).to_excel(w, sheet_name="log", index=False)
        summary.to_excel(w, sheet_name="summary_univariate", index=False)
        pls_df.to_excel(w, sheet_name="summary_PLSDA", index=False)
        meta_out.sort_values("inj_order").to_excel(w, sheet_name="samples")
        feat.to_excel(w, sheet_name="features_QC")
        for name, rr in results.items():
            rr[rr.p_MW < 0.05].to_excel(w, sheet_name=name[:31])
    L.T.join(feat[["rt_min", "mz", "annotation", "annotation_conf"]]).to_csv(out / "processed_log2_matrix.csv")
    pd.concat(results.values()).to_csv(out / "univariate_all.csv")
    (out / "log.txt").write_text("\n".join(log), encoding="utf-8")
    say(f"\nГотово → {out}")


if __name__ == "__main__":
    main()

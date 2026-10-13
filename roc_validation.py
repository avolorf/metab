#!/usr/bin/env python3
"""Корректная оценка диагностической способности (ROC) без утечки информации.

Проблема «наивной» ROC: признаки отобраны по p-value на тех же образцах, на которых затем считается AUC,
поэтому AUC завышен (selection bias). Здесь:
  1. AUC отдельных курированных соединений с 95 % бутстрэп-ДИ: на всех образцах и только у мужчин
     (описательно, остаётся оптимистичным из-за отбора);
  2. базовая модель «только пол» — какой AUC даёт сам конфаундер;
  3. вложенная (nested) повторная CV: отбор top-k признаков (Манн-Уитни) ВНУТРИ обучающего фолда,
     логистическая регрессия, прогноз на тестовом фолде; 20 × 5-fold; перестановочный тест;
  4. PLS-DA на всех признаках (без отбора) — средняя ROC по повторам CV;
  5. для PC — то же без ряда контаминантов (ПАВ).

Запуск после process.py и annotate_visualize.py:
    python roc_validation.py --res results
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold

from process import plsda_cv_predictions

N_REP = 20
N_PERM = 100
TOP_K = 5
N_BOOT = 2000
SEED = 1
FPR_GRID = np.linspace(0, 1, 101)


def boot_auc(y, s, rng):
    a = roc_auc_score(y, s)
    sign = 1 if a >= 0.5 else -1
    idx1, idx0 = np.where(y == 1)[0], np.where(y == 0)[0]
    b = []
    for _ in range(N_BOOT):
        i = np.r_[rng.choice(idx1, len(idx1)), rng.choice(idx0, len(idx0))]
        b.append(roc_auc_score(y[i], sign * s[i]))
    a = max(a, 1 - a)
    return a, np.percentile(b, 2.5), np.percentile(b, 97.5)


def nested_cv(X, y, rng, n_rep=N_REP, k=TOP_K):
    """Отбор признаков внутри фолда -> логистическая регрессия. Возвращает AUC по повторам и средние ROC."""
    aucs, tprs, chosen = [], [], []
    for _ in range(n_rep):
        skf = StratifiedKFold(5, shuffle=True, random_state=rng.randint(1e6))
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            p = stats.mannwhitneyu(X[tr][y[tr] == 1], X[tr][y[tr] == 0], axis=0).pvalue
            sel = np.argsort(p)[:k]
            chosen.extend(sel)
            mu, sd = X[tr][:, sel].mean(0), X[tr][:, sel].std(0) + 1e-9
            m = LogisticRegression(C=1.0, max_iter=1000).fit((X[tr][:, sel] - mu) / sd, y[tr])
            pred[te] = m.predict_proba((X[te][:, sel] - mu) / sd)[:, 1]
        aucs.append(roc_auc_score(y, pred))
        fpr, tpr, _ = roc_curve(y, pred)
        tprs.append(np.interp(FPR_GRID, fpr, tpr))
    return np.array(aucs), np.array(tprs), np.array(chosen)


def plsda_rep(X, y, rng, n_rep=N_REP):
    aucs, tprs = [], []
    for _ in range(n_rep):
        pred = plsda_cv_predictions(X, y, seed=rng.randint(1e6))
        aucs.append(roc_auc_score(y, pred))
        fpr, tpr, _ = roc_curve(y, pred)
        tprs.append(np.interp(FPR_GRID, fpr, tpr))
    return np.array(aucs), np.array(tprs)


def perm_test(X, y, observed, rng):
    null = []
    for _ in range(N_PERM):
        a, _, _ = nested_cv(X, rng.permutation(y), rng, n_rep=1)
        null.append(a[0])
    return (np.sum(np.array(null) >= observed) + 1) / (N_PERM + 1), np.array(null)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results")
    a = ap.parse_args()
    res = Path(a.res)
    rng = np.random.RandomState(SEED)

    M = pd.read_csv(res / "processed_log2_matrix.csv", index_col=0)
    info = M[["rt_min", "mz", "annotation", "annotation_conf"]]
    L = M.drop(columns=info.columns).T
    L.index = L.index.astype(str)
    meta = pd.read_excel(res / "metabolomics_results.xlsx", sheet_name="samples", index_col=0)
    meta.index = meta.index.astype(str)
    cur = pd.read_excel(res / "compounds_BC_PC.xlsx", sheet_name="curated_candidates")
    contam = set(cur.loc[cur.category.str.startswith("контаминант"), "feature"])

    single_rows, model_rows = [], []
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, case in zip(axes, ["BC", "PC"]):
        smp = meta.index[meta["class"].isin([case, "Control"])]
        y = (meta.loc[smp, "class"] == case).astype(int).to_numpy()
        sex = meta.loc[smp, "sex"]
        males = (sex == "м").to_numpy()

        # 1. отдельные соединения
        c = cur[(cur.focus == case) & ~cur.category.str.startswith(("неоднозн", "контаминант"))].drop_duplicates("feature")
        for _, r in c.iterrows():
            s = L.loc[smp, r.feature].to_numpy()
            auc, lo, hi = boot_auc(y, s, rng)
            aucm, lom, him = boot_auc(y[males], s[males], rng)
            single_rows.append(dict(comparison=f"{case}_vs_Control", name=r["name"], category=r.category,
                                    AUC=auc, CI_low=lo, CI_high=hi, AUC_males=aucm, CI_low_males=lom,
                                    CI_high_males=him))

        # 2. только пол
        known = sex.notna().to_numpy()
        auc_sex = roc_auc_score(y[known], males[known].astype(float))
        model_rows.append(dict(comparison=f"{case}_vs_Control", model="только пол (конфаундер)",
                               AUC_mean=auc_sex, AUC_sd=0, perm_p=np.nan))
        ax.plot([0, 1 - (1 - males[known][y[known] == 0]).mean(), 1],
                [0, males[known][y[known] == 1].mean(), 1], ls="--", c="grey", label=f"только пол: {auc_sex:.2f}")

        sets = {"все признаки": L.loc[smp].to_numpy()}
        if case == "PC":
            keep = [f for f in L.columns if f not in contam]
            sets["без ряда ПАВ"] = L.loc[smp, keep].to_numpy()
        colors = iter(["#d7301f", "#2b7bba", "#1a9850", "#7b3294"])
        for nm, X in sets.items():
            # 3. nested CV
            aucs, tprs, chosen = nested_cv(X, y, rng)
            p, _ = perm_test(X, y, aucs.mean(), rng)
            col = next(colors)
            ax.plot(FPR_GRID, tprs.mean(0), c=col, lw=2,
                    label=f"nested CV top-{TOP_K}+LR, {nm}: {aucs.mean():.2f}±{aucs.std():.2f} (perm p={p:.3f})")
            ax.fill_between(FPR_GRID, np.percentile(tprs, 10, 0), np.percentile(tprs, 90, 0), color=col, alpha=0.12)
            feats = pd.Series(chosen).value_counts().head(8)
            cols = L.columns if nm == "все признаки" else pd.Index(keep)
            top_names = "; ".join(f"{info.loc[cols[i], 'annotation'][:30]} ({n / (N_REP * 5):.0%})" for i, n in feats.items())
            model_rows.append(dict(comparison=f"{case}_vs_Control", model=f"nested CV top-{TOP_K} + logistic, {nm}",
                                   AUC_mean=aucs.mean(), AUC_sd=aucs.std(), perm_p=p,
                                   most_selected=top_names))
            # 4. PLS-DA все признаки
            aucs2, tprs2 = plsda_rep(X, y, rng)
            col = next(colors)
            ax.plot(FPR_GRID, tprs2.mean(0), c=col, lw=2, ls="-.",
                    label=f"PLS-DA, {nm}, CV: {aucs2.mean():.2f}±{aucs2.std():.2f}")
            model_rows.append(dict(comparison=f"{case}_vs_Control", model=f"PLS-DA 2 комп., {nm}",
                                   AUC_mean=aucs2.mean(), AUC_sd=aucs2.std(), perm_p=np.nan))

        # то же только у мужчин (снимает конфаундинг по полу)
        Xm = L.loc[smp[males]].to_numpy()
        if case == "PC":
            Xm = L.loc[smp[males], keep].to_numpy()
        aucs_m, _, _ = nested_cv(Xm, y[males], rng)
        model_rows.append(dict(comparison=f"{case}_vs_Control", model=f"nested CV top-{TOP_K} + logistic, только мужчины"
                               + (", без ряда ПАВ" if case == "PC" else ""),
                               AUC_mean=aucs_m.mean(), AUC_sd=aucs_m.std(), perm_p=np.nan,
                               n=f"{y[males].sum()} vs {(1 - y[males]).sum()}"))

        ax.plot([0, 1], [0, 1], ls=":", c="k", lw=0.6)
        ax.set_title(f"{case} vs Control (n={y.sum()} / {len(y) - y.sum()}): ROC без утечки")
        ax.set_xlabel("1 − специфичность")
        ax.set_ylabel("чувствительность")
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(res / "figures/14_roc_validated.png", dpi=150)

    single = pd.DataFrame(single_rows).sort_values(["comparison", "AUC"], ascending=[True, False])
    models = pd.DataFrame(model_rows)
    with pd.ExcelWriter(res / "roc_validation.xlsx") as w:
        models.to_excel(w, sheet_name="models_nested_CV", index=False)
        single.to_excel(w, sheet_name="single_compounds_bootstrap", index=False)
    pd.set_option("display.width", 250)
    print(models.drop(columns=["most_selected"], errors="ignore").round(3).to_string(index=False))
    print(single.round(2).to_string(index=False))


if __name__ == "__main__":
    main()

"""Рисунки: PNG (300 dpi) + PDF."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from . import config as C
from .multivariate import FPR

GC = {"QC": "#7f7f7f", "Control": "#2b7bba", "BC": "#d7301f", "PC": "#fc8d59"}
SEXC = {"M": "#404040", "F": "#f4a6d7", "": "#ffffff"}
CATC = {"contaminant": "#e7298a", "drug": "#7570b3", "food": "#d95f02", "tobacco": "#8c510a",
        "endogenous": "#1b9e77", "ambiguous": "#999999", "annotated (high score)": "#a6d854",
        "putative (low score)": "#e5f5c0", "unknown": "#d9d9d9"}


def lab(g):
    return C.GROUP_LABELS.get(g, g)


def save(fig, path):
    fig.savefig(f"{path}.png", dpi=C.FIG_DPI, bbox_inches="tight")
    fig.savefig(f"{path}.pdf", bbox_inches="tight")
    plt.close(fig)


def catcolor(c):
    c = str(c)
    for k, v in CATC.items():
        if c.startswith(k) or c.startswith(f"possible {k}"):
            return v
    return "#d9d9d9"


def qc_overview(X, meta, order, rsd_raw, factors, path):
    fig, ax = plt.subplots(1, 3, figsize=(18, 4.5), gridspec_kw={"width_ratios": [2, 1, 1]})
    tic = X.sum(axis=1)
    for g, col in GC.items():
        m = (meta.group == g) & (meta.study | (meta.kind == "QC"))
        ax[0].scatter(order[m], tic[m] / 1e9, c=col, s=45 if g == "QC" else 18, marker="s" if g == "QC" else "o",
                      edgecolor="k", lw=0.3, label=lab(g))
    ax[0].set_yscale("log")
    ax[0].set(xlabel="порядок инъекции", ylabel="суммарная интенсивность ×1e9", title="Суммарный сигнал (сырые данные)")
    ax[0].legend(fontsize=8)
    ax[1].hist(rsd_raw.clip(upper=1), np.linspace(0, 1, 51), color="#7f7f7f")
    ax[1].axvline(C.MAX_QC_RSD, ls="--", c="k")
    ax[1].set(xlabel="RSD в QC (без коррекции дрейфа)", ylabel="признаков",
              title=f"медиана {rsd_raw[rsd_raw <= C.MAX_QC_RSD].median():.1%} (после фильтра)")
    s = meta.index[meta.study]
    ax[2].scatter(factors.loc[s, "PQN"], factors.loc[s, "MSTUS"], c=meta.loc[s, "group"].map(GC), s=20, edgecolor="k", lw=0.3)
    low = s[factors.loc[s, "PQN"] < C.LOW_DILUTION_FACTOR]
    for i in low:
        ax[2].annotate(i, (factors.loc[i, "PQN"], factors.loc[i, "MSTUS"]), fontsize=7)
    ax[2].axvline(C.LOW_DILUTION_FACTOR, ls=":", c="k")
    ax[2].set(xscale="log", yscale="log", xlabel="PQN-фактор", ylabel="MSTUS-фактор", title="Факторы разведения")
    save(fig, path)


def pca_panels(S, var, meta, order, path):
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    for g, col in GC.items():
        m = meta.loc[S.index, "group"] == g
        ax[0].scatter(S.loc[m, 0], S.loc[m, 1], c=col, s=50 if g == "QC" else 22, marker="s" if g == "QC" else "o",
                      edgecolor="k", lw=0.3, label=f"{lab(g)} (n={m.sum()})")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Группа")
    smp = S.index[meta.loc[S.index, "study"]]
    sc = ax[1].scatter(S.loc[smp, 0], S.loc[smp, 1], c=meta.loc[smp, "contam_index"], cmap="magma_r", s=24,
                       edgecolor="k", lw=0.3)
    plt.colorbar(sc, ax=ax[1], label="индекс загрязнения")
    ax[1].set_title("Индекс загрязнения (ряд диэтаноламинов)")
    sc = ax[2].scatter(S[0], S[1], c=order.loc[S.index], cmap="viridis", s=24, edgecolor="k", lw=0.3)
    plt.colorbar(sc, ax=ax[2], label="порядок инъекции")
    ax[2].set_title("Порядок инъекции")
    for a in ax:
        a.set(xlabel=f"PC1 ({var[0]:.1%})", ylabel=f"PC2 ({var[1]:.1%})")
    save(fig, path)


def contamination(meta, series, path):
    s = meta[meta.study]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for i, g in enumerate(C.STUDY_GROUPS):
        v = s.contam_index[s.group == g]
        ax[0].boxplot(v, positions=[i], widths=0.6, showfliers=False)
        ax[0].scatter(np.random.RandomState(i).normal(i, 0.07, len(v)), v, c=GC[g], s=18, edgecolor="k", lw=0.3)
    ax[0].axhline(s.contam_index.median(), ls=":", c="k", lw=0.8)
    ax[0].set_xticks(range(3), [f"{lab(g)}\n(n={(s.group == g).sum()})" for g in C.STUDY_GROUPS])
    ax[0].set(ylabel="индекс загрязнения (средний z)", title="Индекс загрязнения; пунктир — граница high/low")
    for inx, mk, nm in [(1, "o", "ряд диэтаноламинов (индекс)"), (0, "^", "вне индекса (C24H30O6, аминоспирт)")]:
        d = series[series.in_index == inx]
        ax[1].scatter(d.mz, d.rt_min, marker=mk, s=45, c="#e7298a" if inx else "#555555", label=nm)
        for _, r in d.iterrows():
            ax[1].annotate(r.formula, (r.mz, r.rt_min), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax[1].set(xlabel="m/z", ylabel="RT, мин", title="Гомологический ряд: шаги C2H4 / O / C2H4O")
    ax[1].legend(fontsize=8)
    save(fig, path)


def volcano(r, name, names, path):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    x, y = r.beta_model, -np.log10(r.p_model)
    cat = r.compound_category.astype(str)
    expo = cat.isin(["contaminant", "drug", "food", "tobacco"])
    ax.scatter(x[~expo], y[~expo], s=5, c="#c8c8c8", lw=0)
    ax.scatter(x[expo], y[expo], s=12, c=[catcolor(c) for c in cat[expo]], marker="x", lw=0.8)
    f = r.finding
    ax.scatter(x[f], y[f], s=28, facecolor="none", edgecolor="k", lw=1)
    for i in r[f].sort_values("p_model").head(12).index:
        ax.annotate(names[i][:30], (x[i], y[i]), fontsize=6)
    q = r.p_model[r.q_model < C.FDR]
    if len(q):
        ax.axhline(-np.log10(q.max()), ls="--", c="k", lw=0.6)
    ax.set(xlabel="log2 FC (коэффициент модели с ковариатами)", ylabel="−log10 p (модель)",
           title=f"{name}: {f.sum()} находок (кружки); × — экспозиции; пунктир — q={C.FDR}")
    handles = [Line2D([], [], marker="x", ls="", color=CATC[k], label=k) for k in ["contaminant", "drug", "food"]]
    handles.append(Line2D([], [], marker="o", ls="", mfc="none", mec="k", label="находка"))
    ax.legend(handles=handles, fontsize=7)
    save(fig, path)


def heatmap(L, meta, r, names, comp_cat, case, ref, path, title, top=40):
    sel = r[~r.compound_category.isin(["contaminant"])].sort_values(["finding", "p_model"], ascending=[False, True])
    sel = sel[(sel.q_model < C.FDR) | sel.finding].head(top)
    if len(sel) < 2:
        return
    smp = list(meta.index[meta.study & meta.group.isin([case, ref])])
    Z = L.loc[smp, sel.index]
    Z = ((Z - Z.mean()) / Z.std()).clip(-3, 3).T
    Z.index = [("★ " if sel.finding[i] else "") + names[i] for i in sel.index]
    age_norm = Normalize(meta.age.min(), meta.age.max())
    ci_norm = Normalize(meta.contam_index[meta.study].quantile(0.02), meta.contam_index[meta.study].quantile(0.98))
    cmap_age, cmap_ci = plt.get_cmap("Greens"), plt.get_cmap("magma_r")
    cc = pd.DataFrame({
        "группа": meta.loc[smp, "group"].map(GC),
        "пол": meta.loc[smp, "sex"].map(SEXC).fillna("#ffffff"),
        "возраст": [matplotlib.colors.to_hex(cmap_age(age_norm(a))) if pd.notna(a) else "#ffffff" for a in meta.loc[smp, "age"]],
        "индекс загр.": [matplotlib.colors.to_hex(cmap_ci(ci_norm(v))) for v in meta.loc[smp, "contam_index"]],
    }, index=smp)
    rc = pd.Series([catcolor(comp_cat.get(i, "")) for i in sel.index], index=Z.index, name="категория")
    g = sns.clustermap(Z, cmap="RdBu_r", center=0, vmin=-3, vmax=3, method="ward", metric="euclidean",
                       col_colors=cc, row_colors=rc, xticklabels=False, yticklabels=True,
                       figsize=(15, 0.26 * len(sel) + 4.5), dendrogram_ratio=(0.07, 0.08),
                       cbar_pos=(0.005, 0.03, 0.012, 0.1), cbar_kws={"label": "z-score"})
    g.ax_heatmap.tick_params(axis="y", labelsize=7)
    g.fig.suptitle(title + "  (★ — находка; одна строка = одно соединение; * — low-score аннотация, см. таблицу)", y=1.01)
    h = [Patch(color=GC[x], label=lab(x)) for x in [case, ref]]
    h += [Patch(color=SEXC["M"], label="мужчины"), Patch(color=SEXC["F"], label="женщины"),
          Patch(color="#ffffff", label="нет данных")]
    h += [Patch(color=catcolor(k), label=k) for k in ["endogenous", "annotated (high score)", "putative (low score)",
                                                       "unknown", "food", "drug"]]
    g.ax_heatmap.legend(handles=h, loc="upper left", bbox_to_anchor=(1.25, 1.25), fontsize=7, frameon=False)
    for i, (cm, nm, nr) in enumerate([(cmap_age, "возраст, лет", age_norm), (cmap_ci, "индекс загрязнения", ci_norm)]):
        cax = g.fig.add_axes([0.93, 0.80 - i * 0.1, 0.06, 0.012])
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=nr, cmap=cm), cax=cax, orientation="horizontal")
        cb.set_label(nm, fontsize=7)
        cb.ax.tick_params(labelsize=6)
    g.savefig(f"{path}.png", dpi=C.FIG_DPI)
    g.savefig(f"{path}.pdf")
    plt.close(g.fig)


def boxplots(L, meta, feats, names, path, title):
    if not len(feats):
        return
    n = len(feats)
    cols = 4
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.2 * rows), squeeze=False)
    for ax, f in zip(axes.ravel(), feats):
        data = [L.loc[meta.index[meta.study & (meta.group == g)], f] for g in C.STUDY_GROUPS]
        ax.boxplot(data, tick_labels=[lab(g) for g in C.STUDY_GROUPS], showfliers=False)
        for i, d in enumerate(data):
            ax.scatter(np.random.RandomState(i).normal(i + 1, 0.06, len(d)), d, s=8, c=GC[C.STUDY_GROUPS[i]], alpha=0.7)
        ax.set_title(names[f][:45], fontsize=7)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(title)
    save(fig, path)


def roc(curves, singles, name, path):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    style = {("all", "panel"): ("#d7301f", "-", 2.5), ("matched", "panel"): ("#d7301f", "--", 2),
             ("all", "plsda"): ("#2b7bba", "-", 1.8), ("matched", "plsda"): ("#2b7bba", "--", 1.5),
             ("all", "sex"): ("#999999", ":", 1.2), ("all", "contam_index"): ("#e7298a", ":", 1.2)}
    names = {"panel": f"панель top-{C.PANEL_K} (CV)", "plsda": "PLS-DA (CV)", "sex": "только пол",
             "contam_index": "только индекс загрязнения"}
    for (st, m), (tpr, auc) in curves.items():
        if (st, m) in style:
            c, ls, lw = style[(st, m)]
            ax.plot(FPR, tpr, c=c, ls=ls, lw=lw,
                    label=f"{names[m]}, {'подбор по возрасту' if st == 'matched' else 'все'}: {auc:.2f}")
    pal = plt.get_cmap("tab10")
    for k, (lbl, fpr, tpr, a, lo, hi) in enumerate(singles[:6]):
        ax.plot(fpr, tpr, lw=0.9, c=pal(k + 2), alpha=0.8, label=f"{lbl[:28]}: {a:.2f} ({lo:.2f}–{hi:.2f})*")
    ax.plot([0, 1], [0, 1], c="k", lw=0.5, ls=":")
    ax.set(xlabel="1 − специфичность", ylabel="чувствительность",
           title=f"{name}: модели с CV и эндогенные находки\n* одиночные соединения отобраны на тех же данных — AUC оптимистичен")
    ax.legend(fontsize=7, loc="lower right")
    save(fig, path)


def norm_concordance(res_by_norm, path):
    comps = list(next(iter(res_by_norm.values())).keys())
    fig, axes = plt.subplots(1, len(comps), figsize=(4.5 * len(comps), 4.3), squeeze=False)
    a, b = list(res_by_norm)
    for ax, c in zip(axes[0], comps):
        x, y = res_by_norm[a][c].beta_model, res_by_norm[b][c].beta_model.reindex(res_by_norm[a][c].index)
        ax.scatter(x, y, s=3, c="#999999", lw=0)
        sig = (res_by_norm[a][c].q_model < C.FDR) | (res_by_norm[b][c].q_model.reindex(x.index) < C.FDR)
        ax.scatter(x[sig], y[sig], s=6, c="#d7301f", lw=0)
        ax.plot([-4, 4], [-4, 4], c="k", lw=0.5)
        ax.set(xlabel=f"log2FC, {a}", ylabel=f"log2FC, {b}", title=f"{c}: ρ={x.corr(y, method='spearman'):.2f}")
    save(fig, path)


def pathways(pw, path, title, top=15):
    pw = pw.sort_values("p-value").head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(pw) + 1.3))
    ax.barh(pw["pathway"], -np.log10(pw["p-value"]),
            color=["#d7301f" if q < 0.05 else ("#fdae6b" if p < 0.05 else "#cccccc") for p, q in zip(pw["p-value"], pw.q_BH)])
    ax.axvline(-np.log10(0.05), ls="--", c="k", lw=0.6)
    ax.set(xlabel="−log10 p (mummichog)", title=title + "\nкрасный — q_BH<0.05; оранжевый — только p<0.05")
    ax.tick_params(axis="y", labelsize=8)
    save(fig, path)

#!/usr/bin/env python3
"""Публикационные рисунки и таблицы (англ.) для журнала — из уже посчитанных results/.

    python make_figures.py --res results        # после run_all.py (вызывается им автоматически)

Выход: results/publication/
    Figure1..Figure6, FigureS1..S3  — PDF (вектор, редактируемый текст) + TIFF (LZW) и PNG, 600 dpi
    Table1_demographics.csv/.xlsx, Table2_findings_BC_vs_Control.csv/.xlsx, TableS1_classification.csv/.xlsx
    figure_captions.md               — черновики подписей
Ширина: 17 см (полная полоса) / 8.5 см (колонка). Шрифт Liberation Sans (метрически = Arial), 7 pt.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize, to_hex
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy import stats

from pipeline import config as C
from pipeline import features, io, multivariate
from pipeline import stats as st

CM = 1 / 2.54
FULL, HALF = 17 * CM, 8.5 * CM
DPI = 600
# валидированная категориальная палитра (первые 3 слота: проходит CVD и normal-vision для всех пар)
COL = {"Control": "#2a78d6", "BC": "#eb6834", "PC": "#1baf7a", "QC": "#6f6e69"}
MRK = {"Control": "o", "BC": "s", "PC": "^", "QC": "D"}      # вторичное кодирование формы
LAB = {"Control": "Control", "BC": "BC", "PC": "PC", "QC": "QC"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
EXPO_COL = {"contaminant": "#6250d6", "drug": "#e34948", "food": "#eda100", "tobacco": "#52514e"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": DPI, "figure.dpi": 150,
})


def panel(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom", ha="left")


def save(fig, out, name):
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(out / f"{name}_tmp.png", dpi=DPI, bbox_inches="tight")
    im = Image.open(out / f"{name}_tmp.png").convert("RGB")
    im.save(out / f"{name}.tif", compression="tiff_lzw", dpi=(DPI, DPI))
    (out / f"{name}_tmp.png").unlink()
    plt.close(fig)


def scatter_groups(ax, x, y, groups, size=10, order=("QC", "Control", "BC", "PC"), legend=True, meta=None):
    for g in order:
        m = groups == g
        if m.any():
            ax.scatter(x[m], y[m], s=size * (1.6 if g == "QC" else 1), c=COL[g], marker=MRK[g],
                       edgecolor="white", lw=0.4, label=f"{LAB[g]} (n = {m.sum()})", zorder=3)
    if legend:
        ax.legend(loc="best", handletextpad=0.2, borderaxespad=0.2)


def pfmt(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


# ------------------------------------------------------------------ загрузка результатов
def load(res):
    t = res / "tables"
    meta = pd.read_csv(t / "samples.csv", index_col=0, dtype={0: str})
    meta.index = meta.index.astype(str)
    meta["study"] = meta["study"].astype(bool)
    M = pd.read_csv(t / "processed_log2_matrix.csv", index_col=0)
    info = M[["rt_min", "mz", "annotation", "annotation_conf"]]
    L = M.drop(columns=info.columns).T
    L.index = L.index.astype(str)
    feat = pd.read_csv(t / "features.csv", index_col=0)
    stats_ = {n: pd.read_csv(t / f"stats_{n}.csv", index_col=0) for n in C.COMPARISONS}
    mv = pd.read_csv(t / "multivariate_CV.csv")
    series = pd.read_csv(t / "contaminant_series.csv", index_col=0)
    cur = pd.read_csv(t / "candidates_curated.csv")
    pairs = pd.read_csv(t / "age_matched_pairs.csv", dtype={"case": str, "ref": str})
    demo = pd.read_csv(t / "demographics.csv")
    demo_t = pd.read_csv(t / "demographics_tests.csv")
    pw = pd.read_csv(t / "pathways_appendix.csv") if (t / "pathways_appendix.csv").exists() else None
    return meta, L, info, feat, stats_, mv, series, cur, pairs, demo, demo_t, pw


def names_en(feat, cur):
    """Подпись соединения: курированное англ. название (MSI) > уверенная аннотация > unknown m/z@RT."""
    out = {}
    cn = cur.dropna(subset=["feature"]).set_index("feature")
    for f, r in feat.iterrows():
        tag = f"m/z {r.mz:.4f} @ {r.rt_min:.2f} min"
        if f in cn.index and isinstance(cn.loc[f, "name_en"], str) and cn.loc[f, "name_en"]:
            out[f] = f"{cn.loc[f, 'name_en']} (MSI {int(cn.loc[f, 'MSI_level'])}), {tag}"
        elif r.annotation_conf == "high":
            out[f] = f"{r.annotation[:34]} (putative), {tag}"
        else:
            out[f] = f"unknown, {tag}"
    return pd.Series(out)


# ------------------------------------------------------------------ рисунки
def fig1_quality(out, meta, L, feat, X_all):
    fig, axs = plt.subplots(2, 2, figsize=(FULL, 12 * CM), gridspec_kw={"hspace": 0.45, "wspace": 0.3})
    ana = meta.index[meta.study | (meta.kind == "QC")]
    grp = meta.loc[ana, "group"]
    tic = np.log10(X_all.loc[ana].sum(axis=1))
    ax = axs[0, 0]
    scatter_groups(ax, meta.loc[ana, "inj_order"].to_numpy(), tic.to_numpy(), grp.to_numpy(), legend=False)
    q = meta.index[meta.kind == "QC"]
    ax.plot(meta.loc[q, "inj_order"], tic[q], c=COL["QC"], lw=0.8, zorder=2)
    ax.set(xlabel="Injection order", ylabel="log$_{10}$ total ion intensity")
    ax.legend(handles=[Line2D([], [], marker=MRK[g], ls="", mfc=COL[g], mec="white", ms=5, label=LAB[g])
                       for g in ["QC", "Control", "BC", "PC"]], ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel(ax, "a")
    ax = axs[0, 1]
    r = feat.loc[feat.pass_detection, "RSD_QC_raw"].clip(upper=1) * 100
    ax.hist(r, np.linspace(0, 100, 51), color=COL["QC"], edgecolor="white", lw=0.3)
    ax.axvline(C.MAX_QC_RSD * 100, c=INK, ls="--", lw=0.7)
    kept = feat.kept.sum()
    ax.text(C.MAX_QC_RSD * 100 + 2, ax.get_ylim()[1] * 0.92,
            f"RSD ≤ {C.MAX_QC_RSD:.0%}: {kept} features\nmedian RSD {r[r <= 30].median():.1f}%", fontsize=6.5, va="top")
    ax.set(xlabel="RSD in pooled QC (%)", ylabel="Number of features")
    panel(ax, "b")
    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    S, var = multivariate.pca_scores(L.loc[ana, reps])
    ax = axs[1, 0]
    scatter_groups(ax, S[0].to_numpy(), S[1].to_numpy(), grp.to_numpy(), legend=False)
    ax.set(xlabel=f"PC1 ({var[0]:.1%})", ylabel=f"PC2 ({var[1]:.1%})")
    panel(ax, "c")
    ax = axs[1, 1]
    s = meta.index[meta.study]
    sc = ax.scatter(S.loc[s, 0], S.loc[s, 1], c=meta.loc[s, "contam_index"], cmap="magma_r", s=10, edgecolor="white", lw=0.3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Contamination index", fontsize=6.5)
    cb.outline.set_linewidth(0.4)
    ax.set(xlabel=f"PC1 ({var[0]:.1%})", ylabel=f"PC2 ({var[1]:.1%})")
    panel(ax, "d")
    save(fig, out, "Figure1_data_quality")


def fig2_contamination(out, meta, series):
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 6.5 * CM), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.3})
    ax = axs[0]
    d = series[series.in_index == 1].sort_values("mz")
    o = series[series.in_index == 0]
    ax.scatter(d.mz, d.rt_min, s=22, c=EXPO_COL["contaminant"], edgecolor="white", lw=0.4, zorder=3,
               label="Diethanolamine series (index)")
    ax.scatter(o.mz, o.rt_min, s=18, c="white", edgecolor=INK2, lw=0.7, marker="D", zorder=3, label="Co-varying, not in index")
    from adjustText import adjust_text
    seen, texts = set(), []
    allp = pd.concat([d, o])
    for _, r in allp.sort_values("rt_min").iterrows():
        key = (r.formula, round(r.mz, 1))
        if key in seen:
            continue
        seen.add(key)
        texts.append(ax.text(r.mz, r.rt_min, r.formula.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")),
                             fontsize=5.5, color=INK2))
    ax.set_xlim(allp.mz.min() - 15, allp.mz.max() + 25)
    adjust_text(texts, x=allp.mz.to_numpy(), y=allp.rt_min.to_numpy(), ax=ax, expand=(1.3, 1.6),
                arrowprops=dict(arrowstyle="-", color=GRID, lw=0.4))
    ax.set(xlabel="m/z ([M+H]$^+$)", ylabel="Retention time (min)")
    ax.legend(loc="lower right")
    ax.set_title("Mass increments: +C$_2$H$_4$ 28.031, +O 15.995, +C$_2$H$_4$O 44.026", fontsize=6.5, color=INK2)
    panel(ax, "a", x=-0.12)
    ax = axs[1]
    s = meta[meta.study]
    for i, g in enumerate(C.STUDY_GROUPS):
        v = s.contam_index[s.group == g]
        ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color=INK, lw=0.9),
                   boxprops=dict(lw=0.6, color=INK2), whiskerprops=dict(lw=0.6, color=INK2), capprops=dict(lw=0.6, color=INK2))
        ax.scatter(np.random.RandomState(i).normal(i, 0.08, len(v)), v, s=9, c=COL[g], marker=MRK[g], edgecolor="white",
                   lw=0.3, zorder=3)
    ax.set_xticks(range(3), [f"{LAB[g]}\n(n = {(s.group == g).sum()})" for g in C.STUDY_GROUPS])
    ax.set_ylabel("Contamination index (mean z-score)")
    ymax = s.contam_index.max()
    for k, (a, b) in enumerate([("Control", "BC"), ("BC", "PC"), ("Control", "PC")]):
        p = stats.mannwhitneyu(s.contam_index[s.group == a], s.contam_index[s.group == b]).pvalue
        i, j = C.STUDY_GROUPS.index(a), C.STUDY_GROUPS.index(b)
        yy = ymax + 0.3 + (0 if k < 2 else 0.75)
        ax.plot([i, i, j, j], [yy - 0.08, yy, yy, yy - 0.08], c=INK2, lw=0.6)
        ax.text((i + j) / 2, yy + 0.06, pfmt(p), ha="center", va="bottom", fontsize=6)
    ax.set_ylim(top=ymax + 1.5)
    panel(ax, "b", x=-0.2)
    save(fig, out, "Figure2_contamination_index")


def volcano_ax(ax, r, nm, title, label_top=6):
    x, y = r.beta_model, -np.log10(r.p_model)
    cat = r.compound_category.astype(str)
    expo = cat.isin(list(EXPO_COL))
    ax.scatter(x[~expo], y[~expo], s=3, c="#bdbcb8", lw=0, rasterized=True)
    for c, col in EXPO_COL.items():
        m = cat == c
        if m.any():
            ax.scatter(x[m], y[m], s=9, c=col, marker="x", lw=0.6, label=c.capitalize(), zorder=3)
    f = r.finding.astype(bool)
    ax.scatter(x[f], y[f], s=14, facecolor=COL["BC"], edgecolor="white", lw=0.4, zorder=4, label="Finding (all criteria met)")
    title = f"{title} ({f.sum()} findings)"
    q = r.p_model[r.q_model < C.FDR]
    if len(q):
        ax.axhline(-np.log10(q.max()), c=INK2, ls="--", lw=0.6)
        ax.text(x.max(), -np.log10(q.max()) - 0.08, f"FDR {C.FDR:.0%}", fontsize=5.5, color=INK2, ha="right", va="top")
    for i in r[f & r.compound_category.isin(["endogenous", "annotated (high score)"])].sort_values("p_model").index[:label_top]:
        ax.annotate(nm[i].split(" (")[0][:24], (x[i], y[i]), fontsize=5.3, xytext=(3, 2), textcoords="offset points")
    ax.axvline(0, c=GRID, lw=0.6, zorder=0)
    ax.set(xlabel="log$_2$ fold change (adjusted)", ylabel="−log$_{10}$ p (adjusted model)", title=title)


def fig3_volcano(out, stats_, nm):
    comps = [("BC_vs_Control", "BC vs Control"), ("BC_vs_PC", "BC vs PC"), ("PC_vs_Control", "PC vs Control")]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 6 * CM), gridspec_kw={"wspace": 0.32})
    for k, (ax, (c, t)) in enumerate(zip(axs, comps)):
        volcano_ax(ax, stats_[c], nm, t)
        panel(ax, "abc"[k], x=-0.22)
    h, l = axs[0].get_legend_handles_labels()
    for a in axs[1:]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh)
                l.append(ll)
    fig.legend(h, l, loc="upper center", ncol=5, bbox_to_anchor=(0.5, -0.02))
    save(fig, out, "Figure3_volcano")


def fig4_heatmap(out, L, meta, r, nm, case="BC", ref="Control"):
    sel = r[r.finding.astype(bool)].sort_values("beta_model")
    if len(sel) < 2:
        return
    smp = list(meta.index[meta.study & meta.group.isin([case, ref])])
    Z = L.loc[smp, sel.index]
    Z = ((Z - Z.mean()) / Z.std()).clip(-2.5, 2.5).T
    Z.index = [nm[i] for i in sel.index]
    age_n = Normalize(meta.age.min(), meta.age.max())
    ci = meta.loc[meta.study, "contam_index"]
    ci_n = Normalize(ci.quantile(0.02), ci.quantile(0.98))
    cmA, cmC = plt.get_cmap("Greys"), plt.get_cmap("magma_r")
    cc = pd.DataFrame({
        "Group": meta.loc[smp, "group"].map(COL),
        "Sex": meta.loc[smp, "sex"].map({"M": "#3a3a38", "F": "#d2d1cc"}).fillna("#ffffff"),
        "Age": [to_hex(cmA(age_n(a))) if pd.notna(a) else "#ffffff" for a in meta.loc[smp, "age"]],
        "Contamination": [to_hex(cmC(ci_n(v))) for v in meta.loc[smp, "contam_index"]],
    }, index=smp)
    order = meta.loc[smp].sort_values(["group", "contam_index"]).index
    g = sns.clustermap(Z[order], cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, row_cluster=False, col_cluster=False,
                       col_colors=cc.loc[order], xticklabels=False, yticklabels=True, figsize=(FULL, 0.33 * CM * len(sel) + 6 * CM),
                       dendrogram_ratio=(0.01, 0.01), colors_ratio=(0.02, 0.025), cbar_pos=None, linewidths=0)
    g.ax_heatmap.tick_params(axis="y", labelsize=5.5, length=0)
    g.ax_heatmap.set_xlabel("")
    g.ax_col_colors.tick_params(axis="y", labelsize=6)
    nb = (meta.loc[order, "group"] == case).sum()
    g.ax_heatmap.axvline(nb if meta.loc[order[0], "group"] == case else len(order) - nb, c="white", lw=1.2)
    h = [Patch(color=COL[x], label=LAB[x]) for x in [case, ref]] + \
        [Patch(color="#3a3a38", label="Male"), Patch(color="#d2d1cc", label="Female"),
         Patch(facecolor="white", edgecolor=INK2, lw=0.4, label="Not available")]
    g.ax_heatmap.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.13), ncol=5, fontsize=6)
    for i, (cm, lab, nr) in enumerate([(plt.get_cmap("RdBu_r"), "z-score", Normalize(-2.5, 2.5)),
                                       (cmA, "Age (years)", age_n), (cmC, "Contamination index", ci_n)]):
        pos = g.ax_heatmap.get_position()
        cax = g.fig.add_axes([pos.x0 + i * 0.22 * pos.width / 0.66, 0.02, 0.15, 0.008])
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=nr, cmap=cm), cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=6)
        cb.ax.tick_params(labelsize=5.5, length=1.5)
        cb.outline.set_linewidth(0.3)
    g.savefig(out / "Figure4_heatmap_BC_vs_Control.pdf", bbox_inches="tight")
    g.savefig(out / "Figure4_heatmap_BC_vs_Control.png", dpi=DPI, bbox_inches="tight")
    Image.open(out / "Figure4_heatmap_BC_vs_Control.png").convert("RGB").save(
        out / "Figure4_heatmap_BC_vs_Control.tif", compression="tiff_lzw", dpi=(DPI, DPI))
    plt.close(g.fig)


def fig5_roc(out, L, meta, feat, mv, pairs):
    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    cat = feat.loc[reps, "compound_category"].fillna("")
    rho = feat.loc[reps, "rho_contam_index"].abs()
    pool = list(reps[~features.panel_excluded(cat) & ~(rho >= C.PANEL_MAX_ABS_RHO_CONTAM)])   # как в run_all.py
    comps = [("BC_vs_Control", "BC vs Control"), ("BC_vs_PC", "BC vs PC"), ("PC_vs_Control", "PC vs Control")]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 6.4 * CM), gridspec_kw={"wspace": 0.3})
    keys = list(C.COMPARISONS)
    for n, (ax, (c, t)) in enumerate(zip(axs, comps)):
        k = keys.index(c)
        case, ref = st.comparison_sets(meta, C.COMPARISONS[c])
        pr = pairs[pairs.comparison == c]
        runs = [("all", "panel", case, ref, "-", "All"), ("all", "panel_adj", case, ref, "-.", "Contamination-adjusted")]
        if len(pr) >= C.MIN_PAIRS_CV:
            runs.append(("matched", "panel", list(pr.case), list(pr.ref), "--", f"Age-matched, {len(pr)} pairs"))
        else:
            ax.text(0.97, 0.03, f"Age-matched: {len(pr)} pairs\n(< {C.MIN_PAIRS_CV}, not estimated)", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=5.5, color=INK2)
        col = COL[C.COMPARISONS[c]["case"]]
        for sname, method, cs, rf, ls, lab in runs:
            y = np.r_[np.ones(len(cs)), np.zeros(len(rf))].astype(int)
            X = L.loc[cs + rf, pool].to_numpy()
            z = meta.loc[cs + rf, "contam_index"].to_numpy(float)
            aucs, tprs, _, _ = multivariate.repeated_cv(X, y, method, C.SEED + k, z)
            key = "скорректированы" if method == "panel_adj" else "логистическая"
            row = mv[(mv.comparison == c) & (mv.set == sname) & mv.model.astype(str).str.contains(key)]
            p = row.perm_p.iloc[0] if len(row) else np.nan
            if sname == "all" and method == "panel":
                ax.fill_between(multivariate.FPR, np.percentile(tprs, 10, 0), np.percentile(tprs, 90, 0), color=col, alpha=0.15, lw=0)
            ax.plot(multivariate.FPR, tprs.mean(0), c=col if method == "panel" else INK, ls=ls, lw=1.4 if method == "panel" else 1.0,
                    label=f"{lab}: AUC {aucs.mean():.2f} ± {aucs.std():.2f}, {pfmt(p)}")
        y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
        for cov, colr, nmc in [("contam_index", EXPO_COL["contaminant"], "Contamination index only"),
                               ("sex", INK2, "Sex only")]:
            v = meta.loc[case + ref, cov]
            v = v.map({"M": 1.0, "F": 0.0}) if cov == "sex" else v.astype(float)
            m = v.notna().to_numpy()
            if v[m].std() > 0:
                from sklearn.metrics import roc_auc_score, roc_curve
                a = roc_auc_score(y[m], v[m])
                fpr, tpr, _ = roc_curve(y[m], v[m] if a >= 0.5 else -v[m])
                ax.plot(fpr, tpr, c=colr, ls=":", lw=1, label=f"{nmc}: AUC {max(a, 1 - a):.2f}")
        
        ax.plot([0, 1], [0, 1], c=GRID, lw=0.6)
        ax.set(xlabel="1 − specificity", ylabel="Sensitivity", title=t, xlim=(-0.01, 1.01), ylim=(-0.01, 1.01), aspect="equal")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.3, handlelength=1.8)
        panel(ax, "abc"[n], x=-0.25)
    save(fig, out, "Figure5_ROC_cross_validated")


def fig6_box(out, L, meta, r, nm):
    sel = r[r.finding.astype(bool) & r.compound_category.isin(["endogenous", "annotated (high score)"])]
    if len(sel) < 8:
        sel = pd.concat([sel, r[r.finding.astype(bool) & ~r.index.isin(sel.index)].sort_values("p_model")])
    sel = sel.sort_values("p_model").head(8)
    if not len(sel):
        return
    fig, axs = plt.subplots(2, 4, figsize=(FULL, 10 * CM), gridspec_kw={"hspace": 0.9, "wspace": 0.45})
    for k, (ax, f) in enumerate(zip(axs.ravel(), sel.index)):
        for i, g in enumerate(C.STUDY_GROUPS):
            v = L.loc[meta.index[meta.study & (meta.group == g)], f]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color=INK, lw=0.9),
                       boxprops=dict(lw=0.6, color=INK2), whiskerprops=dict(lw=0.6, color=INK2), capprops=dict(lw=0.6, color=INK2))
            ax.scatter(np.random.RandomState(i).normal(i, 0.08, len(v)), v, s=5, c=COL[g], marker=MRK[g], edgecolor="none",
                       alpha=0.85, zorder=3)
        ax.set_xticks(range(3), [LAB[g] for g in C.STUDY_GROUPS])
        name, tag = nm[f].split(", m/z")[0], "m/z" + nm[f].split(", m/z")[1] if ", m/z" in nm[f] else ""
        ax.set_title(f"{name[:30]}\n{tag}\nq = {r.loc[f, 'q_model']:.1e}", fontsize=6)
        if k % 4 == 0:
            ax.set_ylabel("log$_2$ normalized intensity")
        panel(ax, "abcdefgh"[k], x=-0.3, y=1.2)
    for ax in axs.ravel()[len(sel):]:
        ax.axis("off")
    save(fig, out, "Figure6_top_findings_BC")


def figS1_age(out, meta, pairs):
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 6 * CM), gridspec_kw={"wspace": 0.3})
    s = meta[meta.study]
    ax = axs[0]
    for i, g in enumerate(C.STUDY_GROUPS):
        for sx, off, mk in [("M", -0.15, MRK[g]), ("F", 0.15, MRK[g])]:
            v = s.age[(s.group == g) & (s.sex == sx)].dropna()
            ax.scatter(np.random.RandomState(i).normal(i + off, 0.04, len(v)), v, s=10, c=COL[g] if sx == "M" else "white",
                       edgecolor=COL[g], lw=0.7, marker=mk, zorder=3)
        na = s.age[s.group == g].isna().sum()
        ax.text(i, s.age.min() - 6, f"age n/a: {na}", ha="center", fontsize=5.5, color=INK2)
    ax.set_xticks(range(3), [LAB[g] for g in C.STUDY_GROUPS])
    ax.set_ylabel("Age (years)")
    ax.set_ylim(s.age.min() - 9, s.age.max() + 3)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=INK2, mec=INK2, ms=4, label="Male"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, ms=4, label="Female")], loc="upper left")
    panel(ax, "a")
    ax = axs[1]
    for k, c in enumerate(["BC_vs_Control", "BC_vs_PC", "PC_vs_Control"]):
        pr = pairs[pairs.comparison == c]
        for _, r in pr.iterrows():
            ax.plot([k - 0.15, k + 0.15], [r.age_case, r.age_ref], c=GRID, lw=0.6, zorder=1)
        ca, re = C.COMPARISONS[c]["case"], C.COMPARISONS[c]["ref"]
        ax.scatter(np.full(len(pr), k - 0.15), pr.age_case, s=9, c=COL[ca], marker=MRK[ca], zorder=3)
        ax.scatter(np.full(len(pr), k + 0.15), pr.age_ref, s=9, c=COL[re], marker=MRK[re], zorder=3)
        ax.text(k, ax.get_ylim()[0] if False else 22, f"{len(pr)} pairs", ha="center", fontsize=6, color=INK2)
    ax.set_xticks(range(3), ["BC vs Control", "BC vs PC", "PC vs Control"])
    ax.set_ylabel("Age (years)")
    ax.set_ylim(18, 92)
    ax.set_title(f"1:1 nearest-age matching within sex (±{C.MATCH_CALIPER_YEARS} years)")
    panel(ax, "b")
    save(fig, out, "FigureS1_age_and_matching")


def figS2_norm(out, meta, res_dir):
    s = meta[meta.study]
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 6 * CM), gridspec_kw={"wspace": 0.3})
    ax = axs[0]
    scatter_groups(ax, s.dilution_PQN.to_numpy(), s.dilution_MSTUS.to_numpy(), s.group.to_numpy(), order=C.STUDY_GROUPS)
    ax.set(xscale="log", yscale="log", xlabel="PQN dilution factor", ylabel="MSTUS dilution factor")
    for i in s.index[s.dilution_PQN < C.LOW_DILUTION_FACTOR]:
        ax.annotate(i, (s.dilution_PQN[i], s.dilution_MSTUS[i]), fontsize=5.5, xytext=(2, 2), textcoords="offset points")
    ax.axvline(C.LOW_DILUTION_FACTOR, c=INK2, ls=":", lw=0.6)
    rho = stats.spearmanr(s.dilution_PQN, s.dilution_MSTUS).statistic
    ax.text(0.02, 0.98, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6)
    panel(ax, "a")
    ax = axs[1]
    ns = pd.read_csv(res_dir / "tables" / "normalization_sensitivity.csv")
    x = np.arange(len(ns))
    w = 0.25
    for j, (col, lab, c) in enumerate([("q<FDR_PQN", "PQN", INK2), ("q<FDR_MSTUS", "MSTUS", "#bdbcb8"),
                                       ("overlap", "Overlap", COL["Control"])]):
        ax.bar(x + (j - 1) * w, ns[col], w * 0.9, color=c, label=lab)
    ax.set_xticks(x, [c.replace("_vs_", " vs ").replace("maleControl", "male Control") for c in ns.comparison], rotation=15)
    ax.set_ylabel(f"Compounds with q < {C.FDR}")
    ax.legend()
    panel(ax, "b")
    save(fig, out, "FigureS2_normalization")


def figS3_pathways(out, pw):
    if pw is None:
        return
    comps = list(pw.comparison.unique())
    fig, axs = plt.subplots(len(comps), 1, figsize=(HALF * 1.45, 6.2 * CM * len(comps)), gridspec_kw={"hspace": 0.45})
    for k, (ax, c) in enumerate(zip(np.atleast_1d(axs), comps)):
        d = pw[pw.comparison == c].sort_values("p-value").head(10).iloc[::-1]
        ax.barh(d.pathway, -np.log10(d["p-value"]), color=[COL["BC"] if q < 0.05 else "#bdbcb8" for q in d.q_BH], height=0.7)
        for yy, q in enumerate(d.q_BH):
            ax.text(-np.log10(d["p-value"].iloc[yy]) + 0.03, yy, f"q = {q:.2f}", va="center", fontsize=5.3, color=INK2)
        ax.axvline(-np.log10(0.05), c=INK2, ls="--", lw=0.6)
        ax.set(xlabel="−log$_{10}$ p (mummichog)", title=c.replace("_vs_", " vs "))
        ax.tick_params(axis="y", labelsize=5.5)
        panel(ax, "abc"[k], x=-0.75, y=1.02)
    save(fig, out, "FigureS3_pathways")


# ------------------------------------------------------------------ таблицы
def tables(out, demo, demo_t, stats_, nm, cur, mv, feat):
    t1 = demo.rename(columns={"group": "Group"})
    t1["Group"] = t1.Group.map(LAB)
    t1["Sex, M/F (n/a)"] = t1.male.astype(str) + "/" + t1.female.astype(str) + " (" + t1.sex_unknown.astype(str) + ")"
    t1["Age, mean ± SD"] = t1.age_mean.round(1).astype(str) + " ± " + t1.age_sd.round(1).astype(str)
    t1["Age, median (range)"] = t1.age_median.astype(int).astype(str) + " (" + t1.age_min.astype(int).astype(str) + "–" + t1.age_max.astype(int).astype(str) + ")"
    t1["Age available, n"] = t1.age_n
    t1 = t1[["Group", "n", "Sex, M/F (n/a)", "Age, mean ± SD", "Age, median (range)", "Age available, n"]]
    tt = demo_t.copy()
    tt["test"] = tt.test.str.replace("возраст", "age").str.replace("пол", "sex").str.replace(", Манн–Уитни", ", Mann–Whitney U") \
        .str.replace(", точный тест Фишера", ", Fisher's exact").str.replace(", Краскел–Уоллис (3 группы)", ", Kruskal–Wallis (3 groups)") \
        .str.replace("мужчины-контроли", "male controls")
    tt["p"] = tt.p.map(lambda p: f"{p:.2g}")
    r = stats_["BC_vs_Control"]
    f = r[r.finding.astype(bool)].copy()
    cn = cur.dropna(subset=["feature"]).set_index("feature")
    t2 = pd.DataFrame({
        "Compound (annotation)": [nm[i].split(", m/z")[0] for i in f.index],
        "m/z": f.mz.round(4), "RT (min)": f.rt_min.round(2),
        "Putative library match (low score)": [feat.loc[i, "annotation"] if feat.loc[i, "annotation_conf"] == "low" else "" for i in f.index],
        "Ions in group": f.n_ions,
        "MSI level": [int(cn.loc[i, "MSI_level"]) if i in cn.index else (3 if feat.loc[i, "annotation_conf"] == "high" else 4) for i in f.index],
        "log2FC (adjusted)": f.beta_model.round(2), "q (adjusted model)": f.q_model.map(lambda x: f"{x:.1e}"),
        "AUC (all)": f.AUC.round(2), "p (age-matched)": f.p_matched.map(lambda x: f"{x:.3f}"),
        "p (low-contamination controls)": f.p_lowindex.map(lambda x: f"{x:.3f}"),
        "Detected, % Control/BC/PC": [f"{100 * f.det_Control[i]:.0f}/{100 * f.det_BC[i]:.0f}/{100 * f.det_PC[i]:.0f}" for i in f.index],
    }).sort_values("q (adjusted model)")
    m = mv[mv.model.str.startswith(("панель", "PLS-DA"))].copy()
    m["model"] = m.model.str.replace(r"панель top-(\d+) \+ логистическая регрессия \(отбор внутри CV\)", r"Top-\1 panel + logistic regression (in-fold selection)", regex=True) \
        .str.replace(r"панель top-(\d+), признаки скорректированы на индекс загрязнения внутри фолда", r"Top-\1 panel, features adjusted for contamination index in-fold", regex=True) \
        .str.replace(r"PLS-DA, 2 компоненты \(все признаки без экспозиций\)", "PLS-DA, 2 components (all non-exposure compounds)", regex=True)
    ts1 = m[["comparison", "set", "model", "n_case", "n_ref", "AUC_mean", "AUC_sd", "Q2", "perm_p"]].round(3)
    with pd.ExcelWriter(out / "Tables.xlsx") as w:
        t1.to_excel(w, sheet_name="Table1_demographics", index=False)
        tt.to_excel(w, sheet_name="Table1_tests", index=False)
        t2.to_excel(w, sheet_name="Table2_findings_BC", index=False)
        ts1.to_excel(w, sheet_name="TableS1_classification", index=False)
    t1.to_csv(out / "Table1_demographics.csv", index=False)
    tt.to_csv(out / "Table1_tests.csv", index=False)
    t2.to_csv(out / "Table2_findings_BC_vs_Control.csv", index=False)
    ts1.to_csv(out / "TableS1_classification.csv", index=False)


CAPTIONS = """# Figure captions (draft)

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the {rsd:.0%} acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound). (d) The same scores coloured by the contamination index.

**Figure 2.** Plastic-derived surfactant contamination. (a) Homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+) used to build the contamination index; open diamonds, co-varying ions not included in the index. (b) Contamination index by group (Mann–Whitney U test).

**Figure 3.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + contamination index). Filled circles, findings that passed all robustness criteria; crosses, exposure-related compounds (contaminants, food, drugs). Dashed line, FDR {fdr:.0%}.

**Figure 4.** Heatmap of bladder-cancer findings (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and contamination index.

**Figure 5.** Cross-validated ROC curves ({rep}×{fold}-fold CV; top-{k} features selected by Mann–Whitney U inside each training fold, logistic regression). Candidate pool defined before CV: exposure-related compounds (contaminants, food, drugs, tobacco) and compounds correlated with the contamination index (|ρ| ≥ {rho}) excluded. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot, features residualized on the contamination index inside each training fold; dashed, 1:1 age-matched pairs (shown only with ≥ {minp} pairs). p, permutation test ({perm} permutations). Dotted lines, single-covariate baselines (sex; contamination index).

**Figure 6.** Top bladder-cancer findings across groups (log2 PQN-normalized intensities). q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure S1.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±{cal} years).

**Figure S2.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < {low}. (b) Number of compounds with q < {fdr} under PQN and MSTUS normalization and their overlap.

**Figure S3.** Pathway enrichment (mummichog, positive mode, {ppm:.0f} ppm; exposure-related features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default=C.OUT)
    ap.add_argument("--data", default=C.DATA)
    a = ap.parse_args()
    res = Path(a.res)
    out = res / "publication"
    out.mkdir(parents=True, exist_ok=True)
    meta, L, info, feat, stats_, mv, series, cur, pairs, demo, demo_t, pw = load(res)
    X_all, _, _, _ = io.load_features(a.data)
    nm = names_en(feat, cur)
    fig1_quality(out, meta, L, feat, X_all)
    fig2_contamination(out, meta, series)
    fig3_volcano(out, stats_, nm)
    fig4_heatmap(out, L, meta, stats_["BC_vs_Control"], nm)
    fig5_roc(out, L, meta, feat, mv, pairs)
    fig6_box(out, L, meta, stats_["BC_vs_Control"], nm)
    figS1_age(out, meta, pairs)
    figS2_norm(out, meta, res)
    figS3_pathways(out, pw)
    tables(out, demo, demo_t, stats_, nm, cur, mv, feat)
    (out / "figure_captions.md").write_text(CAPTIONS.format(
        rsd=C.MAX_QC_RSD, fdr=C.FDR, rep=C.CV_REPEATS, fold=C.CV_FOLDS, k=C.PANEL_K, perm=C.N_PERM,
        cal=C.MATCH_CALIPER_YEARS, rho=C.PANEL_MAX_ABS_RHO_CONTAM, minp=C.MIN_PAIRS_CV, low=C.LOW_DILUTION_FACTOR, ppm=C.PPM), encoding="utf-8")
    print(f"Публикационные рисунки и таблицы → {out}")


if __name__ == "__main__":
    main()

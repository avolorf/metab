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
import re
import textwrap
import shutil
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


def names_en(feat, cur, msi=None):
    """Подпись соединения: курированное англ. название (MSI) > уверенная аннотация > unknown m/z@RT."""
    out = {}
    cn = cur.dropna(subset=["feature"]).set_index("feature")
    for f, r in feat.iterrows():
        tag = f"m/z {r.mz:.4f} @ {r.rt_min:.2f} min"
        if f in cn.index and isinstance(cn.loc[f, "name_en"], str) and cn.loc[f, "name_en"]:
            lv = int(msi[f]) if msi is not None else int(cn.loc[f, "MSI_level"])
            out[f] = f"{cn.loc[f, 'name_en']} (MSI {lv}), {tag}"
        elif r.annotation_conf == "high":
            out[f] = f"{r.annotation} (putative), {tag}"
        else:
            out[f] = f"unknown, {tag}"
    return pd.Series(out)


def msi_levels(feat, res_dir):
    """Уровень MSI для каждого признака (см. правило в pipeline/annotate.msi_level), на уровне соединения:
    совпадение MS2 засчитывается, если оно есть у любого иона соединения."""
    t = res_dir / "tables"
    ok = pd.read_csv(t / "candidates_curated.csv")
    bad = pd.read_csv(t / "candidates_mass_mismatch.csv") if (t / "candidates_mass_mismatch.csv").exists() else pd.DataFrame(columns=["feature"])
    cfile = pd.read_csv(C.CURATION, dtype=str).fillna("")
    amb = set(cfile.loc[cfile.get("isomer_ambiguous", "").isin(["1", "yes", "да", "true"]), "name"]) if "isomer_ambiguous" in cfile else set()
    ms2f = set(pd.read_csv(C.MS2_IDENTIFIED).feature) if Path(C.MS2_IDENTIFIED).exists() else set()
    comp = feat["compound"]
    ms2_comp = set(comp[comp.index.isin(ms2f)].dropna())
    has_ms2 = comp.isin(ms2_comp)
    lvl = pd.Series(4, index=feat.index)
    lvl[has_ms2] = 3
    for _, r in ok.dropna(subset=["feature"]).iterrows():
        if r.feature in lvl.index:
            ms2 = bool(has_ms2.get(r.feature, False))
            lvl[r.feature] = 1 if str(r.get("standard_confirmed", "")).strip().lower() in {"1", "yes", "да", "true"} else \
                (2 if ms2 and r["name"] not in amb else 3)
    for f in bad.feature.dropna():
        if f in lvl.index:
            lvl[f] = 4
    return lvl, has_ms2


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


def table2_order(r, msi):
    """Соединения Table 2 (находки MSI 2–3) в порядке возрастания q — общая нумерация для Table 2 и Figure 3."""
    sel = r[r.finding.astype(bool) & (msi.reindex(r.index) <= 3)].sort_values("q_model")
    return {f: k for k, f in enumerate(sel.index, 1)}


def volcano_ax(ax, r, nm, title, numbered=None):
    from adjustText import adjust_text
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
        ax.axhline(-np.log10(q.max()), c=INK2, ls="--", lw=0.6, label=f"FDR {C.FDR:.0%}")
    if numbered:
        idx = sorted([i for i in numbered if i in r.index], key=lambda i: -y[i])
        ax.scatter(x[idx], y[idx], s=30, facecolor="none", edgecolor=INK, lw=0.8, zorder=5,
                   label="Annotated finding (Table 2 no.)")
        # номера — отдельным столбцом справа от облака, каждый с выносной линией к своей точке
        xr = x.max() + 0.25 * (x.max() - x.min())
        ys = np.linspace(y.max(), y.max() - 0.42 * (y.max() - y.min()), len(idx))
        for (i, yl) in zip(idx, ys):
            ax.plot([x[i], xr - 0.05 * (x.max() - x.min())], [y[i], yl], c=INK, lw=0.4, zorder=5)
            ax.text(xr, yl, str(numbered[i]), fontsize=6, fontweight="bold", va="center", ha="left", zorder=6)
        ax.set_xlim(right=xr + 0.12 * (x.max() - x.min()))
    ax.axvline(0, c=GRID, lw=0.6, zorder=0)
    ax.set(xlabel="log$_2$ fold change (adjusted)", ylabel="−log$_{10}$ p (adjusted model)", title=title)


def fig3_volcano(out, stats_, nm, msi):
    comps = [("BC_vs_Control", "BC vs Control"), ("BC_vs_PC", "BC vs PC"), ("PC_vs_Control", "PC vs Control")]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 6 * CM), gridspec_kw={"wspace": 0.32})
    for k, (ax, (c, t)) in enumerate(zip(axs, comps)):
        volcano_ax(ax, stats_[c], nm, t, table2_order(stats_[c], msi) if c == "BC_vs_Control" else None)
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
    sel = r[r.finding.astype(bool)].sort_values("q_model").head(50).sort_values("beta_model")   # не более 50 строк
    if len(sel) < 2:
        return
    smp = list(meta.index[meta.study & meta.group.isin([case, ref])])
    Z = L.loc[smp, sel.index]
    Z = ((Z - Z.mean()) / Z.std()).clip(-2.5, 2.5).T
    def short_lab(i):
        name, _, tag = nm[i].partition(", m/z")
        name = name.replace(" (putatively 11-oxoetiocholanolone)", "")
        return f"{name[:40]}, m/z{tag}"
    Z.index = [short_lab(i) for i in sel.index]
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
        cax = g.fig.add_axes([pos.x0 + i * pos.width / 3, pos.y0 - 0.035, pos.width / 3.6, 0.007])
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=nr, cmap=cm), cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=6)
        cb.ax.tick_params(labelsize=5.5, length=1.5)
        cb.outline.set_linewidth(0.3)
    g.savefig(out / "FigureS1_heatmap_BC_vs_Control.pdf", bbox_inches="tight")
    g.savefig(out / "FigureS1_heatmap_BC_vs_Control.png", dpi=DPI, bbox_inches="tight")
    Image.open(out / "FigureS1_heatmap_BC_vs_Control.png").convert("RGB").save(
        out / "FigureS1_heatmap_BC_vs_Control.tif", compression="tiff_lzw", dpi=(DPI, DPI))
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
        runs = [("all", "panel", case, ref, "-", "All"), ("all", "panel_adj", case, ref, "-.", "Contam.-adjusted")]
        if len(pr) >= C.MIN_PAIRS_CV:
            runs.append(("matched", "panel", list(pr.case), list(pr.ref), "--", f"Age-matched ({len(pr)} pairs)"))
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
                    label=f"{lab}: {aucs.mean():.2f} ± {aucs.std():.2f}, {pfmt(p)}")
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
                ax.plot(fpr, tpr, c=colr, ls=":", lw=1, label=f"{nmc}: {max(a, 1 - a):.2f}")
        
        ax.plot([0, 1], [0, 1], c=GRID, lw=0.6)
        ax.set(xlabel="1 − specificity", ylabel="Sensitivity", title=f"{t} (AUC)", xlim=(-0.01, 1.01), ylim=(-0.01, 1.01), aspect="equal")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.3, handlelength=1.8)
        panel(ax, "abc"[n], x=-0.25)
    save(fig, out, "Figure5_ROC_cross_validated")


def fig6_box(out, L, meta, r, nm, msi):
    """Figure 4: те же соединения, что в Table 2 (находки MSI 2–3), по возрастанию q."""
    sel = r[r.finding.astype(bool) & (msi.reindex(r.index) <= 3)]
    sel = sel.sort_values("q_model").head(8)
    if not len(sel):
        return
    fig, axs = plt.subplots(2, 4, figsize=(FULL, 11 * CM), gridspec_kw={"hspace": 1.0, "wspace": 0.5})
    for k, (ax, f) in enumerate(zip(axs.ravel(), sel.index)):
        for i, g in enumerate(C.STUDY_GROUPS):
            v = L.loc[meta.index[meta.study & (meta.group == g)], f]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color=INK, lw=0.9),
                       boxprops=dict(lw=0.6, color=INK2), whiskerprops=dict(lw=0.6, color=INK2), capprops=dict(lw=0.6, color=INK2))
            ax.scatter(np.random.RandomState(i).normal(i, 0.08, len(v)), v, s=5, c=COL[g], marker=MRK[g], edgecolor="none",
                       alpha=0.85, zorder=3)
        ax.set_xticks(range(3), [LAB[g] for g in C.STUDY_GROUPS])
        name, tag = nm[f].split(", m/z")[0], "m/z" + nm[f].split(", m/z")[1] if ", m/z" in nm[f] else ""
        name = re.sub(r" \((putatively[^)]*|MSI \d|putative)\)", "", name)
        name = "\n".join(textwrap.wrap(name, 24, break_on_hyphens=True)[:2])
        mz, rt = tag.replace("m/z ", "").replace(" min", "").split(" @ ")
        ax.set_title(f"{name}\nMSI {int(msi[f])} · m/z {mz} · {rt} min\nq = {r.loc[f, 'q_model']:.1e}", fontsize=5.6)
        if k % 4 == 0:
            ax.set_ylabel("log$_2$ normalized intensity")
        panel(ax, "abcdefgh"[k], x=-0.3, y=1.2)
    for ax in axs.ravel()[len(sel):]:
        ax.axis("off")
    save(fig, out, "Figure4_top_findings_BC")


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
    save(fig, out, "FigureS2_age_and_matching")


def figS2_norm(out, meta, res_dir):
    s = meta[meta.study]
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 6 * CM), gridspec_kw={"wspace": 0.3})
    ax = axs[0]
    scatter_groups(ax, s.dilution_PQN.to_numpy(), s.dilution_MSTUS.to_numpy(), s.group.to_numpy(), order=C.STUDY_GROUPS)
    ax.set(xscale="log", yscale="log", xlabel="PQN dilution factor", ylabel="MSTUS dilution factor")
    from adjustText import adjust_text
    low = s.index[s.dilution_PQN < C.LOW_DILUTION_FACTOR]
    texts = [ax.text(s.dilution_PQN[i], s.dilution_MSTUS[i], i, fontsize=5.5) for i in low]
    adjust_text(texts, x=s.dilution_PQN.to_numpy(), y=s.dilution_MSTUS.to_numpy(), ax=ax, expand=(1.4, 1.6),
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.4))
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
    save(fig, out, "FigureS3_normalization")


def figS3_pathways(out, pw):
    if pw is None:
        return
    comps = list(pw.comparison.unique())
    fig, axs = plt.subplots(len(comps), 1, figsize=(HALF * 1.45, 6.2 * CM * len(comps)), gridspec_kw={"hspace": 0.45})
    for k, (ax, c) in enumerate(zip(np.atleast_1d(axs), comps)):
        d = pw[pw.comparison == c].sort_values("p-value").head(10).iloc[::-1]
        ax.barh(d.pathway, -np.log10(d["p-value"]), color=[COL["BC"] if q < 0.05 else "#bdbcb8" for q in d.q_BH], height=0.7)
        xmax = (-np.log10(d["p-value"])).max() * 1.05
        ax.set_xlim(0, xmax)
        for yy, q in enumerate(d.q_BH):     # q — отдельным столбцом справа от осей
            ax.text(1.02, yy, f"{q:.2f}", va="center", fontsize=5.3, color=INK2, transform=ax.get_yaxis_transform(),
                    clip_on=False)
        ax.text(1.02, len(d) - 0.3, "q (BH)", fontsize=5.5, color=INK2, transform=ax.get_yaxis_transform(), clip_on=False)
        ax.axvline(-np.log10(0.05), c=INK2, ls="--", lw=0.6)
        ax.set(xlabel="−log$_{10}$ p (mummichog)", title=c.replace("_vs_", " vs "))
        ax.tick_params(axis="y", labelsize=5.5)
        panel(ax, "abc"[k], x=-0.75, y=1.02)
    save(fig, out, "FigureS4_pathways")


# ------------------------------------------------------------------ таблицы
def tables(out, demo, demo_t, stats_, nm, cur, mv, feat, msi=None, has_ms2=None):
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
        "MS2 library match": ["yes" if has_ms2 is not None and has_ms2.get(i, False) else "no" for i in f.index],
        "MSI level": [int(msi[i]) for i in f.index],
        "log2FC (adjusted)": f.beta_model.round(2), "q (adjusted model)": f.q_model.map(lambda x: f"{x:.1e}"),
        "AUC (all)": f.AUC.round(2), "p (age-matched)": f.p_matched.map(lambda x: f"{x:.3f}"),
        "p (low-contamination controls)": f.p_lowindex.map(lambda x: f"{x:.3f}"),
        "Detected, % Control/BC/PC": [f"{100 * f.det_Control[i]:.0f}/{100 * f.det_BC[i]:.0f}/{100 * f.det_PC[i]:.0f}" for i in f.index],
    })
    t2 = t2.loc[f.sort_values("q_model").index]          # сортировка по числовому q
    m = mv[mv.model.str.startswith(("панель", "PLS-DA"))].copy()
    m["model"] = m.model.str.replace(r"панель top-(\d+) \+ логистическая регрессия \(отбор внутри CV\)", r"Top-\1 panel + logistic regression (in-fold selection)", regex=True) \
        .str.replace(r"панель top-(\d+), признаки скорректированы на индекс загрязнения внутри фолда", r"Top-\1 panel, features adjusted for contamination index in-fold", regex=True) \
        .str.replace(r"PLS-DA, 2 компоненты \(все признаки без экспозиций\)", "PLS-DA, 2 components (all non-exposure compounds)", regex=True)
    ts1 = m[["comparison", "set", "model", "n_case", "n_ref", "AUC_mean", "AUC_sd", "Q2", "perm_p"]].round(3)
    # Table 2 (основной текст): только аннотированные находки (MSI ≤ 3); полный список — Table S1
    t2_main = t2[t2["MSI level"] <= 3].copy()
    t2_main["Compound (annotation)"] = t2_main["Compound (annotation)"].str.replace(r" \((MSI \d|putative)\)$", "", regex=True)
    t2_main.insert(1, "Library match (MS-DIAL)", feat.loc[t2_main.index, "annotation"].to_numpy())
    t2_main = t2_main.drop(columns=["Putative library match (low score)"])
    order = table2_order(r, msi)
    t2_main.insert(0, "No.", [order[i] for i in t2_main.index])
    t2_main = t2_main.sort_values("No.")
    with pd.ExcelWriter(out / "Tables_main.xlsx") as w:
        t1.to_excel(w, sheet_name="Table1_demographics", index=False)
        tt.to_excel(w, sheet_name="Table1_tests", index=False)
        t2_main.to_excel(w, sheet_name="Table2_annotated_findings_BC", index=False)
    t1.to_csv(out / "Table1_demographics.csv", index=False)
    tt.to_csv(out / "Table1_tests.csv", index=False)
    t2_main.to_csv(out / "Table2_annotated_findings_BC.csv", index=False)
    t2.to_csv(out / "TableS1_all_findings_BC.csv", index=False)
    ts1.to_csv(out / "TableS2_classification.csv", index=False)
    return t2, ts1


CAPTIONS = """# Figure and table captions (draft)

## Main text

**Table 1.** Clinical and demographic characteristics of the independent cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the {rsd:.0%} acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound). (d) The same scores coloured by the contamination index.

**Figure 2.** Plastic-derived surfactant contamination. (a) Homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+) used to build the contamination index; open diamonds, co-varying ions not included in the index. (b) Contamination index by group (Mann–Whitney U test).

**Figure 3.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + contamination index). Filled circles, findings that passed all robustness criteria; crosses, exposure-related compounds (contaminants, food, drugs); dashed line, FDR {fdr:.0%}. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3). log2FC and q, adjusted model (sex, age, contamination index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-contamination controls), BC vs controls with a contamination index below the median. MSI 2, precursor within {ppm:.0f} ppm of the curated formula and MS/MS library match without isomeric ambiguity; MSI 3, formula-consistent without isomer-specific MS/MS. For compounds without a curated annotation, the MS-DIAL library match is given for transparency; ethyl 2-amino-5-methyl-3-thiophenecarboxylate is a synthetic building block and this annotation is considered implausible. All findings, including MSI level 4, are listed in Table S1.

**Figure 4.** Annotated bladder-cancer findings (MSI levels 2–3; the compounds of Table 2) across groups (log2 PQN-normalized intensities). q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 5.** Cross-validated ROC curves ({rep}×{fold}-fold CV; top-{k} features selected by Mann–Whitney U inside each training fold, logistic regression). Candidate pool defined before CV: exposure-related compounds (contaminants, food, drugs, tobacco) and compounds correlated with the contamination index (|ρ| ≥ {rho}) excluded. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot, features residualized on the contamination index inside each training fold; dashed, 1:1 age-matched pairs (shown only with ≥ {minp} pairs). p, permutation test ({perm} permutations). Dotted lines, single-covariate baselines (sex; contamination index).

**Figure 6.** Stage-1 bladder-cancer markers (Frolova et al., Anal Bioanal Chem 2026, Table 5) in the independent cohort. Features were matched by m/z (±{ppm:.0f} ppm) and retention time (±{rttol} min); all isomeric candidates are shown. Points, unadjusted log2 fold change (BC vs Control) of quantified features; coloured, same direction and Mann–Whitney p < 0.05; text gives AUC, unadjusted p and FDR-adjusted q of the adjusted model. Markers not detected or removed by the detection/QC filters are labelled. Exact values are given in Table S3.

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and contamination index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±{cal} years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < {low}. (b) Number of compounds with q < {fdr} under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, {ppm:.0f} ppm; exposure-related features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, Q², permutation p).

**Table S3.** Stage-1 markers: all candidate features within the m/z and retention-time windows, detection/filter status and statistics in the independent cohort.

**Table S4.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, supporting ions, MSI level, detection and imputation by group, correlation with the contamination index and the dilution factor; annotations with mass error > {ppm:.0f} ppm are listed separately.

**Table S5.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S6.** Compound-level statistics for all comparisons (all models and robustness checks).
"""


def replication_outputs(out, res_dir, feat, stats_, info_all):
    """Table 3 и Figure 7: маркеры первого этапа (ABC, Table 5) в независимой выборке s266."""
    from pipeline import replication
    tab, summ = replication.run(C.PRIOR_MARKERS, info_all, feat, stats_, lambda *_: None)
    if tab is None:
        return None
    st_lab = {"not detected": "Not detected", "below detection filter": "Detected, < detection filter",
              "failed QC RSD filter": "Detected, QC RSD > 30%", "quantified": "Quantified"}
    t3 = pd.DataFrame({
        "Marker (stage 1)": tab.marker, "Formula": tab.formula, "Ion": tab.ion,
        "m/z (stage 1)": tab.mz, "RT, min (stage 1)": tab.rt_min,
        "Candidate": [f"{int(c)}/{int(n)}" if n else "" for c, n in zip(tab.candidate, tab.n_candidates)],
        "m/z (this study)": tab.get("mz_obs").round(4), "RT, min (this study)": tab.get("rt_obs").round(2),
        "ΔRT, min": tab.get("delta_rt"), "Status (this study)": tab.status.map(st_lab),
        "log2FC BC/Control": tab.get("log2FC").round(2), "AUC": tab.get("AUC").round(2),
        "p (Mann–Whitney)": tab.get("p_MW").map(lambda x: f"{x:.2g}" if pd.notna(x) else ""),
        "q (adjusted model)": tab.get("q_model").map(lambda x: f"{x:.2g}" if pd.notna(x) else ""),
        "Same direction": tab.get("same_direction").map({True: "yes", False: "no"}),
        "Note on stage-1 entry": tab.comment.str.extract(r";\s*(.*)$")[0].fillna(""),
    })
    t3.to_csv(out / "TableS3_replication_stage1_markers.csv", index=False)
    summ.to_csv(out / "TableS3_replication_summary.csv", index=False)
    # Figure 7
    d = tab.reset_index(drop=True)
    d["srank"] = d.status.map({"quantified": 0, "failed QC RSD filter": 1, "below detection filter": 2, "not detected": 3})
    order = d.sort_values(["srank", "marker", "candidate"], ascending=[False, True, False]).index
    labs, hts = [], []
    for i in order:     # название переносится, координаты (m/z, RT) — отдельной строкой; высота строки по числу строк
        r = d.loc[i]
        cand = f" [{int(r.candidate)}/{int(r.n_candidates)}, {r.rt_obs:.2f} min]" if r.n_candidates > 1 else ""
        lines = textwrap.wrap(str(r.marker), 40, break_on_hyphens=True)
        coord = f"({r.mz:.4f}, {r.rt_min:.2f}){cand}"
        lines = lines[:-1] + [f"{lines[-1]} {coord}"] if len(lines) == 1 and len(lines[-1]) + len(coord) <= 48 else lines + [coord]
        labs.append("\n".join(lines))
        hts.append(1 + 0.45 * (len(lines) - 1))
    ys = np.cumsum([0] + [(a + b) / 2 for a, b in zip(hts[:-1], hts[1:])])
    fig, ax = plt.subplots(figsize=(FULL * 0.55, 0.42 * CM * sum(hts) + 2.4 * CM))
    for y, i, lab in zip(ys, order, labs):
        r = d.loc[i]
        ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6, linespacing=1.05)
        if r.status == "quantified":
            col = COL["BC"] if r.nominal else "#bdbcb8"
            ax.scatter(r.log2FC, y, s=26, c=col, marker="o", edgecolor=INK2, lw=0.5, zorder=3)
            ax.text(r.log2FC + 0.08, y + 0.18, f"AUC {r.AUC:.2f}, p {r.p_MW:.2g}, q {r.q_model:.2g}", fontsize=5.3, color=INK2)
        else:
            ax.text(0.02, y, st_lab[r.status], fontsize=5.8, color=INK2, va="center", style="italic",
                    transform=ax.get_yaxis_transform())
    ax.axvline(0, c=GRID, lw=0.8)
    ax.set_yticks([])
    ax.set_ylim(-0.7, ys[-1] + hts[-1] / 2 + 0.2)
    ax.set_xlabel("log$_2$ fold change BC/Control in the independent cohort (unadjusted)")
    ax.spines["left"].set_visible(False)
    q = summ.iloc[0]
    ax.set_title(f"Stage-1 markers (n = {int(q.markers)}): quantified {int(q.quantified)} ({int(q.ambiguous_isomers)} with isomers),\n"
                 f"p < 0.05 in the same direction {int(q.nominal_p05_any_candidate)}, q < {C.FDR} (adjusted) {int(q.replicated_q05_any_candidate)}",
                 fontsize=6.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=COL["BC"], mec=INK2, label="same direction, p < 0.05"),
                       Line2D([], [], marker="o", ls="", mfc="#bdbcb8", mec=INK2, label="not significant")],
              loc="lower right", fontsize=6)
    save(fig, out, "Figure6_replication_stage1_markers")
    return t3, summ


def supplementary_tables(out, res_dir, t2_full, ts_class, rep):
    """Supplementary_Tables.xlsx: S1 все находки РМП, S2 классификация, S3 сверка с первым этапом,
    S4 курированные кандидаты (+ несходящаяся масса), S5 признаки «есть/нет», S6 статистика по соединениям."""
    t = res_dir / "tables"
    with pd.ExcelWriter(out / "Supplementary_Tables.xlsx") as w:
        t2_full.to_excel(w, sheet_name="S1_all_findings_BC", index=False)
        ts_class.to_excel(w, sheet_name="S2_classification", index=False)
        if rep is not None:
            rep[0].to_excel(w, sheet_name="S3_stage1_markers", index=False)
            rep[1].to_excel(w, sheet_name="S3_summary", index=False)
        for sh, fn_ in {"S4_curated_candidates": "candidates_curated.csv", "S4b_mass_mismatch": "candidates_mass_mismatch.csv",
                        "S5_presence_absence": "presence_absence.csv"}.items():
            if (t / fn_).exists():
                pd.read_csv(t / fn_).to_excel(w, sheet_name=sh, index=False)
        for n in C.COMPARISONS:
            if (t / f"stats_{n}.csv").exists():
                pd.read_csv(t / f"stats_{n}.csv").to_excel(w, sheet_name=f"S6_{n}"[:31], index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default=C.OUT)
    ap.add_argument("--data", default=C.DATA)
    a = ap.parse_args()
    res = Path(a.res)
    out = res / "publication"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    meta, L, info, feat, stats_, mv, series, cur, pairs, demo, demo_t, pw = load(res)
    X_all, _, info_all, _ = io.load_features(a.data)
    msi, has_ms2 = msi_levels(feat, res)
    nm = names_en(feat, cur, msi)
    fig1_quality(out, meta, L, feat, X_all)
    fig2_contamination(out, meta, series)
    fig3_volcano(out, stats_, nm, msi)
    fig4_heatmap(out, L, meta, stats_["BC_vs_Control"], nm)
    fig5_roc(out, L, meta, feat, mv, pairs)
    fig6_box(out, L, meta, stats_["BC_vs_Control"], nm, msi)
    figS1_age(out, meta, pairs)
    figS2_norm(out, meta, res)
    figS3_pathways(out, pw)
    t2_full, ts_class = tables(out, demo, demo_t, stats_, nm, cur, mv, feat, msi, has_ms2)
    rep = replication_outputs(out, res, feat, stats_, info_all)
    supplementary_tables(out, res, t2_full, ts_class, rep)
    (out / "figure_captions.md").write_text(CAPTIONS.format(
        rsd=C.MAX_QC_RSD, fdr=C.FDR, rep=C.CV_REPEATS, fold=C.CV_FOLDS, k=C.PANEL_K, perm=C.N_PERM,
        cal=C.MATCH_CALIPER_YEARS, rho=C.PANEL_MAX_ABS_RHO_CONTAM, minp=C.MIN_PAIRS_CV, low=C.LOW_DILUTION_FACTOR, ppm=C.PPM, rttol=C.PRIOR_RT_TOL), encoding="utf-8")
    print(f"Публикационные рисунки и таблицы → {out}")


if __name__ == "__main__":
    main()

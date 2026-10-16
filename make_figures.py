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


# вариант оформления: "full" (results/publication) | "v4" (publication_v4: без отдельного раздела о ряде,
# «exogenous-signal index», рисунки 3–6 -> 2–5, ряд -> Figure S6, без панели 1d, на ROC только базовая линия «пол»)
VARIANT = {"name": "full"}
TERM = {"index": "Contamination index", "adj": "Contam.-adjusted", "series": "Diethanolamine series (index)",
        "contaminant": "Contaminant", "prior": "Stage-1 markers"}
RENAME = {}


def v4():
    return VARIANT["name"] in ("v4", "v5")


def v5():      # пакет для подачи: v4 + без кривой «index-adjusted» на 4c, без Figure 5/Table S3, Tables S4–S8 -> S3–S7
    return VARIANT["name"] == "v5"


def panel(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom", ha="left")


def save(fig, out, name):
    name = RENAME.get(name, name)
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
    src = pd.read_csv(C.CURATION, dtype=str).fillna("")     # курированное совпадение MS/MS берётся из исходного файла курации
    if "ms2_library" in src:
        key = lambda d_: d_.rt_min.astype(float).round(3).astype(str) + "_" + d_.mz.astype(float).round(5).astype(str)
        cur["ms2_library"] = key(cur).map(dict(zip(key(src), src.ms2_library))).fillna("")
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
    """Уровень MSI каждого признака — единый источник для Table 2, Figure 3/4 и Table S4.
    1 — эталонный стандарт; 2 — курированная формула (≤ PPM) и совпадение MS2 с библиотекой (флаг MS-DIAL) у самого
    признака, без неразрешённых изомеров; 3 — формула без MS2, изомеры неразличимы, или совпадение MS2 без курации;
    4 — всё остальное, несходящаяся масса и неправдоподобные библиотечные совпадения (curation/annotation_exclusions.csv)."""
    t = res_dir / "tables"
    ok = pd.read_csv(t / "candidates_curated.csv", dtype={"isomer_ambiguous": str, "standard_confirmed": str})
    bad = pd.read_csv(t / "candidates_mass_mismatch.csv") if (t / "candidates_mass_mismatch.csv").exists() else pd.DataFrame(columns=["feature"])
    yes = {"1", "yes", "true"}
    ms2f = set(pd.read_csv(C.MS2_IDENTIFIED).feature) if Path(C.MS2_IDENTIFIED).exists() else set()
    has_ms2 = pd.Series(feat.index.isin(ms2f), index=feat.index)
    lvl = pd.Series(4, index=feat.index)
    lvl[has_ms2] = 3
    for _, r in ok.dropna(subset=["feature"]).iterrows():
        if r.feature in lvl.index:
            std = str(r.get("standard_confirmed", "")).strip().lower() in yes
            amb = str(r.get("isomer_ambiguous", "")).strip().lower() in yes
            lvl[r.feature] = 1 if std else (2 if has_ms2[r.feature] and not amb else 3)
    for f in bad.feature.dropna():
        if f in lvl.index:
            lvl[f] = 4
    for f in excluded_annotations(feat):
        lvl[f] = 4
    return lvl, has_ms2


def excluded_annotations(feat):
    """Признаки из curation/annotation_exclusions.csv (сопоставление по m/z ± 0.002 и RT ± 0.02 мин) -> {feature: reason}."""
    if not Path(C.ANNOTATION_EXCLUSIONS).exists():
        return {}
    out = {}
    for _, r in pd.read_csv(C.ANNOTATION_EXCLUSIONS).iterrows():
        m = feat.index[((feat.mz - r.mz).abs() < 0.002) & ((feat.rt_min - r.rt_min).abs() < 0.02)]
        out.update({f: r.reason for f in m})
    return out


# ------------------------------------------------------------------ рисунки
def fig1_quality(out, meta, L, feat, X_all):
    if v4():     # три панели в ряд, без окраски PCA по индексу
        fig, a3 = plt.subplots(1, 3, figsize=(FULL, 6.2 * CM), gridspec_kw={"wspace": 0.38})
        axs = np.array([[a3[0], a3[1]], [a3[2], None]], dtype=object)
    else:
        fig, axs = plt.subplots(2, 2, figsize=(FULL, 12 * CM), gridspec_kw={"hspace": 0.45, "wspace": 0.3})
    ana = meta.index[meta.study | (meta.kind == "QC")]
    grp = meta.loc[ana, "group"]
    tic = np.log10(X_all.loc[ana].sum(axis=1))
    ax = axs[0, 0]
    scatter_groups(ax, meta.loc[ana, "inj_order"].to_numpy(), tic.to_numpy(), grp.to_numpy(), legend=False)
    q = meta.index[meta.kind == "QC"]
    ax.plot(meta.loc[q, "inj_order"], tic[q], c=COL["QC"], lw=0.8, zorder=2)
    ax.set(xlabel="Injection order", ylabel="log$_{10}$ total ion intensity")
    hl = [Line2D([], [], marker=MRK[g], ls="", mfc=COL[g], mec="white", ms=5, label=LAB[g]) for g in ["QC", "Control", "BC", "PC"]]
    if v4():
        fig.legend(handles=hl, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.02))
    else:
        ax.legend(handles=hl, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
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
    # PCA только по соединениям, обнаруженным во всех пробах: импутированные значения (½ минимума) есть почти только
    # в отдельных пробах, а не в пулированном QC, и сдвигают QC к краю облака
    full = [f for f in reps if feat.loc[f, [f"det_{g}" for g in C.STUDY_GROUPS]].min() >= 1.0]
    S, var = multivariate.pca_scores(L.loc[ana, full])
    ax = axs[1, 0]
    scatter_groups(ax, S[0].to_numpy(), S[1].to_numpy(), grp.to_numpy(), legend=False)
    ax.set(xlabel=f"PC1 ({var[0]:.1%})", ylabel=f"PC2 ({var[1]:.1%})")
    panel(ax, "c")
    if v4():
        save(fig, out, "Figure1_data_quality")
        return
    ax = axs[1, 1]
    s = meta.index[meta.study]
    sc = ax.scatter(S.loc[s, 0], S.loc[s, 1], c=meta.loc[s, "contam_index"], cmap="magma_r", s=10, edgecolor="white", lw=0.3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label(TERM["index"], fontsize=6.5)
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
               label=TERM["series"])
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
    ax.set_ylabel(f"{TERM['index']} (mean z-score)")
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


T2_FIRST = set()     # v5: соединения, которые идут первыми в Table 2 / Figure 3 (глюкурониды стероидов)


def table2_order(r, msi):
    """Соединения Table 2 (находки MSI 2–3): общая нумерация для Table 2, Figure 2 и Figure 3.
    По возрастанию q; в v5 сначала глюкурониды стероидов (как в тексте), затем остальные."""
    sel = r[r.finding.astype(bool) & (msi.reindex(r.index) <= 3)].copy()
    sel["grp"] = [0 if f in T2_FIRST else 1 for f in sel.index]
    sel = sel.sort_values(["grp", "q_model"])
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
            ax.scatter(x[m], y[m], s=9, c=col, marker="x", lw=0.6, label=TERM["contaminant"] if c == "contaminant" else c.capitalize(), zorder=3)
    f = r.finding.astype(bool)
    ax.scatter(x[f], y[f], s=14, facecolor=COL["BC"], edgecolor="white", lw=0.4, zorder=4, label="Finding (all criteria met)")
    title = f"{title} ({f.sum()} findings)"
    q = r.p_model[r.q_model < C.FDR]
    if len(q):
        ax.axhline(-np.log10(q.max()), c=INK2, ls="--", lw=0.6, label=f"FDR {C.FDR:.0%}")
    if numbered:
        idx = sorted([i for i in numbered if i in r.index], key=lambda i: numbered[i])   # номера сверху вниз по порядку
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
        ("Index" if v4() else "Contamination"): [to_hex(cmC(ci_n(v))) for v in meta.loc[smp, "contam_index"]],
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
                                       (cmA, "Age (years)", age_n), (cmC, TERM["index"], ci_n)]):
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


def fig5_roc(out, L, meta, feat, mv, pairs, zpairs=None):
    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    cat = feat.loc[reps, "compound_category"].fillna("")
    rho = feat.loc[reps, "rho_contam_index"].abs()
    pool = list(reps[~features.panel_excluded(cat) & ~(rho >= C.PANEL_MAX_ABS_RHO_CONTAM)])   # как в run_all.py
    comps = [("BC_vs_Control", "BC vs Control"), ("BC_vs_PC", "BC vs PC"), ("PC_vs_Control", "PC vs Control")]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 7.2 * CM), gridspec_kw={"wspace": 0.3})
    keys = list(C.COMPARISONS)

    def perm_p(c, sname, key):
        row = mv[(mv.comparison == c) & (mv.set == sname) & mv.model.astype(str).str.contains(key, regex=False)]
        return row.perm_p.iloc[0] if len(row) else np.nan

    for n, (ax, (c, t)) in enumerate(zip(axs, comps)):
        k = keys.index(c)
        case, ref = st.comparison_sets(meta, C.COMPARISONS[c])
        pr = pairs[pairs.comparison == c]
        zp = zpairs[zpairs.comparison == c] if zpairs is not None else pr.iloc[:0]
        runs = [("all", "panel", case, ref, "-", "All")]
        if not (v5() and c == "PC_vs_Control"):     # v5: скорректированная модель для РПЖ vs контроль — только в Table S2
            runs.append(("all", "panel_adj", case, ref, "-.", TERM["adj"]))
        notes = []
        msets = [("matched", pr, "--", "Age-matched")] + \
            ([] if v4() else [("contam_matched", zp, (0, (3, 1, 1, 1, 1, 1)), "Contam.-matched")])
        for sname, pp, ls, lab in msets:
            if len(pp) >= C.MIN_PAIRS_CV:
                runs.append((sname, "panel", list(pp.case), list(pp.ref), ls, f"{lab} ({len(pp)} pairs)"))
            else:
                notes.append(f"{lab}: {len(pp)} pairs (< {C.MIN_PAIRS_CV}, n.e.)")
        if notes:
            ax.text(0.97, 0.03, "\n".join(notes), transform=ax.transAxes, ha="right", va="bottom", fontsize=5.3, color=INK2)
        col = COL[C.COMPARISONS[c]["case"]]
        for sname, method, cs, rf, ls, lab in runs:
            y = np.r_[np.ones(len(cs)), np.zeros(len(rf))].astype(int)
            X = L.loc[cs + rf, pool].to_numpy()
            z = meta.loc[cs + rf, "contam_index"].to_numpy(float)
            aucs, tprs, _, _ = multivariate.repeated_cv(X, y, method, C.SEED + k, z)
            p = perm_p(c, sname, "скорректированы" if method == "panel_adj" else "логистическая")
            if sname == "all" and method == "panel":
                ax.fill_between(multivariate.FPR, np.percentile(tprs, 10, 0), np.percentile(tprs, 90, 0), color=col, alpha=0.15, lw=0)
            ax.plot(multivariate.FPR, tprs.mean(0), c=col if method == "panel" else INK, ls=ls,
                    lw=1.4 if (method == "panel" and sname == "all") else 1.0, label=f"{lab}: {aucs.mean():.2f} ± {aucs.std():.2f}, {pfmt(p)}")
        # базовые модели — та же повторная CV (логистическая регрессия на одной ковариате), см. multivariate.evaluate
        bases = [("sex", INK2, "Sex only")] if v4() else [("contam_index", EXPO_COL["contaminant"], "Contamination index only"), ("sex", INK2, "Sex only")]
        for cov, colr, nmc in bases:
            v = meta.loc[case + ref, cov]
            v = (v.map({"M": 1.0, "F": 0.0}) if cov == "sex" else v.astype(float)).to_numpy(float)
            y = np.r_[np.ones(len(case)), np.zeros(len(ref))].astype(int)
            m = ~np.isnan(v)
            if np.nanstd(v) > 0 and min(y[m].sum(), (1 - y[m]).sum()) >= C.CV_FOLDS:
                aucs, tprs, _, _ = multivariate.repeated_cv(v[m][:, None], y[m], "covariate", C.SEED + k)
                ax.plot(multivariate.FPR, tprs.mean(0), c=colr, ls=":", lw=1,
                        label=f"{nmc}: {aucs.mean():.2f} ± {aucs.std():.2f}, {pfmt(perm_p(c, 'all', f'только {cov}'))}")
        ax.plot([0, 1], [0, 1], c=GRID, lw=0.6)
        ax.set(xlabel="1 − specificity", ylabel="Sensitivity", title=f"{t} (AUC)", xlim=(-0.01, 1.01), ylim=(-0.01, 1.01), aspect="equal")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.1, handlelength=2.2)
        panel(ax, "abc"[n], x=-0.25)
    save(fig, out, "Figure5_ROC_cross_validated")


def fig6_box(out, L, meta, r, nm, msi, X_raw=None):
    """Figure 4: те же соединения, что в Table 2 (находки MSI 2–3), по возрастанию q."""
    sel = r[r.finding.astype(bool) & (msi.reindex(r.index) <= 3)]
    order = table2_order(r, msi)
    sel = sel.loc[sorted(sel.index, key=lambda i: order.get(i, 99))].head(8)
    if not len(sel):
        return
    ncol = 3 if len(sel) in (5, 6) else 4       # 5–6 панелей: сетка 3 + 2 / 3 + 3
    fig, axs = plt.subplots(2, ncol, figsize=(FULL * (0.8 if ncol == 3 else 1), (10 if ncol == 3 else 11) * CM),
                            gridspec_kw={"hspace": 0.85 if ncol == 3 else 1.0, "wspace": 0.5})
    for k, (ax, f) in enumerate(zip(axs.ravel(), sel.index)):
        for i, g in enumerate(C.STUDY_GROUPS):
            v = L.loc[meta.index[meta.study & (meta.group == g)], f]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color=INK, lw=0.9),
                       boxprops=dict(lw=0.6, color=INK2), whiskerprops=dict(lw=0.6, color=INK2), capprops=dict(lw=0.6, color=INK2))
            xj = np.random.RandomState(i).normal(i, 0.08, len(v))
            # импутированные значения (в исходной таблице MS-DIAL сигнала нет) — открытые символы
            raw = X_raw.reindex(v.index)[f] if X_raw is not None and f in X_raw.columns else pd.Series(1.0, index=v.index)
            imp = (raw.isna() | (raw <= 0)).to_numpy()
            ax.scatter(xj[~imp], v[~imp], s=5, c=COL[g], marker=MRK[g], edgecolor="none", alpha=0.85, zorder=3)
            if imp.any():
                ax.scatter(xj[imp], v[imp], s=6, facecolor="white", edgecolor=COL[g], marker=MRK[g], lw=0.5, zorder=3)
        ax.set_xticks(range(3), [LAB[g] for g in C.STUDY_GROUPS])
        name, tag = nm[f].split(", m/z")[0], "m/z" + nm[f].split(", m/z")[1] if ", m/z" in nm[f] else ""
        name = re.sub(r" \((putatively[^)]*|MSI \d|putative)\)", "", name)
        name = "\n".join(textwrap.wrap(name, 24, break_on_hyphens=True)[:2])
        mz, rt = tag.replace("m/z ", "").replace(" min", "").split(" @ ")
        ax.set_title(f"{name}\nMSI {int(msi[f])} · m/z {mz} · {rt} min\nq = {r.loc[f, 'q_model']:.1e}", fontsize=5.6)
        if k % ncol == 0:
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
def table1_note(meta, demo):
    """Сноска к Table 1: недостающие возраст/пол и как они учтены."""
    sdf = meta[meta.study]
    na_age = sdf.age.isna()
    na_sex = ~sdf.sex.isin(["M", "F"])
    def by_group(mask):
        return ", ".join(f"{int(mask[sdf.group == g].sum())} {LAB[g]}" for g in C.STUDY_GROUPS if mask[sdf.group == g].sum())
    cc = int((~na_age & ~na_sex).sum())
    parts = []
    if na_age.sum():
        parts.append(f"age was not available for {int(na_age.sum())} participants ({by_group(na_age)})")
    if na_sex.sum():
        parts.append(f"sex was not available for {int(na_sex.sum())} ({by_group(na_sex)})")
    if (na_age & na_sex).sum():
        parts.append(f"{int((na_age & na_sex).sum())} lacked both")
    if not parts:
        return "Age and sex were available for all participants."
    txt = "; ".join(parts)
    return (txt[0].upper() + txt[1:] + ". These participants were retained in the unadjusted analyses (Mann–Whitney U tests, "
            f"contamination index, cross-validated models on all samples) and excluded from the adjusted linear model "
            f"(complete-case analysis, n = {cc} of {len(sdf)}) and from age-matched pairs; missing values were not imputed.")


def figS5_blanks(out, res, meta):
    """Figure S5: ряд диэтаноламинов в бланках растворителя, QC и пробах с низким/высоким индексом (сырые площади)."""
    d0 = res / "blank_check"
    if not (d0 / "blank_check_areas.csv").exists():
        return
    a = pd.read_csv(d0 / "blank_check_areas.csv", dtype={"sample": str})
    summ = pd.read_csv(d0 / "blank_check_summary.csv")
    a = a[a.kind == "series"]
    a["grp"] = [("Blank" if t.endswith("blank") else "QC" if t == "QC" else meta.group.get(str(sm), "sample"))
                for t, sm in zip(a["type"], a["sample"])]
    order = summ[summ.kind == "series"].drop_duplicates("label")
    order = order[order.in_index == 1].sort_values("solvent_blank_frac_of_QC").label.tolist() + \
        [l for l in summ[(summ.kind == "series") & (summ.in_index == 0)].label.drop_duplicates()]
    frac = summ.drop_duplicates("label").set_index("label").solvent_blank_frac_of_QC
    fig = plt.figure(figsize=(FULL, 9.5 * CM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], hspace=0.55, wspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    sty = {"Blank": dict(c="white", ec=INK, m="o"), "QC": dict(c=COL["QC"], ec="none", m=MRK["QC"]),
           "Control": dict(c=COL["Control"], ec="none", m=MRK["Control"]), "PC": dict(c=COL["PC"], ec="none", m=MRK["PC"])}
    offs = {"Blank": -0.27, "QC": -0.09, "Control": 0.09, "PC": 0.27}
    for k, lab in enumerate(order):
        for g, o in offs.items():
            v = a[(a.label == lab) & (a.grp == g)].area
            if len(v):
                st_ = sty[g]
                ax.scatter(np.log10(v.clip(lower=1e3)), np.full(len(v), k + o), s=9, c=st_["c"], edgecolors=st_["ec"],
                           marker=st_["m"], lw=0.5, zorder=3)
    labs = [f"{l}{'' if pd.isna(frac.get(l)) else ('' if frac[l] < 0.10 else ' †')}" for l in order]
    def sub(l):    # нижние индексы в формулах (C12-DEA — длина цепи, без индекса)
        head, *tail = l.split("+", 1)
        f = lambda x: re.sub(r"([A-Za-z\)])(\d+)", r"\1$_{\2}$", x)
        return (f(head) if not head.startswith("C1") and not head.startswith("C2") or head.startswith("C24") else head) + \
            ("+" + f(tail[0]) if tail else "")
    ax.set_yticks(range(len(order)), [sub(l) for l in labs], fontsize=6)
    ax.invert_yaxis()
    ax.axvline(3, c=GRID, lw=0.6)
    ax.set_xlim(2.8, None)
    ax.text(2.82, len(order) - 0.4, "n.d.", fontsize=5.5, color=INK2, va="bottom")
    ax.set_xlabel("log$_{10}$ peak area (raw, XIC ± 5 ppm)")
    ax.legend(handles=[Line2D([], [], marker=sty[g]["m"], ls="", mfc=sty[g]["c"], mec=INK if g == "Blank" else sty[g]["c"], ms=4,
                              label={"Blank": "Solvent blank (methanol)", "QC": "Pooled QC", "Control": "Control, lowest index",
                                     "PC": "PC, highest index"}[g]) for g in offs],
              loc="upper center", bbox_to_anchor=(0.45, -0.1), ncol=2, fontsize=5.5)
    panel(ax, "a", x=-0.62, y=1.01)
    if (d0 / "blank_check_xic.csv").exists():
        x = pd.read_csv(d0 / "blank_check_xic.csv", dtype={"file": str})
        inj = a.drop_duplicates("file").set_index("file")
        pick = {"Blank": inj.index[inj.grp == "Blank"][0], "QC": inj.index[inj.grp == "QC"][0]}
        for g in ["Control", "PC"]:
            ff = inj.index[inj.grp == g]
            if len(ff):
                pick[g] = ff[0]
        for j, lab in enumerate(["C12-DEA", "C14-DEA+O"]):
            ax2 = fig.add_subplot(gs[j, 1])
            for g, f in pick.items():
                q = x[(x.file == f) & (x.label == lab)]
                nm = {"Blank": "Solvent blank", "QC": "QC"}.get(g, f"{g}, sample {inj.loc[f, 'sample']}")
                ax2.plot(q.rt, q.intensity / 1e6, lw=0.8, c=INK if g == "Blank" else COL.get(g, INK2),
                         ls="--" if g == "Blank" else "-", label=nm)
            fr = frac.get(lab, np.nan)
            ax2.set_title(f"{lab}: solvent blank = {fr:.0%} of QC", fontsize=6.5)
            ax2.set(xlabel="Retention time (min)" if j == 1 else "", ylabel="Intensity (×10$^6$)")
            if j == 0:
                ax2.legend(fontsize=5.2, loc="upper right")
            panel(ax2, "bc"[j], x=-0.3, y=1.06)
    save(fig, out, "FigureS5_blank_check")


def tables(out, demo, demo_t, stats_, nm, cur, mv, feat, msi=None, has_ms2=None, meta=None):
    t1 = demo.rename(columns={"group": "Group"})
    t1["Group"] = t1.Group.map(LAB)
    t1["Sex, M/F (n/a)"] = t1.male.astype(str) + "/" + t1.female.astype(str) + " (" + t1.sex_unknown.astype(str) + ")"
    t1["Age, mean ± SD"] = t1.age_mean.round(1).astype(str) + " ± " + t1.age_sd.round(1).astype(str)
    t1["Age, median (range)"] = t1.age_median.astype(int).astype(str) + " (" + t1.age_min.astype(int).astype(str) + "–" + t1.age_max.astype(int).astype(str) + ")"
    t1["Age available, n"] = t1.age_n
    t1 = t1[["Group", "n", "Sex, M/F (n/a)", "Age, mean ± SD", "Age, median (range)", "Age available, n"]]
    t1_note = table1_note(meta, demo) if meta is not None else ""
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
    m = mv[mv.model.astype(str).str.startswith(("панель", "PLS-DA", "только"))].copy()
    m["model"] = m.model.str.replace(r"панель top-(\d+) \+ логистическая регрессия \(отбор внутри CV\)", r"Top-\1 panel + logistic regression (in-fold selection)", regex=True) \
        .str.replace(r"панель top-(\d+), признаки скорректированы на индекс загрязнения внутри фолда", r"Top-\1 panel, features adjusted for contamination index in-fold", regex=True) \
        .str.replace(r"PLS-DA, 2 компоненты \(все признаки без экспозиций\)", "PLS-DA, 2 components (all non-exposure compounds)", regex=True) \
        .str.replace(r"только sex \(базовая модель, та же CV\)", "Sex only (baseline, same CV)", regex=True) \
        .str.replace(r"только contam_index \(базовая модель, та же CV\)", "Contamination index only (baseline, same CV)", regex=True)
    m["set"] = m.set.map({"all": "all samples", "matched": "age-matched pairs", "contam_matched": "contamination-matched pairs"})
    m["comparison"] = m.comparison.str.replace("_vs_", " vs ").str.replace("maleControl", "male Control")
    ts1 = m[["comparison", "set", "model", "n_case", "n_ref", "AUC_mean", "AUC_sd", "AUC_p5", "AUC_p95", "Q2", "perm_p"]].round(3)
    ts1.columns = ["Comparison", "Sample set", "Model", "n case", "n reference", "AUC, mean", "AUC, SD", "AUC, 5th pct",
                   "AUC, 95th pct", "Q²", "Permutation p"]
    # Table 2 (основной текст): только аннотированные находки (MSI ≤ 3); полный список — Table S1
    t2_main = t2[t2["MSI level"] <= 3].copy()
    t2_main["Compound (annotation)"] = t2_main["Compound (annotation)"].str.replace(r" \((MSI \d|putative)\)$", "", regex=True)
    lib = cn["ms2_library"] if "ms2_library" in cn else pd.Series(dtype=str)   # курированное совпадение (напр. по агликону)
    t2_main.insert(1, "Library match (MS-DIAL)", [lib.get(i) if isinstance(lib.get(i), str) and lib.get(i) else feat.loc[i, "annotation"]
                                                   for i in t2_main.index])
    t2_main = t2_main.drop(columns=["Putative library match (low score)"])
    iso = cn["isomers_note"].fillna("") if "isomers_note" in cn else pd.Series(dtype=str)
    t2_main.insert(t2_main.columns.get_loc("MSI level") + 1, "Isomers not excluded",
                   [iso.get(i, "") if isinstance(iso.get(i, ""), str) else "" for i in t2_main.index])
    order = table2_order(r, msi)
    t2_main.insert(0, "No.", [order[i] for i in t2_main.index])
    t2_main = t2_main.sort_values("No.")
    exc = excluded_annotations(feat)
    t2.insert(4, "Note", [f"library match judged implausible: {exc[i]}" if i in exc else "" for i in t2.index])
    t2_note = ("Isomers: " + "; ".join(f"{a}: {b}" for a, b in zip(t2_main["Compound (annotation)"], t2_main["Isomers not excluded"]) if b)) \
        if t2_main["Isomers not excluded"].astype(bool).any() else ""
    with pd.ExcelWriter(out / "Tables_main.xlsx") as w:
        t1.to_excel(w, sheet_name="Table1_demographics", index=False)
        pd.DataFrame({"Note": [t1_note]}).to_excel(w, sheet_name="Table1_demographics", index=False, startrow=len(t1) + 2)
        tt.to_excel(w, sheet_name="Table1_tests", index=False)
        t2_main.to_excel(w, sheet_name="Table2_annotated_findings_BC", index=False)
        pd.DataFrame({"Note": [t2_note]}).to_excel(w, sheet_name="Table2_annotated_findings_BC", index=False, startrow=len(t2_main) + 2)
    (out / "table_notes.md").write_text(f"**Table 1 note.** {t1_note}\n\n**Table 2 note.** {t2_note}\n", encoding="utf-8")
    t1.to_csv(out / "Table1_demographics.csv", index=False)
    tt.to_csv(out / "Table1_tests.csv", index=False)
    t2_main.to_csv(out / "Table2_annotated_findings_BC.csv", index=False)
    t2.to_csv(out / "TableS1_all_findings_BC.csv", index=False)
    ts1.to_csv(out / "TableS2_classification.csv", index=False)
    return t2, ts1


CAPTIONS = """# Figure and table captions (draft)

## Main text

**Table 1.** Clinical and demographic characteristics of the independent cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). The footnote on missing age and sex is given in table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the {rsd:.0%} acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound; the {nfull} compounds detected in all study samples, i.e., without imputed values). (d) The same scores coloured by the contamination index.

**Figure 2.** Plastic-derived surfactant contamination. (a) Homologous series of alkyl diethanolamines and their oxygenated/ethoxylated analogues ([M+H]+) used to build the contamination index; open diamonds, co-varying ions not included in the index. (b) Contamination index by group (Mann–Whitney U test).

**Figure 3.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + contamination index). Filled circles, findings that passed all robustness criteria; crosses, exposure-related compounds (contaminants, food, drugs); dashed line, FDR {fdr:.0%}. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3). log2FC and q, adjusted model (sex, age, contamination index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-contamination controls), BC vs controls with a contamination index below the median. MSI 2, precursor within {ppm:.0f} ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the feature itself; MSI 3, formula-consistent annotation without MS/MS confirmation, MS/MS match to several isomers, or MS/MS match without a curated formula. "Isomers not excluded" lists isomers that the MS/MS spectrum does not distinguish (e.g., Leu-Val, Val-Ile and Ile-Val for Val-Leu; sn-1/sn-2 regioisomers for LysoPC). For compounds without a curated annotation, the MS-DIAL library match is given for transparency. Library matches judged implausible in human urine (ethyl 2-amino-5-methyl-3-thiophenecarboxylate, a synthetic building block matched to an in-source fragment) were assigned MSI level 4. All findings, including MSI level 4, are listed in Table S1.

**Figure 4.** Annotated bladder-cancer findings (MSI levels 2–3; the compounds of Table 2) across groups (log2 PQN-normalized intensities). Open symbols, values below the detection limit, imputed as half of the feature minimum. q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 5.** Cross-validated ROC curves ({rep}×{fold}-fold CV; top-{k} features selected by Mann–Whitney U inside each training fold, logistic regression). Candidate pool defined before CV: exposure-related compounds (contaminants, food, drugs, tobacco) and compounds correlated with the contamination index (|ρ| ≥ {rho}) excluded. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot, features residualized on the contamination index inside each training fold; dashed, 1:1 age-matched pairs; dash-dot-dot, 1:1 contamination-matched pairs (within sex, caliper {zcal} SD of the index; sensitivity analysis); matched sets shown only with ≥ {minp} pairs. p, permutation test ({perm} permutations). Dotted lines, single-covariate baselines (sex only; contamination index only) evaluated with the same repeated CV (logistic regression on one covariate) and permutation test. AUC, mean ± SD over repeats. In (c) the contamination index alone nearly separates PC from controls; the contamination-adjusted curve is then an extrapolation beyond the overlap of the two groups and is not interpretable as a contamination-independent signal (see Section 3.6).

**Figure 6.** Stage-1 bladder-cancer markers (Frolova et al., Anal Bioanal Chem 2026, Table 5) in the independent cohort. Features were matched by m/z (±{ppm:.0f} ppm) and retention time (±{rttol} min); all isomeric candidates are shown. Points, unadjusted log2 fold change (BC vs Control) of quantified features; coloured, same direction and Mann–Whitney p < 0.05; text gives AUC, unadjusted p and FDR-adjusted q of the adjusted model. Markers not detected or removed by the detection/QC filters are labelled. Exact values are given in Table S3.

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and contamination index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±{cal} years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < {low}. (b) Number of compounds with q < {fdr} under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, {ppm:.0f} ppm; exposure-related features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Figure S5.** Alkyl diethanolamine series in solvent blanks, pooled QC and urine samples (batch 1; raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas of each series member in methanol solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest contamination index; † members with solvent-blank area ≥ 10% of the QC median (system/solvent background); n.d., not detected. (b, c) Extracted-ion chromatograms of C12-DEA (present in blanks) and its oxygenated analogue C14-DEA+O (absent from blanks).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, Q², permutation p).

**Table S3.** Stage-1 markers: all candidate features within ±{ppm:.0f} ppm of the reported or the theoretical m/z (formula and ion) and ±{rttol} min of the reported retention time, detection/filter status and statistics in the independent cohort. For Leu-Leu-Tyr and Leu-Thr-Ser the formula of the named peptide was used to compute the theoretical m/z.

**Table S4.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (assigned by the same rule as in Table 2), supporting ions, detection by group, and correlation with the contamination index and the dilution factor; annotations with mass error > {ppm:.0f} ppm are listed separately (S4b). No reference standards were analysed (no MSI level 1). MS-DIAL match scores were not part of the exported table; the MS/MS match is the MS-DIAL identification flag.

**Table S5.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S6.** Compound-level statistics for all comparisons (all models and robustness checks).

**Table S7.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < 0.2, with a contamination index built only from the series members that were absent from solvent blanks, restricted to compounds without imputed values, and with imputed values treated as missing; overlap with the main analysis and cross-validated panel AUC.

**Table S8.** Solvent-blank check of the alkyl diethanolamine series (batch 1): median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.
"""


CAPTIONS_V4 = """# Figure and table captions — v4 ("Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study")

## Main text

**Table 1.** Clinical and demographic characteristics of the study cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). Footnote: see table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the {rsd:.0%} acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound; the {nfull} compounds detected in all study samples, i.e., without imputed values).

**Figure 2.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + exogenous-signal index). Filled circles, findings that passed all robustness criteria; crosses, exogenous compounds excluded from marker selection (drugs, food- and tobacco-related compounds, background ions); dashed line, FDR {fdr:.0%}. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3). log2FC and q, adjusted model (sex, age, exogenous-signal index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-index controls), BC vs controls with an exogenous-signal index below the median. MSI 2, precursor within {ppm:.0f} ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the feature itself; MSI 3, formula-consistent annotation without MS/MS confirmation, MS/MS match to several isomers, or MS/MS match without a curated formula. "Isomers not excluded" lists isomers that the MS/MS spectrum does not distinguish. All findings, including MSI level 4, are listed in Table S1.

**Figure 3.** Annotated bladder-cancer findings (the compounds of Table 2) across groups (log2 PQN-normalized intensities). Open symbols, values below the detection limit, imputed as half of the feature minimum. q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 4.** Cross-validated ROC curves ({rep}×{fold}-fold CV; top-{k} compounds selected by Mann–Whitney U inside each training fold, logistic regression). Exogenous compounds and compounds correlated with the exogenous-signal index (|ρ| ≥ {rho}) were excluded from the candidate pool before CV. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot, compounds adjusted for the exogenous-signal index inside each training fold; dashed, 1:1 age- and sex-matched pairs (shown only with ≥ {minp} pairs); dotted, sex-only baseline evaluated with the same CV. AUC, mean ± SD over repeats; p, permutation test ({perm} permutations). In (c) the exogenous-signal index nearly separates PC from controls (Table S2), so the index-adjusted estimate extrapolates beyond the overlap of the two groups and is not interpreted as a PC-specific signal.

**Figure 5** (optional; omit together with Section 3.6 if the comparison with the previous study is not included). Previously reported bladder-cancer markers in the present cohort. Features were matched by m/z (± {ppm:.0f} ppm of the reported or theoretical value) and retention time (± {rttol} min); all isomeric candidates are shown. Points, unadjusted log2 fold change (BC vs Control); coloured, same direction and Mann–Whitney p < 0.05; text gives AUC, unadjusted p and FDR-adjusted q of the adjusted model. Exact values in Table S3.

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and exogenous-signal index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±{cal} years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < {low}. (b) Number of compounds with q < {fdr} under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, {ppm:.0f} ppm; exogenous features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Figure S5.** Members of the homologous exogenous series in methanol solvent blanks, pooled QC and urine samples (raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas in solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest index; † members with solvent-blank area ≥ 10% of the QC median; n.d., not detected. (b, c) Extracted-ion chromatograms of a member present in the blanks (C12) and of one absent from the blanks (C14+O).

**Figure S6.** Homologous exogenous series used for the exogenous-signal index. (a) Members ([M+H]+; mass increments C2H4, O and C2H4O) in the m/z–retention-time plane; open diamonds, co-varying ions not included in the index. (b) Exogenous-signal index (mean z-score of the ten members) by group (Mann–Whitney U test).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, 5th–95th percentiles, Q², permutation p), including the single-covariate baselines (sex only; exogenous-signal index only).

**Table S3.** Previously reported bladder-cancer markers: all candidate features within ± {ppm:.0f} ppm of the reported or theoretical m/z and ± {rttol} min of the reported retention time, detection/filter status and statistics in the present cohort.

**Table S4.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (same rule as in Table 2), supporting ions, detection by group, and correlation with the exogenous-signal index and the dilution factor; annotations with mass error > {ppm:.0f} ppm are listed separately (S4b). No reference standards were analysed (no MSI level 1).

**Table S5.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S6.** Compound-level statistics for all comparisons (all models and robustness checks).

**Table S7.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < {low}, with the exogenous-signal index built only from series members absent from solvent blanks, restricted to compounds without imputed values, and with imputed values treated as missing; overlap with the main analysis and cross-validated panel AUC.

**Table S8.** Solvent-blank check of the homologous exogenous series: median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.
"""

CAPTIONS_V5 = """# Figure and table captions — submission version (v5) ("Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study")

## Main text

**Table 1.** Clinical and demographic characteristics of the study cohort. Age is given as mean ± SD and median (range); p-values: Kruskal–Wallis test (three groups), Mann–Whitney U test (pairwise age) and Fisher's exact test (sex). Footnote: see table_notes.md (Table 1 note).

**Figure 1.** Data quality. (a) Total ion intensity of pooled QC samples and study samples versus injection order. (b) Distribution of feature-wise relative standard deviation (RSD) in pooled QC samples; dashed line, the {rsd:.0%} acceptance threshold. (c) PCA score plot of QC and study samples (log2, Pareto scaling; one representative feature per compound; the {nfull} compounds detected in all study samples, i.e., without imputed values).

**Figure 2.** Volcano plots from the adjusted linear model (log2 intensity ~ group + sex + age + exogenous-signal index). Filled circles, findings that passed all robustness criteria; crosses, exogenous compounds excluded from marker selection (drugs, food- and tobacco-related compounds, the homologous exogenous series); dashed line, FDR {fdr:.0%}. In (a), circled numbers mark the annotated findings listed in Table 2 (No.).

**Table 2.** Annotated bladder-cancer findings (MSI levels 2–3), steroid glucuronides first, then other compounds, each ordered by q. log2FC and q, adjusted model (sex, age, exogenous-signal index; Benjamini–Hochberg); AUC, all samples; p (age-matched), Mann–Whitney U test in 1:1 age- and sex-matched pairs; p (low-index controls), BC vs controls with an exogenous-signal index below the median. Findings met all of the following: q < 0.05 in the adjusted model; p < 0.05 in age-matched pairs and against low-index controls, with the same direction of change; and no difference between controls with high and low index (p ≥ 0.05), which excludes compounds whose level itself depends on the index. MSI 2, precursor within {ppm:.0f} ppm of the curated formula and an MS/MS library match flagged by MS-DIAL (NIST 23 HRMS / MS-DIAL MSP) for the feature itself; MSI 3, formula-consistent annotation without MS/MS confirmation, MS/MS match to several isomers, or MS/MS match without a curated formula. "Isomers not excluded" lists isomers that the MS/MS spectrum does not distinguish. All findings, including MSI level 4, are listed in Table S1.

**Figure 3.** Annotated bladder-cancer findings (the compounds of Table 2) across groups (log2 PQN-normalized intensities). Open symbols, values below the detection limit, imputed as half of the feature minimum. q, FDR-adjusted p of the adjusted model (BC vs Control).

**Figure 4.** Cross-validated ROC curves ({rep}×{fold}-fold CV; top-{k} compounds selected by Mann–Whitney U inside each training fold, logistic regression). Exogenous compounds and compounds correlated with the exogenous-signal index (|ρ| ≥ {rho}) were excluded from the candidate pool before CV. Solid, all samples (shaded band, 10th–90th percentile across repeats); dash-dot (a, b), compounds adjusted for the exogenous-signal index inside each training fold (for PC vs Control, see Table S2 and its note); dashed, 1:1 age- and sex-matched pairs (shown only with ≥ {minp} pairs); dotted, sex-only baseline evaluated with the same CV. AUC, mean ± SD over repeats; p, permutation test ({perm} permutations).

## Supplementary Materials

**Figure S1.** Heatmap of the 50 bladder-cancer findings with the smallest q (one row per compound; z-scored log2 intensities). Column bars: group, sex, age and exogenous-signal index.

**Figure S2.** (a) Age by group and sex (filled, male; open, female). (b) 1:1 nearest-age matching within sex (caliper ±{cal} years).

**Figure S3.** (a) PQN versus MSTUS dilution factors; labelled samples have PQN factor < {low}. (b) Number of compounds with q < {fdr} under PQN and MSTUS normalization and their overlap.

**Figure S4.** Pathway enrichment (mummichog, positive mode, {ppm:.0f} ppm; exogenous features excluded). Coloured bars, Benjamini–Hochberg q < 0.05.

**Figure S5.** Members of the homologous exogenous series in methanol solvent blanks, pooled QC and urine samples (raw peak areas from extracted-ion chromatograms, ± 5 ppm). (a) Peak areas in solvent blanks (n = 7, including two injected before the first sample), pooled QC (n = 6), the two controls with the lowest and the three PC samples with the highest index; † members with solvent-blank area ≥ 10% of the QC median; n.d., not detected. (b, c) Extracted-ion chromatograms of a member present in the blanks (C12) and of one absent from the blanks (C14+O).

**Figure S6.** Homologous exogenous series used for the exogenous-signal index. (a) Members ([M+H]+; mass increments C2H4, O and C2H4O) in the m/z–retention-time plane; open diamonds, co-varying ions not included in the index. (b) Exogenous-signal index (mean z-score of the ten members) by group (Mann–Whitney U test).

**Table S1.** All bladder-cancer findings (BC vs Control), including MSI level 4 features; low-score library matches are listed for transparency only.

**Table S2.** Cross-validated classification for all comparisons and models (AUC mean ± SD over repeats, 5th–95th percentiles, Q², permutation p), including the single-covariate baselines (sex only; exogenous-signal index only).

**Table S3.** Curated candidate annotations: formula, ion, mass error, HMDB/KEGG identifiers, MS-DIAL MS/MS library entry, isomers not excluded, MSI level (same rule as in Table 2), supporting ions, detection by group, and correlation with the exogenous-signal index and the dilution factor; annotations with mass error > {ppm:.0f} ppm are listed separately (S3b). No reference standards were analysed (no MSI level 1).

**Table S4.** Presence/absence features (detection < 50% in one group and ≥ 80% in another) with Fisher's exact tests.

**Table S5.** Compound-level statistics for all comparisons (all models and robustness checks).

**Table S6.** Sensitivity analyses: number of compounds with q < 0.05 and of findings under MSTUS instead of PQN normalization, after exclusion of samples with a PQN dilution factor < {low}, with the exogenous-signal index built only from series members absent from solvent blanks, restricted to compounds without imputed values, and with imputed values treated as missing; overlap with the main analysis and cross-validated panel AUC.

**Table S7.** Solvent-blank check of the homologous exogenous series: median raw peak areas in methanol solvent blanks (n = 7), pooled QC (n = 6) and five study samples, and the blank-to-QC ratio.

**Table S8.** Keyword list used to flag drug-, food- and tobacco-related compounds: regular expressions (case-insensitive) matched against MS-DIAL library names (NIST 23 HRMS, MS-DIAL MSP), compiled before the statistical analysis, and the number of compounds flagged by each expression (high-scoring library matches; compounds with low-scoring drug- or tobacco-related matches were excluded only from the candidate pool of the classification models).
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
        "m/z reported (stage 1)": tab.mz, "m/z theoretical": tab.theo_mz,
        "Δ reported vs theoretical": [f"{p_:+.1f} ppm" if abs(p_) < 1000 else f"{r_ - t_:+.4f} Da"
                                      for p_, r_, t_ in zip(tab.reported_ppm_vs_theo, tab.mz, tab.theo_mz)], "RT, min (stage 1)": tab.rt_min,
        "Candidate": [f"{int(c)}/{int(n)}" if n else "" for c, n in zip(tab.candidate, tab.n_candidates)],
        "m/z (this study)": tab.get("mz_obs").round(4), "Δ vs theoretical, ppm": tab.get("ppm_obs_vs_theo"),
        "RT, min (this study)": tab.get("rt_obs").round(2),
        "ΔRT, min": tab.get("delta_rt"), "Status (this study)": tab.status.map(st_lab),
        "log2FC BC/Control": tab.get("log2FC").round(2), "AUC": tab.get("AUC").round(2),
        "p (Mann–Whitney)": tab.get("p_MW").map(lambda x: f"{x:.2g}" if pd.notna(x) else ""),
        "q (adjusted model)": tab.get("q_model").map(lambda x: f"{x:.2g}" if pd.notna(x) else ""),
        "Same direction": tab.get("same_direction").map({True: "yes", False: "no"}),
    })
    summ_en = summ.rename(columns={"source": "Source", "comparison": "Comparison", "markers": "Markers",
                                "not_detected": "Not detected", "below_detection_filter": "Below detection filter",
                                "failed_rsd": "Failed QC RSD filter", "quantified": "Quantified",
                                "ambiguous_isomers": "Quantified with isomeric candidates",
                                "same_direction_all_candidates": "Same direction (all candidates)",
                                "same_direction_any_candidate": "Same direction (any candidate)",
                                "nominal_p05_any_candidate": "p < 0.05, same direction (any candidate)",
                                "replicated_q05_any_candidate": "q < 0.05, same direction (any candidate)"})
    t3.to_csv(out / "TableS3_replication_stage1_markers.csv", index=False)
    summ_en.to_csv(out / "TableS3_replication_summary.csv", index=False)
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
    ax.set_xlabel("log$_2$ fold change BC/Control in this cohort (unadjusted)")
    ax.spines["left"].set_visible(False)
    q = summ.iloc[0]
    ax.set_title(f"{TERM['prior']} (n = {int(q.markers)}): quantified {int(q.quantified)} ({int(q.ambiguous_isomers)} with isomers),\n"
                 f"p < 0.05 in the same direction {int(q.nominal_p05_any_candidate)}, q < {C.FDR} (adjusted) {int(q.replicated_q05_any_candidate)}",
                 fontsize=6.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=COL["BC"], mec=INK2, label="same direction, p < 0.05"),
                       Line2D([], [], marker="o", ls="", mfc="#bdbcb8", mec=INK2, label="not significant")],
              loc="lower right", fontsize=6)
    save(fig, out, "Figure6_replication_stage1_markers")
    return t3, summ_en


def s4_table(df, feat, msi, has_ms2):
    """Table S4/S4b: курированные кандидаты — английские столбцы; MSI и MS2 из того же источника, что Table 2."""
    if not len(df):
        return df
    f = df.feature
    def g(col):
        return df[col] if col in df else pd.Series("", index=df.index)
    ms2 = [bool(has_ms2.get(x, False)) if isinstance(x, str) else False for x in f]
    out = pd.DataFrame({
        "Compound": [a if isinstance(a, str) and a else b for a, b in zip(g("name_en"), g("name"))],
        "Feature RT (min)": df.rt_min.round(3), "Feature m/z": df.mz.round(5), "Formula": df.formula, "Ion": df.ion,
        "Theoretical m/z": df.theo_mz, "Mass error (ppm)": df.ppm, "HMDB": g("hmdb"), "KEGG": g("kegg"),
        "Category": g("category"), "Comment": g("comment"),
        "MS2 library match (MS-DIAL)": ["yes" if m else "no" for m in ms2],
        "Library entry (MS-DIAL)": [feat.loc[x, "annotation"] if m and x in feat.index else "" for x, m in zip(f, ms2)],
        "Isomers not excluded": g("isomers_note").fillna(""),
        "MSI level": [int(msi[x]) if isinstance(x, str) and x in msi.index else 4 for x in f],
        "Supporting ions": g("confirming_ions"), "HMDB alternatives within 3 ppm": g("hmdb_alternatives_3ppm"),
    })
    for grp in C.STUDY_GROUPS:
        if f"det_{grp}" in df:
            out[f"Detected, % {LAB[grp]}"] = (100 * df[f"det_{grp}"].astype(float)).round(0)
    extra = {"kept_after_filters": "Passed detection and RSD filters", "compound": "Compound group (ions)",
             "is_representative": "Representative ion",
             "rho_contam": "ρ with contamination index", "p_contam": "p (ρ contamination)",
             "rho_contam_controls": "ρ with contamination index (controls)", "p_contam_controls": "p (ρ contamination, controls)",
             "rho_dilution": "ρ with PQN dilution factor", "p_dilution": "p (ρ dilution)",
             "rho_dilution_controls": "ρ with PQN dilution factor (controls)", "p_dilution_controls": "p (ρ dilution, controls)"}
    for k, v in extra.items():
        if k in df:
            out[v] = df[k].round(3) if df[k].dtype.kind == "f" else df[k]
    lab = {"beta_model": "log2FC (adjusted)", "q_model": "q (adjusted)", "p_matched": "p (age-matched)",
           "p_lowindex": "p (low-contamination controls)", "p_ctrl_hi_vs_lo": "p (controls, high vs low index)", "finding": "finding"}
    for c in df.columns:
        if ":" in c:
            comp, stat = c.split(":", 1)
            out[f"{comp.replace('_vs_', ' vs ').replace('maleControl', 'male Control')}: {lab.get(stat, stat)}"] = df[c]
    return out


def s5_table(d):
    d = d.rename(columns={d.columns[0]: "Feature"})
    cols = {"Feature": "Feature", "rt_min": "RT (min)", "mz": "m/z", "annotation": "MS-DIAL annotation",
            "annotation_conf": "Annotation score (MS-DIAL)", "compound_category": "Category",
            "det_QC": "Detected, fraction QC", "det_Control": "Detected, fraction Control", "det_BC": "Detected, fraction BC",
            "det_PC": "Detected, fraction PC", "RSD_QC_raw": "RSD in QC"}
    cols.update({c: c.replace("p_fisher_", "p (Fisher), ").replace("_vs_", " vs ").replace("maleControl", "male Control")
                 for c in d.columns if c.startswith("p_fisher_")})
    return d[[c for c in cols if c in d]].rename(columns=cols)


def s6_table(d):
    cols = {"feature": "Representative feature", "name": "Label", "n_ions": "Ions in group", "rt_min": "RT (min)", "mz": "m/z",
            "annotation": "MS-DIAL annotation", "annotation_conf": "Annotation score (MS-DIAL)", "compound_category": "Category",
            "log2FC": "log2FC (unadjusted)", "AUC": "AUC", "p_MW": "p (Mann–Whitney)", "q_MW": "q (Mann–Whitney)",
            "beta_model": "log2FC (adjusted)", "p_model": "p (adjusted model)", "q_model": "q (adjusted model)",
            "log2FC_matched": "log2FC (age-matched)", "AUC_matched": "AUC (age-matched)", "p_matched": "p (age-matched)",
            "q_matched": "q (age-matched)", "log2FC_lowindex": "log2FC (low-contamination controls)",
            "p_lowindex": "p (low-contamination controls)", "log2FC_ctrl_hi_vs_lo": "log2FC (controls, high vs low index)",
            "p_ctrl_hi_vs_lo": "p (controls, high vs low index)", "finding": "Finding",
            "det_Control": "Detected, fraction Control", "det_BC": "Detected, fraction BC", "det_PC": "Detected, fraction PC"}
    return d[[c for c in cols if c in d]].rename(columns=cols)


def supplementary_tables(out, res_dir, t2_full, ts_class, rep, feat, msi, has_ms2):
    """Supplementary_Tables.xlsx: S1 все находки РМП, S2 классификация, S3 сверка с первым этапом,
    S4 курированные кандидаты (+ несходящаяся масса), S5 признаки «есть/нет», S6 статистика по соединениям."""
    t = res_dir / "tables"
    s4 = s4_table(pd.read_csv(t / "candidates_curated.csv"), feat, msi, has_ms2)
    s4b = s4_table(pd.read_csv(t / "candidates_mass_mismatch.csv"), feat, msi, has_ms2) \
        if (t / "candidates_mass_mismatch.csv").exists() else None
    with pd.ExcelWriter(out / "Supplementary_Tables.xlsx") as w:
        t2_full.to_excel(w, sheet_name="S1_all_findings_BC", index=False)
        ts_class.to_excel(w, sheet_name="S2_classification", index=False)
        if rep is not None:
            rep[0].to_excel(w, sheet_name="S3_stage1_markers", index=False)
            rep[1].to_excel(w, sheet_name="S3_summary", index=False)
        s4.to_excel(w, sheet_name="S4_curated_candidates", index=False)
        if s4b is not None and len(s4b):
            s4b.to_excel(w, sheet_name="S4b_mass_mismatch", index=False)
        if (t / "presence_absence.csv").exists():
            s5_table(pd.read_csv(t / "presence_absence.csv")).to_excel(w, sheet_name="S5_presence_absence", index=False)
        for n in C.COMPARISONS:
            if (t / f"stats_{n}.csv").exists():
                s6_table(pd.read_csv(t / f"stats_{n}.csv")).to_excel(w, sheet_name=f"S6_{n}"[:31], index=False)
    s4.to_csv(out / "TableS4_curated_candidates.csv", index=False)
    s7, s8 = sensitivity_table(t), blank_table(res_dir)
    if s7 is not None or s8 is not None:
        with pd.ExcelWriter(out / "Supplementary_Tables.xlsx", mode="a", engine="openpyxl") as w:
            if s7 is not None:
                s7.to_excel(w, sheet_name="S7_sensitivity", index=False)
            if s8 is not None:
                s8.to_excel(w, sheet_name="S8_blank_check", index=False)
        for x, n in [(s7, "TableS7_sensitivity.csv"), (s8, "TableS8_blank_check.csv")]:
            if x is not None:
                x.to_csv(out / n, index=False)


def sensitivity_table(t):
    """Table S7: число значимых соединений и находок при альтернативной обработке (относительно основного анализа)."""
    rows = []
    if (t / "normalization_sensitivity.csv").exists():
        for _, r in pd.read_csv(t / "normalization_sensitivity.csv").iterrows():
            c = r.comparison.replace("_vs_", " vs ").replace("maleControl", "male Control")
            rows.append({"Comparison": c, "Analysis": "Main (PQN, all samples, 10-member index)", "q < 0.05": r["q<FDR_PQN"],
                         "Findings": r.findings_PQN, "Shared with main (q < 0.05)": r["q<FDR_PQN"], "Shared findings": r.findings_PQN,
                         "Panel AUC (CV)": r.get("panel_AUC_PQN")})
            rows.append({"Comparison": c, "Analysis": "MSTUS normalization", "q < 0.05": r["q<FDR_MSTUS"], "Findings": r.findings_MSTUS,
                         "Shared with main (q < 0.05)": r.overlap, "Shared findings": r.findings_overlap,
                         "Panel AUC (CV)": r.get("panel_AUC_MSTUS")})
            rows.append({"Comparison": c, "Analysis": f"Without samples with PQN factor < {C.LOW_DILUTION_FACTOR} (n = {r.n_without_lowdilution})",
                         "q < 0.05": r.q_without_lowdilution, "Findings": r.findings_without_lowdilution,
                         "Shared with main (q < 0.05)": r.overlap_with_main, "Shared findings": r.findings_overlap_main,
                         "Panel AUC (CV)": r.get("panel_AUC_without_lowdilution")})
    if (t / "sensitivity_contam_index.csv").exists():
        for _, r in pd.read_csv(t / "sensitivity_contam_index.csv").iterrows():
            c = r.comparison.replace("_vs_", " vs ").replace("maleControl", "male Control")
            rows.append({"Comparison": c, "Analysis": f"Contamination index from {int(r.members_clean_index)} members absent from solvent blanks",
                         "q < 0.05": r.q05_clean, "Findings": r.findings_clean, "Shared with main (q < 0.05)": r.q05_overlap,
                         "Shared findings": r.findings_overlap, "Panel AUC (CV)": r.get("AUC_panel_adj"),
                         "Note": "panel AUC: features adjusted for this index in-fold" if pd.notna(r.get("AUC_panel_adj")) else ""})
    if (t / "sensitivity_imputation.csv").exists():
        for _, r in pd.read_csv(t / "sensitivity_imputation.csv").iterrows():
            c = r.comparison.replace("_vs_", " vs ").replace("maleControl", "male Control")
            rows.append({"Comparison": c, "Analysis": f"Only the {int(r.compounds_without_imputation):,} compounds detected in all samples (no imputed values; FDR within this set)",
                         "q < 0.05": r.A_q05, "Findings": r.A_findings, "Shared findings": r.A_shared,
                         "Note": f"{int(r.main_findings_without_imputation)} of the {int(r.main_findings)} main findings belong to this set"})
            rows.append({"Comparison": c, "Analysis": "Imputed values treated as missing (each compound analysed on its detected values)",
                         "q < 0.05": r.B_q05, "Findings": r.B_findings, "Shared findings": r.B_shared})
    if not rows:
        return None
    d = pd.DataFrame(rows)
    d["Panel AUC (CV)"] = d["Panel AUC (CV)"].astype(float).round(2)
    return d.sort_values("Comparison", kind="stable")


def blank_table(res_dir):
    """Table S8: проверка бланков — площади членов ряда в бланках растворителя, QC и пробах."""
    p = res_dir / "blank_check" / "blank_check_summary.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    d = d[d.kind == "series"]
    cols = {"label": "Series member", "mz": "m/z [M+H]+ (theoretical)", "rt_expected": "RT (min)",
            "in_index": "In contamination index", "solvent_blank": "Solvent blank, median area", "QC": "QC, median area",
            "sample": "Samples (2 lowest-, 3 highest-index), median area", "solvent_blank_frac_of_QC": "Solvent blank / QC"}
    d = d[[c for c in cols if c in d]].rename(columns=cols)
    d["In contamination index"] = d["In contamination index"].map({1: "yes", 0: "no"})
    d["Solvent blank / QC"] = d["Solvent blank / QC"].round(3)
    d["Background in solvent blank (≥ 10% of QC)"] = np.where(d["Solvent blank / QC"] >= 0.10, "yes", "no")
    for c in d.columns:
        if "median area" in c:
            d[c] = d[c].round(-2)
    return d


V4_TABLE_TERMS = [   # замены в таблицах и сносках варианта v4 (порядок важен: длинные фразы раньше коротких)
    ("features adjusted for contamination index in-fold", "features adjusted for exogenous-signal index in-fold"),
    ("Contamination index only", "Exogenous-signal index only"),
    ("contamination-matched pairs", "index-matched pairs"),
    ("low-contamination controls", "low-index controls"),
    ("Contamination index from", "Exogenous-signal index from"),
    ("In contamination index", "In exogenous-signal index"),
    ("ρ with contamination index", "ρ with exogenous-signal index"),
    ("p (ρ contamination", "p (ρ index"),
    ("10-member index", "10-member exogenous-signal index"),
    ("contamination index", "exogenous-signal index"),
    ("Contamination index", "Exogenous-signal index"),
    ("age- and contamination-matched pairs", "age-matched pairs"),
    ("Marker (stage 1)", "Marker (previous study)"),
    ("(stage 1)", "(previous study)"),
]
V4_CELL_VALUES = {"contaminant": "background ion"}


V5_SHEETS = {"S4_curated_candidates": "S3_curated_candidates", "S4b_mass_mismatch": "S3b_mass_mismatch",
             "S5_presence_absence": "S4_presence_absence", "S6_BC_vs_Control": "S5_BC_vs_Control", "S6_BC_vs_PC": "S5_BC_vs_PC",
             "S6_PC_vs_Control": "S5_PC_vs_Control", "S6_PC_vs_maleControl": "S5_PC_vs_maleControl",
             "S7_sensitivity": "S6_sensitivity", "S8_blank_check": "S7_blank_check"}
V5_CSV = {"TableS4_curated_candidates.csv": "TableS3_curated_candidates.csv", "TableS7_sensitivity.csv": "TableS6_sensitivity.csv",
          "TableS8_blank_check.csv": "TableS7_blank_check.csv"}


KEYWORD_CLASSES = ["Paracetamol and related analgesics", "NSAIDs, metamizole and salicylates", "Proton-pump inhibitors and H2-receptor antagonists",
                   "Antidiabetic drugs", "Lipid-lowering drugs", "Calcium-channel blockers, ACE inhibitors and sartans", "Beta-blockers",
                   "Diuretics", "Anticoagulants and antiplatelet drugs",
                   "Urological drugs (alpha-blockers, 5-alpha-reductase inhibitors, antimuscarinics, beta3-agonists, PDE5 inhibitors)",
                   "Hormonal antitumour drugs", "Cytotoxic and targeted antitumour drugs", "Antiepileptic drugs",
                   "Antidepressants, antipsychotics, anxiolytics and hypnotics", "Opioids and local anaesthetics", "Antibacterial and antifungal drugs",
                   "Corticosteroids", "Other drugs", "Dietary methylxanthines and caffeine metabolites",
                   "Other food constituents (polyphenols, alkaloids, terpenoids, betaines and others)", "Tobacco alkaloids"]


def exposure_keywords_table(res=C.OUT):
    """Table S8: список ключевых слов (curation/exposure_keywords.csv) и число соединений, помеченных каждым выражением."""
    import re
    kw = pd.read_csv(C.EXPOSURE_KEYWORDS)
    assert len(kw) == len(KEYWORD_CLASSES)
    ft = pd.read_csv(Path(res) / "tables" / "features.csv", index_col=0)
    ann = ft.annotation.fillna("").astype(str)
    first = pd.Series(-1, index=ft.index)
    for k, r in kw.iterrows():
        first[(first < 0) & ann.str.contains(r.pattern, flags=re.I, regex=True)] = k
    hit = (ft.annotation_conf == "high") & (ft.category == ft.compound_category)
    n = [ft.loc[hit & (first == k) & (ft.compound_category == r.category), "compound"].nunique() for k, r in kw.iterrows()]
    return pd.DataFrame({"No.": range(1, len(kw) + 1), "Category": kw.category.str.capitalize(), "Class": KEYWORD_CLASSES,
                         "Regular expression (case-insensitive)": kw.pattern, "Compounds flagged": n})


def v5_package(out, mv):
    """Пакет для подачи: перенумерация таблиц приложения, сноска к Table S2, раскладка по папкам submission_analytica/."""
    from openpyxl import load_workbook
    root = out.parent
    f = out / "Supplementary_Tables.xlsx"
    wb = load_workbook(f)
    for ws in wb.worksheets:
        ws.title = V5_SHEETS.get(ws.title, ws.title)
    ws = wb["S2_classification"]
    pc = mv[(mv.comparison == "PC_vs_Control") & (mv.set == "all") & mv.model.astype(str).str.contains("скорректированы")]
    note = ("Note: for PC vs Control the exogenous-signal index alone nearly separates the groups (index-only AUC in this table); "
            "the index-adjusted panel AUC for this comparison is therefore an extrapolation beyond the overlap of the two groups "
            "and is not interpreted as a PC-specific signal; it is not shown in Figure 4c."
            + (f" (index-adjusted AUC {pc.AUC_mean.iloc[0]:.2f})" if len(pc) else ""))
    ws.cell(row=ws.max_row + 2, column=1, value=note)
    s8 = exposure_keywords_table()
    s8.to_csv(out / "TableS8_exposure_keywords.csv", index=False)
    ws8 = wb.create_sheet("S8_exposure_keywords")
    ws8.append(list(s8.columns))
    for row in s8.itertuples(index=False):
        ws8.append(list(row))
    wb.save(f)
    for a_, b_ in V5_CSV.items():
        if (out / a_).exists():
            (out / a_).rename(out / b_)
    with open(out / "table_notes.md", "a", encoding="utf-8") as fh:
        fh.write("\n**Table S2 note.** " + note + "\n")
    dirs = {"figures": root / "1_main_figures", "tables": root / "2_main_tables", "supp": root / "3_supplementary"}
    for d in dirs.values():
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    for p_ in out.iterdir():
        n_ = p_.name
        if n_.startswith("Figure") and not n_.startswith("FigureS"):
            shutil.move(str(p_), dirs["figures"] / n_)
        elif n_.startswith(("Table1", "Table2", "Tables_main")):
            shutil.move(str(p_), dirs["tables"] / n_)
        elif n_.startswith(("FigureS", "TableS", "Supplementary_Tables")):
            shutil.move(str(p_), dirs["supp"] / n_)
        else:
            shutil.move(str(p_), root / n_)
    out.rmdir()


def v4_retitle(out):
    """Вариант v4: замена терминов в CSV/XLSX/MD каталога out (заголовки и текстовые ячейки)."""
    from openpyxl import load_workbook

    def fix(x):
        if not isinstance(x, str):
            return x
        if x in V4_CELL_VALUES:
            return V4_CELL_VALUES[x]
        for a_, b_ in V4_TABLE_TERMS:
            x = x.replace(a_, b_)
        return x
    for f in out.glob("*.xlsx"):
        wb = load_workbook(f)
        for ws in wb.worksheets:
            ws.title = ws.title.replace("S3_stage1_markers", "S3_previous_markers")
            for row in ws.iter_rows():
                for c in row:
                    if isinstance(c.value, str):
                        c.value = fix(c.value)
        wb.save(f)
    for f in out.glob("*.csv"):
        d = pd.read_csv(f, dtype=str, keep_default_na=False)
        d.columns = [fix(c) for c in d.columns]
        d = d.map(fix)
        d.to_csv(f, index=False)
    for f in out.glob("table_notes.md"):
        t = f.read_text(encoding="utf-8")
        f.write_text(fix(t), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default=C.OUT)
    ap.add_argument("--data", default=C.DATA)
    ap.add_argument("--variant", default="full", choices=["full", "v4", "v5"])
    a = ap.parse_args()
    res = Path(a.res)
    VARIANT["name"] = a.variant
    if v4():
        TERM.update({"index": "Exogenous-signal index", "adj": "Index-adjusted", "series": "Homologous series (index)",
                     "contaminant": "Background ion", "prior": "Previously reported BC markers"})
        RENAME.update({"Figure2_contamination_index": "FigureS6_exogenous_signal_series", "Figure3_volcano": "Figure2_volcano",
                       "Figure4_top_findings_BC": "Figure3_annotated_findings_BC",
                       "Figure5_ROC_cross_validated": "Figure4_ROC_cross_validated",
                       "Figure6_replication_stage1_markers": "Figure5_previous_BC_markers",
                       "FigureS5_blank_check": "FigureS5_solvent_blanks"})
        if v5():
            TERM["contaminant"] = "Homologous exogenous series"
            V4_CELL_VALUES["contaminant"] = "homologous exogenous series"
        out = Path("submission_analytica") / "_build" if v5() else Path("publication_v4") / "figures_tables"
    else:
        out = res / "publication"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    meta, L, info, feat, stats_, mv, series, cur, pairs, demo, demo_t, pw = load(res)
    X_all, _, info_all, _ = io.load_features(a.data)
    msi, has_ms2 = msi_levels(feat, res)
    nm = names_en(feat, cur, msi)
    if v5():
        cn_ = cur.dropna(subset=["feature"])
        glc = set(cn_.feature[cn_.name_en.fillna("").str.contains("glucuronide", case=False)])
        T2_FIRST.update(feat.index[feat.compound.isin(feat.loc[feat.index.intersection(list(glc)), "compound"])])
    nfull = f"{sum(feat.loc[feat.is_representative.fillna(False).astype(bool), [f'det_{g}' for g in C.STUDY_GROUPS]].min(axis=1) >= 1.0):,}"
    fig1_quality(out, meta, L, feat, X_all)
    fig2_contamination(out, meta, series)
    fig3_volcano(out, stats_, nm, msi)
    fig4_heatmap(out, L, meta, stats_["BC_vs_Control"], nm)
    zp = pd.read_csv(res / "tables" / "contam_matched_pairs.csv", dtype={"case": str, "ref": str}) \
        if (res / "tables" / "contam_matched_pairs.csv").exists() else None
    fig5_roc(out, L, meta, feat, mv, pairs, zp)
    fig6_box(out, L, meta, stats_["BC_vs_Control"], nm, msi, X_all)
    figS1_age(out, meta, pairs)
    figS2_norm(out, meta, res)
    figS3_pathways(out, pw)
    figS5_blanks(out, res, meta)
    t2_full, ts_class = tables(out, demo, demo_t, stats_, nm, cur, mv, feat, msi, has_ms2, meta)
    rep = None if v5() else replication_outputs(out, res, feat, stats_, info_all)
    supplementary_tables(out, res, t2_full, ts_class, rep, feat, msi, has_ms2)
    (out / "figure_captions.md").write_text((CAPTIONS_V5 if v5() else CAPTIONS_V4 if v4() else CAPTIONS).format(
        rsd=C.MAX_QC_RSD, fdr=C.FDR, rep=C.CV_REPEATS, fold=C.CV_FOLDS, k=C.PANEL_K, perm=C.N_PERM,
        cal=C.MATCH_CALIPER_YEARS, rho=C.PANEL_MAX_ABS_RHO_CONTAM, minp=C.MIN_PAIRS_CV, low=C.LOW_DILUTION_FACTOR, ppm=C.PPM, rttol=C.PRIOR_RT_TOL,
        zcal=C.CONTAM_MATCH_CALIPER_SD, nfull=nfull), encoding="utf-8")
    if v4():
        v4_retitle(out)
        for f in out.glob("TableS3_replication_*.csv"):
            f.rename(out / f.name.replace("replication_stage1_markers", "previous_BC_markers").replace("replication_summary", "previous_BC_markers_summary"))
    if v5():
        v5_package(out, mv)
    print(f"Публикационные рисунки и таблицы → {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Выборка 107 проб (РМП и контроль, этап 1 / ABC): рисунки и сводные таблицы отчёта.

    METAB_DATASET=s107 python run_all.py      # обработка (results_s107/)
    python report_s107.py                     # рисунки и таблицы -> report_s107/

Таблицы: демография, 30 маркеров отчёта (Табл. 20 / ABC) в полной выборке, воспроизведение находок выборки s266
(Table 2) в s107, сверка маркеров РПЖ старой обработки (Табл. 22/28/29 отчёта) с обработкой s266, находки s107.
"""
import shutil
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

RES = Path("results_s107/tables")
S266 = Path("results/tables")
OUT = Path("report_s107")
COL = {"Control": "#2a78d6", "BC": "#e8673a", "QC": "#7a7a7a"}
MRK = {"Control": "o", "BC": "s", "QC": "D"}
LAB = {"Control": "Control", "BC": "BC", "QC": "Pooled QC"}
AA = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
PEP = re.compile(rf"\b(?:{AA})(?:-(?:{AA}))+\b|alanyl|glycyl|prolyl|leucyl|valyl|isoleucyl|glutamyl", re.I)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
                     "pdf.fonttype": 42, "savefig.dpi": 300, "legend.frameon": False})
CM = 1 / 2.54

# соединения Table 2 (s266) и их признаки в s107 (те же ионы; сдвиг RT ≈ −0.3 мин, проверено по группам ионов)
S266_T2 = [("11-Oxoetiocholanolone glucuronide", 498.2697, 11.05, 10.73),
           ("Tetrahydrocortisone glucuronide, isomer I", 541.2642, 11.06, 10.76),
           ("Tetrahydrocortisone glucuronide, isomer II", 558.2914, 10.76, 10.43),
           ("Biliverdin", 583.2553, 13.48, 13.13),
           ("LysoPC(16:0)", 496.3397, 16.71, None)]


VERSION, VERSION_DATE = "v2", "2026-11-08"


def load():
    meta = pd.read_csv(RES / "samples.csv", index_col=0)
    meta.index = meta.index.astype(str)
    L = pd.read_csv(RES / "processed_log2_matrix.csv", index_col=0)
    L = L[[c for c in L.columns if c in meta.index]].T
    L.index = L.index.astype(str)
    feat = pd.read_csv(RES / "features.csv", index_col=0)
    r = pd.read_csv(RES / "stats_BC_vs_Control.csv", index_col=0)
    return meta, L, feat, r


def rep_of(feat, f):
    c = feat.loc[f, "compound"]
    reps = feat.index[(feat.compound == c) & feat.is_representative.fillna(False).astype(bool)]
    return reps[0] if len(reps) else f


def find(feat, mz, rt, tol_rt=0.05, tol_ppm=5):
    d = feat[((feat.mz - mz).abs() / mz * 1e6 <= tol_ppm) & ((feat.rt_min - rt).abs() <= tol_rt)]
    return None if not len(d) else (d.rt_min - rt).abs().idxmin()


def fig_batch(meta, feat):
    """Межсерийный сдвиг: PCA проб и QC по признакам, обнаруженным во всех инъекциях, до и после коррекции по QC."""
    raw = pd.read_csv("data/s107/s107_table.csv", low_memory=False).set_index("Sample")
    raw.index = raw.index.astype(str)
    X = raw.drop(columns="Group").astype(float).replace(0, np.nan)
    kept = feat.index[feat.kept.fillna(False).astype(bool)]
    ana = meta.index[meta.study | (meta.kind == "QC")]
    X = X.loc[ana, kept]
    X = X.loc[:, X.notna().all()]
    qc = meta.index[meta.kind == "QC"]
    b = meta.loc[ana, "batch"]
    Y = X.copy()
    ref = X.loc[qc].median()
    for k in sorted(b.unique()):
        Y.loc[b.index[b == k]] = X.loc[b.index[b == k]] * (ref / X.loc[qc[meta.loc[qc, "batch"] == k]].median())
    fig, axs = plt.subplots(1, 2, figsize=(17 * CM, 6.5 * CM), gridspec_kw={"wspace": 0.3})
    bc = {1: "#2a78d6", 2: "#e8673a", 3: "#2ca58d"}
    for ax, (T, ttl) in zip(axs, [(X, "before correction"), (Y, "after QC-based batch correction")]):
        Lg = np.log2(T.div(T.median(axis=1), axis=0))
        mu, sd = Lg.mean(), np.sqrt(Lg.std())
        p = PCA(2).fit((Lg - mu) / sd)
        S = pd.DataFrame(p.transform((Lg - mu) / sd), index=Lg.index)
        for k in (1, 2, 3):
            m = (b == k)
            smp = m & (meta.loc[ana, "kind"] == "sample")
            q = m & (meta.loc[ana, "kind"] == "QC")
            ax.scatter(S.loc[smp, 0], S.loc[smp, 1], s=8, c=bc[k], alpha=0.55, lw=0, label=f"Batch {k}, samples")
            ax.scatter(S.loc[q, 0], S.loc[q, 1], s=26, c=bc[k], marker="D", edgecolor="k", lw=0.5, label=f"Batch {k}, QC")
        v = p.explained_variance_ratio_
        ax.set(xlabel=f"PC1 ({v[0]:.1%})", ylabel=f"PC2 ({v[1]:.1%})", title=f"{ttl.capitalize()} ({X.shape[1]} features)")
        ax.title.set_fontsize(6.5)
    axs[0].legend(fontsize=5.5, ncol=2, loc="best")
    for a, t in zip(axs, "ab"):
        a.text(-0.16, 1.04, f"({t})", transform=a.transAxes, fontweight="bold", fontsize=9)
    fig.savefig(OUT / "FigureS107_0_batches.png", bbox_inches="tight")
    fig.savefig(OUT / "FigureS107_0_batches.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_pca(meta, L, feat):
    reps = feat.index[feat.is_representative.fillna(False).astype(bool)]
    full = [f for f in reps if feat.loc[f, ["det_Control", "det_BC"]].min() >= 1.0]
    ana = meta.index[meta.study | (meta.kind == "QC")]
    X = L.loc[ana, full]
    mu, sd = X.mean(), np.sqrt(X.std())
    p = PCA(3).fit((X - mu) / sd)
    S = pd.DataFrame(p.transform((X - mu) / sd), index=ana)
    grp = meta.loc[ana, "group"]
    fig, ax = plt.subplots(1, 2, figsize=(17 * CM, 6.5 * CM), gridspec_kw={"wspace": 0.35})
    for g in ["Control", "BC", "QC"]:
        m = grp == g
        ax[0].scatter(S.loc[m, 0], S.loc[m, 1], s=12, c=COL[g], marker=MRK[g], edgecolor="white", lw=0.3,
                      label=f"{LAB[g]} (n = {int(m.sum())})")
        ax[1].scatter(meta.loc[S.index[m], "inj_order"], S.loc[m, 0], s=12, c=COL[g], marker=MRK[g], edgecolor="white", lw=0.3)
    v = p.explained_variance_ratio_
    ax[0].set(xlabel=f"PC1 ({v[0]:.1%})", ylabel=f"PC2 ({v[1]:.1%})", title=f"PCA, {len(full)} compounds detected in all samples (Pareto)")
    ax[0].legend(loc="best", fontsize=6)
    ax[1].set(xlabel="Injection order (3 batches, QC every 6 samples)", ylabel=f"PC1 ({v[0]:.1%})")
    st = meta.index[meta.study]
    auc = roc_auc_score((meta.loc[st, "group"] == "BC").astype(int), S.loc[st, 0])
    rho = stats.spearmanr(S.loc[st, 0], meta.loc[st, "inj_order"]).statistic
    ax[1].set_title(f"PC1: AUC BC vs Control {max(auc, 1 - auc):.2f}; ρ with order {rho:+.2f}", fontsize=6.5)
    for a, t in zip(ax, "ab"):
        a.text(-0.16, 1.04, f"({t})", transform=a.transAxes, fontweight="bold", fontsize=9)
    fig.savefig(OUT / "FigureS107_1_PCA.png", bbox_inches="tight")
    fig.savefig(OUT / "FigureS107_1_PCA.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(full), max(auc, 1 - auc), rho, v


def fig_volcano(r):
    pep = r.annotation.astype(str).str.contains(PEP)
    fig, ax = plt.subplots(figsize=(9 * CM, 7 * CM))
    y = -np.log10(r.p_model)
    ax.scatter(r.beta_model, y, s=3, c="#c9c8c4", lw=0)
    f = r.finding.astype(bool)
    ax.scatter(r.beta_model[f & ~pep], y[f & ~pep], s=6, c=COL["BC"], lw=0, label="Finding")
    ax.scatter(r.beta_model[f & pep], y[f & pep], s=9, c="#5b3fb5", lw=0, label="Finding, peptide library match")
    thr = r.p_model[r.q_model < 0.05].max()
    ax.axhline(-np.log10(thr), ls="--", lw=0.6, c="#555")
    ax.set(xlabel="log$_2$ fold change (adjusted)", ylabel="−log$_{10}$ p (adjusted model)",
           title=f"BC vs Control, 107-sample cohort ({int(f.sum())} findings)")
    ax.legend(fontsize=6, loc="upper left")
    fig.savefig(OUT / "FigureS107_2_volcano.png", bbox_inches="tight")
    fig.savefig(OUT / "FigureS107_2_volcano.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_box(meta, L, feat, items, fname, title_note=""):
    st = meta.index[meta.study]
    n = len(items)
    fig, axs = plt.subplots(1, n, figsize=(3.6 * n * CM, 5.8 * CM), gridspec_kw={"wspace": 0.55})
    for ax, (lab, f, sub) in zip(np.atleast_1d(axs), items):
        for i, g in enumerate(["Control", "BC"]):
            v = L.loc[st[meta.loc[st, "group"] == g], f]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, medianprops=dict(color="k", lw=0.9))
            ax.scatter(np.random.RandomState(i).normal(i, 0.07, len(v)), v, s=4, c=COL[g], marker=MRK[g], lw=0)
        ax.set_xticks([0, 1], ["Control", "BC"])
        ax.set_title(f"{lab}\n{sub}", fontsize=5.8)
    np.atleast_1d(axs)[0].set_ylabel("log$_2$ normalized intensity")
    fig.savefig(OUT / f"{fname}.png", bbox_inches="tight")
    fig.savefig(OUT / f"{fname}.pdf", bbox_inches="tight")
    plt.close(fig)


def table_s266_replication(feat, r):
    s266 = {n: pd.read_csv(S266 / f"stats_{n}.csv", index_col=0) for n in ["BC_vs_Control"]}["BC_vs_Control"]
    rows = []
    for name, mz, rt266, rt107 in S266_T2:
        f266 = s266.index[((s266.mz - mz).abs() < 0.002) & ((s266.rt_min - rt266).abs() < 0.03)]
        a = s266.loc[f266[0]] if len(f266) else None
        f = find(feat, mz, rt107, 0.05) if rt107 else None
        if f is None:
            rows.append(dict(Compound=name, **{"m/z": mz}, **{"RT s266": rt266}, **{"RT s107": ""},
                             **{"log2FC s266": a.beta_model if a is not None else np.nan},
                             **{"q s266": a.q_model if a is not None else np.nan},
                             Status_s107="not detected after filters (no feature ±5 ppm within ±0.6 min)"))
            continue
        b = r.loc[rep_of(feat, f)]
        rows.append({"Compound": name, "m/z": mz, "RT s266": rt266, "RT s107": feat.loc[f, "rt_min"],
                     "log2FC s266": a.beta_model, "q s266": a.q_model, "AUC s266": a.AUC,
                     "log2FC s107": b.beta_model, "q s107": b.q_model, "AUC s107": b.AUC, "p matched s107": b.p_matched,
                     "p low-index ctrl s107": b.p_lowindex, "Finding s107": bool(b.finding),
                     "Status_s107": "same direction, q < 0.05" if (np.sign(b.beta_model) == np.sign(a.beta_model) and b.q_model < 0.05)
                     else ("same direction, n.s." if np.sign(b.beta_model) == np.sign(a.beta_model) else "opposite direction")})
    return pd.DataFrame(rows)


def table_abc(r, feat):
    rep = pd.read_csv(RES / "replication_markers.csv")
    cur = pd.read_csv(RES / "candidates_curated.csv", index_col=0) if (RES / "candidates_curated.csv").exists() else None
    rep = rep.sort_values("q_model").drop_duplicates("marker", keep="first")      # один (наиболее значимый) признак на маркер
    t20 = pd.read_csv("data/s107/table_nastya.csv")
    order = {c: i for i, c in enumerate(t20.Compound)}
    rep = rep.assign(_o=rep.marker.map(order)).sort_values("_o")
    from pipeline.features import formula_mass, ADDUCTS
    def ppm_formula(x):
        if x.ion not in ADDUCTS:
            return np.nan
        th = formula_mass(x.formula) + ADDUCTS[x.ion]
        return (x.mz - th) / th * 1e6
    out = pd.DataFrame({
        "No. (report Table 20)": rep.marker.map(lambda m: int(t20.loc[t20.Compound == m, "No"].iloc[0])),
        "Compound (report)": rep.marker, "Reported direction": rep.direction, "m/z": rep.mz, "RT": rep.rt_min, "Ion": rep.ion,
        "Formula": rep.formula, "Formula vs m/z, ppm": rep.apply(ppm_formula, axis=1).round(1),
        "Category (this analysis)": rep.category, "log2FC (adjusted)": rep.beta_model if "beta_model" in rep else rep.log2FC,
        "q (adjusted)": rep.q_model, "p (matched pairs, Wilcoxon)": rep.p_matched, "AUC": rep.AUC,
        "Same direction (adjusted)": np.sign(rep.beta_model) == rep.direction.map({"up": 1, "down": -1}),
        "Finding (all criteria)": rep.finding_this_study})
    return out


def table_findings(r, feat):
    f = r[r.finding.astype(bool)].copy()
    f["peptide_library_match"] = f.annotation.astype(str).str.contains(PEP)
    cols = ["rt_min", "mz", "annotation", "annotation_conf", "compound_category", "beta_model", "q_model", "p_matched",
            "p_lowindex", "p_ctrl_hi_vs_lo", "finding_via_ctrl_not_detected", "AUC", "det_Control", "det_BC",
            "peptide_library_match"]
    return f[[c for c in cols if c in f]].sort_values("q_model")


def main():
    OUT.mkdir(exist_ok=True)
    meta, L, feat, r = load()
    fig_batch(meta, feat)
    nfull, auc_pc1, rho_pc1, var = fig_pca(meta, L, feat)
    fig_volcano(r)
    # воспроизведение Table 2 (s266) и пептиды
    t2 = table_s266_replication(feat, r)
    items = []
    for name, mz, rt266, rt107 in S266_T2:
        f = find(feat, mz, rt107, 0.05) if rt107 else None
        if f is not None:
            b = r.loc[rep_of(feat, f)]
            items.append((name.replace("Tetrahydrocortisone", "THE").replace(", isomer", " isomer"),
                          rep_of(feat, f), f"m/z {mz:.4f} · q = {b.q_model:.3f}"))
    fig_box(meta, L, feat, items, "FigureS107_3_s266_findings_in_s107")
    pep_items = []
    for lab, mz, rt in [("Val-Leu", 231.1701, 4.145), ("Gln-Pro (report: N-mevalonylcreatinine)", 244.1293, 2.578),
                        ("Glu-Ile-Val", 360.2129, 7.951), ("Thr-Phe", 267.1337, 3.432)]:
        f = find(feat, mz, rt, 0.05)
        if f is not None:
            b = r.loc[rep_of(feat, f)]
            pep_items.append((lab, rep_of(feat, f), f"AUC {b.AUC:.2f} · det. C/BC {b.det_Control:.0%}/{b.det_BC:.0%}"))
    fig_box(meta, L, feat, pep_items, "FigureS107_4_peptides")
    abc = table_abc(r, feat)
    fnd = table_findings(r, feat)
    pc_old = pd.read_csv("results_s107_extra/PC_markers_old_report_vs_s266.csv")
    # пептиды среди находок: обогащение (Fisher)
    pep = r.annotation.astype(str).str.contains(PEP)
    up = (r.q_model < 0.05) & (r.beta_model > 0)
    ct = [[int((pep & up).sum()), int((~pep & up).sum())], [int((pep & ~up).sum()), int((~pep & ~up).sum())]]
    fisher = stats.fisher_exact(ct)
    mv = pd.read_csv(RES / "multivariate_CV.csv")
    version = pd.DataFrame({"field": ["version", "date", "changes vs v1"], "value": [
        VERSION, VERSION_DATE,
        "1) Finding rule: if a compound is not detected in any control, the comparison of high- vs low-index controls is "
        "undefined (p = NaN) and the criterion is considered met (column finding_via_ctrl_not_detected); "
        "2) categories transferred from report Table 20 / ABC annotations are kept only when confirmed by data "
        "(MS-DIAL annotation of the same feature names the same compound); otherwise 'ambiguous' "
        "(incl. 269.0880 @ 2.02 min and 5-hydroxyomeprazole)."]})
    with pd.ExcelWriter(OUT / "Tables_s107.xlsx") as w:
        version.to_excel(w, sheet_name="Version", index=False)
        abc.to_excel(w, sheet_name="T1_report_markers_in_s107", index=False)
        t2.to_excel(w, sheet_name="T2_s266_Table2_in_s107", index=False)
        pc_old.to_excel(w, sheet_name="T3_old_PC_markers_vs_s266", index=False)
        fnd.to_excel(w, sheet_name="S_all_findings_s107")
        mv.to_excel(w, sheet_name="S_classification_s107", index=False)
    summ = dict(n_full_compounds_pca=nfull, pc1_auc=auc_pc1, pc1_rho_order=rho_pc1, pc1_var=var[0],
                n_findings=int(r.finding.sum()), n_q05=int((r.q_model < 0.05).sum()),
                n_up=int(((r.q_model < 0.05) & (r.beta_model > 0)).sum()),
                pep_up=ct[0][0], pep_total=int(pep.sum()), n_total=len(r), fisher_or=fisher.statistic, fisher_p=fisher.pvalue,
                n_find_up=int((r.finding.astype(bool) & (r.beta_model > 0)).sum()),
                n_find_down=int((r.finding.astype(bool) & (r.beta_model < 0)).sum()),
                n_find_pep=int((r.finding.astype(bool) & pep).sum()))
    pd.Series(summ).to_csv(OUT / "summary_numbers.csv")
    shutil.copy(OUT / "Tables_s107.xlsx", OUT / f"Tables_s107_{VERSION}.xlsx")
    print(pd.Series(summ).to_string())


if __name__ == "__main__":
    main()

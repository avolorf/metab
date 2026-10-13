#!/usr/bin/env python3
"""Аннотация, визуализация и экспорт для внешних ресурсов (MetaboAnalyst, HMDB, METLIN, KEGG).

Запускается после process.py:
    python annotate_visualize.py --res results --seq data/sequence.xlsx

Делает:
  * проверку курированного списка кандидатов (curation/candidates_curated.csv): ошибка массы (ppm),
    сверка HMDB/KEGG ID с формулой по офлайн-словарям mummichog (HMDB ~114 тыс., KEGG-модель человека);
  * поиск HMDB-кандидатов по точной массе для всех значимых признаков;
  * mummichog (тот же алгоритм, что в модуле MetaboAnalyst «Functional Analysis / MS Peaks to Pathways»);
  * тепловые карты, ROC-кривые, сравнение BC и PC, ряд контаминантов;
  * файлы для загрузки в MetaboAnalyst, METLIN, HMDB, KEGG Mapper.
"""
import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve

from process import GROUP_COLORS, plsda_cv_predictions  # noqa: E402

PPM = 5.0
PROTON = 1.007276
ADDUCTS = {  # m/z = M + shift
    "[M+H]+": PROTON, "[M+Na]+": 22.989218, "[M+K]+": 38.963158, "[M+NH4]+": 18.033823,
    "[M+H-H2O]+": PROTON - 18.010565, "[M+H-2H2O]+": PROTON - 2 * 18.010565,
}
ELEM = {"C": 12.0, "H": 1.00782503, "N": 14.0030740, "O": 15.9949146, "P": 30.9737620,
        "S": 31.9720707, "Cl": 34.9688527, "Br": 78.9183371, "Na": 22.9897693, "K": 38.9637064}
SOURCE_COLORS = {"ДР": "#e41a1c", "ОД": "#ff7f00", "МНН": "#4daf4a", "RR": "#377eb8", "РМП": "#984ea3",
                 "КГ": "#a6cee3", "num": "#b2df8a", "Вместо": "#e41a1c", "QC": "#7f7f7f"}
CAT_COLORS = {"эндогенный": "#1b9e77", "пища": "#d95f02", "лекарство": "#7570b3",
              "контаминант": "#e7298a", "неоднозначная": "#999999"}


def formula_mass(f):
    return sum(ELEM[e] * (int(n) if n else 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", f))


def short(name, n=38):
    return name if len(name) <= n else name[: n - 1] + "…"


def cat_color(c):
    for k, v in CAT_COLORS.items():
        if str(c).startswith(k):
            return v
    return "#666666"


# ----------------------------------------------------------------------------- загрузка
def load(res, seq_path):
    M = pd.read_csv(res / "processed_log2_matrix.csv", index_col=0)
    info = M[["rt_min", "mz", "annotation", "annotation_conf"]]
    L = M.drop(columns=info.columns).T
    L.index = L.index.astype(str)
    uni = pd.read_csv(res / "univariate_all.csv", index_col=0)
    meta = pd.read_excel(res / "metabolomics_results.xlsx", sheet_name="samples", index_col=0)
    meta.index = meta.index.astype(str)
    seq = pd.read_excel(seq_path)
    seq.index = seq.iloc[:, 0].astype("Int64").astype(str)
    src = seq["Шифр"].astype(str).str.extract(r"^([A-Za-zА-Яа-я]+)")[0].fillna("num")
    meta["source"] = src.reindex(meta.index).fillna("QC")
    return L, info, uni, meta


def match_feature(info, rt, mz):
    d = info[(info.mz - mz).abs() < 0.002]
    d = d.iloc[(d.rt_min - rt).abs().argsort()]
    if len(d) and abs(d.rt_min.iloc[0] - rt) < 0.02:
        return d.index[0]
    return None


# ----------------------------------------------------------------------------- курация
def check_curation(path, info, uni):
    from mummichog.compound_dicts import hmdb_dict, kegg_dict
    cur = pd.read_csv(path)
    cur["name_en"] = cur["name_en"].fillna("")
    rows = []
    for _, r in cur.iterrows():
        f = match_feature(info, r.rt_min, r.mz)
        m = formula_mass(r.formula)
        theo = m + ADDUCTS[r.ion]
        ppm = (r.mz - theo) / theo * 1e6
        ok_mass = abs(ppm) <= 10
        hmdb = r.hmdb if isinstance(r.hmdb, str) else ""
        kegg = r.kegg if isinstance(r.kegg, str) else ""
        hmdb_ok = bool(hmdb) and hmdb in hmdb_dict and abs(hmdb_dict[hmdb][1] - m) < 0.005
        kegg_ok = bool(kegg) and kegg in kegg_dict and abs(float(kegg_dict[kegg][1] or 0) - m) < 0.005
        d = dict(feature=f, **r.to_dict())
        d.update(theo_mz=round(theo, 5), ppm=round(ppm, 2), mass_ok=ok_mass,
                 hmdb=hmdb if (hmdb_ok and ok_mass) else "",
                 kegg=kegg if (kegg_ok and ok_mass) else "")
        if not ok_mass and not str(r.category).startswith("неоднозначная"):
            d["category"] = "неоднозначная аннотация"
        for comp in ["BC_vs_Control", "PC_vs_Control", "PC_vs_maleControl", "BC_vs_PC"]:
            u = uni[uni.comparison == comp]
            if f in u.index:
                d[f"{comp}_log2FC"] = u.loc[f, "log2FC_mean"]
                d[f"{comp}_AUC"] = u.loc[f, "AUC"]
                d[f"{comp}_q"] = u.loc[f, "q_MW"]
                d[f"{comp}_q_sexadj"] = u.loc[f, "q_sex_adj"]
        rows.append(d)
    return pd.DataFrame(rows)


def hmdb_mass_candidates(mzs, ppm=PPM):
    """Для каждого m/z — число HMDB-структур и формулы-кандидаты по основным аддуктам."""
    from mummichog.compound_dicts import hmdb_dict
    ids = np.array(list(hmdb_dict.keys()))
    form = np.array([v[0] for v in hmdb_dict.values()])
    mass = np.array([float(v[1]) for v in hmdb_dict.values()])
    o = np.argsort(mass)
    ids, form, mass = ids[o], form[o], mass[o]
    out = []
    for mz in mzs:
        hits = []
        for ad in ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+H-H2O]+"]:
            M = mz - ADDUCTS[ad]
            tol = M * ppm * 1e-6
            lo, hi = np.searchsorted(mass, [M - tol, M + tol])
            for i in range(lo, hi):
                hits.append((ad, form[i], ids[i], (M - mass[i]) / mass[i] * 1e6))
        h = pd.DataFrame(hits, columns=["adduct", "formula", "hmdb", "ppm"])
        if len(h):
            g = h.groupby(["adduct", "formula"]).agg(n=("hmdb", "size"), ids=("hmdb", lambda s: ";".join(s[:5])),
                                                     ppm=("ppm", "first")).reset_index().sort_values("n", ascending=False)
            txt = " | ".join(f"{a} {f} ({p:+.1f} ppm, {n} шт.: {i})" for a, f, n, i, p in
                             g[["adduct", "formula", "n", "ids", "ppm"]].head(4).itertuples(index=False))
            out.append((len(h), h.formula.nunique(), txt))
        else:
            out.append((0, 0, ""))
    return pd.DataFrame(out, columns=["HMDB_n_hits", "HMDB_n_formulas", "HMDB_top_candidates"], index=mzs.index)


# ----------------------------------------------------------------------------- mummichog
def peak_table(uni, comp, exclude):
    """m/z, RT, p, знаковая статистика. p и знак – из модели с поправкой на пол (если она есть).
    Знаковая статистика = sign(log2FC) * -log10(p): mummichog использует её только для отчёта,
    GSEA в MetaboAnalyst – для ранжирования."""
    d = uni[(uni.comparison == comp) & ~uni.index.isin(exclude)]
    adj = d.p_sex_adj.notna().all()
    p = d.p_sex_adj if adj else d.p_MW
    fc = d.log2FC_sex_adj if adj else d.log2FC_mean
    return pd.DataFrame({"mz": d.mz, "rt_sec": d.rt_min * 60, "p": p, "t": np.sign(fc) * -np.log10(p),
                         "id": d.index}), ("линейная модель с поправкой на пол" if adj else "Манн-Уитни")


def run_mummichog(uni, comp, outdir, exclude, top_frac=0.10, n_perm=500):
    """mummichog 2: значимыми считаются top_frac признаков с наименьшим p (как top 10 % в MetaboAnalyst)."""
    tab, ptype = peak_table(uni, comp, exclude)
    cutoff = float(np.quantile(tab.p, top_frac))
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tab.to_csv(Path(tmp) / "input.txt", sep="\t", index=False)
        subprocess.run(["mummichog", "-f", "input.txt", "-o", comp, "-m", "positive", "-u", str(PPM),
                        "-c", f"{cutoff:.3g}", "-p", str(n_perm)], cwd=tmp, check=True, capture_output=True)
        res = next(p for p in Path(tmp).iterdir() if p.is_dir())
        shutil.rmtree(outdir)
        # только таблицы и PDF-графики (HTML-отчёт тянет ~10 МБ plotly)
        for sub in ["tables", "figures"]:
            shutil.copytree(res / sub, outdir / sub)
    pw = pd.read_csv(outdir / "tables" / f"mcg_pathwayanalysis_{comp}.tsv", sep="\t")
    settings = dict(comparison=comp, p_source=ptype, n_features=len(tab), n_significant=int((tab.p < cutoff).sum()),
                    p_cutoff=cutoff, excluded_contaminants=len(exclude), permutations=n_perm, ppm=PPM,
                    mode="positive", model="human MFN (mummichog 2.7)")
    return pw, settings


def plot_pathways(pw, path, title, top=15):
    pw = pw.sort_values("p-value").head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(pw) + 1.2))
    col = ["#d7301f" if p < 0.05 else "#bbbbbb" for p in pw["p-value"]]
    ax.barh(pw["pathway"], -np.log10(pw["p-value"]), color=col)
    for y, (o, s) in enumerate(zip(pw["overlap_size"], pw["pathway_size"])):
        ax.text(0.02, y, f"{o}/{s}", va="center", fontsize=7)
    ax.axvline(-np.log10(0.05), ls="--", c="k", lw=0.6)
    ax.set_xlabel("-log10 p (mummichog, перестановочный)")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------- графики
def clustermap(L, meta, feats, labels, samples, path, title):
    Z = L.loc[samples, feats]
    Z = ((Z - Z.mean()) / Z.std()).clip(-3, 3).T
    Z.index = labels
    sex_col = {"м": "#1f78b4", "ж": "#fb9a99"}
    col_colors = pd.DataFrame({
        "группа": meta.loc[samples, "class"].map(GROUP_COLORS),
        "пол": meta.loc[samples, "sex"].map(sex_col).fillna("#ffffff"),
        "источник": meta.loc[samples, "source"].map(SOURCE_COLORS).fillna("#ffffff"),
    }, index=samples)
    g = sns.clustermap(Z, cmap="RdBu_r", center=0, vmin=-3, vmax=3, method="average", metric="correlation",
                       col_colors=col_colors, xticklabels=False, yticklabels=True,
                       figsize=(14, 0.22 * len(feats) + 4), dendrogram_ratio=(0.08, 0.08),
                       cbar_pos=(0.005, 0.03, 0.012, 0.1), cbar_kws={"label": "z-score"})
    g.ax_heatmap.tick_params(axis="y", labelsize=7)
    g.fig.suptitle(title, y=1.0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in GROUP_COLORS.values()]
    groups_present = [k for k in GROUP_COLORS if (meta.loc[samples, "class"] == k).any()]
    g.ax_col_dendrogram.legend([handles[list(GROUP_COLORS).index(k)] for k in groups_present], groups_present,
                               loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7, title="группа")
    g.savefig(path, dpi=150)
    plt.close(g.fig)


def group_mean_heatmap(L, meta, cur, path):
    groups = ["Control", "BC", "PC", "KC", "PanC"]
    feats = cur.feature.tolist()
    smp = meta.index[meta["class"].isin(groups)]
    Z = L.loc[smp, feats]
    Z = (Z - Z.mean()) / Z.std()
    G = Z.groupby(meta.loc[smp, "class"]).mean().reindex(groups).T
    G.index = [f"{short(n, 45)}  [{r:.2f}/{m:.3f}]" for n, r, m in zip(cur.name, cur.rt_min, cur.mz)]
    fig, ax = plt.subplots(figsize=(9, 0.28 * len(G) + 1.5))
    sns.heatmap(G, cmap="RdBu_r", center=0, vmin=-1.2, vmax=1.2, ax=ax, cbar_kws={"label": "средний z-score"},
                linewidths=0.3, linecolor="white")
    for i, c in enumerate(cur.category):
        ax.add_patch(plt.Rectangle((-0.35, i), 0.3, 1, color=cat_color(c), clip_on=False, transform=ax.transData))
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlabel("")
    handles = [plt.Rectangle((0, 0), 1, 1, color=v) for v in CAT_COLORS.values()]
    ax.legend(handles, CAT_COLORS.keys(), bbox_to_anchor=(1.25, 1), loc="upper left", fontsize=7, title="категория")
    ax.set_title("Курированные кандидаты: среднее по группам (z-score)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def roc_panel(L, meta, case_cls, cur, path, title):
    smp = meta.index[meta["class"].isin([case_cls, "Control"])]
    y = (meta.loc[smp, "class"] == case_cls).astype(int).to_numpy()
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for _, r in cur.iterrows():
        s = L.loc[smp, r.feature].to_numpy()
        auc = roc_auc_score(y, s)
        if auc < 0.5:
            s, auc = -s, 1 - auc
        fpr, tpr, _ = roc_curve(y, s)
        ax.plot(fpr, tpr, lw=1.2, label=f"{short(r['name'], 30)} ({auc:.2f})")
    pred = plsda_cv_predictions(L.loc[smp].to_numpy(), y)
    fpr, tpr, _ = roc_curve(y, pred)
    ax.plot(fpr, tpr, c="k", lw=2.5, label=f"PLS-DA, все признаки, CV ({roc_auc_score(y, pred):.2f})")
    ax.plot([0, 1], [0, 1], ls=":", c="grey")
    ax.set_xlabel("1 − специфичность")
    ax.set_ylabel("чувствительность")
    ax.set_title(title)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fc_scatter(uni, cur, path):
    a = uni[uni.comparison == "BC_vs_Control"]
    b = uni[uni.comparison == "PC_vs_Control"].reindex(a.index)
    sa, sb = a.q_MW < 0.05, b.q_MW < 0.05
    fig, ax = plt.subplots(figsize=(7, 6.5))
    ax.scatter(a.log2FC_mean[~sa & ~sb], b.log2FC_mean[~sa & ~sb], s=3, c="#dddddd", lw=0)
    ax.scatter(a.log2FC_mean[sa & ~sb], b.log2FC_mean[sa & ~sb], s=6, c="#d7301f", lw=0, label=f"только BC ({(sa & ~sb).sum()})")
    ax.scatter(a.log2FC_mean[~sa & sb], b.log2FC_mean[~sa & sb], s=10, c="#fc8d59", lw=0, label=f"только PC ({(~sa & sb).sum()})")
    ax.scatter(a.log2FC_mean[sa & sb], b.log2FC_mean[sa & sb], s=14, c="k", lw=0, label=f"оба ({(sa & sb).sum()})")
    for _, r in cur[~cur.category.str.startswith("неоднозн")].iterrows():
        if r.feature in a.index:
            ax.annotate(short(r["name"], 22), (a.log2FC_mean[r.feature], b.log2FC_mean[r.feature]), fontsize=6,
                        color=cat_color(r.category))
    ax.axhline(0, c="grey", lw=0.5)
    ax.axvline(0, c="grey", lw=0.5)
    ax.set_xlabel("log2FC BC / Control")
    ax.set_ylabel("log2FC PC / Control")
    ax.set_title("Сравнение изменений при раке мочевого пузыря и простаты (q<0.05)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def contaminant_plot(L, meta, cur, path):
    c = cur[cur.category.str.startswith("контаминант")].drop_duplicates("feature")
    smp = meta.index[meta["class"] != "QC"]
    Z = L.loc[smp, c.feature]
    score = ((Z - Z.mean()) / Z.std()).mean(axis=1)
    df = pd.DataFrame({"score": score, "source": meta.loc[smp, "source"].replace({"Вместо": "ДР"}),
                       "class": meta.loc[smp, "class"]})
    order = ["КГ", "num", "RR", "МНН", "РМП", "ОД", "ДР"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
    for i, s in enumerate(order):
        d = df[df.source == s]
        ax[0].boxplot(d.score, positions=[i], widths=0.6, showfliers=False)
        ax[0].scatter(np.random.normal(i, 0.07, len(d)), d.score, c=d["class"].map(GROUP_COLORS), s=18,
                      edgecolor="k", lw=0.3)
    ax[0].set_xticks(range(len(order)), [f"{s}\n(n={(df.source == s).sum()})" for s in order])
    ax[0].set_xlabel("источник проб (префикс шифра)")
    ax[0].set_ylabel("индекс ряда ПАВ (средний z-score)")
    ax[0].set_title("Ряд алкилдиэтаноламинов/этоксилатов по источникам")
    handles = [plt.Line2D([], [], marker="o", ls="", color=GROUP_COLORS[k]) for k in ["Control", "BC", "PC", "KC", "PanC"]]
    ax[0].legend(handles, ["Control", "BC", "PC", "KC", "PanC"], fontsize=7)
    cs = c.sort_values("mz")
    ax[1].scatter(cs.mz, cs.rt_min, c=cat_color("контаминант"), s=40)
    for _, r in cs.iterrows():
        ax[1].annotate(f"{r.formula}", (r.mz, r.rt_min), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax[1].set_xlabel("m/z")
    ax[1].set_ylabel("RT, мин")
    ax[1].set_title("Гомологический ряд: шаги C2H4 (28.031), O (15.995), C2H4O (44.026)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return df


# ----------------------------------------------------------------------------- экспорт
SRC_LATIN = {"ДР": "DR", "Вместо": "DR", "ОД": "OD", "МНН": "MNN", "РМП": "RMP", "КГ": "KG", "RR": "RR",
             "num": "num", "QC": "QC"}


def export_external(L, info, uni, meta, cur, out, n_iqr=5000):
    out.mkdir(parents=True, exist_ok=True)
    for f in out.iterdir():
        f.unlink()
    names = (info.mz.round(4).astype(str) + "__" + info.rt_min.round(3).astype(str))
    X = 2 ** L  # PQN-нормализованные, импутированные (без log)
    X.columns = names[X.columns]
    smp_all = meta.index[meta["class"] != "QC"]
    # неконтролируемый фильтр по IQR (как в MetaboAnalyst) — не использует метки групп
    iqr = L.loc[smp_all].quantile(0.75) - L.loc[smp_all].quantile(0.25)
    top_iqr = names[iqr.sort_values(ascending=False).index[:n_iqr]]

    # 1. Statistical Analysis [one factor] и Biomarker Analysis
    sets = {"BC_vs_Control": ["BC", "Control"], "PC_vs_Control": ["PC", "Control"],
            "BC_PC_Control": ["BC", "PC", "Control"]}
    for n, groups in sets.items():
        smp = meta.index[meta["class"].isin(groups)]
        t = X.loc[smp].copy()
        t.insert(0, "Label", meta.loc[smp, "class"])
        t.index.name = "Sample"
        t.to_csv(out / f"1_stats_{n}_all{X.shape[1]}.csv", float_format="%.5g")
        t[["Label"] + list(top_iqr)].to_csv(out / f"1_stats_{n}_IQRtop{n_iqr}.csv", float_format="%.5g")

    # 2. Statistical Analysis [metadata table] — поправка на пол (2 образца без пола исключены)
    smp = meta.index[meta["class"].isin(["BC", "PC", "Control"]) & meta["sex"].notna()]
    d = X.loc[smp, top_iqr].copy()
    d.index.name = "Sample"
    d.to_csv(out / f"2_metadata_module_data_IQRtop{n_iqr}.csv", float_format="%.5g")
    md = pd.DataFrame({"Class": meta.loc[smp, "class"],
                       "Sex": meta.loc[smp, "sex"].map({"м": "M", "ж": "F"}),
                       "Source": meta.loc[smp, "source"].map(SRC_LATIN)}, index=smp)
    md.index.name = "Sample"
    md.to_csv(out / "2_metadata_module_metadata.csv")

    # 3. Functional Analysis (MS peaks to pathways): m.z, p.value, t.score, r.t
    contam = set(cur.loc[cur.category.str.startswith("контаминант"), "feature"])
    for comp in ["BC_vs_Control", "PC_vs_Control", "BC_vs_PC"]:
        tab, _ = peak_table(uni, comp, contam)
        tab.rename(columns={"mz": "m.z", "p": "p.value", "t": "t.score", "rt_sec": "r.t"})[
            ["m.z", "p.value", "t.score", "r.t"]].to_csv(out / f"3_peaks_{comp}.txt", sep="\t", index=False)

    # 4–5. Enrichment Analysis / Pathway Analysis — списки соединений, по одному типу ID на файл
    good = cur[~cur.category.str.startswith(("неоднозн", "контаминант"))]
    for focus in ["BC", "PC"]:
        g = good[good.focus == focus]
        (out / f"4_compounds_{focus}_HMDB.txt").write_text("\n".join(g.hmdb[g.hmdb != ""].drop_duplicates()))
        (out / f"4_compounds_{focus}_names.txt").write_text("\n".join(g.name_en[g.name_en != ""].drop_duplicates()))
        k = g[g.kegg != ""]
        col = np.where(k[f"{focus}_vs_Control_log2FC"] > 0, "red", "blue")
        (out / f"6_kegg_mapper_{focus}.txt").write_text("\n".join(f"{a} {c}" for a, c in zip(k.kegg, col)))
    # 7. METLIN / HMDB LC-MS search
    sig = uni[(uni.q_MW < 0.05) & uni.comparison.isin(["BC_vs_Control", "PC_vs_Control"])]
    sig = sig[~sig.index.duplicated() & ~sig.index.isin(contam)]
    (out / "7_metlin_hmdb_mz.txt").write_text("\n".join(f"{m:.5f}" for m in sig.mz.sort_values()))
    sig[["mz", "rt_min", "annotation"]].sort_values("mz").to_csv(out / "7_metlin_hmdb_mz_rt.csv")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results")
    ap.add_argument("--seq", default="data/sequence.xlsx")
    ap.add_argument("--curation", default="curation/candidates_curated.csv")
    ap.add_argument("--no-mummichog", action="store_true")
    a = ap.parse_args()
    res = Path(a.res)
    fig = res / "figures"
    L, info, uni, meta = load(res, a.seq)

    cur = check_curation(a.curation, info, uni)
    print(f"Курированных записей: {len(cur)}, не найдено в матрице: {cur.feature.isna().sum()}, "
          f"масса не сходится (>10 ppm): {(~cur.mass_ok).sum()}")
    cur = cur.dropna(subset=["feature"])

    # таблица значимых признаков с HMDB-кандидатами по массе
    sig_tabs = {}
    for comp in ["BC_vs_Control", "PC_vs_Control", "BC_vs_PC"]:
        u = uni[(uni.comparison == comp) & (uni.q_MW < 0.05)].copy()
        u = u.join(hmdb_mass_candidates(u.mz))
        u["curated_name"] = u.index.map(cur.drop_duplicates("feature").set_index("feature")["name"])
        u["curated_category"] = u.index.map(cur.drop_duplicates("feature").set_index("feature")["category"])
        sig_tabs[comp] = u

    # mummichog
    pws, mcg_settings = {}, []
    contam = set(cur.loc[cur.category.str.startswith("контаминант"), "feature"])
    if not a.no_mummichog:
        for comp in ["BC_vs_Control", "PC_vs_Control", "BC_vs_PC"]:
            pw, st = run_mummichog(uni, comp, res / "mummichog" / comp, contam)
            pws[comp] = pw
            mcg_settings.append(st)
            plot_pathways(pw, fig / f"10_pathways_{comp}.png",
                          f"mummichog: {comp} (топ-10 %: p<{st['p_cutoff']:.2g}, {st['n_significant']} пиков)")
            print(f"mummichog {comp}: {st}; путей с p<0.05 — {(pw['p-value'] < 0.05).sum()}")

    # тепловые карты
    for comp, case, ctrl in [("BC_vs_Control", ["BC"], ["Control"]), ("PC_vs_Control", ["PC"], ["Control"]),
                             ("BC_vs_PC", ["BC"], ["PC"])]:
        u = uni[uni.comparison == comp]
        u = u[(u.q_MW < 0.05) & ((u.q_sex_adj < 0.05) | u.q_sex_adj.isna())].sort_values("p_MW").head(60)
        cn = cur.drop_duplicates("feature").set_index("feature")["name"]
        labels = [f"{short(cn.get(f, r.annotation), 40)} [{r.rt_min:.2f}/{r.mz:.3f}]"
                  + ("" if r.annotation_conf == "high" or f in cn.index else " *") for f, r in u.iterrows()]
        smp = meta.index[meta["class"].isin(case + ctrl)]
        clustermap(L, meta, u.index.tolist(), labels, smp, fig / f"08_heatmap_{comp}.png",
                   f"{comp}: топ-{len(u)} признаков (q<0.05, устойчивы к поправке на пол); * — low score/unknown")

    good = cur[~cur.category.str.startswith("неоднозн")].drop_duplicates("feature")
    group_mean_heatmap(L, meta, good.sort_values(["focus", "category"]), fig / "09_heatmap_curated_groupmeans.png")

    for focus in ["BC", "PC"]:
        c = cur[(cur.focus == focus) & ~cur.category.str.startswith(("неоднозн",))].drop_duplicates("feature")
        c = c.assign(auc=(c[f"{focus}_vs_Control_AUC"] - 0.5).abs()).sort_values("auc", ascending=False).head(8)
        roc_panel(L, meta, focus, c, fig / f"11_roc_{focus}.png", f"ROC: {focus} против контроля")
    fc_scatter(uni, cur, fig / "12_fc_BC_vs_PC.png")
    cont = contaminant_plot(L, meta, cur, fig / "13_contaminant_series.png")
    print("Индекс ряда ПАВ по источникам:\n", cont.groupby("source").score.agg(["count", "mean"]).round(2))

    export_external(L, info, uni, meta, cur, res / "external")

    with pd.ExcelWriter(res / "compounds_BC_PC.xlsx") as w:
        cols = ["focus", "category", "name", "name_en", "formula", "ion", "rt_min", "mz", "theo_mz", "ppm", "hmdb", "kegg",
                "BC_vs_Control_log2FC", "BC_vs_Control_AUC", "BC_vs_Control_q", "BC_vs_Control_q_sexadj",
                "PC_vs_Control_log2FC", "PC_vs_Control_AUC", "PC_vs_Control_q", "PC_vs_Control_q_sexadj",
                "PC_vs_maleControl_q", "BC_vs_PC_log2FC", "BC_vs_PC_q", "comment", "feature"]
        cur[cols].to_excel(w, sheet_name="curated_candidates", index=False)
        for comp, u in sig_tabs.items():
            u.drop(columns=["comparison"]).to_excel(w, sheet_name=f"sig_{comp}"[:31])
        if mcg_settings:
            pd.DataFrame(mcg_settings).to_excel(w, sheet_name="mummichog_settings", index=False)
        for comp, pw in pws.items():
            pw.to_excel(w, sheet_name=f"pathways_{comp}"[:31], index=False)
        cont.to_excel(w, sheet_name="contaminant_index")
    print("Готово.")


if __name__ == "__main__":
    main()

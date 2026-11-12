"""Проверка 14 соединений табл. 5 статьи ABC в обработке s107 (BC против контролей).

    METAB_DATASET=s107 python report_s107/abc_table5_check.py

Для каждого соединения ищется признак (m/z ±5 ppm, RT ±0.1 мин) и берутся результаты из results_s107/tables.
Признаки, удалённые фильтрами (например, RSD QC > 30 %), пересчитываются отдельно той же цепочкой
(коррекция серий по QC → PQN с факторами основной обработки → ½ минимума → log2 → run_comparison);
воспроизводимость цепочки проверяется на признаках, прошедших фильтры (результат совпадает с основной таблицей).
Выход: report_s107/abc_table5/ — T1_ABC_markers_in_s107_checked.xlsx, рисунок ABC_Table5_boxplots.(png|pdf|tiff).
"""
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import config as C, io, preprocess, stats as st  # noqa: E402

assert C.DATASET == "s107", "запускать с METAB_DATASET=s107"
T = Path(C.OUT) / "tables"
OUT = Path(__file__).parent / "abc_table5"
OUT.mkdir(exist_ok=True)
SRC = Path("data/abc/T1_ABC_markers_in_s107.xlsx")

# ------------------------------------------------------------------ данные основной обработки
feat = pd.read_csv(T / "features.csv", index_col=0)
res = pd.read_csv(T / "stats_BC_vs_Control.csv", index_col=0)
smp = pd.read_csv(T / "samples.csv", index_col=0)
Lmain = pd.read_csv(T / "processed_log2_matrix.csv", index_col=0).T
Lmain.index = Lmain.index.astype(str)

X_all, kind, info, _ = io.load_features(C.DATA)
meta = io.load_metadata(C.METADATA, kind)
meta["inj_order"] = io.injection_order(meta)
meta["batch"] = io.batch_of(meta)
for c in ["contam_index", "contam_level", "dilution_PQN"]:
    meta[c] = smp[c].reindex(meta.index)
use = meta.index[meta.study | meta.kind.isin(["QC", "blank"])]
ana = meta.index[meta.study | (meta.kind == "QC")]
qc = meta.index[meta.kind == "QC"]
study = meta.index[meta.study]


def recompute(cols):
    """Та же цепочка, что в run_all, для произвольных признаков (без фильтров)."""
    det = preprocess.detection(X_all.loc[use], meta.loc[use])
    keep = pd.Series(True, index=X_all.columns)
    Xb = preprocess.batch_correct(X_all.loc[use], meta.loc[use], keep & X_all.columns.isin(cols), lambda s: None)
    Xn = Xb.loc[ana, cols].div(meta.loc[ana, "dilution_PQN"], axis=0)
    L = preprocess.impute_log(Xn, ana)
    cat = pd.Series("endogenous", index=cols)
    spec = C.COMPARISONS["BC_vs_Control"]
    r, _, _ = st.run_comparison(L.loc[study], meta, "BC_vs_Control", spec, cat)
    r["RSD_QC"] = preprocess.rsd(Xb.loc[qc, cols])
    for g in C.STUDY_GROUPS:
        r[f"det_{g}"] = det.loc[cols, f"det_{g}"]
    r["det_QC"] = det.loc[cols, "det_QC"]
    return r, L


def find(mz, rt, ppm=5, drt=0.1):
    d = info[((info.mz - mz).abs() / mz * 1e6 <= ppm) & ((info.rt_min - rt).abs() <= drt)]
    return list(d.index)


# ------------------------------------------------------------------ проверка воспроизводимости цепочки
chk = [f for f in res.index if f in X_all.columns][:5]
r_chk, _ = recompute(chk)
dev = (r_chk[["beta_model", "p_matched", "AUC"]] - res.loc[chk, ["beta_model", "p_matched", "AUC"]]).abs().max().max()
print(f"Проверка цепочки на {len(chk)} признаках: макс. отклонение beta/p_matched/AUC = {dev:.2e}")

# ------------------------------------------------------------------ сопоставление
src = pd.read_excel(SRC, sheet_name=None)
t1 = src["T1_ABC_markers_in_s107"].copy()
# опечатка в табл. 5: признак N-ацетилтриптофана в исходной таблице res_annotated_2 — m/z 245.0923 при 6.34 мин
LOOKUP = {8: (245.0923, 6.34)}

rows, EXTRA_L = [], {}
for _, r in t1.iterrows():
    no = int(r.iloc[0])
    mz, rt = LOOKUP.get(no, (r["m/z"], r["RT"]))
    cand = find(mz, rt)
    out = dict(no=no, feature=None, rt_s107=np.nan, mz_s107=np.nan, ppm=np.nan, kept=None, rep=None)
    if cand:
        f = min(cand, key=lambda c: abs(info.rt_min[c] - rt))
        out.update(feature=f, rt_s107=info.rt_min[f], mz_s107=info.mz[f], ppm=(info.mz[f] - mz) / mz * 1e6,
                   kept=bool(feat.kept.get(f, False)))
        if out["kept"]:
            rep = feat.representative.get(f, f)
            out["rep"] = rep
            s = res.loc[rep]
            out.update(src="основная обработка", category=s.compound_category, n_ions=s.n_ions,
                       det_Control=feat.det_Control[f], det_BC=feat.det_BC[f], RSD_QC=feat.RSD_QC_raw[f])
        else:
            s, Lx = recompute([f])
            s = s.loc[f].copy()
            # q: BH по всем соединениям основной обработки + этот признак
            allp = pd.concat([res.p_model, pd.Series({f: s.p_model})])
            s["q_model"] = st.bh(allp.to_numpy())[-1]
            s["finding"] = False                     # признак не прошёл фильтры — находкой не считается
            EXTRA_L[f] = Lx.loc[study, f]
            out.update(src="пересчёт (признак удалён фильтром)", category="не классифицирован (вне фильтров)",
                       n_ions=1, det_Control=s.det_Control, det_BC=s.det_BC, RSD_QC=s.RSD_QC)
        for k in ["beta_model", "q_model", "p_model", "log2FC_matched", "p_matched", "n_pairs_nonzero", "AUC",
                  "p_lowindex", "log2FC_lowindex", "p_ctrl_hi_vs_lo", "log2FC_ctrl_hi_vs_lo", "finding"]:
            out[k] = s.get(k, np.nan)
    rows.append(out)
R = pd.DataFrame(rows).set_index("no")
print(R.drop(columns=["feature", "rep"]).to_string())

# ------------------------------------------------------------------ соседние изомеры и изотопы
def corr(f1, f2):
    a, b = Lmain.loc[study, f1], Lmain.loc[study, f2]
    return a.corr(b, method="spearman")


extra = []
for lab, (mz, rt) in {"Leu-Val (MS-DIAL), 231.1705 / 3.01": (231.1705, 3.008),
                      "Asn-Ala-Arg (MS-DIAL), 360.2129 / 7.95": (360.2129, 7.951)}.items():
    f = find(mz, rt, drt=0.02)[0]
    rep = feat.representative.get(f, f)
    s = res.loc[rep]
    extra.append(dict(label=lab, feature=f, beta_model=s.beta_model, q_model=s.q_model, p_matched=s.p_matched,
                      AUC=s.AUC, finding=s.finding))
E = pd.DataFrame(extra)
f_vl, f_lv = R.loc[1, "rep"], E.feature[0]
f_eiv, f_nar = R.loc[2, "rep"], E.feature[1]
rho = {"Val-Leu 4.15 vs Leu-Val 3.01": corr(f_vl, feat.representative.get(f_lv, f_lv)),
       "Glu-Ile-Val 2.57 vs 360.2129 7.95": corr(f_eiv, feat.representative.get(f_nar, f_nar))}
print(E.to_string()); print(rho)
iso = info[(info.rt_min - 2.016).abs() <= 0.03]
iso = iso[iso.mz.between(269.0, 272.2)]
print("Ионы при 2.016 мин в диапазоне m/z 269–272:\n", iso.to_string())

# ------------------------------------------------------------------ итоговая таблица
RU = {1: "Val-Leu", 2: "Glu-Ile-Val", 3: "Phe-Ile-Arg", 4: "Leu-Leu-Arg", 5: "Leu-Leu-Tyr", 6: "Глутаминилпролин",
      7: "Leu-Thr-Ser", 8: "N-Ацетилтриптофан", 9: "LysoPE(20:4/0:0)", 10: "Изомер инозина",
      11: "5-Гидроксиомепразол", 12: "2-Иминобиотин", 13: "(3S)-3-[(2S)-2-Амино-3-гидроксибутил]пирролидин-2-он",
      14: "2-Амино-3-{[2-(пиридин-4-ил)этил]сульфанил}пропановая кислота"}


def reason(no, x):
    """Почему соединение не стало находкой (по фактическим значениям критериев)."""
    if no == 8:
        return (f"Признак удалён фильтром RSD QC ≤ 30 % (RSD {x.RSD_QC:.1%}); статистика рассчитана отдельно той же "
                "цепочкой; 245.0923 не является [M+H]+ C13H14N2O3 (247.1077) -> MSI 4")
    if x.finding is True or x.finding == 1:
        if pd.isna(x.p_ctrl_hi_vs_lo):
            return "Находка; не обнаружен ни в одном контроле — критерий «контроли high vs low» выполнен по правилу ctrl_not_detected"
        return ""
    why = []
    if x.category in ("contaminant", "drug", "tobacco", "food"):
        why.append(f"категория «{x.category}» (экзогенное)")
    if not x.q_model < 0.05:
        why.append(f"q = {x.q_model:.2g}")
    if not x.p_matched < 0.05:
        why.append(f"p_matched = {x.p_matched:.2g}")
    if not x.p_lowindex < 0.05:
        why.append(f"p против контролей с низким индексом = {x.p_lowindex:.2g}")
    if x.p_ctrl_hi_vs_lo < 0.05:
        why.append(f"различие между контролями с высоким и низким индексом загрязнения (p = {x.p_ctrl_hi_vs_lo:.3f}, "
                   f"log2FC = {x.log2FC_ctrl_hi_vs_lo:.2f})")
    return "Не находка: " + "; ".join(why)


EXTRA_NOTE = {
    10: "Категория: в v1 «пищевое» (перенесена из аннотации отчёта, C10H20O4S2); по данным не подтверждена -> "
        "«неоднозначное». Проверка спектров в FreeStyle (пробы 4, 56, QC17, Blank15): в холостой пробе нет; фрагмента "
        "m/z 137.0458 нет ни в пробе, ни в QC; изотопная картина противоречива (заряд пика определяется как z = 16) — "
        "ни C10H12N4O5 (ABC), ни C10H20O4S2 (отчёт) не подтверждены",
    11: "Категория: в v1 «лекарственное» (аннотация ABC); по данным не подтверждена -> «неоднозначное». m/z на 21.8 ppm "
        "от [M+H]+ C17H19N3O4S; MS-DIAL — безафибрат (low score). Проверка в FreeStyle (пробы 4, 56, QC17, Blank15): "
        "пик слабый, изотопа хлора (364.1219) на этом времени нет, MS/MS для прекурсора 362.12 нет (362.1512 при "
        "8.2–8.3 мин — другое соединение)",
}
MSI = pd.read_csv(T / "candidates_curated.csv")
MSI = pd.concat([MSI, pd.read_csv(T / "candidates_mass_mismatch.csv")])


def msi_level(no, x):
    if no == 8:
        return 4
    m = MSI[((MSI.mz - x.mz_s107).abs() < 0.002) & ((MSI.rt_min - x.rt_s107).abs() < 0.02)]
    return int(m.MSI_level.iloc[0]) if len(m) and pd.notna(m.MSI_level.iloc[0]) else 4


# исправленные/дополненные столбцы M–T исходного файла
# N-ацетилтриптофан: m/z исправлен на 245.0923 (опечатка в табл. 5), ppm пересчитан от опубликованной формулы
from pipeline.features import formula_mass, ADDUCTS
t1["m/z"] = t1["m/z"].astype(float)
t1.loc[t1.index[7], "m/z"] = 245.0923
_th = formula_mass("C13H14N2O3") + ADDUCTS["[M+H]+"]
t1.loc[t1.index[7], "Formula vs m/z, ppm"] = round((245.0923 - _th) / _th * 1e6, 1)
t1.insert(t1.columns.get_loc("Category (this analysis)") + 1, "MSI (this analysis)", [msi_level(int(r.iloc[0]), R.loc[int(r.iloc[0])]) for _, r in t1.iterrows()])
cols = {"Category (this analysis)": "category", "log2FC (adjusted)": "beta_model", "q (adjusted)": "q_model",
        "p (matched pairs)": "p_matched", "AUC": "AUC"}
chk_rows = []
for i, (_, r) in enumerate(t1.iterrows()):
    no = int(r.iloc[0])
    x = R.loc[no]
    diffs = []
    for c, k in cols.items():
        old, new = r[c], x[k]
        if isinstance(new, str) or pd.isna(old):
            changed = str(old) != str(new) and not (pd.isna(old) and pd.isna(new))
        else:
            changed = not np.isclose(float(old), float(new), rtol=1e-6)
        if changed:
            diffs.append(c)
        t1.loc[t1.index[i], c] = new
    same = bool(np.sign(x.beta_model) > 0 and x.log2FC_matched > 0 and x.p_matched < 0.05 and x.q_model < 0.05)
    t1.loc[t1.index[i], "Same direction (adjusted)"] = "yes" if same else "no"
    t1.loc[t1.index[i], "Finding (all criteria)"] = "yes" if x.finding is True or x.finding == 1 else "no"
    chk_rows.append(dict(no=no, prefilled="подтверждено" if not diffs else
                         ("заполнено" if pd.isna(r["log2FC (adjusted)"]) else "исправлено: " + ", ".join(diffs))))
    note = "; ".join(v for v in [reason(no, x), EXTRA_NOTE.get(no, "")] if v)
    if no == 8:
        r = r.copy(); r["Note"] = ("m/z corrected: 249.0923 in Table 5 is a typo for 245.0923 (source table res_annotated_2); "
                                "ppm recalculated for 245.0923")
    old_note = r["Note"] if pd.notna(r["Note"]) else ""
    t1.loc[t1.index[i], "Note"] = "; ".join(s for s in [old_note, note] if s)
chk = pd.DataFrame(chk_rows).set_index("no")

det = R[["feature", "rt_s107", "mz_s107", "ppm", "kept", "src", "category", "n_ions", "det_Control", "det_BC", "RSD_QC",
         "beta_model", "q_model", "p_model", "log2FC_matched", "p_matched", "n_pairs_nonzero", "AUC",
         "log2FC_lowindex", "p_lowindex", "log2FC_ctrl_hi_vs_lo", "p_ctrl_hi_vs_lo", "finding"]].copy()
det.insert(0, "compound", [RU[n] for n in det.index])
det.insert(1, "prefilled_values", chk.prefilled)
det = det.rename(columns={"beta_model": "log2FC_adjusted (beta_model)", "q_model": "q_adjusted (BH)",
                          "log2FC_matched": "log2FC_matched (median paired diff.)",
                          "p_matched": "p_matched (Wilcoxon signed-rank, paired)", "ppm": "Δm/z, ppm"})
E2 = E.copy()
E2["spearman_rho_with_ABC_feature"] = list(rho.values())

notes = [
    "Источник значений: обработка s107 (results_s107, коррекция серий по QC, PQN, ½ минимума, log2); BC n = 64, контроль n = 40.",
    "log2FC (adjusted) — коэффициент группы в модели log2 ~ группа + пол + возраст + индекс загрязнения; q — BH по всем соединениям.",
    "p (matched pairs) — критерий знаковых рангов Уилкоксона для связанных выборок (пары по возрасту внутри пола, 1:1, "
    "оптимальное назначение; нулевые разности исключаются; точный вариант при отсутствии совпадений). Манн–Уитни не "
    "использовался; пересчёт не требуется. Минимально возможное точное p при 24 парах = 2·2^-24 ≈ 1.2e-7.",
    "AUC — по всем пробам (64 против 40).",
    "Критерии находки: q < 0.05; p_matched < 0.05 в том же направлении; p против контролей с низким индексом < 0.05 "
    "в том же направлении; нет различия между контролями с высоким и низким индексом (p ≥ 0.05; если соединение не "
    "обнаружено ни в одном контроле, сравнение не определено и критерий считается выполненным — столбец "
    "finding_via_ctrl_not_detected); соединение не экзогенное (контаминант, лекарство, табак, пища).",
    "Leu-Val (3.01 мин) и 360.2129 (7.95 мин) — отдельные хроматографические пики, а не сдвиг RT: признаки при RT "
    "статьи (4.15 и 2.57 мин) присутствуют в той же таблице; разница RT 1.1 и 5.4 мин при воспроизводимости RT "
    "в серии < 0.02 мин. См. лист Isomers.",
    "269.0880 / 2.02 мин и 362.1248 / 8.79 мин: категории из отчёта/ABC по данным не подтверждены (проверка спектров "
    "в FreeStyle, см. Note) -> «неоднозначное».",
    "Когорта: 64 РМП (54 М / 10 Ж) и 40 контролей (16 М / 24 Ж); исключены № 32 (лейкоплакия), № 54 (опухоль "
    "мочеточника), № 105 (гнёзда Брунна); пол № 35 и № 86 исправлен на Ж по первичной таблице.",
    "N-ацетилтриптофан: статистика признака 245.0923 / 6.34 мин рассчитана вне фильтров (RSD QC 31 % > 30 %); "
    "q — BH по всем соединениям основной обработки плюс этот признак. Результат следует приводить с оговоркой "
    "о воспроизводимости в QC.",
]
with pd.ExcelWriter(OUT / "T1_ABC_markers_in_s107_checked.xlsx", engine="openpyxl") as w:
    t1.to_excel(w, sheet_name="T1_ABC_markers_in_s107", index=False)
    det.to_excel(w, sheet_name="Details_s107")
    E2.to_excel(w, sheet_name="Isomers", index=False)
    pd.DataFrame({"Notes": notes}).to_excel(w, sheet_name="Notes", index=False)
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
wb = load_workbook(OUT / "T1_ABC_markers_in_s107_checked.xlsx")
for ws in wb.worksheets:
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    for col in ws.columns:
        L_ = max(len(str(c.value or "")) for c in col)
        ws.column_dimensions[col[0].column_letter].width = min(max(10, L_ * 0.9), 60)
ws = wb["T1_ABC_markers_in_s107"]
yel, grn = PatternFill("solid", fgColor="FFF2CC"), PatternFill("solid", fgColor="E2EFDA")
for k, no in enumerate(R.index, start=2):
    fill = grn if t1.iloc[k - 2]["Finding (all criteria)"] == "yes" else yel
    for c in ws[k]:
        c.alignment = Alignment(wrap_text=True, vertical="top")
        if c.column > 11:
            c.fill = fill
    for c in ws[k][14:16]:
        c.number_format = "0.0E+00"
wb.save(OUT / "T1_ABC_markers_in_s107_checked.xlsx")

# ------------------------------------------------------------------ рисунки
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

COL = {"BC": "#2a5cb8", "Control": "#b5262b"}
DARK = {"BC": "#1d3f80", "Control": "#7e1a1e"}
BOX = {"BC": "#dbe6f5", "Control": "#f5dcda"}
LAB = {"BC": "РМП", "Control": "Контроль"}
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
                     "mathtext.it": "Liberation Sans:italic"})
raw = X_all.loc[study]
grp = meta.loc[study, "group"]
ION = {14: "[M+Na]$^+$", 8: "", 10: ""}      # № 8 и № 10: ион не установлен (MSI 4) — только m/z и RT


def panel(ax, no, letter):
    x = R.loc[no]
    f = x.feature
    y = EXTRA_L[f] if f in EXTRA_L else Lmain.loc[study, x.rep]
    imp = raw[f].isna() | (raw[f] <= 0)
    for i, g in enumerate(["BC", "Control"]):
        idx = grp.index[grp == g]
        v = y[idx].to_numpy()
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, patch_artist=True, whis=(5, 95),
                   boxprops=dict(facecolor=BOX[g], edgecolor=COL[g], lw=0.8),
                   medianprops=dict(color=DARK[g], lw=1.6), whiskerprops=dict(color=COL[g], lw=0.8),
                   capprops=dict(color=COL[g], lw=0.8))
        xj = np.random.RandomState(no * 10 + i).normal(i, 0.08, len(v))
        m = imp[idx].to_numpy()
        ax.scatter(xj[~m], v[~m], s=9, c=COL[g], alpha=0.75, lw=0.3, edgecolors="white", zorder=3)
        ax.scatter(xj[m], v[m], s=9, facecolors="none", edgecolors=COL[g], lw=0.6, alpha=0.8, zorder=3)
    ax.set_xticks([0, 1], [f"{LAB[g]}\n(n = {(grp == g).sum()})" for g in ["BC", "Control"]], fontsize=7)
    ax.set_xlim(-0.6, 1.6)
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_, hi_ + 0.12 * (hi_ - lo_))
    ax.grid(axis="y", color="#efefef", lw=0.5)
    ax.set_axisbelow(True)
    ion = ION.get(no, "[M+H]$^+$")
    ax.set_title(f"№ {no}\n$\\mathit{{m/z}}$ {x.mz_s107:.4f}{' ' + ion if ion else ''} · {x.rt_s107:.2f} мин",
                 fontsize=7, loc="left", pad=4)
    ax.text(0.98, 0.98, f"AUC {x.AUC:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.3,
            color="#555555")
    ax.text(-0.32, 1.12, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def figure(nos, name, ncol=3):
    nrow = int(np.ceil(len(nos) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(16 / 2.54, (6.2 * nrow + 0.5) / 2.54),
                            gridspec_kw={"hspace": 0.85, "wspace": 0.42}, squeeze=False)
    for k, (ax, no) in enumerate(zip(axs.ravel(), nos)):
        panel(ax, no, "абвгдежзиклмн"[k])
        if k % ncol == 0:
            ax.set_ylabel("log$_2$ интенсивность (PQN)")
    for ax in axs.ravel()[len(nos):]:
        ax.axis("off")
    for ext in ("pdf", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(600, 600))


fnd = [n for n in R.index if R.loc[n, "finding"] is True or R.loc[n, "finding"] == 1]
rest = [n for n in R.index if n not in fnd]
figure(fnd, "Рисунок_ABC_табл5_находки")
figure(rest, "Рисунок_ABC_табл5_не_находки")


def caption(nos, head):
    parts = []
    for k, no in enumerate(nos):
        x = R.loc[no]
        extra = " (признак удалён фильтром RSD QC; m/z исправлен с 249.0923, MSI 4)" if no == 8 else ""
        parts.append(f"{'абвгдежзиклмн'[k]} — № {no}, {RU[no]}{extra}")
    return (head + " Номера соответствуют табл. 5 статьи ABC: " + "; ".join(parts) + ". Показаны log2-интенсивности "
            "после коррекции серий по QC и нормализации PQN; открытые символы — значения, не обнаруженные MS-DIAL и "
            "заполненные ½ минимума признака; ящик — медиана и квартили, усы — 5–95-й процентили; РМП n = 64, "
            "контроль n = 40. В правом верхнем углу — AUC по всем пробам.")


with open(OUT / "Подписи_к_рисункам.txt", "w", encoding="utf-8") as fh:
    fh.write("Рисунок_ABC_табл5_находки\n" + caption(fnd, "Соединения табл. 5 статьи ABC, подтверждённые в выборке "
             "107 проб как находки (все критерии выполнены).") + "\n\n")
    fh.write("Рисунок_ABC_табл5_не_находки\n" + caption(rest, "Соединения табл. 5 статьи ABC, обнаруженные в выборке "
             "107 проб, но не вошедшие в находки (причины — в табл. T1).") + "\n")
print("находки:", fnd, "остальные:", rest, "->", OUT)

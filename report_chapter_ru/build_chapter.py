#!/usr/bin/env python3
"""Сборка главы отчёта (рус.) о фреймворке обработки метаболомных данных: текст, таблицы и рисунки берутся из
файлов результатов, поэтому все числа в главе совпадают с данными.

    python report_chapter_ru/build_chapter.py      # после run_all.py (s266), METAB_DATASET=s107 run_all.py и report_s107.py

Результат: report_chapter_ru/Глава_фреймворк_биомаркеры.docx (+ .md), рисунки — report_chapter_ru/figures/,
приложения — report_chapter_ru/appendix/ (таблицы xlsx, архив кода).
"""
import shutil
import subprocess
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report_chapter_ru"
FIG = OUT / "figures"
APP = OUT / "appendix"
SUB = ROOT / "submission_analytica"
R266 = ROOT / "results" / "tables"
R107 = ROOT / "results_s107" / "tables"
REP107 = ROOT / "report_s107"

MODEL_RU = {
    "Top-5 panel + logistic regression (in-fold selection)": "Панель из 5 соединений (отбор внутри фолда)",
    "Top-5 panel, features adjusted for exogenous-signal index in-fold": "Панель, поправка на индекс внутри фолда",
    "PLS-DA, 2 components (all non-exposure compounds)": "PLS-DA, 2 компоненты",
    "Sex only (baseline, same CV)": "Только пол (базовая модель)",
    "Exogenous-signal index only (baseline, same CV)": "Только экзогенный индекс (базовая модель)",
    "панель top-5 + логистическая регрессия (отбор внутри CV)": "Панель из 5 соединений (отбор внутри фолда)",
    "панель top-5, признаки скорректированы на индекс загрязнения внутри фолда": "Панель, поправка на индекс внутри фолда",
    "PLS-DA, 2 компоненты (все признаки без экспозиций)": "PLS-DA, 2 компоненты",
    "только sex (базовая модель, та же CV)": "Только пол (базовая модель)",
    "только contam_index (базовая модель, та же CV)": "Только экзогенный индекс (базовая модель)",
}
SET_RU = {"all samples": "все пробы", "age-matched pairs": "пары по возрасту и полу", "index-matched pairs": "пары по индексу",
          "all": "все пробы", "matched": "пары по возрасту и полу", "contam_matched": "пары по индексу"}
CMP_RU = {"BC vs Control": "РМП – контроль", "BC vs PC": "РМП – РПЖ", "PC vs Control": "РПЖ – контроль",
          "PC vs male Control": "РПЖ – контроль (мужчины)", "BC_vs_Control": "РМП – контроль"}


def md_table(df, floatfmt=None):
    floatfmt = floatfmt or {}
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                v = "" if np.isnan(v) else (floatfmt.get(c, "{:.3g}").format(v))
            cells.append(str(v).replace("|", "/"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def p_fmt(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


# ------------------------------------------------------------------ схема фреймворка
def workflow_figure():
    steps = [
        ("1. Входные данные", "таблица MS-DIAL (площади пиков),\nметаданные: группа, пол, возраст,\nпорядок инъекций"),
        ("2. Предобработка", "фильтры детекции и RSD в QC;\nпроверка дрейфа (leave-one-QC-out);\nPQN, импутация ½ min, log₂"),
        ("3. Соединения", "группировка ионов (ΔRT ≤ 0.02 мин,\nρ ≥ 0.8): аддукты, изотопы,\nфрагменты → одно соединение"),
        ("4. Происхождение сигнала", "категории: лекарства, пища, табак\n(курация + ключевые слова);\nгомологический ряд → экзогенный индекс"),
        ("5. Одномерный анализ", "модель log₂ ~ группа + пол + возраст +\nиндекс, FDR; пары (Уилкоксон);\nконтроли с низким индексом"),
        ("6. Правило «находки»", "все проверки согласованы\nпо знаку; экзогенные\nсоединения исключены"),
        ("7. Классификация", "повторная CV 20×5, отбор внутри фолда,\n1000 перестановок, базовые модели,\nподобранные подмножества"),
        ("8. Аннотация и проверки", "MSI 2–4, курация формул и ионов;\nанализы чувствительности;\nсверка с прежними маркерами"),
    ]
    fig, ax = plt.subplots(figsize=(17 / 2.54, 11.5 / 2.54))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 4)
    ax.axis("off")
    pos = [(0, 3), (1, 3), (0, 2), (1, 2), (0, 1), (1, 1), (0, 0), (1, 0)]
    for (t, d), (x, y) in zip(steps, pos):
        ax.add_patch(FancyBboxPatch((x + 0.04, y + 0.08), 0.92, 0.84, boxstyle="round,pad=0.01,rounding_size=0.04",
                                    fc="#eef3fb" if y % 2 else "#fbf1ec", ec="#555", lw=0.6))
        ax.text(x + 0.5, y + 0.78, t, ha="center", va="center", fontsize=7.5, fontweight="bold")
        ax.text(x + 0.5, y + 0.4, d, ha="center", va="center", fontsize=6.2, linespacing=1.3)
    arrows = [((0.96, 3.5), (1.04, 3.5)), ((1.5, 3.08), (0.5, 2.92)), ((0.96, 2.5), (1.04, 2.5)), ((1.5, 2.08), (0.5, 1.92)),
              ((0.96, 1.5), (1.04, 1.5)), ((1.5, 1.08), (0.5, 0.92)), ((0.96, 0.5), (1.04, 0.5))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.8))
    fig.savefig(FIG / "Рис_1_схема_фреймворка.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ таблицы
def tables():
    T = {}
    d = pd.read_csv(SUB / "2_main_tables" / "Table1_demographics.csv")
    d["Group"] = d.Group.map({"Control": "Контроль", "BC": "РМП", "PC": "РПЖ"})
    d.columns = ["Группа", "n", "Пол, М/Ж (н/д)", "Возраст, среднее ± SD", "Возраст, медиана (диапазон)", "Возраст известен, n"]
    T["demo266"] = d
    t = pd.read_csv(SUB / "2_main_tables" / "Table1_tests.csv")
    T["demo266_tests"] = t
    t2 = pd.read_csv(SUB / "2_main_tables" / "Table2_compact.csv")
    t2.columns = ["№", "Соединение (уровень MSI)", "m/z (ион)", "RT, мин", "log₂FC", "q", "AUC", "p в парах"]
    T["t2"] = t2
    s2 = pd.read_excel(SUB / "3_supplementary" / "Supplementary_Tables.xlsx", sheet_name="S2_classification").dropna(subset=["Model"])
    s2 = s2[s2.Comparison.isin(CMP_RU)]
    T["s2"] = pd.DataFrame({"Сравнение": s2.Comparison.map(CMP_RU), "Набор проб": s2["Sample set"].map(SET_RU),
                            "Модель": s2.Model.map(MODEL_RU).fillna(s2.Model), "n (случаи/референс)": [f"{int(a)}/{int(b)}" for a, b in zip(s2["n case"], s2["n reference"])],
                            "AUC": [f"{a:.2f} ± {b:.2f}" for a, b in zip(s2["AUC, mean"], s2["AUC, SD"])],
                            "p (перестановки)": s2["Permutation p"].map(p_fmt)})
    s6 = pd.read_excel(SUB / "3_supplementary" / "Supplementary_Tables.xlsx", sheet_name="S6_sensitivity").dropna(subset=["Analysis"])
    s6 = s6[s6.Comparison.isin(["BC vs Control", "BC vs PC"])]
    an = {"Main (PQN": "Основной анализ (PQN, индекс из 10 членов ряда)", "MSTUS normalization": "Нормализация MSTUS",
          "Without samples with PQN factor": "Без проб с PQN-фактором < 0.2",
          "Exogenous-signal index from 5 members": "Индекс из 5 членов ряда, отсутствующих в бланках",
          "Only the 1,366 compounds": "Только 1366 соединений без импутированных значений",
          "Imputed values treated as missing": "Импутированные значения как пропуски"}
    def a_ru(x):
        for k, v in an.items():
            if str(x).startswith(k):
                return v
        return x
    T["s6"] = pd.DataFrame({"Сравнение": s6.Comparison.map(CMP_RU), "Анализ": s6.Analysis.map(a_ru), "q < 0.05": s6["q < 0.05"].astype(int),
                            "Находки": s6.Findings.astype(int), "Общие находки с основным": s6["Shared findings"].astype(int),
                            "AUC панели": s6["Panel AUC (CV)"].map(lambda v: "" if pd.isna(v) else f"{v:.2f}")})
    s2b = pd.read_csv(SUB / "3_supplementary" / "TableS2b_panel_composition.csv")
    s2b = s2b[(s2b.Comparison == "BC vs Control") & (s2b.Model == "Top-5 panel")]
    T["s2b"] = pd.DataFrame({"m/z": s2b["m/z"].map(lambda v: f"{v:.4f}"), "RT, мин": s2b["RT (min)"],
                             "Отобрано в фолдах (из 100)": s2b["Selected in training folds (of 100)"],
                             "Ранг в панели на всех пробах": s2b["All-sample panel (rank)"].fillna("").map(lambda v: "" if v == "" else str(int(v))),
                             "Одномерная находка": s2b["Univariate finding"].map({"yes": "да", "no": "нет"})})
    s7 = pd.read_csv(SUB / "3_supplementary" / "TableS7_blank_check.csv")
    T["s7"] = pd.DataFrame({"Член ряда": s7.iloc[:, 0], "m/z [M+H]⁺": s7.iloc[:, 1].map(lambda v: f"{v:.4f}"), "RT, мин": s7.iloc[:, 2].round(2),
                            "В индексе": s7.iloc[:, 3].map({"yes": "да", "no": "нет"}),
                            "Бланк / QC": s7["Solvent blank / QC"].round(3)})
    # s107
    d7 = pd.read_csv(R107 / "demographics.csv")
    T["demo107"] = pd.DataFrame({"Группа": d7.group.map({"Control": "Контроль", "BC": "РМП"}), "n": d7.n,
                                 "Пол, М/Ж": [f"{a}/{b}" for a, b in zip(d7.male, d7.female)],
                                 "Возраст, среднее ± SD": [f"{a:.1f} ± {b:.1f}" for a, b in zip(d7.age_mean, d7.age_sd)],
                                 "Возраст, медиана (диапазон)": [f"{a:.0f} ({b:.0f}–{c:.0f})" for a, b, c in zip(d7.age_median, d7.age_min, d7.age_max)]})
    m7 = pd.read_csv(R107 / "multivariate_CV.csv")
    m7 = m7[m7.model.astype(str).isin(MODEL_RU)]
    T["mv107"] = pd.DataFrame({"Набор проб": m7.set.map(SET_RU), "Модель": m7.model.map(MODEL_RU),
                               "n (РМП/контроль)": [f"{int(a)}/{int(b)}" for a, b in zip(m7.n_case, m7.n_ref)],
                               "AUC": [f"{a:.3f} ± {b:.3f}" for a, b in zip(m7.AUC_mean, m7.AUC_sd)], "p (перестановки)": m7.perm_p.map(p_fmt)})
    x = pd.read_excel(REP107 / "Tables_s107.xlsx", sheet_name=None)
    t2r = x["T2_s266_Table2_in_s107"]
    T["rep"] = pd.DataFrame({"Соединение": t2r.Compound, "m/z": t2r["m/z"].map(lambda v: f"{v:.4f}"),
                             "log₂FC (s266)": t2r["log2FC s266"].round(2), "q (s266)": t2r["q s266"].map(p_fmt),
                             "log₂FC (s107)": t2r["log2FC s107"].round(2), "q (s107)": t2r["q s107"].map(lambda v: "" if pd.isna(v) else p_fmt(v)),
                             "AUC (s107)": t2r["AUC s107"].round(2), "p в парах (s107)": t2r["p matched s107"].map(lambda v: "" if pd.isna(v) else p_fmt(v)),
                             "Находка в s107": t2r["Finding s107"].map(lambda v: "" if pd.isna(v) else ("да" if v else "нет"))})
    t1 = x["T1_report_markers_in_s107"]
    cat = {"food": "пища", "drug": "лекарство", "tobacco": "табак", "endogenous": "эндогенное", "putative (low score)": "не установлено"}
    T["abc"] = pd.DataFrame({"№": t1["No. (report Table 20)"], "Соединение (отчёт)": t1["Compound (report)"],
                             "Заявлено": t1["Reported direction"].map({"up": "↑", "down": "↓"}), "m/z": t1["m/z"].map(lambda v: f"{v:.4f}"),
                             "Формула vs m/z, ppm": t1["Formula vs m/z, ppm"].map(lambda v: "—" if pd.isna(v) else f"{v:.1f}"),
                             "Категория": t1["Category (this analysis)"].map(lambda v: cat.get(v, v)),
                             "log₂FC": t1["log2FC (adjusted)"].round(2), "q": t1["q (adjusted)"].map(p_fmt), "AUC": t1["AUC"].round(2),
                             "Находка": t1["Finding (all criteria)"].map({True: "да", False: "нет"})})
    pc = pd.read_csv(ROOT / "results_s107_extra" / "PC_markers_old_report_vs_s266.csv")
    T["pc"] = pd.DataFrame({"№": pc.No, "RT": pc.rt, "m/z": pc.mz.map(lambda v: f"{v:.3f}"), "Аннотация (отчёт)": pc.name.fillna(""),
                            "Член экзогенного ряда": pc.series_label.fillna("").replace("", "—"),
                            "ρ с индексом": pc.rho_index.map(lambda v: "" if pd.isna(v) else f"{v:.2f}"),
                            "log₂FC (без поправок)": pc.log2FC_unadj.round(2), "q (без поправок)": pc.q_MW.map(lambda v: f"{v:.3f}"),
                            "q (модель с поправками)": pc.q_adj.map(lambda v: f"{v:.2f}"),
                            "q (мужчины-контроли)": pc.q_adj_male.map(lambda v: f"{v:.2f}")})
    return T


def copy_figures():
    m = {SUB / "1_main_figures" / "Figure1_data_quality.png": "Рис_2_качество_данных.png",
         SUB / "3_supplementary" / "FigureS6_exogenous_signal_series.png": "Рис_3_экзогенный_ряд.png",
         SUB / "3_supplementary" / "FigureS5_solvent_blanks.png": "Рис_4_бланки.png",
         SUB / "3_supplementary" / "FigureS2_age_and_matching.png": "Рис_5_возраст_и_пары.png",
         SUB / "1_main_figures" / "Figure2_volcano.png": "Рис_6_вулкан_s266.png",
         SUB / "3_supplementary" / "FigureS1_heatmap_BC_vs_Control.png": "Рис_7_тепловая_карта.png",
         SUB / "1_main_figures" / "Figure3_annotated_findings_BC.png": "Рис_8_аннотированные_находки.png",
         SUB / "1_main_figures" / "Figure4_ROC_cross_validated.png": "Рис_9_ROC.png",
         SUB / "3_supplementary" / "FigureS3_normalization.png": "Рис_10_нормализация.png",
         SUB / "3_supplementary" / "FigureS4_pathways.png": "Рис_11_пути.png",
         REP107 / "FigureS107_1_PCA.png": "Рис_12_s107_PCA.png",
         REP107 / "FigureS107_2_volcano.png": "Рис_13_s107_вулкан.png",
         REP107 / "FigureS107_4_peptides.png": "Рис_14_s107_пептиды.png",
         REP107 / "FigureS107_3_s266_findings_in_s107.png": "Рис_15_s107_воспроизведение.png"}
    for a, b in m.items():
        shutil.copy(a, FIG / b)


def appendix_files():
    shutil.copy(SUB / "3_supplementary" / "Supplementary_Tables.xlsx", APP / "Прил_В1_таблицы_s266.xlsx")
    shutil.copy(SUB / "2_main_tables" / "Tables_main.xlsx", APP / "Прил_В2_основные_таблицы_s266.xlsx")
    shutil.copy(REP107 / "Tables_s107.xlsx", APP / "Прил_В3_таблицы_s107_и_проверка_маркеров.xlsx")
    files = ["run_all.py", "make_figures.py", "report_s107.py", "blank_check.py", "sensitivity_index.py",
             "sensitivity_imputation.py", "requirements.txt", "README.md", "report_chapter_ru/build_chapter.py"]
    dirs = ["pipeline", "curation", "curation_s107", "metadata"]
    with zipfile.ZipFile(APP / "Прил_Г_код_фреймворка.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in files:
            if (ROOT / f).exists():
                z.write(ROOT / f, f)
        for d in dirs:
            for p in sorted((ROOT / d).rglob("*")):
                if p.is_file() and "__pycache__" not in p.parts:
                    z.write(p, p.relative_to(ROOT))


def format_docx(path):
    """Оформление: Times New Roman 12 (таблицы 8), поля 2 см, таблицы с сеткой и шириной столбцов по содержимому."""
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt
    doc = Document(path)
    for sec in doc.sections:
        sec.left_margin = sec.right_margin = Cm(2)
        sec.top_margin = sec.bottom_margin = Cm(2)
    for st in doc.styles:
        try:
            if st.type == 1:
                st.font.name = "Times New Roman"
                st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman") if st.element.rPr is not None and st.element.rPr.rFonts is not None else None
        except Exception:
            pass
    doc.styles["Normal"].font.size = Pt(12)
    width = 17.0
    for t in doc.tables:
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        tblPr = t._tbl.tblPr
        borders = OxmlElement("w:tblBorders")
        for e in ("top", "left", "bottom", "right", "insideH", "insideV"):
            b = OxmlElement(f"w:{e}")
            b.set(qn("w:val"), "single"); b.set(qn("w:sz"), "4"); b.set(qn("w:color"), "808080")
            borders.append(b)
        tblPr.append(borders)
        ncol = len(t.columns)
        lens = []
        for j in range(ncol):
            L = [max(max((len(w) for w in c.text.split()), default=1), min(len(c.text), 28) * 0.7) for c in t.columns[j].cells]
            lens.append(max(L) if L else 1)
        tot = sum(lens)
        ws = [max(1.35, width * l / tot) for l in lens]
        k = width / sum(ws)
        t.autofit = False
        for el in tblPr.findall(qn("w:tblW")):
            tblPr.remove(el)
        tw = OxmlElement("w:tblW"); tw.set(qn("w:w"), str(int(width * 567))); tw.set(qn("w:type"), "dxa"); tblPr.append(tw)
        lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblPr.append(lay)
        grid = t._tbl.tblGrid
        for j, gc in enumerate(grid.findall(qn("w:gridCol"))):
            gc.set(qn("w:w"), str(int(ws[j] * k * 567)))
        for j, col in enumerate(t.columns):
            for c in col.cells:
                c.width = Cm(ws[j] * k)
                for par in c.paragraphs:
                    par.paragraph_format.space_after = Pt(0)
                    for run in par.runs:
                        run.font.size = Pt(8)
                        run.font.name = "Times New Roman"
    doc.save(path)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    APP.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Liberation Sans"})
    workflow_figure()
    copy_figures()
    appendix_files()
    T = tables()
    tpl = (OUT / "chapter_template.md").read_text(encoding="utf-8")
    fmt = {k: md_table(v) for k, v in T.items()}
    md = tpl.format(**fmt)
    (OUT / "Глава_фреймворк_биомаркеры.md").write_text(md, encoding="utf-8")
    subprocess.run(["pandoc", "Глава_фреймворк_биомаркеры.md", "-o", "Глава_фреймворк_биомаркеры.docx",
                    "--resource-path", ".", "--toc", "--toc-depth=2", "-M", "toc-title=Содержание"], cwd=OUT, check=True)
    format_docx(OUT / "Глава_фреймворк_биомаркеры.docx")
    print("готово:", OUT / "Глава_фреймворк_биомаркеры.docx")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Распределение 30 маркеров первого этапа (таблица 20 отчёта) в группах РМП и контроль: 84 пробы этапа 1
(№ 1–84 таблицы 107 проб: 53 РМП, 31 контроль), исходные площади пиков из экспорта MS-DIAL.
    python figures_markers/plot_markers.py   -> figures_markers/out/Рисунок_markers_1..3 (PNG 600 dpi, TIFF, PDF, SVG)
Признак каждого маркера — ближайший по RT к заявленному (results_s107/tables/replication_markers.csv)."""
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).parent / "out"
RU = ["Непетазид", "5-Метилмиосмин", "Ацетилцистеин", "4-Гидрокси-3-(2-гидроксиэтил)ацетофенон", "Изоферуловая кислота",
      "Кмериол", "N-Мевалонилкреатинин", "N-Ацетил-S-(N-аллилтиокарбамоил)цистеин", "Треонил-фенилаланин",
      "3,6-Бис(2-(метилтио)этил)-1,4-диоксан-2,5-диол", "D-Гликозид циклогексанона", "Венлафаксин",
      "Аспарагинил-фенилаланин", "3-Гидроксиоктаноилкарнитин", "4-Метилциннамилкарнитин",
      "(2S)-2-((2-Морфолиноацетил)амино)-4-фенилбутановая кислота", "N,N-Диметилсфингозин",
      "Глутамил-изолейцил-валин", "D-Гликозид ванилинового гликоля", "Фраксин",
      "D-Глюкуронид гликоля перилловой кислоты", "Самурациклиб", "D-Глюкуронид 3-(гидроксиметил)антипирина",
      "Пептид EFTPAVAH", "Грейанотоксин I", "9-β-D-Глюкуронид 12Z-9,10-дигидроксиоктадеценовой кислоты",
      "D-Глюкуронид γ-CEHC", "12-D-Глюкуронид 9Z-12,13-дигидроксиоктадеценовой кислоты",
      "Тауринил-(9Z)-3-гидроксидодеценоилкарнитин", "N-Ацетилцистеинил-2,6-диизопропилфенол"]
COL = {"BC": "#2a5cb8", "Control": "#b5262b"}
DARK = {"BC": "#1d3f80", "Control": "#7e1a1e"}
BOX = {"BC": "#dbe6f5", "Control": "#f5dcda"}
# уточнённые аннотации (проверка формулы и иона, альтернативные аннотации того же признака)
FIX = {2: "5-Метилмиосмин (не подтверждено: Δm 63 ppm)",
       7: "Gln-Pro (в отчёте — N-мевалонилкреатинин)",
       24: "Пептид (в отчёте — EFTPAVAH; формула не соответствует m/z)",
       25: "Не установлено (в отчёте — грейанотоксин I; формула не соответствует m/z)"}
LAB = {"BC": "РМП", "Control": "Контроль"}
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic"})


def load():
    t = pd.read_csv(ROOT / "data/s107/s107_table.csv", low_memory=False).set_index("Sample")
    t.index = t.index.astype(str)
    smp = [s for s in t.index if s.isdigit() and int(s) <= 84]
    rep = pd.read_csv(ROOT / "results_s107/tables/replication_markers.csv")
    rep = rep.assign(ad=rep.delta_rt.abs()).sort_values("ad").drop_duplicates("marker")
    t20 = pd.read_csv(ROOT / "data/s107/table_nastya.csv")
    rows = []
    for _, r in t20.iterrows():
        x = rep[rep.marker == r.Compound].iloc[0]
        mz, ion = r.Characteristic_ion_mz.split(" ", 1)
        rows.append(dict(no=int(r.No), name=RU[int(r.No) - 1], feature=x.feature, mz=float(mz), ion=ion,
                         rt=float(r.Retention_time_min)))
    return t.loc[smp], rows


def ion_tex(ion):
    return ion.replace("[M+H]+", "[M+H]$^+$").replace("[M+NH4]+", "[M+NH$_4$]$^+$").replace("[M+H]2+", "[M+2H]$^{2+}$")


def panel(ax, t, m):
    """Логарифмическая шкала; не обнаружено (площадь 0) — открытые символы на уровне «н.о.» под минимумом."""
    y = t[m["feature"]].astype(float).fillna(0)
    grp = t["Group"]
    pos = y[y > 0]
    lo, hi = np.log10(pos.min()), np.log10(pos.max())
    nd = lo - 0.45                                         # уровень «не обнаружено»
    for i, g in enumerate(["BC", "Control"]):
        v = y[grp == g].to_numpy()
        lv = np.where(v > 0, np.log10(np.where(v > 0, v, 1)), nd)
        ax.boxplot(lv, positions=[i], widths=0.5, showfliers=False, patch_artist=True, whis=(5, 95),
                   boxprops=dict(facecolor=BOX[g], edgecolor=COL[g], lw=0.8),
                   medianprops=dict(color=DARK[g], lw=1.6), whiskerprops=dict(color=COL[g], lw=0.8),
                   capprops=dict(color=COL[g], lw=0.8))
        xj = np.random.RandomState(m["no"] * 10 + i).normal(i, 0.08, len(v))
        det = v > 0
        ax.scatter(xj[det], lv[det], s=9, c=COL[g], alpha=0.75, lw=0.3, edgecolors="white", zorder=3)
        ax.scatter(xj[~det], lv[~det], s=9, facecolors="none", edgecolors=COL[g], lw=0.6, alpha=0.8, zorder=3)
    ax.set_xticks([0, 1], [f"{LAB[g]} (n = {(grp == g).sum()})" for g in ["BC", "Control"]], fontsize=7)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(nd - 0.35, hi + 0.35)
    ticks = list(range(int(np.ceil(lo)), int(np.floor(hi)) + 1))
    ax.set_yticks([nd] + ticks, ["н.о."] + [f"$10^{{{k}}}$" for k in ticks], fontsize=6.8)
    ax.axhline(nd + 0.22, color="#cfcfcf", lw=0.5, ls=":")
    ax.grid(axis="y", color="#efefef", lw=0.5)
    ax.set_axisbelow(True)
    title = FIX.get(m["no"], m["name"])
    name = "\n".join(textwrap.wrap(title, 52)[:2])
    ax.set_title(f"{m['no']}. {name}\n$\\mathit{{m/z}}$ {m['mz']:.4f} {ion_tex(m['ion'])} · {m['rt']:.1f} мин",
                 fontsize=7, loc="left", pad=4)


def main():
    OUT.mkdir(exist_ok=True)
    t, rows = load()
    for k in range(3):
        fig, axs = plt.subplots(5, 2, figsize=(16 / 2.54, 23 / 2.54), gridspec_kw={"hspace": 0.95, "wspace": 0.32})
        for ax, m in zip(axs.ravel(), rows[k * 10:(k + 1) * 10]):
            panel(ax, t, m)
        for ax in axs[:, 0]:
            ax.set_ylabel("Площадь пика, отн. ед.")
        name = f"Рисунок_маркеры_{k * 10 + 1}-{k * 10 + 10}"
        for ext in ("pdf", "svg"):
            fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
        fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
        plt.close(fig)
        Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(600, 600))
    print("готово:", OUT)


if __name__ == "__main__":
    main()

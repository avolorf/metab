#!/usr/bin/env python3
"""Диаграммы Парето, нормальные графики стандартизованных эффектов и графики взаимодействий (ДОЭ пробоподготовки).
    python figures_doe/plot_doe.py   -> figures_doe/out/*.png (600 dpi), *.tif, *.pdf, *.svg"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from PIL import Image
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent))
from doe_data import F, PB, FF, VLEV, V_ORG, V_EVAP, ORG_EVAP, V

BAR, SIG, NS, CRIT, INK = "#6f9fd8", "#b5262b", "#2a5cb8", "#c4402b", "#222222"
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8.5, "axes.linewidth": 0.7,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
OUT = Path(__file__).parent / "out"


def save(fig, name):
    for ext in ("pdf", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(600, 600))


def pareto(ax, items, crit, resp):
    items = sorted(items, key=lambda x: abs(x[1]))
    y = np.arange(len(items))
    ax.barh(y, [abs(e) for _, e in items], color=[SIG if abs(e) > crit else BAR for _, e in items],
            edgecolor="#3d5f8f", lw=0.5, height=0.75)
    ax.set_yticks(y, [n for n, _ in items])
    ax.axvline(crit, color=CRIT, ls="--", lw=1)
    ax.text(crit, len(items) - 0.35, f"{crit:.3f}", color=CRIT, ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax.set_xlabel("Стандартизованный эффект (|t|)")
    ax.set_ylabel("Фактор")
    ax.set_xlim(0, max(max(abs(e) for _, e in items), crit) * 1.1)
    ax.set_ylim(-0.6, len(items) - 0.2)
    ax.set_title(f"Диаграмма Парето эффектов\n(отклик: {resp}; α = 0.05)", fontsize=9)
    ax.grid(axis="x", color="#e3e3e3", lw=0.5)
    ax.set_axisbelow(True)


def normal_plot(ax, eff, crit, resp):
    names = list(eff)
    e = np.array([eff[k] for k in names])
    o = np.argsort(e)
    p = (np.arange(1, len(e) + 1) - 0.5) / len(e)
    xs = e[o]
    sig = np.abs(xs) > crit
    zq = norm.ppf(p)
    lim = max(6.5, np.abs(xs).max() + 0.5)
    xx = np.linspace(-lim, lim, 200)
    ax.plot(xx, norm.cdf(xx) * 100, color="#d08a7e", lw=1)            # нулевая линия: эффекты ~ N(0, 1)
    ax.scatter(xs[~sig], p[~sig] * 100, s=22, c=NS, zorder=3)
    ax.scatter(xs[sig], p[sig] * 100, s=26, c=SIG, marker="s", zorder=3)
    for i, k in enumerate(np.array(names)[o]):
        if sig[i] or k == "V":
            ax.annotate(F[k], (xs[i], p[i] * 100), xytext=(6, 2), textcoords="offset points", fontsize=8)
    ax.set_yscale("function", functions=(lambda v: norm.ppf(np.clip(v, 0.01, 99.99) / 100),
                                         lambda v: norm.cdf(v) * 100))
    ticks = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99]
    ax.set_yticks(ticks, [str(t) for t in ticks])
    ax.set_ylim(1, 99)
    ax.set_xlim(-lim, 3.5)
    ax.set_xlabel("Стандартизованный эффект")
    ax.set_ylabel("Накопленная вероятность, %")
    ax.set_title(f"График нормального распределения стандартизованных эффектов\n(отклик: {resp}; α = 0.05)", fontsize=9)
    ax.grid(color="#ececec", lw=0.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=NS, label="не значимый"),
                       Line2D([], [], marker="s", ls="", color=SIG, label="значимый")],
              title="Тип эффекта", loc="lower right", fontsize=7.5, title_fontsize=7.5, frameon=True)


def letter(ax, t):
    ax.text(-0.02, 1.02, f"({t})", transform=ax.transAxes, ha="right", va="bottom", fontsize=11, fontweight="bold")


def interactions():
    C = ["#2a5cb8", "#b5262b", "#2e9a4f", "#8e44ad"]
    M = ["o", "s", "D", "^"]
    fig, axs = plt.subplots(3, 3, figsize=(17 / 2.54, 15 / 2.54), sharey=True, gridspec_kw={"hspace": 0.55, "wspace": 0.35})
    xV = np.arange(4)
    def line(ax, x, ys, col, m, lab, ls):
        ax.plot(x, ys, color=col, marker=m, ms=4.5, lw=1.2, ls=ls, label=lab)
    ax = axs[0, 1]   # % орг. р-ля, линии — V
    for k, v in enumerate(VLEV):
        line(ax, [0, 1], [V_ORG[50][k], V_ORG[100][k]], C[k], M[k], f"{v}", "-" if k == 0 else "--")
    ax.set_xticks([0, 1], ["50", "100"]); ax.set_title("% орг. р-ля × " + V, fontsize=7.8)
    ax = axs[0, 2]   # упаривание, линии — V
    for k, v in enumerate(VLEV):
        line(ax, [0, 1], [V_EVAP["+"][k], V_EVAP["−"][k]], C[k], M[k], f"{v}", "-" if k == 0 else "--")
    ax.set_xticks([0, 1], ["+", "−"]); ax.set_title("упаривание × " + V, fontsize=7.8)
    ax = axs[1, 0]   # V, линии — % орг.
    for k, o in enumerate([50, 100]):
        line(ax, xV, V_ORG[o], C[k], M[k], f"{o}", "-" if k == 0 else "--")
    ax.set_xticks(xV, [str(v) for v in VLEV]); ax.set_title(V + " × % орг. р-ля", fontsize=7.8)
    ax = axs[1, 2]   # упаривание, линии — % орг.
    for k, o in enumerate([50, 100]):
        line(ax, [0, 1], [ORG_EVAP["+"][k], ORG_EVAP["−"][k]], C[k], M[k], f"{o}", "-" if k == 0 else "--")
    ax.set_xticks([0, 1], ["+", "−"]); ax.set_title("упаривание × % орг. р-ля", fontsize=7.8)
    ax = axs[2, 0]   # V, линии — упаривание
    for k, s in enumerate(["+", "−"]):
        line(ax, xV, V_EVAP[s], C[k], M[k], s, "-" if k == 0 else "--")
    ax.set_xticks(xV, [str(v) for v in VLEV]); ax.set_title(V + " × упаривание", fontsize=7.8)
    ax = axs[2, 1]   # % орг., линии — упаривание
    for k, s in enumerate(["+", "−"]):
        line(ax, [0, 1], ORG_EVAP[s], C[k], M[k], s, "-" if k == 0 else "--")
    ax.set_xticks([0, 1], ["50", "100"]); ax.set_title("% орг. р-ля × упаривание", fontsize=7.8)
    for i, (src, ttl) in enumerate([(axs[0, 1], V), (axs[1, 0], "% орг. р-ля"), (axs[2, 0], "упаривание")]):
        h, l = src.get_legend_handles_labels()
        axs[i, i].axis("off")
        axs[i, i].legend(h, l, title="Линии ряда:\n" + ttl, loc="center", fontsize=7.5, title_fontsize=7.8, frameon=False)
    for ax in (axs[0, 1], axs[1, 2]):
        ax.tick_params(labelleft=True)
        ax.set_ylabel("Число пиков")
    for ax in axs.ravel():
        ax.set_ylim(2000, 9200)
        ax.grid(axis="y", color="#ececec", lw=0.5)
        ax.margins(x=0.15)
    axs[2, 0].set_xlabel(V)
    axs[2, 1].set_xlabel("% орг. р-ля")
    axs[1, 2].set_xlabel("упаривание")
    for i in range(3):
        axs[i, 0].set_ylabel("Число пиков")
    fig.suptitle("Графики взаимодействия эффектов\n(отклик: число пиков; средние значения)", fontsize=10, y=1.0)
    save(fig, "Рисунок_5")


def main():
    OUT.mkdir(exist_ok=True)
    for name, (resp, crit, eff) in PB.items():
        fig, axs = plt.subplots(2, 1, figsize=(14 / 2.54, 17 / 2.54), gridspec_kw={"hspace": 0.55})
        pareto(axs[0], [(F[k], v) for k, v in eff.items()], crit, resp)
        normal_plot(axs[1], eff, crit, resp)
        letter(axs[0], "а"); letter(axs[1], "б")
        save(fig, name)
    resp, crit, items = FF
    fig, ax = plt.subplots(figsize=(14 / 2.54, 8.5 / 2.54))
    pareto(ax, items, crit, resp)
    save(fig, "Рисунок_4")
    interactions()
    print("готово:", OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Зеркальные спектры MS/MS (эксперимент — вверх, библиотека NIST — вниз), подписи на русском.
    python figures_ms2/plot_mirror.py      -> figures_ms2/out/*.png (600 dpi), *.tif, *.pdf, *.svg
Значения пиков — figures_ms2/spectra.py."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from spectra import FIGS

EXP, LIB, INK = "#2a5cb8", "#c4402b", "#222222"
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.7,
                     "mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic", "mathtext.rm": "Liberation Sans",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "xtick.major.width": 0.6, "ytick.major.width": 0.6})
OUT = Path(__file__).parent / "out"


def labels(ax, peaks, sign, xr):
    """Подписи m/z над (под) пиками; при наложении соседние подписи сдвигаются по вертикали."""
    placed = []
    for mz, h, lab in sorted(peaks, key=lambda p: -p[1]):
        if not lab:
            continue
        y = sign * (h + 3)
        for (x0, y0) in placed:
            if abs(x0 - mz) < (xr[1] - xr[0]) * 0.055 and abs(y0 - y) < 9:
                y = y0 + sign * 9
        placed.append((mz, y))
        ax.text(mz, y, f"{mz:.4f}", ha="center", va="bottom" if sign > 0 else "top", fontsize=6.3, color=INK)


def panel(ax, d, letter):
    mzmax = max(m for m, _, _ in d["exp"] + d["lib"])
    xr = (d["xlim"][0], max(d["xlim"][1], mzmax + 0.06 * (d["xlim"][1] - d["xlim"][0])))
    for mz, h, _ in d["exp"]:
        ax.vlines(mz, 0, h, color=EXP, lw=1.0)
    for mz, h, _ in d["lib"]:
        ax.vlines(mz, 0, -h, color=LIB, lw=1.0)
    labels(ax, d["exp"], 1, xr)
    labels(ax, d["lib"], -1, xr)
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xlim(*xr)
    ax.set_ylim(-118, 118)
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(25))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.075, 1.04, f"({letter})", transform=ax.transAxes, ha="left", va="bottom", fontsize=10, fontweight="bold")
    ax.text(0.015, 0.95, "Эксперимент", transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color=EXP)
    ax.text(0.015, 0.05, "Библиотека NIST", transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, color=LIB)


def main():
    OUT.mkdir(exist_ok=True)
    for name, panels in FIGS.items():
        fig, axs = plt.subplots(2, 1, figsize=(16 / 2.54, 15 / 2.54), gridspec_kw={"hspace": 0.28})
        for ax, (letter, d) in zip(axs, panels.items()):
            panel(ax, d, letter)
            ax.set_ylabel("")
        axs[1].set_xlabel("$\\mathit{m/z}$", fontsize=9)
        fig.supylabel("Относительная интенсивность, %", fontsize=9, x=0.035)
        for ext in ("pdf", "svg"):
            fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
        fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
        plt.close(fig)
        Image.open(OUT / f"{name}.png").convert("RGB").save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(600, 600))
    print("готово:", OUT)


if __name__ == "__main__":
    main()

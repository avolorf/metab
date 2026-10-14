#!/usr/bin/env python3
"""Проверка бланков на алкилдиэтаноламины (ряд индекса загрязнения) по сырым данным.

    python blank_check.py                         # mzML и список инъекций из data/raw_check/
    python blank_check.py --dir путь/к/mzml --ppm 5

Вход (см. docs/BLANK_CHECK_GUIDE.md):
  data/raw_check/*.mzML           — MS1, центроиды; можно урезанные по m/z 250–450 и RT 9–19 мин
  data/raw_check/injections.csv   — file, type, batch, order[, sample]
      type: solvent_blank | process_blank | container_blank | QC | sample
Мишени: curation/contaminant_series.csv (m/z — по формуле и иону) + эндогенные контрольные ионы ниже
(в настоящем бланке их нет; если есть — перенос с предыдущей инъекции).

Для каждой мишени в каждой инъекции: XIC ± PPM, поиск апекса в окне RT ± RT_SEARCH около ожидаемого
(медианный сдвиг RT по QC учитывается), площадь в окне апекс ± RT_INT. Выход — results/blank_check/.

Правило (задано заранее): соединение «присутствует в бланке», если площадь в бланке ≥ 10 % медианы QC;
«происходит из расходников», если в process/container-бланке ≥ 10 % QC, а в solvent-бланке < 10 % QC.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.features import ADDUCTS, formula_mass

RT_SEARCH = 0.5      # мин; поиск апекса около ожидаемого RT
RT_INT = 0.10        # мин; интегрирование апекс ± RT_INT
MAX_GAP_MIN = 35     # мин; запуски ближе — соседние инъекции (градиент 23 мин + уравновешивание)
BLANK_FRAC = 0.10    # доля медианы QC, начиная с которой соединение считается присутствующим в бланке
TYPES = ["solvent_blank", "process_blank", "container_blank", "QC", "sample"]
COL = {"solvent_blank": "#6f6e69", "process_blank": "#7a5af5", "container_blank": "#c2185b", "QC": "#2a78d6",
       "sample": "#eb6834"}
# эндогенные контрольные ионы внутри окна m/z 250–450 / RT 9–19 мин (должны отсутствовать в бланках)
REFERENCE = [("ref: C19H28O3 steroid", "C19H28O3", "[M+H-H2O]+", 11.056),
             ("ref: tetrahydrocortisone", "C21H32O5", "[M+H-H2O]+", 10.482)]


def targets(series_path):
    s = pd.read_csv(series_path)
    rows = [dict(label=r.label, formula=r.formula, ion=r.ion, rt=r.rt_min, in_index=int(r.in_index), kind="series")
            for _, r in s.iterrows()]
    rows += [dict(label=l, formula=f, ion=i, rt=rt, in_index=0, kind="endogenous reference") for l, f, i, rt in REFERENCE]
    t = pd.DataFrame(rows)
    t["mz"] = [formula_mass(f) + ADDUCTS[i] for f, i in zip(t.formula, t.ion)]
    return t


def read_ms1(path):
    """Возвращает список (rt_min, mz[], int[]) для MS1-спектров mzML."""
    from pyteomics import mzml
    out = []
    with mzml.MzML(str(path)) as f:
        for sp in f:
            if sp.get("ms level", 1) != 1:
                continue
            sc = sp["scanList"]["scan"][0]
            rt = float(sc["scan start time"])
            unit = getattr(sc["scan start time"], "unit_info", "minute")
            if unit and "second" in str(unit):
                rt /= 60
            out.append((rt, np.asarray(sp["m/z array"]), np.asarray(sp["intensity array"])))
    return out


def start_time(path):
    """Время запуска инъекции (атрибут startTimeStamp элемента run в mzML; msconvert его сохраняет)."""
    import re
    with open(path, "rb") as f:
        head = f.read(200000).decode("utf-8", "ignore")
    m = re.search(r'<run[^>]*startTimeStamp="([^"]+)"', head)
    return pd.to_datetime(m.group(1), utc=True) if m else pd.NaT


def xic(spectra, mz, ppm):
    tol = mz * ppm * 1e-6
    rt = np.array([s[0] for s in spectra])
    y = np.array([s[2][(s[1] >= mz - tol) & (s[1] <= mz + tol)].sum() for s in spectra])
    return rt, y


def peak(rt, y, rt0, shift=0.0):
    m = np.abs(rt - (rt0 + shift)) <= RT_SEARCH
    if not m.any() or y[m].max() <= 0:
        return np.nan, 0.0, 0.0
    i = np.where(m)[0][np.argmax(y[m])]
    w = np.abs(rt - rt[i]) <= RT_INT
    area = float(np.trapezoid(y[w], rt[w] * 60)) if w.sum() > 1 else float(y[i])
    return float(rt[i]), float(y[i]), area


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="data/raw_check")
    ap.add_argument("--list", default=None, help="по умолчанию <dir>/injections.csv")
    ap.add_argument("--series", default="curation/contaminant_series.csv")
    ap.add_argument("--ppm", type=float, default=5.0)
    ap.add_argument("--out", default="results/blank_check")
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    inj = pd.read_csv(a.list or d / "injections.csv", dtype={"file": str, "sample": str})
    for c in ["batch", "order", "sample"]:
        if c not in inj:
            inj[c] = np.nan if c == "order" else ""
    inj["type"] = inj["type"].str.strip()
    bad = set(inj["type"]) - set(TYPES)
    if bad:
        raise SystemExit(f"Неизвестные type в injections.csv: {bad}; допустимо: {TYPES}")
    tg = targets(a.series)
    rows, traces = [], {}
    for _, r in inj.iterrows():
        p = d / r.file
        if not p.exists():
            cand = list(d.glob(Path(r.file).stem + ".mzML")) + list(d.glob(Path(r.file).stem + ".mzml"))
            if not cand:
                print(f"нет файла: {p}")
                continue
            p = cand[0]
        sp = read_ms1(p)
        print(f"{p.name}: {len(sp)} MS1-спектров, RT {sp[0][0]:.2f}–{sp[-1][0]:.2f} мин")
        for _, t in tg.iterrows():
            rt, y = xic(sp, t.mz, a.ppm)
            traces[(r.file, t.label, t.rt)] = (rt, y)
            rows.append(dict(file=r.file, type=r["type"], batch=r.get("batch", ""), order=r.get("order", np.nan),
                             sample=r.get("sample", ""), label=t.label, kind=t.kind, in_index=t.in_index,
                             mz=round(t.mz, 4), rt_expected=t.rt))
    long = pd.DataFrame(rows)
    # сдвиг RT: медиана по QC для каждой мишени (по апексам без сдвига), затем повторный поиск со сдвигом
    tmp = [peak(*traces[(f, l, rt0)], rt0) for f, l, rt0 in zip(long.file, long.label, long.rt_expected)]
    long["rt_apex0"] = [x[0] for x in tmp]
    qc = long[long["type"] == "QC"]
    shift = (qc.rt_apex0 - qc.rt_expected).groupby([qc.label, qc.rt_expected]).median()
    res = []
    for _, x in long.iterrows():
        sh = shift.get((x.label, x.rt_expected), 0.0)
        res.append(peak(*traces[(x.file, x.label, x.rt_expected)], x.rt_expected, 0.0 if pd.isna(sh) else sh))
    long["rt_apex"], long["height"], long["area"] = zip(*res)
    long = long.drop(columns="rt_apex0")
    long.to_csv(out / "blank_check_areas.csv", index=False)

    # сводка: медиана площади по типу инъекции и доля от медианы QC
    key = ["label", "kind", "in_index", "mz", "rt_expected"]
    summ = long.pivot_table(index=key, columns="type", values="area", aggfunc="median").reindex(columns=TYPES)
    for t in TYPES:
        if t != "QC" and t in summ:
            summ[f"{t}_frac_of_QC"] = summ[t] / summ["QC"]
    flags = []
    for _, s in summ.iterrows():
        fr = {t: s.get(f"{t}_frac_of_QC", np.nan) for t in ["solvent_blank", "process_blank", "container_blank"]}
        present = [t for t, v in fr.items() if pd.notna(v) and v >= BLANK_FRAC]
        if not present:
            flags.append("not in blanks (< 10% of QC)")
        elif "solvent_blank" in present:
            flags.append("in solvent blank: instrument/solvent background or carryover")
        else:
            flags.append("in process/container blank only: consumables (tubes/containers)")
    summ["interpretation"] = flags
    summ.reset_index().to_csv(out / "blank_check_summary.csv", index=False)

    # перенос: бланк сразу после пробы/QC. Порядок — столбец order (внутри батча), иначе время запуска из mzML;
    # по времени пара считается соседней, если между запусками < MAX_GAP_MIN (иначе между ними были неконвертированные инъекции)
    co = []
    if "order" in inj and inj.order.notna().all():
        o = inj.sort_values(["batch", "order"]).reset_index(drop=True)
        o["adjacent"] = True
    else:
        o = inj.assign(t=inj.file.map(lambda f: start_time(d / f) if (d / f).exists() else pd.NaT)).dropna(subset=["t"])
        o = o.sort_values(["batch", "t"]).reset_index(drop=True)
        o["gap_min"] = o.groupby("batch").t.diff().dt.total_seconds() / 60
        o["adjacent"] = o.gap_min < MAX_GAP_MIN
        o.to_csv(out / "injection_order_from_timestamps.csv", index=False)
    for i in range(1, len(o)):
        if o.loc[i, "type"].endswith("blank") and o.loc[i - 1, "type"] in ("QC", "sample") \
                and o.loc[i, "batch"] == o.loc[i - 1, "batch"] and bool(o.loc[i, "adjacent"]):
            for lab in tg.label.unique():
                a_ = long[(long.file == o.loc[i, "file"]) & (long.label == lab)].area.sum()
                b_ = long[(long.file == o.loc[i - 1, "file"]) & (long.label == lab)].area.sum()
                co.append(dict(blank=o.loc[i, "file"], previous=o.loc[i - 1, "file"], label=lab,
                               carryover_pct=100 * a_ / b_ if b_ > 0 else np.nan))
    pd.DataFrame(co, columns=["blank", "previous", "label", "carryover_pct"]).to_csv(out / "blank_check_carryover.csv", index=False)

    # рисунок: площади по типам инъекций (log10) + XIC главного члена ряда
    labs = list(tg.label.drop_duplicates())
    fig, axs = plt.subplots(1, 2, figsize=(17 / 2.54, 0.42 / 2.54 * len(labs) + 3), gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axs[0]
    present_types = [t for t in TYPES if t in set(long["type"])]
    for k, lab in enumerate(labs):
        for j, t in enumerate(present_types):
            v = long[(long.label == lab) & (long["type"] == t)].area
            yy = np.full(len(v), k) + (j - (len(present_types) - 1) / 2) * 0.14
            ax.scatter(np.log10(v.clip(lower=1)), yy, s=8, c=COL[t], label=t if k == 0 else None, lw=0)
    ax.set_yticks(range(len(labs)), labs, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("log$_{10}$ peak area (XIC ± %.0f ppm)" % a.ppm, fontsize=7)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
    ax.tick_params(labelsize=6)
    ax = axs[1]
    main_lab = tg[tg.in_index == 1].label.iloc[0]
    rt0 = tg[tg.label == main_lab].rt.iloc[0]
    for t in present_types:
        f = inj[inj["type"] == t].file.iloc[0]
        if (f, main_lab, rt0) in traces:
            rt, y = traces[(f, main_lab, rt0)]
            m = np.abs(rt - rt0) <= 1.5
            ax.plot(rt[m], y[m], c=COL[t], lw=0.8, label=f"{t} ({Path(f).stem})")
    ax.set(xlabel="Retention time (min)", ylabel="Intensity", title=f"XIC {main_lab} m/z {tg[tg.label == main_lab].mz.iloc[0]:.4f}")
    ax.title.set_fontsize(7)
    ax.xaxis.label.set_fontsize(7)
    ax.yaxis.label.set_fontsize(7)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=5.5, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "blank_check.png", dpi=300)
    fig.savefig(out / "blank_check.pdf")
    print("\n" + summ.reset_index()[["label", "interpretation"] + [c for c in summ.columns if c.endswith("frac_of_QC")]].round(3).to_string(index=False))
    print(f"\nГотово → {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Экспорт выравнивания MS-DIAL (Alignment result, Raw data matrix: Area) -> таблица для фреймворка.

    python tools/msdial_to_table.py export.txt --out data/s300/s300_table.csv
    python tools/msdial_to_table.py export.txt --out data/s300/s300_table.csv --sample-pattern "(\\d+)$"

Вход — текстовый файл MS-DIAL (разделитель — табуляция), как он сохраняется командой Export → Alignment result:
над строкой заголовков («Alignment ID», «Average Rt(min)», «Average Mz», «Metabolite name», …) MS-DIAL пишет служебные
строки «Class», «File type», «Injection order», «Batch ID» для столбцов инъекций; после столбцов инъекций могут идти
средние и SD по классам («Average», «Stdev») — они отбрасываются.

Выход (рядом с --out):
  <out>                         таблица: Sample, Group, затем по столбцу на признак «RT_mz_название» (как в s266/s107);
                                значения — площади, 0 = не обнаружено (заполнение пропусков в MS-DIAL не выполняется);
  <out>_injections.csv          file, sample, type, class, batch, injection_order — журнал инъекций из проекта MS-DIAL;
  <out>_metadata_template.csv   заготовка metadata/*.csv: sample, group, sex, age, include, exclude_reason, injection_order,
                                batch, comment (пол и возраст заполняются вручную);
  <out>_ms2.csv                 feature, Adduct type, Formula, MS/MS assigned, MS/MS matched, Total score, Dot product, …
                                (признаки с MS/MS matched = True -> ms2_identified.csv фреймворка для уровня MSI 2).

Имя признака собирается из исходных строк «Average Rt(min)», «Average Mz» и «Metabolite name» без переформатирования
чисел; пустое название -> «Unknown». Префиксы MS-DIAL («low score:», «w/o MS2:») сохраняются.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

INFO_ROWS = {"class": "class", "file type": "type", "injection order": "injection_order", "batch id": "batch",
             "batch": "batch", "analytical order": "injection_order"}
META_COLS = ["Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type", "Post curation result",
             "Fill %", "MS/MS assigned", "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY", "SMILES",
             "Annotation tag (VS1.0)", "RT matched", "m/z matched", "MS/MS matched", "Comment",
             "Manually modified for quantification", "Manually modified for annotation", "Isotope tracking parent ID",
             "Isotope tracking weight number", "RT similarity", "m/z similarity", "Simple dot product",
             "Weighted dot product", "Reverse dot product", "Matched peaks count", "Matched peaks percentage",
             "Total score", "Dot product", "Fragment presence %", "S/N average", "Spectrum reference file name",
             "MS1 isotopic spectrum", "MS/MS spectrum"]
MS2_COLS = ["Adduct type", "Formula", "Ontology", "INCHIKEY", "MS/MS assigned", "MS/MS matched", "Total score",
            "Dot product", "Simple dot product", "Weighted dot product", "Reverse dot product", "Matched peaks count",
            "Matched peaks percentage", "Spectrum reference file name"]
SAMPLE_TYPES = {"sample", "qc", "blank", "standard"}


def read_export(path):
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    hdr = next((i for i in range(min(len(raw), 50))
                if {"Alignment ID", "Average Mz"} <= set(raw.iloc[i].str.strip())), None)
    if hdr is None:
        sys.exit("не найдена строка заголовков с «Alignment ID» и «Average Mz» — это экспорт выравнивания MS-DIAL?")
    cols = [c.strip() for c in raw.iloc[hdr]]
    info = {}
    for i in range(hdr):
        row = raw.iloc[i].str.strip()
        lab = next((INFO_ROWS[v.lower()] for v in row if v.lower() in INFO_ROWS), None)
        if lab:
            info[lab] = list(row)
    body = raw.iloc[hdr + 1:].reset_index(drop=True)
    body.columns = cols
    body = body[body["Alignment ID"].str.strip() != ""]
    return cols, info, body


def injection_columns(cols, info):
    """Столбцы инъекций: после последнего служебного столбца; по строке File type (если есть) — только Sample/QC/Blank."""
    last_meta = max(i for i, c in enumerate(cols) if c in META_COLS)
    idx = list(range(last_meta + 1, len(cols)))
    if "type" in info:
        idx = [i for i in idx if info["type"][i].strip().lower() in SAMPLE_TYPES]
    else:                                                  # без служебных строк: до первого блока Average/Stdev
        stop = next((i for i in idx if cols[i].lower() in {"average", "stdev", "std", "sd"}), None)
        idx = [i for i in idx if stop is None or i < stop]
    names = [cols[i] for i in idx]
    rep = {c for c in names if names.count(c) > 1}         # имена классов повторяются в блоках Average и Stdev
    out = []
    for i in idx:
        if cols[i] in rep:
            break
        out.append(i)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", help="файл экспорта выравнивания MS-DIAL (.txt, табуляция)")
    ap.add_argument("--out", required=True, help="путь к таблице для фреймворка (.csv)")
    ap.add_argument("--sample-pattern", default=None,
                    help="регулярное выражение с одной группой: номер образца из имени файла (напр. «(\\d+)$»); "
                         "QC и бланки сохраняют имя файла")
    ap.add_argument("--qc-names", default=r"^qc", help="регулярное выражение для QC по имени, если нет строки File type")
    ap.add_argument("--blank-names", default=r"^(blank|blk|бланк)", help="то же для бланков")
    a = ap.parse_args()

    cols, info, body = read_export(a.export)
    inj = injection_columns(cols, info)
    if not inj:
        sys.exit("не найдены столбцы инъекций")
    files = [cols[i] for i in inj]

    def typ(i, name):
        if "type" in info and info["type"][i]:
            t = info["type"][i].strip().lower()
            return {"qc": "QC", "blank": "Blank", "sample": "Sample", "standard": "Standard"}.get(t, t)
        if re.search(a.qc_names, name, re.I):
            return "QC"
        if re.search(a.blank_names, name, re.I):
            return "Blank"
        return "Sample"

    rows = []
    for i, f in zip(inj, files):
        t = typ(i, f)
        cls = info["class"][i].strip() if "class" in info else ""
        sample = f
        if t == "Sample" and a.sample_pattern:
            m = re.search(a.sample_pattern, f)
            sample = str(int(m.group(1))) if m and m.group(1).isdigit() else (m.group(1) if m else f)
        rows.append(dict(file=f, sample=sample, type=t,
                         **{"class": cls or {"QC": "QC", "Blank": "Blank"}.get(t, "")},
                         batch=info.get("batch", [""] * len(cols))[i].strip(),
                         injection_order=info.get("injection_order", [""] * len(cols))[i].strip()))
    inj_tab = pd.DataFrame(rows)
    if inj_tab["sample"].duplicated().any():
        dup = inj_tab.loc[inj_tab["sample"].duplicated(keep=False), "file"].tolist()
        sys.exit(f"повторяющиеся имена образцов после --sample-pattern: {dup}")

    name = body["Metabolite name"].str.strip().replace("", "Unknown")
    feat = (body["Average Rt(min)"].str.strip() + "_" + body["Average Mz"].str.strip() + "_" + name).tolist()
    seen = {}
    for k, f in enumerate(feat):                           # одинаковые RT_mz_название (редко) — суффикс #2, #3
        seen[f] = seen.get(f, 0) + 1
        if seen[f] > 1:
            feat[k] = f"{f} #{seen[f]}"
    vals = body[files].replace("", "0").apply(pd.to_numeric, errors="coerce").fillna(0)
    vals.index = feat
    tab = vals.T
    tab.index = inj_tab["sample"].tolist()
    grp = [r["class"] if r["type"] == "Sample" else r["type"] for _, r in inj_tab.iterrows()]
    tab.insert(0, "Group", grp)
    tab.index.name = "Sample"

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out)
    stem = out.with_suffix("")
    inj_tab.to_csv(f"{stem}_injections.csv", index=False)
    smp = inj_tab[inj_tab.type == "Sample"]
    pd.DataFrame({"sample": smp["sample"], "group": smp["class"], "sex": "", "age": "", "include": 1, "exclude_reason": "",
                  "injection_order": smp.injection_order, "batch": smp.batch, "comment": ""}).to_csv(
        f"{stem}_metadata_template.csv", index=False)
    ms2 = body[[c for c in MS2_COLS if c in body.columns]].copy()
    ms2.insert(0, "feature", feat)
    ms2.to_csv(f"{stem}_ms2.csv", index=False)
    nmatched = int(ms2["MS/MS matched"].str.lower().eq("true").sum()) if "MS/MS matched" in ms2 else None
    print(f"{len(feat)} признаков; инъекций {len(files)}: " +
          ", ".join(f"{k} {v}" for k, v in inj_tab.type.value_counts().items()) +
          (f"; MS/MS matched: {nmatched}" if nmatched is not None else "") + f"\n-> {out}")


if __name__ == "__main__":
    main()

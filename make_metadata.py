#!/usr/bin/env python3
"""Однократное создание файла метаданных metadata/metadata.csv из листа последовательности.

После создания этот CSV — единственный источник групп, пола, возраста и исключений:
его правят вручную (например, уточнённые возрасты), а run_all.py читает только его.
Шифры проб не переносятся (внутренние коды; содержат служебные пометки).

    python make_metadata.py --seq data/sequence.xlsx --out metadata/metadata.csv
"""
import argparse
from pathlib import Path

import pandas as pd

# Правила (из комментариев к обработке):
#   контроль = «КГ...»; РМП -> BC; РПЖ -> PC;
#   71 (доброкачественное образование МП + рак простаты) -> PC;
#   42 (РПЖ + РМП + рак кожи) -> исключить; рак почки и поджелудочной -> исключить.
MANUAL = {
    "71": dict(group="PC", comment="доброкачественное образование мочевого пузыря + РПЖ -> группа РПЖ"),
    "42": dict(include=0, exclude_reason="множественные опухоли (РПЖ + РМП + рак кожи)"),
}


def classify(text):
    t = str(text).strip().lower()
    if t.startswith("кг"):
        return "Control", 1, ""
    if "почк" in t:
        return "KC", 0, "рак почки (вне дизайна)"
    if "поджелуд" in t:
        return "PanC", 0, "рак поджелудочной железы (вне дизайна)"
    if "рпж" in t or "рак предст" in t:
        return "PC", 1, ""
    if "рмп" in t:
        return "BC", 1, ""
    return "?", 0, "диагноз не распознан"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq", default="data/sequence.xlsx")
    ap.add_argument("--out", default="metadata/metadata.csv")
    a = ap.parse_args()
    seq = pd.read_excel(a.seq).dropna(subset=[pd.read_excel(a.seq).columns[0]])
    rows = []
    for _, r in seq.iterrows():
        s = str(int(r.iloc[0]))
        g, inc, why = classify(r["Диагноз / группа"])
        sex = str(r["Пол"]).strip()
        d = dict(sample=s, group=g, sex={"м": "M", "ж": "F"}.get(sex, ""),
                 age=r["Возраст"] if pd.notna(r["Возраст"]) else "", include=inc, exclude_reason=why, comment="")
        d.update(MANUAL.get(s, {}))
        rows.append(d)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df["age"] = pd.to_numeric(df["age"], errors="coerce").astype("Int64")
    df.to_csv(out, index=False)
    print(df.groupby(["group", "include"]).size())
    print(f"возраст отсутствует: {df.age.isna().sum()}; пол отсутствует: {(df.sex == '').sum()} -> {out}")


if __name__ == "__main__":
    main()

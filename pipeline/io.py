"""Чтение таблицы признаков и метаданных, порядок инъекций."""
import re
import warnings

import numpy as np
import pandas as pd

from . import config as C


def parse_features(cols):
    rows = []
    for c in cols:
        rt, mz, name = c.split("_", 2)
        low = name.startswith("low score:")
        clean = re.sub(r"^low score:\s*", "", name)
        conf = "unknown" if clean == "Unknown" else ("low" if low else "high")
        rows.append(dict(feature=c, rt_min=float(rt), mz=float(mz), annotation=clean, annotation_conf=conf))
    return pd.DataFrame(rows).set_index("feature")


def load_features(path=C.DATA):
    raw = pd.read_csv(path, encoding="utf-8-sig", dtype={"Sample": str}, low_memory=False)
    meta_cols = [c for c in ["Sample", "Group", "Diag"] if c in raw.columns]
    feat_cols = [c for c in raw.columns if c not in meta_cols]
    n_dup = raw.duplicated("Sample").sum()
    raw = raw.drop_duplicates("Sample")  # пробы с двойным диагнозом записаны дважды с одинаковыми данными
    X = raw.set_index("Sample")[feat_cols].astype(float).replace(0, np.nan)
    grp = raw.set_index("Sample").get("Group", pd.Series("", index=X.index)).astype(str)
    kind = pd.Series("sample", index=X.index)
    kind[X.index.str.upper().str.startswith(C.QC_PREFIX)] = "QC"
    is_blank = X.index.str.contains(C.BLANK_REGEX, flags=re.I, regex=True) | grp.str.contains(C.BLANK_REGEX, flags=re.I, regex=True)
    kind[is_blank] = "blank"
    return X, kind, parse_features(feat_cols), int(n_dup)


def load_metadata(path, kind):
    md = pd.read_csv(path, dtype={"sample": str})
    need = {"sample", "group", "sex", "age", "include"}
    miss = need - set(md.columns)
    if miss:
        raise ValueError(f"в {path} нет столбцов: {miss}")
    md = md.set_index("sample")
    md["sex"] = md["sex"].fillna("").astype(str).str.strip().str.upper().replace({"М": "M", "Ж": "F"})
    md["age"] = pd.to_numeric(md["age"], errors="coerce")
    md["include"] = md["include"].fillna(0).astype(int)
    samples = kind.index[kind == "sample"]
    not_in_md = samples.difference(md.index)
    if len(not_in_md):
        warnings.warn(f"пробы без метаданных (исключены): {list(not_in_md)}")
    not_in_data = md.index[md.include == 1].difference(samples)
    if len(not_in_data):
        warnings.warn(f"в метаданных, но нет в данных: {list(not_in_data)}")
    meta = pd.DataFrame(index=kind.index)
    meta["kind"] = kind
    meta = meta.join(md[["group", "sex", "age", "include"] + [c for c in ["exclude_reason", "comment"] if c in md]])
    meta.loc[meta.kind == "QC", "group"] = "QC"
    meta.loc[meta.kind == "blank", "group"] = "Blank"
    meta["include"] = meta["include"].fillna(0).astype(int)
    meta.loc[meta.kind != "sample", "include"] = 0
    bad = meta[(meta.include == 1) & ~meta.group.isin(C.STUDY_GROUPS)]
    if len(bad):
        warnings.warn(f"include=1, но группа вне {C.STUDY_GROUPS}: {list(bad.index)} -> исключены")
        meta.loc[bad.index, "include"] = 0
    meta["study"] = (meta.kind == "sample") & (meta.include == 1)
    return meta


def batch_of(meta):
    """Номер серии для проб и QC (C.BATCHES); QC распределяются по сериям по порядку: в каждой серии
    QC в начале, после каждых QC_EVERY проб и в конце серии (если последний блок неполный)."""
    if not C.BATCHES:
        return pd.Series(1, index=meta.index)
    b = pd.Series(np.nan, index=meta.index)
    for k, (lo, hi) in enumerate(C.BATCHES, 1):
        for s in meta.index[meta.kind == "sample"]:
            n = int(re.sub(r"\D", "", s) or 0)
            if lo <= n <= hi:
                b[s] = k
    qcs = sorted(meta.index[meta.kind == "QC"], key=lambda s: int(re.sub(r"\D", "", s) or 0))
    it = iter(qcs)
    for k, (lo, hi) in enumerate(C.BATCHES, 1):
        nq = 1 + (hi - lo + 1) // C.QC_EVERY + (1 if (hi - lo + 1) % C.QC_EVERY else 0)
        for _ in range(nq):
            q = next(it, None)
            if q is not None:
                b[q] = k
    return b


def injection_order(meta):
    if C.BATCHES:
        return _injection_order_batches(meta)
    return _injection_order_single(meta)


def _injection_order_batches(meta):
    """Серии подряд; в каждой: QC, затем блоки по QC_EVERY проб, после каждого блока — QC."""
    b = batch_of(meta)
    order, pos = {}, 0
    for k, (lo, hi) in enumerate(C.BATCHES, 1):
        qcs = iter(sorted(meta.index[(meta.kind == "QC") & (b == k)], key=lambda s: int(re.sub(r"\D", "", s) or 0)))
        smp = sorted(meta.index[(meta.kind == "sample") & (b == k)], key=lambda s: int(re.sub(r"\D", "", s) or 0))
        q = next(qcs, None)
        if q:
            order[q] = pos; pos += 1
        for i, s in enumerate(smp, 1):
            order[s] = pos; pos += 1
            if i % C.QC_EVERY == 0 or i == len(smp):
                q = next(qcs, None)
                if q:
                    order[q] = pos; pos += 1
    for x in meta.index[meta.kind == "blank"]:
        order[x] = np.nan
    return pd.Series(order).reindex(meta.index)


def _injection_order_single(meta):
    """Порядок инъекций: образцы по номеру, QC1 – первым, QC(k+1) – после каждых QC_EVERY образцов.
    Если в metadata.csv есть столбец injection_order, используется он (для образцов)."""
    qcs = sorted(meta.index[meta.kind == "QC"], key=lambda s: int(re.sub(r"\D", "", s) or 0))
    smp = sorted(meta.index[meta.kind == "sample"], key=lambda s: (len(s), s))
    order, pos = {}, 0
    qi = iter(qcs)
    q = next(qi, None)
    if q:
        order[q] = pos
    for i, s in enumerate(smp, 1):
        pos += 1
        order[s] = pos
        if i % C.QC_EVERY == 0:
            q = next(qi, None)
            if q:
                pos += 1
                order[q] = pos
    for b in meta.index[meta.kind == "blank"]:
        order[b] = np.nan
    return pd.Series(order).reindex(meta.index)

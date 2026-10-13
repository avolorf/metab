"""Курированные кандидаты: проверка массы и ID, уровень MSI, подтверждающие ионы, детекция, связи с конфаундерами."""
import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .features import ADDUCTS, confirming_ions, formula_mass, ppm


def hmdb_alternatives(mz, tol=3.0, top=6):
    from mummichog.compound_dicts import hmdb_dict
    out = []
    for ion in ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+H-H2O]+"]:
        M = mz - ADDUCTS[ion]
        forms = {}
        for k, (f, m) in hmdb_dict.items():
            if m and abs(M - float(m)) / M * 1e6 <= tol:
                forms.setdefault(f, [float(m), []])[1].append(k)
        for f, (m, ids) in forms.items():
            out.append(f"{ion} {f} ({ppm(M, m):+.1f} ppm; {len(ids)} HMDB: {','.join(ids[:3])})")
    return " | ".join(out[:top])


def msi_level(r):
    if str(r.get("standard_confirmed", "")).strip().lower() in {"1", "yes", "да", "true"}:
        return 1
    if not r["mass_ok"]:
        return 4
    try:
        if float(r.get("ms2_score")) >= 0.7:
            return 2
    except (TypeError, ValueError):
        pass
    return 3


def curate(path, info_all, info, Lraw_all, L, meta, groups, feat_table):
    from mummichog.compound_dicts import hmdb_dict, kegg_dict
    cur = pd.read_csv(path, dtype=str).fillna("")
    cur["rt_min"] = cur.rt_min.astype(float)
    cur["mz"] = cur.mz.astype(float)
    study = meta.index[meta.study]
    rows = []
    for _, r in cur.iterrows():
        d = info_all[((info_all.mz - r.mz).abs() < 0.002) & ((info_all.rt_min - r.rt_min).abs() < 0.02)]
        f = d.index[0] if len(d) else None
        M = formula_mass(r.formula)
        theo = M + ADDUCTS[r.ion]
        e = ppm(r.mz, theo)
        x = r.to_dict()
        x.update(feature=f, theo_mz=round(theo, 5), ppm=round(e, 2), mass_ok=abs(e) <= C.PPM)
        x["hmdb"] = r.hmdb if r.hmdb and r.hmdb in hmdb_dict and abs(float(hmdb_dict[r.hmdb][1]) - M) < 0.005 and x["mass_ok"] else ""
        x["kegg"] = r.kegg if r.kegg and r.kegg in kegg_dict and abs(float(kegg_dict[r.kegg][1] or 0) - M) < 0.005 and x["mass_ok"] else ""
        x["MSI_level"] = msi_level(x)
        if f is not None:
            x["confirming_ions"] = confirming_ions(f, r.formula, info_all, Lraw_all) if x["mass_ok"] else ""
            x["hmdb_alternatives_3ppm"] = hmdb_alternatives(r.mz) if abs(e) > 3 or not x["mass_ok"] else ""
            if f in feat_table.index:
                ft = feat_table.loc[f]
                for g in C.STUDY_GROUPS:
                    x[f"det_{g}"] = ft[f"det_{g}"]
                    x[f"imputed_{g}"] = 1 - ft[f"det_{g}"]
                x["kept_after_filters"] = bool(ft.kept)
            if f in L.columns:
                x["compound"] = groups.loc[f, "compound"]
                x["is_representative"] = bool(groups.loc[f, "is_representative"])
                v = L.loc[study, f]
                for cov, lab in [("contam_index", "contam"), ("dilution_PQN", "dilution")]:
                    rho = stats.spearmanr(v, meta.loc[study, cov])
                    x[f"rho_{lab}"], x[f"p_{lab}"] = rho.statistic, rho.pvalue
                    ctrl = meta.index[meta.study & (meta.group == "Control")]
                    rho = stats.spearmanr(L.loc[ctrl, f], meta.loc[ctrl, cov])
                    x[f"rho_{lab}_controls"], x[f"p_{lab}_controls"] = rho.statistic, rho.pvalue
        rows.append(x)
    df = pd.DataFrame(rows)
    ok = df[df.mass_ok]
    bad = df[~df.mass_ok]
    return ok, bad

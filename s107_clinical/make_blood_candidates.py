"""Кандидаты в «кровяную» оценку (для согласования до расчёта): соединения основного анализа s107 (РМП vs контроль),
вероятно поступающие с кровью, — пептиды, повышенные при РМП, плюс биливердин и LysoPC.

Пептид: библиотечная аннотация представителя — ди-/трипептид (то же правило PEP, что в report_s107.py) или пептид по
курации (curation_s107: Leu-Thr-Ser, Glu-Ile-Val 2,57 мин, Leu-Leu-Tyr, Gln-Pro / N-мевалонилкреатинин 244,129).
Повышен при РМП: log2FC > 0 и q < 0,05 в основной модели. Предлагаемый список (proposed = True) — пептиды, которые
в основном анализе выполнили все критерии находки (устойчивы к поправке на возраст, пол, индекс и к подбору пар).
"""
import re
import sys
from pathlib import Path

CWD = Path.cwd()
sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_matrix  # noqa: E402

import pandas as pd  # noqa: E402

AA = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
PEP = re.compile(rf"\b(?:{AA})(?:-(?:{AA}))+\b|alanyl|glycyl|prolyl|leucyl|valyl|isoleucyl|glutamyl", re.I)
CURATED_PEPTIDES = [("Gln-Pro (отчёт: N-mevalonylcreatinine)", 244.1293, 2.578), ("Glu-Ile-Val (2,57 мин)", 360.2132, 2.569),
                    ("Leu-Thr-Ser", 320.1823, 2.683), ("Leu-Leu-Tyr (neurotensin 11–13)", 408.2497, 9.225)]

d = prepare_matrix.build(lambda s: None)
info, g, main, names = d["info"], d["groups"], d["main"], d["names"]
r = main.join(info[["rt_min", "mz", "annotation", "annotation_conf"]])
rep_of = lambda mz, rt: g.loc[info.index[((info.mz - mz).abs() < 0.003) & ((info.rt_min - rt).abs() < 0.05)][0], "representative"]
cur = {rep_of(mz, rt): lab for lab, mz, rt in CURATED_PEPTIDES}
pep = r.annotation.astype(str).str.contains(PEP) | r.index.isin(list(cur))
up = (r.q_model < 0.05) & (r.beta_model > 0)
# стероиды (глюкурониды блока А и соединения с аннотацией стероида у любого иона) в оценку не входят:
# так, представитель глюкуронида тетрагидрокортизона II (m/z 558,29) аннотирован MS-DIAL как «…valyl…valine»
from run_clinical import STEROIDS, STEROID_RE  # noqa: E402
ster_comp = {g.loc[rep_of(mz, rt), "compound"] for _, mz, rt in STEROIDS}
ster_comp |= set(g.loc[info.loc[g.index, "annotation"].astype(str).str.contains(STEROID_RE)].compound)
not_steroid = ~g.loc[r.index, "compound"].isin(ster_comp).to_numpy()
t = r[pep & up & not_steroid].copy()
t["kind"] = "пептид"
t["label"] = [cur.get(i, a) for i, a in zip(t.index, t.annotation)]
extra = []
bil = rep_of(583.2553, 13.134)
b = r.loc[[bil]].assign(kind="желчный пигмент", label="Biliverdin")
lyso = info[info.annotation.str.contains(r"lysoPC|lysophosphatidylcholine|LPC\b|PC\(\d+:\d+/0:0\)|glycerophosphocholine", case=False)]
out = pd.concat([t, b])
out["proposed"] = out.finding.astype(bool) | (out.kind == "желчный пигмент")
cols = ["kind", "label", "rt_min", "mz", "annotation", "annotation_conf", "beta_model", "q_model", "AUC", "finding", "proposed"]
out = out[cols].rename(columns={"beta_model": "log2FC_BC_vs_Control", "q_model": "q_BC_vs_Control", "AUC": "AUC_BC_vs_Control",
                                "finding": "finding_main_analysis"}).sort_values(["kind", "q_BC_vs_Control"])
dest = CWD / (sys.argv[1] if len(sys.argv) > 1 else "../results_s107_clinical/blood_score_candidates.csv")
dest.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(dest, encoding="utf-8-sig", index_label="compound")
print(f"исключено как стероиды: {int((pep & up & ~not_steroid).sum())}; пептидов, повышенных при РМП (q < 0,05): {len(t)}; из них находок: {int(t.finding.sum())}; "
      f"биливердин: q = {b.q_model.iloc[0]:.2f}; LysoPC после фильтров: {len(lyso)}")
print(f"предложено: {int(out.proposed.sum())} → {dest}")

#!/usr/bin/env python3
"""Кандидаты для курации в выборке 107 проб: 30 маркеров отчёта (Табл. 20 / ABC) и альтернативные аннотации
(res_annotated_2.csv). Каждый кандидат привязывается к ближайшему признаку таблицы (±5 ppm, ±0.3 мин);
формула, ион, категория — из источника. Категория: endogenous / food / drug / tobacco (экзогенные не считаются находками)."""
import re
import pandas as pd
from pipeline.features import formula_mass, ADDUCTS, ppm

CAT = {"Nepetazide": "food: vegetable glycoside", "5-Methylmyosmine": "tobacco", "Acetylcysteine": "drug",
       "4-hydroxy-3-(2-hydroxyethyl)acetophenone": "food: plant phenolic", "Isoferulic acid": "food: phenolic acid",
       "Kmeriol": "food: plant phenylpropanoid", "N-Mevalonylcreatinine": "endogenous: creatinine adduct",
       "N-acetyl-S-(N-allylthiocarbamoyl)cysteine": "food: allyl isothiocyanate (cruciferous) mercapturate",
       "Threonyl-phenylalanine (Thr-Phe)": "endogenous: dipeptide",
       "3,6-bis(2-(methylthio)ethyl)-1,4-dioxane-2,5-diol": "food: methionine supplement product",
       "D-cyclohexanone glycoside": "food: plant glycoside", "Venlafaxine": "drug",
       "Asparaginyl-phenylalanine (Asn-Phe)": "endogenous: dipeptide", "3-Hydroxyoctanoylcarnitine": "endogenous: acylcarnitine",
       "4-Methylcinnamylcarnitine": "endogenous: acylcarnitine",
       "(2S)-2-((2-morpholinoacetyl)amino)-4-phenylbutanoic acid": "drug: carfilzomib metabolite",
       "N, N-Dimethylsphingosine": "endogenous: sphingolipid", "Glutamyl-isoleucyl-valine (Glu-Ile-Val)": "endogenous: tripeptide",
       "D-Vanillyl glycol glycoside": "food: plant glycoside", "Fraxin": "food: coumarin glycoside",
       "D-glucuronide of perillyl glycolic acid": "food: monoterpene metabolite", "Samuraciclib": "drug",
       "D-glucuronide of 3-(hydroxymethyl)antipyrine": "drug: antipyrine metabolite",
       "Grayanotoxin I": "food: honey phytotoxin",
       "9-β-D-glucuronide of 12Z-9, 10-dihydroxyoctadecenoic acid": "endogenous: oxylipin glucuronide",
       "D-glucuronide of γ-CEHC": "endogenous: vitamin E metabolite",
       "12-D-glucuronide of 9Z-12, 13-dihydroxyoctadecenoic acid": "endogenous: oxylipin glucuronide",
       "Taurinyl-(9Z)-3-hydroxydodecenoylcarnitine": "endogenous: acylcarnitine",
       "N-acetylcysteinyl-2, 6-diisopropylphenol": "drug: propofol–N-acetylcysteine adduct"}


# Категория из отчёта сохраняется, только если она подтверждается данными: аннотация MS-DIAL того же признака
# (в т. ч. low score) называет то же соединение (синонимы указаны). Остальные категории отчёта не подтверждены
# (MS-DIAL предлагает другое соединение или масса не соответствует формуле) -> «ambiguous»; исходная категория — в comment.
CONFIRMED_BY_DATA = {"Acetylcysteine": "N-Acetyl-L-cysteine",
                     "4-hydroxy-3-(2-hydroxyethyl)acetophenone": "4-Hydroxy-3-(2-hydroxyethyl)acetophenone",
                     "Isoferulic acid": "3-Hydroxy-4-methoxycinnamic acid", "Venlafaxine": "Venlafaxine (low score)",
                     "Threonyl-phenylalanine (Thr-Phe)": "Thr-Phe", "Asparaginyl-phenylalanine (Asn-Phe)": "Asn-Phe"}


def snap(info, mz, rt, rt_tol=0.3):
    d = info[(ppm(info.mz, mz).abs() <= 5) & ((info.rt_min - rt).abs() <= rt_tol)]
    if not len(d):
        return None
    return (d.rt_min - rt).abs().idxmin()


def main():
    cols = list(pd.read_csv("data/s107/s107_table.csv", nrows=0).columns[2:])
    info = pd.DataFrame([dict(feature=c, rt_min=float(c.split("_")[0]), mz=float(c.split("_")[1])) for c in cols]).set_index("feature")
    rows, log = [], []
    t = pd.read_csv("data/s107/table_nastya.csv")
    for _, r in t.iterrows():
        mz, ion = r.Characteristic_ion_mz.split(" ", 1)
        mz = float(mz)
        if ion not in ADDUCTS:
            log.append(f"ABC {r.No} {r.Compound}: ион {ion} не поддерживается, пропущен")
            continue
        f = snap(info, mz, float(r.Retention_time_min))
        if f is None:
            log.append(f"ABC {r.No} {r.Compound}: нет признака ±5 ppm/±0.3 мин")
            continue
        cat = CAT.get(r.Compound, "ambiguous")
        msd = f.split("_", 2)[2]
        if r.Compound in CONFIRMED_BY_DATA:
            note = f"category confirmed by data: MS-DIAL annotation of this feature = {CONFIRMED_BY_DATA[r.Compound]}"
        else:
            note = f"report category '{cat}' not confirmed by data (MS-DIAL annotation of this feature: {msd}) -> ambiguous"
            if r.No == 10:
                note += ("; spectra checked in FreeStyle (samples 4, 56, QC17, Blank15): absent in blank, no fragment "
                         "m/z 137.0458 in sample or QC, inconsistent isotope pattern (charge assigned z = 16) -> "
                         "annotation not confirmed (MSI 4)")
            cat = "ambiguous"
        rows.append(dict(rt_min=info.loc[f, "rt_min"], mz=info.loc[f, "mz"], focus="BC", name=f"ABC{r.No}: {r.Compound}",
                         name_en=r.Compound, formula=r.Molecular_formula, ion=ion, hmdb="", kegg="",
                         category=cat, comment=f"Report Table 20 / ABC No. {r.No}: {r.Notes}; {note}",
                         msi_override=4 if r.No == 10 else "",
                         ms2_library="", standard_confirmed="", isomer_ambiguous="", isomers_note=""))
    a = pd.read_csv("data/s107/res_annotated_2.csv", index_col=0)
    for _, r in a.iterrows():
        rt, mz = map(float, re.match(r"X([\d.]+)_([\d.]+)", r.x).groups())
        f = snap(info, mz, rt, 0.02)
        name = str(r["Putative structure"]).strip()
        if f is None:
            log.append(f"res_annotated {name}: нет признака")
            continue
        if any(abs(x["mz"] - info.loc[f, "mz"]) < 0.002 and abs(x["rt_min"] - info.loc[f, "rt_min"]) < 0.02 for x in rows):
            for x in rows:
                if abs(x["mz"] - info.loc[f, "mz"]) < 0.002 and abs(x["rt_min"] - info.loc[f, "rt_min"]) < 0.02:
                    x["comment"] += f"; alternative annotation (res_annotated_2): {name}, {r['elemental composition']}"
            continue
        # пептиды — эндогенные; 5-гидроксиомепразол (m/z на 22 ppm от формулы, MS-DIAL — безафибрат) и
        # N-ацетилтриптофан (m/z не соответствует формуле) по данным не подтверждены -> ambiguous
        cat = "endogenous: peptide" if re.search(r"^[A-Z][a-z]{2}(-[A-Za-z]{3})+|proline", name, re.I) else "ambiguous"
        extra = ("; spectra checked in FreeStyle (samples 4, 56, QC17, Blank15): weak peak, no Cl isotope (364.1219) "
                 "at this RT, no MS/MS for precursor 362.12 (362.1512 at 8.2-8.3 min is another compound) -> ambiguous"
                 if "omeprazole" in name.lower() else "")
        rows.append(dict(rt_min=info.loc[f, "rt_min"], mz=info.loc[f, "mz"], focus="BC", name=f"res_annotated: {name}",
                         name_en=name, formula=r["elemental composition"], ion="[M+H]+", hmdb="", kegg="", category=cat,
                         comment=f"res_annotated_2.csv; MS-DIAL: {r.x}{extra}", msi_override="", ms2_library="", standard_confirmed="",
                         isomer_ambiguous="", isomers_note=""))
    pd.DataFrame(rows).to_csv("curation_s107/candidates_curated.csv", index=False)
    print(f"{len(rows)} кандидатов"); print("\n".join(log))


if __name__ == "__main__":
    main()

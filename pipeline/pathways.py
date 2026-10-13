"""mummichog (приложение): без экспозиций, p из основной модели, поправка BH по путям."""
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C
from .stats import bh


def run(stats_feat, info, comp, outdir, exclude, seed):
    """stats_feat: признаки (не соединения) с p_model / beta_model для данного сравнения."""
    d = stats_feat[~stats_feat.index.isin(exclude)].dropna(subset=["p_model"])
    tab = pd.DataFrame({"mz": info.loc[d.index, "mz"], "rt": info.loc[d.index, "rt_min"] * 60, "p": d.p_model,
                        "t": np.sign(d.beta_model) * -np.log10(d.p_model), "id": d.index})
    cutoff = float(np.quantile(tab.p, C.MUMMICHOG_TOP_FRAC))
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tab.to_csv(Path(tmp) / "input.txt", sep="\t", index=False)
        np.random.seed(seed)
        subprocess.run(["mummichog", "-f", "input.txt", "-o", comp, "-m", "positive", "-u", str(C.PPM),
                        "-c", f"{cutoff:.3g}", "-p", str(C.MUMMICHOG_PERM)], cwd=tmp, check=True,
                       capture_output=True)
        res = next(p for p in Path(tmp).iterdir() if p.is_dir())
        shutil.rmtree(outdir)
        shutil.copytree(res / "tables", outdir)
    pw = pd.read_csv(outdir / f"mcg_pathwayanalysis_{comp}.tsv", sep="\t")
    pw["q_BH"] = bh(pw["p-value"])
    pw.insert(0, "comparison", comp)
    settings = dict(comparison=comp, n_features=len(tab), n_significant=int((tab.p < cutoff).sum()),
                    p_cutoff=cutoff, excluded_exposure_features=len(exclude), permutations=C.MUMMICHOG_PERM,
                    n_pathways_p05=int((pw["p-value"] < 0.05).sum()), n_pathways_q05=int((pw.q_BH < 0.05).sum()))
    return pw, settings

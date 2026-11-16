"""Python reference for the end-to-end test of the R package: the one-call path of threadfin 4.6.1 on the user's
files of Tutorial 5 (GSE287123), with the arguments of the tutorial's Step 9, from raw counts and AIRR contigs.

usage: python python_reference.py <folder with gse287123_rbd.h5ad and gse287123_rbd_contigs.tsv> <n_jobs>
Writes <folder>/python_reference.pkl (every table and text of the result) and python_summary.txt.
"""
import json
import os
import sys
import time
from importlib.metadata import version

import anndata as ad
import pandas as pd
import threadfin as tf

folder, n_jobs = sys.argv[1], int(sys.argv[2])
h5ad = os.path.join(folder, "gse287123_rbd.h5ad")
tsv = os.path.join(folder, "gse287123_rbd_contigs.tsv")
adata = ad.read_h5ad(h5ad)
t0 = time.time()
res = tf.run(adata, bcr=tsv, donor_key="donor", sample_key="donor", strata_key="sample",
             test=["division_gate", "rbd_bait", "zone_gate"], n_jobs=n_jobs)
seconds = time.time() - t0
clone_key = res.settings["clone_key"]
cols = [c for c in adata.obs.columns
        if c in (clone_key, "clone_programme", "clone_programme_assigned") or str(c).startswith("bcr_")]
cells = adata.obs.loc[:, cols].copy()
cells.index = pd.Index(adata.obs_names.astype(str), name="cell", dtype=object)
ref = {
    "coherence": res.coherence, "programmes": res.programmes, "clones": res.clones, "composition": res.composition,
    "tests": res.tests, "profile_tests": res.profile_tests, "memory": res.memory, "settings": res.settings,
    "skipped": res.skipped, "summary": res.summary(), "methods": res.methods(), "cells": cells,
}
pd.to_pickle(ref, os.path.join(folder, "python_reference.pkl"))
with open(os.path.join(folder, "python_summary.txt"), "w") as fh:
    fh.write(res.summary() + "\n")
info = {"seconds": round(seconds, 1), "n_jobs": n_jobs, "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        "versions": {d: version(d) for d in ("threadfin", "sc-dandelion", "pandas", "numpy", "anndata", "scanpy")},
        "clone_calling": {k: v for k, v in (res.settings.get("clone_calling") or {}).items()
                          if k in ("dandelion_hash_seed", "dandelion_subprocess", "n_clones", "n_clones_ge2")},
        "icc": res.coherence["icc"], "n_profiled": int(res.clones.shape[0])}
with open(os.path.join(folder, "python_reference_info.json"), "w") as fh:
    json.dump(info, fh, indent=1, default=str)
print(json.dumps(info, indent=1, default=str))
print(res.summary())

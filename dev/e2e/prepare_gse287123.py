"""The user's files for the end-to-end test, made as Tutorial 5 (threadfin 4.6.1, docs/tutorials/05_end_to_end_tutorial.ipynb)
makes them: steps 1-4 and 6 (download GSE287123 with MD5 checks, read each tube, hashtag demultiplexing, one table of
mice and gates, the AIRR contigs of the singlets). Writes the expression data with raw counts (.h5ad) and the contig
table (.tsv) that a user would pass to threadfin.run.

usage: python prepare_gse287123.py <download folder> <output folder>

These files feed the opt-in end-to-end test of the R package (tests/testthat/test-e2e.R). From the repository folder,
with E2E a folder outside it (about 600 MB of downloads; each run needs about 15 GB of memory):

    python dev/e2e/prepare_gse287123.py $E2E/GSE287123 $E2E
    python dev/e2e/python_reference.py $E2E 4          # threadfin.run() in Python: python_reference.pkl
    THREADFIN_E2E_DIR=$E2E THREADFIN_E2E_NJOBS=4 Rscript -e 'testthat::test_local(".", filter = "e2e")'

Run the last two with the same Python environment (threadfin 4.6.1 from PyPI) and the same number of BLAS threads
(OMP_NUM_THREADS), on which the automatic clone map depends.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

DATA = Path(sys.argv[1])
OUT = Path(sys.argv[2])
DATA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- Step 1 (as in the tutorial)
TUBES = {  # tube: (GEO sample of the expression library, GEO sample of the BCR library, library name)
    "A": ("GSM8739143", "GSM8739145", "Rebutal_4-5_Tube_A"),
    "B": ("GSM8739146", "GSM8739148", "Rebutal_4-5_Tube_B"),
    "C": ("GSM8739149", "GSM8739151", "Rebutal_4-5_Tube_C"),
    "D": ("GSM8739152", "GSM8739154", "Rebutal_4-5_Tube_D"),
}
MD5 = {
    "GSM8739143_Rebutal_4-5_Tube_A_GEX_HTO_barcodes.tsv.gz": "589eb6503c2ac806de9b2db3fb81e433",
    "GSM8739143_Rebutal_4-5_Tube_A_GEX_HTO_features.tsv.gz": "f89add0fc8f0a6126328592c4f98e2b3",
    "GSM8739143_Rebutal_4-5_Tube_A_GEX_HTO_matrix.mtx.gz": "871c4c27db40fef6e85b22925f669939",
    "GSM8739146_Rebutal_4-5_Tube_B_GEX_HTO_barcodes.tsv.gz": "5ba3969c131f02d8bcb058f5650e3f87",
    "GSM8739146_Rebutal_4-5_Tube_B_GEX_HTO_features.tsv.gz": "fd0d96207c7e96ce218a7a73564c204a",
    "GSM8739146_Rebutal_4-5_Tube_B_GEX_HTO_matrix.mtx.gz": "6e7891bf6cf9a406aa94847ccab4028e",
    "GSM8739149_Rebutal_4-5_Tube_C_GEX_HTO_barcodes.tsv.gz": "a8a6ef4a3016159f776ca60b0c1feb3d",
    "GSM8739149_Rebutal_4-5_Tube_C_GEX_HTO_features.tsv.gz": "653fa5bdfd7707dc82a84e1281e36b69",
    "GSM8739149_Rebutal_4-5_Tube_C_GEX_HTO_matrix.mtx.gz": "fc59fab812357466c2d47a680d0491ef",
    "GSM8739152_Rebutal_4-5_Tube_D_GEX_HTO_barcodes.tsv.gz": "f4e934d72a9b49c3c479a3a54d01d3da",
    "GSM8739152_Rebutal_4-5_Tube_D_GEX_HTO_features.tsv.gz": "475ef5a635671345d145aa0a42782b22",
    "GSM8739152_Rebutal_4-5_Tube_D_GEX_HTO_matrix.mtx.gz": "76987b838c609f0a68887ec63fb8cf2c",
    "GSM8739145_Rebutal_4-5_Tube_A_airr_rearrangement.tsv.gz": "c673685b5ebc2179e1cdbeb54877086a",
    "GSM8739148_Rebutal_4-5_Tube_B_airr_rearrangement.tsv.gz": "e6460e998ff9d2099fe8cff671530f48",
    "GSM8739151_Rebutal_4-5_Tube_C_airr_rearrangement.tsv.gz": "ac1ceacd942cc268b9ce424b10287669",
    "GSM8739154_Rebutal_4-5_Tube_D_airr_rearrangement.tsv.gz": "08bce424fa0e9d1acdf7943062e75882",
}


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fetch(name):
    gsm = name.split("_")[0]
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
    path = DATA / name
    if path.exists() and md5(path) == MD5[name]:
        return f"{name}: already present, ok"
    for attempt in range(3):                       # GEO occasionally drops a connection
        try:
            urllib.request.urlretrieve(url, path)
            break
        except OSError as err:
            print(f"  retrying {name} ({err})")
    if md5(path) != MD5[name]:
        raise RuntimeError(f"checksum mismatch for {name}; delete it and run this cell again")
    return f"{name}: {path.stat().st_size / 1e6:.1f} MB, ok"


for name in MD5:
    print(fetch(name), flush=True)

# ---------------------------------------------------------------- Step 2
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1
tubes = {}
for tube, (gsm_gex, gsm_bcr, library) in TUBES.items():
    a = sc.read_10x_mtx(DATA, prefix=f"{gsm_gex}_{library}_GEX_HTO_", var_names="gene_symbols",
                        gex_only=False, make_unique=True)
    a.obs_names = [f"{tube}|{barcode}" for barcode in a.obs_names]
    a.obs["tube"] = tube
    tubes[tube] = a
    kinds = a.var["feature_types"].value_counts().to_dict()
    print(f"tube {tube}: {a.n_obs:,} barcodes; features {kinds}")

# ---------------------------------------------------------------- Step 3
from sklearn.cluster import KMeans


def demultiplex(hto_counts: pd.DataFrame, min_total: int = 5) -> pd.Series:
    x = hto_counts.to_numpy(dtype=float)
    clr = np.log1p(x) - np.log1p(x).mean(axis=0, keepdims=True)
    positive = np.zeros_like(clr, dtype=bool)
    for j in range(clr.shape[1]):
        km = KMeans(2, n_init=10, random_state=0).fit(clr[:, [j]])
        positive[:, j] = km.labels_ == int(np.argmax(km.cluster_centers_.ravel()))
    n_pos = positive.sum(axis=1)
    call = np.where(n_pos == 1, np.asarray(hto_counts.columns)[positive.argmax(axis=1)], "negative")
    call = np.where(n_pos > 1, "doublet", call)
    call = np.where(x.sum(axis=1) < min_total, "negative", call)
    return pd.Series(call, index=hto_counts.index, name="hashtag")


expression = []
for tube, a in tubes.items():
    is_hto = (a.var["feature_types"] == "Antibody Capture").to_numpy()
    hto = pd.DataFrame(a[:, is_hto].X.toarray(), index=a.obs_names, columns=a.var_names[is_hto])
    gex = a[:, ~is_hto].copy()
    gex.obs["hashtag"] = demultiplex(hto).values
    kind = gex.obs["hashtag"].where(gex.obs["hashtag"].isin(["doublet", "negative"]), "singlet")
    print(f"tube {tube}:", kind.value_counts().to_dict())
    expression.append(gex[kind.eq("singlet").to_numpy()].copy())

adata = ad.concat(expression, join="outer", fill_value=0)
print(f"\n{adata.n_obs:,} singlet cells x {adata.n_vars:,} genes")

# ---------------------------------------------------------------- Step 4
tag = adata.obs["hashtag"].astype(str)
adata.obs["donor"] = tag.str.split("_").str[0]
adata.obs["arm"] = np.where(tag.str.contains("_RBD_"), "RBD protein", "RBD mRNA")
adata.obs["rbd_bait"] = tag.map(lambda s: "RBD+" if "_RBD_pos" in s else ("RBD-" if "_RBD_neg" in s else None))
adata.obs["zone_gate"] = tag.map(lambda s: "LZ" if "_LZ_" in s else ("DZ" if "_DZ_" in s else None))
adata.obs["division_gate"] = np.where(tag.str.contains("mCherry_low"), "mCherry-low", "mCherry-high")
adata.obs["sample"] = (adata.obs["donor"] + "_" + adata.obs["rbd_bait"].fillna(adata.obs["zone_gate"]).astype(str)
                       + "_" + adata.obs["division_gate"])
for col in ["donor", "arm", "rbd_bait", "zone_gate", "division_gate", "sample", "tube"]:
    adata.obs[col] = adata.obs[col].astype("category")
print(adata.obs.groupby("arm", observed=True)["donor"].nunique().rename("mice").to_string(), "\n")

# ---------------------------------------------------------------- Step 6
contigs = []
for tube, (gsm_gex, gsm_bcr, library) in TUBES.items():
    t = pd.read_csv(DATA / f"{gsm_bcr}_{library}_airr_rearrangement.tsv.gz", sep="\t", low_memory=False)
    t["cell_id"] = tube + "|" + t["cell_id"].astype(str)
    contigs.append(t)
contigs = pd.concat(contigs, ignore_index=True)
contigs = contigs[contigs["cell_id"].isin(adata.obs_names)]            # singlets only
contigs["locus"] = contigs["v_call"].astype(str).str[:3]
print(f"{len(contigs):,} contigs from {contigs['cell_id'].nunique():,} singlet cells")

# ---------------------------------------------------------------- the user's two files
h5ad = OUT / "gse287123_rbd.h5ad"
tsv = OUT / "gse287123_rbd_contigs.tsv"
print("X dtype:", adata.X.dtype, type(adata.X).__name__)
adata.write_h5ad(h5ad, compression="gzip")
contigs.to_csv(tsv, sep="\t", index=False)
print(f"wrote {h5ad} ({h5ad.stat().st_size / 1e6:.0f} MB) and {tsv} ({tsv.stat().st_size / 1e6:.0f} MB)")

"""Example figures of the analysis part of the tutorial (3-1 to 3-5), from the three public zebrafish sections of 1-5.

Usage: python make_analysis_figures.py <folder with the .imzml/.IBD files> <output folder, e.g. ../figures> [cache.h5ad]
Environment: msi-analysis. Preprocessing as for make_figures.py (example_data.py): reference-ion tissue masks,
stand-in peak calling, ppm matching, requantification; then a permissive feature filter (detected in >= 5% of the
tissue pixels of at least one section, 2-5) and TIC normalization on these features (sc.pp.normalize_total, 2-6).
The optional cache file keeps the preprocessed data, so that reruns of the analysis take a few minutes.
"""
import os, re, sys, glob, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyimzml.ImzMLParser import ImzMLParser
from example_data import tic_images, tissue_mask, requantified
import analysis_examples
warnings.filterwarnings("ignore")

D, OUT = sys.argv[1], sys.argv[2]
CACHE = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(OUT, exist_ok=True)
PIXEL_UM = 7
acqs = sorted((os.path.basename(f)[:-6] for f in glob.glob(f"{D}/*.imzml")),
              key=lambda a: int(re.search(r"section(\d+)", a).group(1)))
lab = {a: re.search(r"section(\d+)", a).group(0).replace("section", "section ") for a in acqs}

if CACHE and os.path.exists(CACHE):
    T = ad.read_h5ad(CACHE)
else:
    M = {}
    for a in acqs:
        p = ImzMLParser(f"{D}/{a}.imzml", ibd_file=open(f"{D}/{a}.IBD", "rb"))
        M[a] = tissue_mask(tic_images(p)[2], PIXEL_UM, smooth=1)
    A = requantified({a: (f"{D}/{a}.imzml", f"{D}/{a}.IBD") for a in acqs}, M)
    for a in acqs:
        xy = A[a].obsm["spatial"]
        A[a].obs["sample"] = lab[a]; A[a].obs["x"], A[a].obs["y"] = xy[:, 0], xy[:, 1]
        A[a].var = A[a].var[["mz", "window"]]
    T = ad.concat([A[a] for a in acqs], merge="same")
    T.uns["windows"] = A[acqs[0]].uns["windows"]
    if CACHE:
        T.write_h5ad(CACHE, compression="gzip")

# permissive feature filter (2-5) and TIC normalization on the retained features (2-6)
det = pd.DataFrame({s: np.asarray((T[T.obs["sample"] == s].X > 0).mean(axis=0)).ravel() for s in T.obs["sample"].unique()},
                   index=T.var_names)
T = T[:, (det >= 0.05).any(axis=1).to_numpy()].copy()
T.layers["raw"] = T.X.copy()
sc.pp.normalize_total(T, target_sum=None)
print(f"{T.n_obs} tissue pixels x {T.n_vars} features")

num = analysis_examples.run(T, OUT, batch="sample", resolution=0.5)
for k, v in num.items():
    print(f"{k}: {v}")
print("figures written to", OUT)

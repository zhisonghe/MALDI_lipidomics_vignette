"""Reproduces the gap-filling and censoring numbers of tutorial 2-4 and 2-6 from the three zebrafish sections of 1-5.

Usage: python zebrafish_gapfill_censoring.py <folder with the .imzml/.IBD files>
Same preprocessing as the figures (example_data.py): reference-ion tissue masks (2-2), stand-in peak calling
(features in >= 5% of the tissue pixels), ppm matching at 3 ppm and requantification (2-4).
"""
import sys, re, glob, os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from pyimzml.ImzMLParser import ImzMLParser
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from example_data import tic_images, tissue_mask, requantified
from censor_stats import censor_summary

D = sys.argv[1]                                          # folder with the zebrafish .imzml/.IBD files
acqs = sorted(os.path.basename(f)[:-6] for f in glob.glob(f"{D}/*.imzml"))
lab = lambda a: re.search(r'section\d+', a).group(0)
pixel_um = 7

M = {}
for a in acqs:
    p = ImzMLParser(f"{D}/{a}.imzml", ibd_file=open(f"{D}/{a}.IBD", "rb"))
    M[a] = tissue_mask(tic_images(p)[2], pixel_um, smooth=1)   # summed reference ions (2-2)
    print(lab(a), "pixels", len(p.coordinates), "tissue", M[a].sum())
A = requantified({a: (f"{D}/{a}.imzml", f"{D}/{a}.IBD") for a in acqs}, M)

called = pd.DataFrame({a: A[a].var.called.values for a in acqs}, index=A[acqs[0]].var_names)
ncalled = called.sum(axis=1)
print("features:", len(called), "| called in all:", (ncalled == len(acqs)).sum(), "| in some but not all:",
      ((ncalled > 0) & (ncalled < len(acqs))).sum())
dets = []
for a in acqs:
    gap = ~called[a] & (ncalled > 0)
    d = np.asarray((A[a][:, gap.values].X > 0).mean(axis=0)).ravel(); dets.append(d)
    print(f"{lab(a)}: gaps {gap.sum()}, detection after requant: median {np.median(d):.3f}, >=5% {(d>=.05).sum()}, "
          f"1-5% {((d<.05)&(d>=.01)).sum()}, <1% {(d<.01).sum()}")
d = np.concatenate(dets); print(f"all gaps: {len(d)}; >=1% px: {(d>=.01).mean():.2f}; 1-5%: {((d<.05)&(d>=.01)).mean():.2f}; <1%: {(d<.01).mean():.2f}")

censor_summary(A, lab)

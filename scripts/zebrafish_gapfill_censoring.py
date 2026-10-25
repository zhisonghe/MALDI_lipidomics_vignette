"""Reproduces the gap-filling and censoring numbers of tutorial 2-4 and 2-6 from the three zebrafish sections of 1-5.

Usage: python zebrafish_gapfill_censoring.py <folder with the .imzml/.IBD files>
"""
import sys, re, glob, os, itertools, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from pyimzml.ImzMLParser import ImzMLParser
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
D = sys.argv[1]                                          # folder with the zebrafish .imzml/.IBD files
G = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, G); from gapfill import consensus_intervals, requantify

acqs = sorted(os.path.basename(f)[:-6] for f in glob.glob(f"{D}/*.imzml"))
lab = lambda a: re.search(r'section\d+', a).group(0)
pixel_um = 7

def tissue_mask(p):
    xs = np.array([c[0] for c in p.coordinates]) - 1; ys = np.array([c[1] for c in p.coordinates]) - 1
    # Otsu on the full-spectrum TIC fails on these sections (strong matrix signal off tissue, see the 2-2 figure);
    # the summed intensity in the lipid-rich range m/z 700-850, lightly smoothed, separates tissue well
    lip = np.zeros((ys.max() + 1, xs.max() + 1)); mzr = [np.inf, -np.inf]
    for i, (x, y) in enumerate(zip(xs, ys)):
        mz, it = p.getspectrum(i); lip[y, x] = it[(mz >= 700) & (mz <= 850)].sum()
        if len(mz): mzr = [min(mzr[0], mz[0]), max(mzr[1], mz[-1])]
    lt = np.log1p(ndi.gaussian_filter(lip, 1)); m = lt > threshold_otsu(lt[lt > 0])
    lab, _ = ndi.label(m); size = np.bincount(lab.ravel()); size[0] = 0
    m = (size >= int(np.ceil(0.1e6 / pixel_um**2)))[lab]
    lab, _ = ndi.label(~m); size = np.bincount(lab.ravel()); size[0] = 0
    m |= ((size > 0) & (size < int(np.ceil(0.01e6 / pixel_um**2))))[lab]
    return m, mzr

BIN = 0.25e-6
def call_peaks(p, mask, L0, nbin, min_frac=0.05):
    idx = [i for i, (x, y, z) in enumerate(p.coordinates) if mask[y - 1, x - 1]]
    counts = np.zeros(nbin)
    for i in idx:
        mz, _ = p.getspectrum(i)
        b = ((np.log(mz) - L0) / BIN).astype(int); b = b[(b >= 0) & (b < nbin)]
        np.add.at(counts, b, 1)
    sm = gaussian_filter1d(counts, 2); tot = np.convolve(counts, np.ones(21), mode="same")
    pk, _ = find_peaks(sm, distance=4); pk = pk[tot[pk] >= min_frac * len(idx)]
    c = np.exp(L0 + (pk + 0.5) * BIN); lo, hi = c * (1 - 2.5e-6), c * (1 + 2.5e-6)
    mid = (hi[:-1] + lo[1:]) / 2; ov = hi[:-1] > lo[1:]; hi[:-1][ov], lo[1:][ov] = mid[ov], mid[ov]
    return pd.DataFrame({"mz_estimated": c, "min": lo, "max": hi}), len(idx)

P, M = {}, {}; mzr = [np.inf, -np.inf]
for a in acqs:
    P[a] = ImzMLParser(f"{D}/{a}.imzml", ibd_file=open(f"{D}/{a}.IBD", "rb"))
    M[a], r = tissue_mask(P[a]); mzr = [min(mzr[0], r[0]), max(mzr[1], r[1])]
    print(lab(a), "pixels", len(P[a].coordinates), "tissue", M[a].sum(), "m/z", np.round(r, 1))
L0 = np.log(mzr[0] - 1); nbin = int((np.log(mzr[1] + 1) - L0) / BIN) + 1
dfl = []
for a in acqs:
    df, n = call_peaks(P[a], M[a], L0, nbin); dfl.append(df); print(lab(a), len(df), "peaks in >= 5% of", n, "tissue px")

def match_ppm(df_list, acquisitions, ppm=3):
    """Tutorial 2-4: pool peak centres, start a new feature where neighbours are > ppm apart."""
    pk = pd.concat([d.assign(acq=a)[["mz_estimated", "min", "max", "acq"]] for d, a in zip(df_list, acquisitions)])
    pk = pk.sort_values("mz_estimated").reset_index(drop=True)
    gap = np.diff(pk.mz_estimated.values) / pk.mz_estimated.values[1:] * 1e6
    pk["feature"] = np.concatenate([[0], np.cumsum(gap > ppm)])
    n_dup = pk.groupby("feature").acq.apply(lambda a: a.duplicated().any()).sum()
    span = pk.groupby("feature").mz_estimated.agg(lambda m: (m.max() - m.min()) / m.mean() * 1e6)
    print(f"{pk.feature.nunique()} features; {n_dup} contain >1 peak of one acquisition; "
          f"{(span > 2 * ppm).sum()} span more than {2 * ppm} ppm")
    return pk

pk = match_ppm(dfl, acqs, ppm=3)
iv = consensus_intervals(pk[["feature", "acq", "min", "max", "mz_estimated"]], acqs)
A = {}
for a in acqs:
    A[a] = requantify(f"{D}/{a}.imzml", iv[iv.acq == a])
    xy = A[a].obsm["spatial"]; A[a] = A[a][M[a][xy[:, 1] - 1, xy[:, 0] - 1]].copy()

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

from censor_stats import censor_summary
censor_summary(A, lab)

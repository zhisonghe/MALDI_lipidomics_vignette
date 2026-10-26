"""Preprocessing of the example data for the tutorial figures (make_figures.py, make_analysis_figures.py).

Tissue masks as in 2-2 (Otsu on the smoothed summed image of a few abundant tissue ions), a simple frequency-based stand-in for uMAIA peak
calling (features found in >= 5% of the tissue pixels of an acquisition), ppm matching at 3 ppm (2-4) and
requantification of every feature in every acquisition (2-4). Runs on a laptop without the msi-umaia environment.
"""
import os, sys
import numpy as np, pandas as pd
from pyimzml.ImzMLParser import ImzMLParser
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gapfill import consensus_intervals, requantify


# abundant tissue lipids in positive mode with DHB (the reference ions of 2-2): PC 34:1 [M+H]+, [M+Na]+, [M+K]+,
# PC 32:0 [M+H]+, PC 36:1 [M+H]+
MASK_IONS = [760.58508, 782.56702, 798.54096, 734.56943, 788.61638]


def tic_images(p, lipid_range=(700, 850), ions=MASK_IONS, ppm=5.0):
    """Full-spectrum TIC, the summed intensity in a lipid-rich m/z range, and the summed intensity of a few
    abundant tissue ions (nearest centroid within ppm of each)."""
    xs = np.array([c[0] for c in p.coordinates]) - 1; ys = np.array([c[1] for c in p.coordinates]) - 1
    tic = np.zeros((ys.max() + 1, xs.max() + 1)); lip = np.zeros_like(tic); ion = np.zeros_like(tic)
    ions = np.asarray(ions)
    for i, (x, y) in enumerate(zip(xs, ys)):
        mz, it = p.getspectrum(i)
        tic[y, x] = it.sum(); lip[y, x] = it[(mz >= lipid_range[0]) & (mz <= lipid_range[1])].sum()
        if len(mz) > 1:
            j = np.clip(np.searchsorted(mz, ions), 1, len(mz) - 1)
            k = np.where(np.abs(mz[j - 1] - ions) < np.abs(mz[j] - ions), j - 1, j)
            ok = np.abs(mz[k] - ions) / ions <= ppm * 1e-6
            ion[y, x] = it[k[ok]].sum()
    return tic, lip, ion


def tissue_mask(img, pixel_um, smooth=0):                          # tutorial 2-2
    logtic = np.log1p(ndi.gaussian_filter(img, smooth) if smooth else img)
    mask = logtic > threshold_otsu(logtic[logtic > 0])
    lb, _ = ndi.label(mask); size = np.bincount(lb.ravel()); size[0] = 0
    mask = (size >= int(np.ceil(0.1e6 / pixel_um**2)))[lb]
    lb, _ = ndi.label(~mask); size = np.bincount(lb.ravel()); size[0] = 0
    return mask | ((size > 0) & (size < int(np.ceil(0.01e6 / pixel_um**2))))[lb]


BIN = 0.25e-6
L0 = np.log(400 * 0.999); NBIN = int((np.log(1600 * 1.001) - L0) / BIN) + 1


def call_peaks(p, mask, min_frac=0.05):
    """Stand-in for uMAIA: centroid counts over tissue pixels on a 0.25-ppm grid, peaks in >= min_frac of them."""
    idx = [i for i, (x, y, _) in enumerate(p.coordinates) if mask[y - 1, x - 1]]
    counts = np.zeros(NBIN)
    for i in idx:
        mz, it = p.getspectrum(i)
        b = ((np.log(mz[it > 0]) - L0) / BIN).astype(int); b = b[(b >= 0) & (b < NBIN)]
        np.add.at(counts, b, 1)
    sm = gaussian_filter1d(counts, 2); tot = np.convolve(counts, np.ones(21), mode="same")
    pk, _ = find_peaks(sm, distance=4); pk = pk[tot[pk] >= min_frac * len(idx)]
    c = np.exp(L0 + (pk + 0.5) * BIN); lo, hi = c * (1 - 2.5e-6), c * (1 + 2.5e-6)
    mid = (hi[:-1] + lo[1:]) / 2; ov = hi[:-1] > lo[1:]; hi[:-1][ov], lo[1:][ov] = mid[ov], mid[ov]
    return pd.DataFrame({"mz_estimated": c, "min": lo, "max": hi})


def requantified(files, masks, ppm=3.0):
    """files: {acquisition: (imzml, ibd)}; masks: {acquisition: boolean raster}. Returns {acquisition: AnnData of the
    tissue pixels x matched features, raw summed intensities}, all with the same features."""
    acqs = list(files)
    P = {a: ImzMLParser(files[a][0], ibd_file=open(files[a][1], "rb")) for a in acqs}
    dfl = [call_peaks(P[a], masks[a]) for a in acqs]
    pk = pd.concat([d.assign(acq=a) for d, a in zip(dfl, acqs)]).sort_values("mz_estimated").reset_index(drop=True)
    gap = np.diff(pk.mz_estimated.values) / pk.mz_estimated.values[1:] * 1e6
    pk["feature"] = np.concatenate([[0], np.cumsum(gap > ppm)])
    pk = pk.drop_duplicates(["feature", "acq"])
    iv = consensus_intervals(pk[["feature", "acq", "min", "max", "mz_estimated"]], acqs)
    complete = iv.groupby("feature").acq.nunique() == len(acqs)
    iv = iv[iv.feature.isin(complete[complete].index)]
    A = {}
    for a in acqs:
        B = requantify(files[a][0], iv[iv.acq == a], ibd=files[a][1]); xy = B.obsm["spatial"]
        A[a] = B[masks[a][xy[:, 1] - 1, xy[:, 0] - 1]].copy()
    return A

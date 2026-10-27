"""Gap filling (tutorial 2-4) and the low-centroid diagnostic used for censoring sensitivity analysis (2-6)."""
import warnings
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pyimzml.ImzMLParser import ImzMLParser


def consensus_intervals(matches, acquisitions, pad_ppm=1.0):
    """One m/z interval per feature and acquisition.

    matches: one row per (feature, acquisition) in which the feature was called, with columns
             feature, acq, min, max, mz_estimated.
    Called: the acquisition's own interval. Not called: median centre and width of the called
    acquisitions, widened by pad_ppm and shifted by the acquisition's calibration offset.
    Overlapping intervals are cut at the midpoint between the neighbouring feature centres.
    """
    m = (matches.groupby(["feature", "acq"], as_index=False)
                .agg(min=("min", "min"), max=("max", "max"), mz_estimated=("mz_estimated", "median")))
    centre = m.groupby("feature").mz_estimated.median()
    m["dev_ppm"] = (m.mz_estimated / m.feature.map(centre) - 1) * 1e6
    offset = m.groupby("acq").dev_ppm.median()                         # calibration offset per acquisition
    half_ppm = ((m["max"] - m["min"]) / 2 / m.mz_estimated * 1e6).groupby(m.feature).median()

    full = pd.MultiIndex.from_product([centre.index, acquisitions], names=["feature", "acq"]).to_frame(index=False)
    full = full.merge(m[["feature", "acq", "min", "max"]], on=["feature", "acq"], how="left")
    full["called"] = full["min"].notna()
    mid = full.feature.map(centre) * (1 + full.acq.map(offset).fillna(0).values * 1e-6)
    w = (full.feature.map(half_ppm) + pad_ppm) * 1e-6 * mid
    full["min"] = full["min"].fillna(mid - w)
    full["max"] = full["max"].fillna(mid + w)
    full["mz"] = full.feature.map(centre)
    out = []
    for a, g in full.groupby("acq", sort=False):
        g = g.sort_values("mz").copy()                                 # order by feature centre
        lo, hi, c = g["min"].values.copy(), g["max"].values.copy(), g["mz"].values
        ov = hi[:-1] > lo[1:]
        cut = (c[:-1] + c[1:]) / 2
        hi[:-1][ov] = np.minimum(hi[:-1][ov], cut[ov])
        lo[1:][ov] = np.maximum(lo[1:][ov], cut[ov])
        g["min"], g["max"] = lo, hi
        out.append(g)
    iv = pd.concat(out, ignore_index=True)
    bad = ~(iv["min"] < iv["max"])
    if bad.any():                                                      # e.g. an interval entirely beyond its neighbour
        warnings.warn(f"{bad.sum()} intervals became empty after removing overlaps and are dropped")
        iv = iv[~bad]
    return iv


def requantify(imzml, intervals, window_da=25.0, q=0.05, ibd=None):
    """Raw summed intensity of every interval in every pixel of one acquisition (sparse AnnData).

    obs['tic_spectrum']: total intensity of each spectrum (for full-spectrum TIC normalization).
    obsm['low_centroid']: per pixel and m/z window of width window_da, the q-quantile of the centroid
    intensities (raw counts). This describes how low the retained centroids go in that part of the
    spectrum; it is a diagnostic, not a limit of detection. var['window'] maps features to windows.
    The windows are defined from the feature centres, so they are identical for all acquisitions.
    ibd: path of the binary file if it is not next to the imzML with the same name.
    """
    if intervals.empty:
        raise ValueError("no intervals for this acquisition")
    if not intervals["feature"].is_unique:
        raise ValueError("duplicate features among the intervals of one acquisition: resolve the matching first")
    iv = intervals.sort_values("min")
    acq = str(iv["acq"].iat[0])
    lo, hi = iv["min"].values, iv["max"].values
    mz0 = np.floor(iv["mz"].min() / window_da) * window_da
    nwin = int((iv["mz"].max() - mz0) // window_da) + 1
    p = ImzMLParser(imzml, ibd_file=open(ibd, "rb")) if ibd else ImzMLParser(imzml)
    n = len(p.coordinates)
    low = np.full((n, nwin), np.nan, dtype=np.float32)
    tic = np.zeros(n)
    rows, cols, vals = [], [], []
    for i in range(n):
        mz, it = p.getspectrum(i)
        keep = np.isfinite(it) & (it > 0)
        mz, it = mz[keep], it[keep]
        if len(mz) == 0:
            continue
        tic[i] = it.sum()                                              # full-spectrum TIC (2-6)
        if np.any(np.diff(mz) < 0):
            o = np.argsort(mz); mz, it = mz[o], it[o]
        cs = np.concatenate([[0.0], np.cumsum(it)])
        v = cs[np.searchsorted(mz, hi, side="right")] - cs[np.searchsorted(mz, lo, side="left")]
        nz = np.flatnonzero(v > 0)
        rows.append(np.full(len(nz), i)); cols.append(nz); vals.append(v[nz])
        win = ((mz - mz0) // window_da).astype(int)
        ok = (win >= 0) & (win < nwin)
        win, itw = win[ok], it[ok]
        o = np.lexsort((itw, win)); win, itw = win[o], itw[o]
        cnt = np.bincount(win, minlength=nwin)
        start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
        has = cnt > 0
        low[i, has] = itw[start[has] + (q * (cnt[has] - 1)).astype(int)]
    if rows:
        X = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, len(iv)), dtype=np.float32)
    else:
        X = sp.csr_matrix((n, len(iv)), dtype=np.float32)
    extra = [c for c in ("chain", "split") if c in iv]                  # from match() (scripts/matching.py), if present
    var = iv.assign(feature=iv.feature.astype(str)).set_index("feature")[["mz", "min", "max", "called"] + extra]
    var["window"] = ((var["mz"] - mz0) // window_da).astype(int).clip(0, nwin - 1)
    A = ad.AnnData(X, var=var)
    A.obs_names = [f"{acq}:px{i}" for i in range(n)]
    A.obs["acq"] = acq
    A.obs["tic_spectrum"] = tic
    A.obsm["low_centroid"] = low
    A.obsm["spatial"] = np.array([c[:2] for c in p.coordinates])      # imzML x, y (1-based)
    A.uns["windows"] = {"mz0": float(mz0), "width": float(window_da), "q": float(q)}
    return A

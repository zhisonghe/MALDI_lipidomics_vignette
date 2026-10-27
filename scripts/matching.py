"""Matching of peak-calling intervals across acquisitions by ppm clustering (tutorial 2-4).

1. Within each acquisition, peaks closer than merge_ppm are merged: closer than the scatter of one ion between
   acquisitions, they are one peak split by the peak calling.
2. The centres are corrected by each acquisition's calibration offset (median deviation from the feature centres
   of a first matching).
3. Single-linkage ppm clustering of the pooled centres gives chains of peaks.
4. A chain that holds two peaks of one acquisition, or is wider than max_span_ppm, is split into consecutive
   groups so that each group has at most one peak per acquisition and spans at most max_span_ppm. Among the valid
   splits, the one with the fewest groups is taken, and among those the one with the smallest within-group sum of
   squares (in ppm). This is exact (dynamic programming over the sorted peaks).
"""
import numpy as np, pandas as pd


def pool(df_list, acquisitions):
    """One row per peak: acq, mz_estimated, min, max, num_pixels."""
    cols = ["mz_estimated", "min", "max", "num_pixels"]
    return pd.concat([d[cols].assign(acq=a) for d, a in zip(df_list, acquisitions)], ignore_index=True)


def chains(mz, ppm):
    """Single-linkage clusters of sorted m/z values: a new cluster wherever neighbours are > ppm apart."""
    gap = np.diff(mz) / mz[1:] * 1e6
    return np.concatenate([[0], np.cumsum(gap > ppm)])


def merge_close(pk, merge_ppm):
    """Merge peaks of the same acquisition closer than merge_ppm (chains of them, too)."""
    pk = pk.sort_values(["acq", "mz_estimated"]).reset_index(drop=True)
    new = np.ones(len(pk), bool)
    same = (pk.acq.values[1:] == pk.acq.values[:-1])
    close = np.diff(pk.mz_estimated.values) / pk.mz_estimated.values[1:] * 1e6 < merge_ppm
    new[1:] = ~(same & close)
    pk["g"] = np.cumsum(new)
    w = pk.num_pixels.clip(lower=1)
    pk["wmz"] = pk.mz_estimated * w
    m = pk.assign(w=w).groupby("g").agg(acq=("acq", "first"), min=("min", "min"), max=("max", "max"),
                                        wmz=("wmz", "sum"), w=("w", "sum"), num_pixels=("num_pixels", "max"),
                                        n=("acq", "size"))
    m["mz_estimated"] = m.wmz / m.w
    return m[["acq", "mz_estimated", "min", "max", "num_pixels"]].reset_index(drop=True), int((m.n > 1).sum())


def calibration_offsets(pk, ppm):
    """Median deviation (ppm) of each acquisition from the centres of unambiguous features of a first matching."""
    s = pk.sort_values("mz_estimated")
    f = pd.Series(chains(s.mz_estimated.values, ppm), index=s.index)
    s = s.assign(f=f)
    ok = s.groupby("f").acq.transform(lambda a: not a.duplicated().any() and a.size > 1)
    s = s[ok.astype(bool)]
    dev = (s.mz_estimated / s.groupby("f").mz_estimated.transform("median") - 1) * 1e6
    return dev.groupby(s.acq).median()


def segment(u, acq, max_span):
    """Split sorted positions u (ppm) into consecutive groups with at most one peak per acquisition and a span of
    at most max_span: fewest groups, then smallest within-group sum of squares. Returns a group label per peak."""
    n = len(u)
    codes = pd.factorize(np.asarray(acq))[0]
    lo = np.empty(n, int)                       # u[lo[j]..j] is the longest valid group ending at j
    last, start = {}, 0
    for j in range(n):
        c = codes[j]
        if last.get(c, -1) >= start:
            start = last[c] + 1
        last[c] = j
        while u[j] - u[start] > max_span:
            start += 1
        lo[j] = start
    s1 = np.concatenate([[0.0], np.cumsum(u)])
    s2 = np.concatenate([[0.0], np.cumsum(u * u)])
    best = [(0, 0.0)] + [None] * n               # best[j]: (groups, sum of squares) of the first j peaks
    prev = [0] * (n + 1)
    for j in range(1, n + 1):
        for i in range(lo[j - 1], j):           # last group = peaks i..j-1
            m, t = j - i, s1[j] - s1[i]
            cand = (best[i][0] + 1, best[i][1] + (s2[j] - s2[i]) - t * t / m)
            if best[j] is None or cand < best[j]:
                best[j], prev[j] = cand, i
    lab, j, k = np.empty(n, int), n, best[n][0]
    while j > 0:
        k -= 1
        lab[prev[j]:j] = k
        j = prev[j]
    return lab


def match(df_list, acquisitions, ppm, merge_ppm=1.0, max_span_ppm=4.0, split=True):
    """Match peaks across acquisitions. split=False: no merging or splitting; a chain with two peaks of one
    acquisition keeps the one detected in most pixels (the earlier behaviour).
    Returns the peaks (columns acq, mz_estimated, min, max, num_pixels, mz_corr, chain, feature, split; at most
    one peak per feature and acquisition), the calibration offsets (ppm, per acquisition) and a dict of counts."""
    pk = pool(df_list, acquisitions)
    n_peaks = len(pk)
    pk, n_merged = merge_close(pk, merge_ppm) if split else (pk, 0)
    off = calibration_offsets(pk, ppm)
    pk["mz_corr"] = pk.mz_estimated * (1 - pk.acq.map(off).fillna(0).values * 1e-6)
    pk = pk.sort_values("mz_corr").reset_index(drop=True)
    pk["chain"] = chains(pk.mz_corr.values, ppm)
    g = pk.groupby("chain")
    dup = g.acq.transform(lambda a: a.duplicated().any()).astype(bool)
    span = g.mz_corr.transform(lambda m: (m.max() - m.min()) / m.mean() * 1e6)
    stats = {"peaks": n_peaks, "merged_close_peaks": n_merged, "chains": int(pk.chain.nunique()),
             "chains_with_two_peaks_of_one_acquisition": int(pk.chain[dup].nunique()),
             "chains_wider_than_max_span": int(pk.chain[span > max_span_ppm].nunique())}
    pk["group"] = 0
    if split:
        todo = dup | (span > max_span_ppm)
        for c, s in pk[todo].groupby("chain"):
            u = (s.mz_corr.values / s.mz_corr.values[0] - 1) * 1e6
            pk.loc[s.index, "group"] = segment(u, s.acq.values, max_span_ppm)
    else:                                       # keep the peak detected in most pixels
        pk = pk.sort_values("num_pixels", ascending=False).drop_duplicates(["chain", "acq"]).sort_values("mz_corr")
    pk["feature"] = pk.groupby(["chain", "group"], sort=True).ngroup()
    pk["split"] = pk.groupby("chain").feature.transform("nunique") > 1
    stats.update(features=int(pk.feature.nunique()), chains_split=int(pk.chain[pk.split].nunique()),
                 features_from_split_chains=int(pk.feature[pk.split].nunique()))
    return pk.drop(columns="group").reset_index(drop=True), off, stats

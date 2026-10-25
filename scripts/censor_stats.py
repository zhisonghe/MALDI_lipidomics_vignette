"""Summary of the low-centroid censoring diagnostic (tutorial 2-6) for requantified acquisitions."""
import numpy as np, pandas as pd


def censor_summary(A, label=lambda a: a):
    """A: dict acquisition -> AnnData from gapfill.requantify(), restricted to tissue pixels.
    Pixels are TIC-normalized over all features; zero-TIC pixels are ignored."""
    acqs = list(A)
    tics = {a: np.asarray(A[a].X.sum(axis=1)).ravel() for a in acqs}
    ok = {a: tics[a] > 0 for a in acqs}
    med = np.median(np.concatenate([tics[a][ok[a]] for a in acqs]))
    pf = {a: np.where(ok[a], med / np.where(ok[a], tics[a], 1), np.nan) for a in acqs}
    raw = {label(a): float(np.nanmedian(A[a].obsm["low_centroid"][ok[a]])) for a in acqs}
    per_acq = pd.DataFrame({label(a): np.nanmedian(A[a].obsm["low_centroid"][ok[a]] * pf[a][ok[a], None], axis=0)
                            for a in acqs}).T                              # acquisitions x windows
    ratio = per_acq.max() / per_acq.min()
    thr = np.nan_to_num(per_acq.max().to_numpy(), nan=0.0)
    kept = nz = 0
    for a in acqs:
        X = A[a].X[ok[a]].toarray() * pf[a][ok[a], None]
        t = thr[A[a].var["window"].to_numpy()]
        kept += ((X >= t[None, :]) & (X > 0)).sum(); nz += (X > 0).sum()
    print("median raw low-centroid intensity (counts):", {k: round(v, 1) for k, v in raw.items()})
    print("median tissue TIC:", {label(a): round(float(np.median(tics[a][ok[a]]))) for a in acqs})
    print("median low-centroid intensity on the TIC scale:", per_acq.median(axis=1).round(1).to_dict())
    print(f"max/min across acquisitions per m/z window: median {ratio.median():.2f}, "
          f"range {ratio.min():.2f}-{ratio.max():.2f}")
    print("fraction of nonzero values kept in the censored copy:", round(kept / nz, 3))
    return per_acq

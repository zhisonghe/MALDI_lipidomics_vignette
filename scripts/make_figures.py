"""Example figures of the tutorial, computed from the three public zebrafish sections of 1-5 (Zenodo 8410583).

Usage: python make_figures.py <folder with the .imzml/.IBD files> <output folder, e.g. ../figures>
Environment: msi-analysis (pyimzML, numpy, scipy, pandas, anndata, scikit-image, squidpy, matplotlib).
Peak calling uses a simple frequency-based stand-in for uMAIA (as in zebrafish_gapfill_censoring.py), so the
script runs on a laptop in a few minutes without the msi-umaia environment.
"""
import os, re, sys, glob, warnings
import numpy as np, pandas as pd, anndata as ad, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser, getionimage
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gapfill import consensus_intervals, requantify
from lipid_reference import lipid_reference, mono
warnings.filterwarnings("ignore")

D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
acqs = sorted((os.path.basename(f)[:-6] for f in glob.glob(f"{D}/*.imzml")),
              key=lambda a: int(re.search(r"section(\d+)", a).group(1)))
lab = {a: re.search(r"section(\d+)", a).group(0).replace("section", "section ") for a in acqs}
COL = dict(zip(acqs, ["#2a78d6", "#eb6834", "#1baf7a"]))       # reference palette, categorical slots 1-3
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "savefig.dpi": 150, "savefig.bbox": "tight", "image.cmap": "magma"})
PIXEL_UM = 7
P = {a: ImzMLParser(f"{D}/{a}.imzml", ibd_file=open(f"{D}/{a}.IBD", "rb")) for a in acqs}


def tic_images(p, lipid_range=(700, 850)):
    """Full-spectrum TIC and the summed intensity in a lipid-rich m/z range."""
    xs = np.array([c[0] for c in p.coordinates]) - 1; ys = np.array([c[1] for c in p.coordinates]) - 1
    tic = np.zeros((ys.max() + 1, xs.max() + 1)); lip = np.zeros_like(tic)
    for i, (x, y) in enumerate(zip(xs, ys)):
        mz, it = p.getspectrum(i)
        tic[y, x] = it.sum(); lip[y, x] = it[(mz >= lipid_range[0]) & (mz <= lipid_range[1])].sum()
    return tic, lip


def tissue_mask(img, smooth=0):                         # tutorial 2-2
    logtic = np.log1p(ndi.gaussian_filter(img, smooth) if smooth else img)
    mask = logtic > threshold_otsu(logtic[logtic > 0])
    lb, _ = ndi.label(mask); size = np.bincount(lb.ravel()); size[0] = 0
    mask = (size >= int(np.ceil(0.1e6 / PIXEL_UM**2)))[lb]
    lb, _ = ndi.label(~mask); size = np.bincount(lb.ravel()); size[0] = 0
    return mask | ((size > 0) & (size < int(np.ceil(0.01e6 / PIXEL_UM**2))))[lb]


IMG = {a: tic_images(P[a]) for a in acqs}
TIC = {a: IMG[a][0] for a in acqs}
M_tic = {a: tissue_mask(TIC[a]) for a in acqs}                         # fails on these sections (matrix signal)
M = {a: tissue_mask(IMG[a][1], smooth=1) for a in acqs}                 # lipid-range signal: used downstream

# ---- Figure 1 (2-2): TIC, an ion image and the tissue mask of one section
a0 = acqs[0]; z = P[a0].coordinates[0][2]
mz_k = 798.54096                                        # PC 34:1 [M+K]+
ion = getionimage(P[a0], mz_k, tol=mz_k * 5e-6, z=z)
fig, ax = plt.subplots(2, 2, figsize=(11, 4.6))
for a_, (im, t) in zip(ax.ravel(), [(np.log1p(TIC[a0]), "log TIC (full spectrum)"),
                                    (np.log1p(IMG[a0][1]), "log summed intensity, m/z 700-850"),
                                    (M_tic[a0], "mask: Otsu on log TIC (fails: matrix signal off tissue)"),
                                    (M[a0], "mask: Otsu on smoothed m/z 700-850 signal (used)")]):
    a_.imshow(im, cmap="gray" if im.dtype == bool else "magma", interpolation="nearest"); a_.set_title(t, fontsize=9)
    a_.set_xticks([]); a_.set_yticks([]); a_.grid(False)
fig.suptitle(f"Zebrafish {lab[a0]} (7 µm raster)", fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUT}/fig_2-2_qc_mask.png"); plt.close()

# ---- Figure 2 (2-2): mass error of reference ions against m/z and acquisition order
ION = {"H": 1.00727646662, "Na": 22.98921832, "K": 38.96315813}
refs = {"PC 32:0 +H": mono("C40H80NO8P") + ION["H"], "PC 32:0 +Na": mono("C40H80NO8P") + ION["Na"],
        "PC 32:0 +K": mono("C40H80NO8P") + ION["K"], "PC 34:1 +H": mono("C42H82NO8P") + ION["H"],
        "PC 34:1 +K": mono("C42H82NO8P") + ION["K"], "PC 36:1 +H": mono("C44H86NO8P") + ION["H"],
        "PC 36:1 +Na": mono("C44H86NO8P") + ION["Na"], "PC 36:1 +K": mono("C44H86NO8P") + ION["K"],
        "PC 34:1 [M+H-183]+": mono("C37H68O4") + ION["H"]}
rows = []
for a in acqs:
    p = P[a]
    idx = [i for i, (x, y, _) in enumerate(p.coordinates) if M[a][y - 1, x - 1]][::4]
    for i in idx:
        mz, it = p.getspectrum(i)
        if len(mz) == 0:
            continue
        for ion, m in refs.items():
            j = np.argmin(np.abs(mz - m))
            e = (mz[j] - m) / m * 1e6
            if abs(e) <= 10 and it[j] >= 500:
                rows.append((a, ion, m, i, e))
err = pd.DataFrame(rows, columns=["acq", "ion", "mz", "order", "ppm"])
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
for a in acqs:
    d = err[err.acq == a]; med = d.groupby("mz").ppm.median()
    ax[0].plot(med.index, med.values, "o-", color=COL[a], lw=2, ms=6, label=lab[a])
    o = d.assign(b=pd.qcut(d.order, 20, labels=False)).groupby("b").ppm.median()
    ax[1].plot(np.arange(20) * 5 + 2.5, o.values, "o-", color=COL[a], lw=2, ms=5, label=lab[a])
ax[0].set(xlabel="m/z of reference ion", ylabel="median mass error (ppm)", title="Error against m/z")
ax[1].set(xlabel="acquisition order (percentile of pixels)", ylabel="median mass error (ppm)",
          title="Error against acquisition order")
for a_ in ax: a_.axhline(0, color=MUTED, lw=0.8)
ax[1].legend(loc="best")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_2-2_calibration.png"); plt.close()

# ---- peak calling (stand-in), matching and requantification (2-3, 2-4)
BIN = 0.25e-6
mzr = [400.0, 1600.0]
L0 = np.log(mzr[0] * 0.999); nbin = int((np.log(mzr[1] * 1.001) - L0) / BIN) + 1


def call_peaks(p, mask, min_frac=0.05):
    idx = [i for i, (x, y, _) in enumerate(p.coordinates) if mask[y - 1, x - 1]]
    counts = np.zeros(nbin)
    for i in idx:
        mz, it = p.getspectrum(i)
        b = ((np.log(mz[it > 0]) - L0) / BIN).astype(int); b = b[(b >= 0) & (b < nbin)]
        np.add.at(counts, b, 1)
    sm = gaussian_filter1d(counts, 2); tot = np.convolve(counts, np.ones(21), mode="same")
    pk, _ = find_peaks(sm, distance=4); pk = pk[tot[pk] >= min_frac * len(idx)]
    c = np.exp(L0 + (pk + 0.5) * BIN); lo, hi = c * (1 - 2.5e-6), c * (1 + 2.5e-6)
    mid = (hi[:-1] + lo[1:]) / 2; ov = hi[:-1] > lo[1:]; hi[:-1][ov], lo[1:][ov] = mid[ov], mid[ov]
    return pd.DataFrame({"mz_estimated": c, "min": lo, "max": hi})


dfl = [call_peaks(P[a], M[a]) for a in acqs]
pk = pd.concat([d.assign(acq=a) for d, a in zip(dfl, acqs)]).sort_values("mz_estimated").reset_index(drop=True)
gap = np.diff(pk.mz_estimated.values) / pk.mz_estimated.values[1:] * 1e6
pk["feature"] = np.concatenate([[0], np.cumsum(gap > 3)])
pk = pk.drop_duplicates(["feature", "acq"])
iv = consensus_intervals(pk[["feature", "acq", "min", "max", "mz_estimated"]], acqs)
complete = iv.groupby("feature").acq.nunique() == len(acqs)
iv = iv[iv.feature.isin(complete[complete].index)]
A = {}
for a in acqs:
    B = requantify(f"{D}/{a}.imzml", iv[iv.acq == a]); xy = B.obsm["spatial"]
    A[a] = B[M[a][xy[:, 1] - 1, xy[:, 0] - 1]].copy()

# ---- Figure 3 (2-4): detection of gap-filled features after requantification
called = pd.DataFrame({a: A[a].var.called.values for a in acqs}, index=A[acqs[0]].var_names)
nc = called.sum(axis=1)
dets = np.concatenate([np.asarray((A[a][:, (~called[a] & (nc > 0)).values].X > 0).mean(axis=0)).ravel() for a in acqs])
fig, ax = plt.subplots(figsize=(6.5, 3.2))
bins = np.logspace(-4, 0, 41)
ax.hist(np.clip(dets, 1e-4, 1), bins=bins, color=COL[acqs[0]], edgecolor="white", linewidth=0.5)
ax.set_xscale("log"); ax.axvline(0.01, color=MUTED, lw=0.8, ls="--"); ax.axvline(0.05, color=MUTED, lw=0.8, ls="--")
ax.text(0.011, ax.get_ylim()[1] * 0.92, "1%", color=MUTED); ax.text(0.055, ax.get_ylim()[1] * 0.92, "5% (calling threshold)", color=MUTED)
ax.set(xlabel="fraction of tissue pixels with signal after requantification",
       ylabel="acquisition x feature gaps", title=f"{len(dets):,} gaps of features not called in every section")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_2-4_gapfill.png"); plt.close()

# ---- Figure 4 (2-5): Moran's I of one section, with images just above and below 0.4
import squidpy as sq
S = A[a0].copy()
tic = np.asarray(S.X.sum(axis=1)).ravel(); S = S[tic > 0].copy(); tic = tic[tic > 0]
S.X = (S.X.multiply(np.median(tic) / tic[:, None])).tocsr()
S = S[:, np.asarray((S.X > 0).mean(axis=0)).ravel() >= 0.05].copy()
S.obsm["spatial"] = S.obsm["spatial"].astype(float)
sq.gr.spatial_neighbors(S, coord_type="grid", n_neighs=4)
sq.gr.spatial_autocorr(S, mode="moran", genes=S.var_names.tolist(), n_perms=None, show_progress_bar=False)
mi = S.uns["moranI"]["I"].reindex(S.var_names)
detf = pd.Series(np.asarray((S.X > 0).mean(axis=0)).ravel(), index=S.var_names)
well = mi[detf.reindex(mi.index) >= 0.3]
examples = [(well - t).abs().idxmin() for t in (0.2, 0.4, 0.7)]          # well-detected features near I = 0.2, 0.4, 0.7
fig = plt.figure(figsize=(13, 3.0)); gs = fig.add_gridspec(1, 4, width_ratios=[1.3, 1, 1, 1])
ax0 = fig.add_subplot(gs[0]); ax0.hist(mi.dropna(), bins=40, color=COL[a0], edgecolor="white", linewidth=0.5)
ax0.axvline(0.4, color=MUTED, lw=0.8, ls="--"); ax0.set(xlabel="Moran's I (4 neighbours)", ylabel="features",
                                                         title=f"{lab[a0]}: {len(mi)} features detected in >= 5% of pixels")
xy = S.obsm["spatial"].astype(int)
for k, f in enumerate(examples):
    img = np.full((xy[:, 1].max(), xy[:, 0].max()), np.nan)
    img[xy[:, 1] - 1, xy[:, 0] - 1] = S[:, f].X.toarray().ravel()
    a_ = fig.add_subplot(gs[k + 1]); a_.imshow(img, vmax=np.nanpercentile(img, 99), interpolation="nearest")
    a_.set_title(f"m/z {S.var.loc[f, 'mz']:.3f}, I = {mi[f]:.2f}", fontsize=9)
    a_.set_xticks([]); a_.set_yticks([]); a_.grid(False)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_2-5_moran.png"); plt.close()

# ---- Figure 5 (2-6): low-centroid diagnostic per m/z window, raw counts and TIC scale
tics = {a: np.asarray(A[a].X.sum(axis=1)).ravel() for a in acqs}
med = np.median(np.concatenate(list(tics.values())))
w = A[acqs[0]].uns["windows"]; nwin = A[acqs[0]].obsm["low_centroid"].shape[1]
centres = w["mz0"] + w["width"] * (np.arange(nwin) + 0.5)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4), sharex=True)
for a in acqs:
    low = A[a].obsm["low_centroid"]; ok = tics[a] > 0
    ax[0].plot(centres, np.nanmedian(low[ok], axis=0), "o-", color=COL[a], lw=2, ms=4, label=lab[a])
    ax[1].plot(centres, np.nanmedian(low[ok] * (med / tics[a][ok])[:, None], axis=0), "o-", color=COL[a], lw=2, ms=4, label=lab[a])
ax[0].set(xlabel="m/z window (25 Da)", ylabel="median 5th-percentile centroid\nintensity (raw counts)", title="Raw counts")
ax[1].set(xlabel="m/z window (25 Da)", ylabel="same, after TIC normalization", title="TIC scale (the scale samples are compared on)")
ax[1].legend(loc="best")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_2-6_censoring.png"); plt.close()

# ---- Figure 6 (2-8): annotation tolerance diagnostics on the matched features
ref = lipid_reference()
ADD = {"[M+H]+": ION["H"], "[M+Na]+": ION["Na"], "[M+K]+": ION["K"]}
DECOY = {"[M+Ag]+": 106.90450, "[M+Cs]+": 132.90490}
peaks = A[acqs[0]].var["mz"].to_numpy()


def ions(adds):
    return np.sort(np.concatenate([ref.exact_mass.to_numpy() + dm for dm in adds.values()]))


def errors(th, tol):
    out = []
    for f, m in enumerate(peaks):
        lo, hi = np.searchsorted(th, m * (1 - tol * 1e-6)), np.searchsorted(th, m * (1 + tol * 1e-6))
        out += [(f, (m - t) / t * 1e6, round(t, 5)) for t in th[lo:hi]]
    return pd.DataFrame(out, columns=["f", "ppm", "ion"])


th, thd = ions(ADD), ions(DECOY)
wide, wided = errors(th, 10), errors(thd, 10)
curve = []
for tol in [0.5, 1, 1.5, 2, 3, 4, 5, 7.5, 10]:
    a_ = wide[wide.ppm.abs() <= tol]; n = a_.groupby("f").ion.nunique()
    curve.append((tol, n.size / len(peaks), (n > 1).sum() / len(peaks), wided[wided.ppm.abs() <= tol].f.nunique() / len(peaks)))
curve = pd.DataFrame(curve, columns=["ppm", "annotated", "multiple", "decoy"])
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
ax[0].hist(wide.ppm, bins=80, color=COL[acqs[0]], edgecolor="white", linewidth=0.3)
ax[0].set(xlabel="mass error of all candidate ions within 10 ppm (ppm)", ylabel="candidates",
          title="Mass error histogram")
for c, k, nm in [("annotated", 0, "annotated (>= 1 candidate)"), ("multiple", 1, "several candidate ions"), ("decoy", 2, "decoy adducts (Ag+, Cs+)")]:
    ax[1].plot(curve.ppm, curve[c], "o-", color=list(COL.values())[k], lw=2, ms=5, label=nm)
ax[1].set(xlabel="tolerance (ppm)", ylabel="fraction of features", title="Annotation gain against multiple annotation")
ax[1].legend(loc="upper left")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_2-8_tolerance.png"); plt.close()

print(f"features {len(peaks)}; gaps {len(dets)}; Moran's I median {mi.median():.2f}")
print(curve.round(3).to_string(index=False))
print("figures written to", OUT)

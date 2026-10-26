"""Example figures of the analysis part of the tutorial (3-1 to 3-5, and 3-6 when a censored copy is given).

Used by make_analysis_figures.py (public zebrafish sections) and by the organoid working copy. Follows the generic
code of section 3: feature selection, NMF + Harmony, Leiden clustering, restoration from the embedding (and its
circularity), marker lipids. Environment: msi-analysis (scanpy, squidpy, harmonypy, scikit-learn, igraph,
leidenalg, matplotlib).
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import scanpy as sc, squidpy as sq, harmonypy as hm
from sklearn.decomposition import NMF

PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#8a5cd6"]                     # reference palette, categorical slots
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "savefig.dpi": 150, "savefig.bbox": "tight", "image.cmap": "magma"})


def image(S, values):
    x, y = S.obs["x"].to_numpy(int), S.obs["y"].to_numpy(int)
    img = np.full((y.max(), x.max()), np.nan)
    img[y - 1, x - 1] = values
    return img


def noaxis(a):
    a.set_xticks([]); a.set_yticks([]); a.grid(False)


def dense(M):
    return M.toarray() if hasattr(M, "toarray") else np.asarray(M)


def run(adata, out, batch="sample", name=None, resolution=1.0, n_components=16, censored=None, condition=None):
    """adata: tissue pixels x QC'd features; X = pixel-normalized values; obs[batch] (one value per acquisition),
    obs x, y (raster, 1-based); var mz. name: Series feature -> label for plots. censored: layer name of a censored
    copy (enables the 3-6 figure, with obs[condition]). Writes fig_3-*.png to out and returns a dict of numbers."""
    num = {}
    samples = list(dict.fromkeys(adata.obs[batch]))
    col = dict(zip(samples, PAL))
    name = name if name is not None else adata.var["mz"].map(lambda m: f"m/z {m:.4f}")
    X = dense(adata.X).astype(np.float32)

    # ---- 3-1 feature selection: Moran's I per sample (2-5), the LBA heuristic, and a standardized alternative
    morans = pd.DataFrame(index=adata.var_names, columns=samples, dtype=float)
    for s in samples:
        S = adata[adata.obs[batch] == s].copy()
        S.obsm["spatial"] = S.obs[["x", "y"]].to_numpy(float)
        sq.gr.spatial_neighbors(S, coord_type="grid", n_neighs=4)
        sq.gr.spatial_autocorr(S, mode="moran", genes=S.var_names.tolist(), n_perms=None, show_progress_bar=False)
        morans[s] = S.uns["moranI"]["I"].reindex(adata.var_names)
    features = list(adata.var_names)
    df = pd.DataFrame(X, columns=features, index=adata.obs_names).assign(Section=adata.obs[batch].to_numpy())
    # the LBA heuristic as in 3-1, on the data at their own scale and after dividing all values by a constant
    frac_pass = {}
    for c in [1e-4, 1e-3, 1e-2, 1e-1, 1, 10]:
        d = df[features] * c
        sec_var = d.groupby(df["Section"]).var()
        score = sec_var.mean() - sec_var.var() / 2
        frac_pass[c] = float((score > 0).mean())
    det = pd.DataFrame({s: (X[(adata.obs[batch] == s).to_numpy()] > 0).mean(axis=0) for s in samples}, index=features)
    mean_moran = morans.mean(axis=1)
    selected = mean_moran.index[(mean_moran > 0.4) & (det.min(axis=1) >= 0.05)]   # structured and detected everywhere
    # one well-detected, structured lipid is kept out of the embedding, to be restored from it in 3-4
    target = mean_moran[selected][det.loc[selected].min(axis=1) >= 0.3].sort_values().index
    target = target[int(0.75 * (len(target) - 1))]
    selected = selected.drop(target)
    num.update(features=len(features), selected=len(selected), score_pass=frac_pass)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    ax[0].plot(list(frac_pass), list(frac_pass.values()), "o-", color=PAL[0], lw=2, ms=6)
    ax[0].set_xscale("log"); ax[0].set_ylim(-0.03, 1.03)
    ax[0].set(xlabel="all intensities multiplied by", ylabel="fraction of features with score > 0",
              title="LBA score: the cutoff depends on the scale")
    ok = det.min(axis=1) >= 0.05
    ax[1].scatter(det.min(axis=1), mean_moran, s=4, color=MUTED, alpha=0.4, label="all QC'd features")
    ax[1].scatter(det.min(axis=1)[selected], mean_moran[selected], s=6, color=PAL[1], label=f"selected ({len(selected)})")
    ax[1].axhline(0.4, color=MUTED, lw=0.8, ls="--"); ax[1].axvline(0.05, color=MUTED, lw=0.8, ls="--")
    ax[1].set_xscale("log")
    ax[1].set(xlabel="lowest detection fraction over samples", ylabel="mean Moran's I over samples",
              title="Selection on separate criteria"); ax[1].legend(loc="upper left", markerscale=3)
    plt.tight_layout(); plt.savefig(f"{out}/fig_3-1_features.png"); plt.close()

    # ---- 3-2 NMF on the selected features, projected on all pixels; Harmony on the sample
    Xsel = adata[:, selected].X
    Xsel = dense(Xsel)
    Xsel = Xsel / np.percentile(Xsel, 99, axis=0).clip(1e-12)                  # comparable feature scales for NMF
    assert Xsel.min() >= 0, "NMF needs non-negative input"
    nmf = NMF(n_components=n_components, init="nndsvda", random_state=0, max_iter=500)
    nmf.fit(Xsel)
    adata.obsm["X_nmf"] = nmf.transform(Xsel)
    ho = hm.run_harmony(adata.obsm["X_nmf"], adata.obs, [batch], random_state=0)
    Z = np.asarray(ho.Z_corr)
    n_obs, n_comp = adata.obsm["X_nmf"].shape
    if Z.shape == (n_obs, n_comp) and Z.shape != (n_comp, n_obs):
        adata.obsm["X_harmony"] = Z
    elif Z.shape == (n_comp, n_obs) and Z.shape != (n_obs, n_comp):
        adata.obsm["X_harmony"] = Z.T
    else:
        raise ValueError(f"cannot determine the orientation of Z_corr {Z.shape}")
    umaps = {}
    for rep in ["X_nmf", "X_harmony"]:
        sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=0, key_added=f"nb_{rep}")
        sc.tl.umap(adata, neighbors_key=f"nb_{rep}", random_state=0)
        umaps[rep] = adata.obsm["X_umap"].copy()
    s0 = samples[0]; S0 = adata[adata.obs[batch] == s0]
    top = np.argsort(-adata.obsm["X_nmf"].mean(axis=0))[:3]
    fig, ax = plt.subplots(1, 5, figsize=(17, 3.4), gridspec_kw={"width_ratios": [1.2, 1.2, 1, 1, 1]})
    order = np.random.default_rng(0).permutation(adata.n_obs)
    for k, (rep, t) in enumerate([("X_nmf", "NMF, before Harmony"), ("X_harmony", f"after Harmony on '{batch}'")]):
        for s in samples:
            i = order[(adata.obs[batch].to_numpy()[order] == s)]
            ax[k].scatter(umaps[rep][i, 0], umaps[rep][i, 1], s=0.3, color=col[s], alpha=0.5, rasterized=True)
            ax[k].scatter([], [], s=20, color=col[s], label=s)
        ax[k].set(title=f"UMAP: {t}", xlabel="UMAP 1", ylabel="UMAP 2"); noaxis(ax[k]); ax[k].legend(loc="best", fontsize=7)
    for k, f in enumerate(top):
        a_ = ax[k + 2]; a_.imshow(image(S0, S0.obsm["X_nmf"][:, f]), interpolation="nearest"); noaxis(a_)
        lead = selected[np.argmax(nmf.components_[f])]
        a_.set_title(f"{s0}: NMF factor {f + 1}\n(top lipid {name[lead]})", fontsize=8)
    plt.tight_layout(); plt.savefig(f"{out}/fig_3-2_embedding.png"); plt.close()
    mix = {}
    for rep in ["X_nmf", "X_harmony"]:                                       # mixing: same-sample neighbours
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=16).fit(adata.obsm[rep])
        idx = nn.kneighbors(adata.obsm[rep], return_distance=False)[:, 1:]
        lab = adata.obs[batch].to_numpy()
        mix[rep] = float((lab[idx] == lab[:, None]).mean())
    num["same_sample_neighbours"] = mix

    # ---- 3-3 Leiden on the integrated embedding
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30, random_state=0)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, key_added="leiden", random_state=0)
    cl = adata.obs["leiden"].astype(int).to_numpy()
    ncl = cl.max() + 1
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap((list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors))[:max(ncl, 2)])
    num["clusters"] = int(ncl)
    num["cluster_share"] = pd.crosstab(adata.obs["leiden"], adata.obs[batch], normalize="columns").round(3)
    fig, ax = plt.subplots(1, len(samples), figsize=(4.2 * len(samples), 3.6), squeeze=False)
    for k, s in enumerate(samples):
        S = adata[adata.obs[batch] == s]
        img = image(S, S.obs["leiden"].astype(int).to_numpy())
        ax[0, k].imshow(np.ma.masked_invalid(img), cmap=cmap, vmin=-0.5, vmax=cmap.N - 0.5, interpolation="nearest")
        ax[0, k].set_title(f"{s}: Leiden clusters", fontsize=9); noaxis(ax[0, k])
    handles = [plt.Line2D([], [], marker="s", ls="", color=cmap(c), label=str(c)) for c in range(ncl)]
    fig.legend(handles=handles, loc="center right", ncol=1 if ncl <= 16 else 2, fontsize=7, title="cluster")
    plt.tight_layout(rect=(0, 0, 0.93, 1)); plt.savefig(f"{out}/fig_3-3_clusters.png"); plt.close()

    # ---- 3-4 restoration from the embedding: train on all but the last sample, predict the last
    from sklearn.ensemble import HistGradientBoostingRegressor
    held = samples[-1]; train = (adata.obs[batch] != held).to_numpy(); test = ~train
    j = features.index(target)
    y = X[:, j]
    model = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=0)   # boosted trees, as XGBoost
    model.fit(adata.obsm["X_harmony"][train], y[train])
    pred = model.predict(adata.obsm["X_harmony"])
    r = float(np.corrcoef(pred[test], y[test])[0, 1])

    def eta2(v, g):                                                          # variance explained by the clusters
        m = pd.Series(v).groupby(g).transform("mean").to_numpy()
        return float(((m - v.mean()) ** 2).sum() / ((v - v.mean()) ** 2).sum())
    num.update(restored=str(name[target]), restored_r=r, eta2_measured=eta2(y[test], cl[test]),
               eta2_restored=eta2(pred[test], cl[test]))
    St = adata[test]
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.6))
    vmax = np.percentile(y[test], 99)
    for a_, (v, t) in zip(ax[:2], [(y[test], f"measured (held-out {held})"),
                                   (pred[test], f"restored from the embedding, r = {r:.2f}")]):
        a_.imshow(image(St, v), vmin=0, vmax=vmax, interpolation="nearest"); a_.set_title(f"{name[target]}: {t}", fontsize=8); noaxis(a_)
    ax[2].imshow(np.ma.masked_invalid(image(St, cl[test])), cmap=cmap, vmin=-0.5, vmax=cmap.N - 0.5, interpolation="nearest")
    ax[2].set_title(f"Leiden clusters; variance explained by clusters:\nmeasured {num['eta2_measured']:.2f}, "
                    f"restored {num['eta2_restored']:.2f}", fontsize=8); noaxis(ax[2])
    plt.tight_layout(); plt.savefig(f"{out}/fig_3-4_restoration.png"); plt.close()

    # ---- 3-5 marker lipids on 0-1 scaled measured values
    lo, hi = np.percentile(X, 0.5, axis=0), np.percentile(X, 99.5, axis=0)
    adata.layers["measured_scaled"] = np.clip((X - lo) / np.where(hi > lo, hi - lo, 1), 0, 1).astype(np.float32)
    B = adata.copy(); B.var_names = pd.Index(name.astype(str)).where(~pd.Index(name.astype(str)).duplicated(),
                                                                     pd.Index(name.astype(str)) + "_" + pd.Index(adata.var_names))
    sc.tl.rank_genes_groups(B, groupby="leiden", method="wilcoxon", layer="measured_scaled", use_raw=False)
    dp = sc.pl.rank_genes_groups_dotplot(B, n_genes=2, layer="measured_scaled", show=False,
                                        standard_scale="var", figsize=(max(8, 0.32 * 2 * ncl + 3), 0.32 * ncl + 2))
    plt.savefig(f"{out}/fig_3-5_markers.png"); plt.close("all")

    # ---- 3-6 (only with a censored copy): detection and intensity per sample and cluster, uncensored vs censored
    if censored is not None and condition is not None:
        conds = list(dict.fromkeys(adata.obs[condition])); c0, c1 = conds[:2]
        res = {}
        for lay in [None, censored]:
            M = X if lay is None else dense(adata.layers[lay])
            dfp = pd.DataFrame(M, columns=features).assign(cond=adata.obs[condition].to_numpy(), lz=cl)
            g = dfp.groupby(["cond", "lz"])
            detf = g[features].agg(lambda v: (v > 0).mean())
            n_px = g.size()
            logi = g[features].agg(lambda v: np.log2(v[v > 0]).mean() if (v > 0).sum() >= 20 else np.nan)
            k = detf.mul(n_px, axis=0)
            elogit = np.log((k + 0.5) / (n_px.to_numpy()[:, None] - k + 0.5))
            d_det = (elogit.xs(c1) - elogit.xs(c0)).stack()
            d_int = (logi.xs(c1) - logi.xs(c0)).stack()
            res[lay or "uncensored"] = (d_det, d_int)
        du, iu = res["uncensored"]; dc, ic = res[censored]
        both = du.index.intersection(dc.index)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].scatter(du[both], dc[both], s=3, alpha=0.3, color=PAL[0], rasterized=True)
        lim = np.nanpercentile(np.abs(np.r_[du[both], dc[both]]), 99.5)
        ax[0].plot([-lim, lim], [-lim, lim], color=MUTED, lw=0.8); ax[0].set_xlim(-lim, lim); ax[0].set_ylim(-lim, lim)
        ax[0].set(xlabel=f"uncensored: log odds of detection, {c1} - {c0}", ylabel="censored copy: same",
                  title="Detection: depends on censoring")
        bi = iu.index.intersection(ic.index); ok = iu[bi].notna() & ic[bi].notna()
        ax[1].scatter(iu[bi][ok], ic[bi][ok], s=3, alpha=0.3, color=PAL[1], rasterized=True)
        lim2 = np.nanpercentile(np.abs(np.r_[iu[bi][ok], ic[bi][ok]]), 99.5)
        ax[1].plot([-lim2, lim2], [-lim2, lim2], color=MUTED, lw=0.8); ax[1].set_xlim(-lim2, lim2); ax[1].set_ylim(-lim2, lim2)
        ax[1].set(xlabel=f"uncensored: log2 intensity where detected, {c1} - {c0}", ylabel="censored copy: same",
                  title="Intensity where detected")
        plt.tight_layout(); plt.savefig(f"{out}/fig_3-6_two_part.png"); plt.close()
        num.update(two_part_n=int(len(both)), det_median_unc=float(du[both].median()), det_median_cens=float(dc[both].median()),
                   det_sign_flip=float((np.sign(du[both]) != np.sign(dc[both])).mean()),
                   int_median_unc=float(iu[bi][ok].median()), int_median_cens=float(ic[bi][ok].median()))
    return num

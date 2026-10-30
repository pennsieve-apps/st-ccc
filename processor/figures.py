"""
Figures (PNG, 150 dpi, sample id in the title): the dot plot of the top interactions,
the source x target heatmap of interaction counts, and the maps of local scores.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm  # noqa: E402

DPI = 150
SEQUENTIAL = "Blues"        # one hue, light -> dark: magnitude and counts
DIVERGING = "RdBu_r"        # two hues around a neutral zero: local scores that can be negative
SEQUENTIAL_LOCAL = "Reds"   # local scores that are all >= 0 (mask_negatives)


def _clean_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#9a9a9a")
    ax.tick_params(colors="#444444", labelsize=8)


def dotplot(cluster: pd.DataFrame, sample_id: str, label: str, top_n: int, path: Path) -> Path:
    """Top TOP_N interactions of a sample: rows ligand -> receptor, columns source -> target;
    colour = mean ligand/receptor expression (magnitude), size = specificity (-log10 specificity_rank)."""
    top = cluster.head(top_n).copy()
    top["pair"] = top["ligand_complex"].astype(str) + " → " + top["receptor_complex"].astype(str)
    top["st"] = top["source"].astype(str) + " → " + top["target"].astype(str)
    pairs = list(dict.fromkeys(top["pair"]))                     # in rank order
    sts = list(dict.fromkeys(top["st"]))
    spec = -np.log10(np.clip(top["specificity_rank"].astype(float).values, 1e-3, 1.0))
    sizes = 15 + 55 * np.clip(spec, 0, 3)
    fig_w = max(4.5, 1.2 + 0.45 * len(sts) + 1.6)
    fig_h = max(3.5, 0.8 + 0.3 * len(pairs) + 1.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sc = ax.scatter([sts.index(s) for s in top["st"]], [pairs.index(p) for p in top["pair"]], s=sizes,
                    c=top["lr_means"].astype(float).values, cmap=SEQUENTIAL, edgecolors="#444444", linewidths=0.4,
                    zorder=3)
    ax.set_xticks(range(len(sts)))
    ax.set_xticklabels(sts, rotation=45, ha="right")
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(sts) - 0.4)
    ax.grid(True, color="#e6e6e6", linewidth=0.6, zorder=0)
    _clean_axes(ax)
    ax.set_xlabel("source → target", fontsize=9)
    ax.set_title(f"{sample_id} · top {len(top)} ligand–receptor interactions ({label})", fontsize=10)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("mean expression (lr_means)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for v, name in ((1, "0.1"), (2, "0.01"), (3, "≤0.001")):
        ax.scatter([], [], s=15 + 55 * v, c="white", edgecolors="#444444", linewidths=0.4, label=name)
    ax.legend(title="specificity rank", fontsize=7, title_fontsize=7, loc="upper left", bbox_to_anchor=(1.22, 1.0),
              frameon=False, labelspacing=1.1)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def source_target_heatmap(summary: pd.DataFrame, sample_id: str, label: str, column: str, path: Path) -> Path:
    """Source x target matrix of `column` (interaction counts), annotated."""
    labels = sorted(set(summary["source"].astype(str)) | set(summary["target"].astype(str)))
    mat = np.zeros((len(labels), len(labels)))
    for _, r in summary.iterrows():
        mat[labels.index(str(r["source"])), labels.index(str(r["target"]))] = float(r[column])
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * n + 2.5), max(3.5, 0.5 * n + 2)))
    im = ax.imshow(mat, cmap=SEQUENTIAL, vmin=0, vmax=max(1.0, mat.max()))
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels)
    ax.set_xlabel("target", fontsize=9)
    ax.set_ylabel("source", fontsize=9)
    ax.tick_params(labelsize=8, colors="#444444")
    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=7,
                    color="white" if v > 0.6 * max(1.0, mat.max()) else "#333333")
    ax.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    what = {"n_interactions_top": "interactions in the top set", "n_specific": "specific interactions"}.get(column, column)
    ax.set_title(f"{sample_id} · {what}\nper source → target ({label})", fontsize=10)
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label(what, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def _lowres_image(adata, sample_id: str):
    """(image, scalefactor) of the lowres H&E for Visium, else (None, 1)."""
    try:
        entry = adata.uns["spatial"][sample_id]
        img = entry["images"]["lowres"]
        scale = float(entry["scalefactors"]["tissue_lowres_scalef"])
        return np.asarray(img), scale
    except (KeyError, TypeError, ValueError):
        return None, 1.0


def local_map(adata, scores: np.ndarray, sample_id: str, ligand: str, receptor: str, morans: float, pval,
              local_name: str, path: Path, n_perms: int = 0) -> Path:
    """Every observation at obsm['spatial'] (y inverted), coloured by the local score: sequential from zero when
    the scores are all >= 0 (mask_negatives), else diverging around a neutral zero."""
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    img, scale = _lowres_image(adata, sample_id)
    x, y = xy[:, 0] * scale, xy[:, 1] * scale
    n = len(x)
    size = float(np.clip(30000.0 / max(n, 1), 0.3, 12.0))
    vmax = float(np.nanmax(np.abs(scores))) if np.isfinite(scores).any() else 1.0
    vmax = max(vmax, 1e-6)
    span_x, span_y = np.ptp(x) or 1.0, np.ptp(y) or 1.0
    fig_h = 5.0
    fig_w = float(np.clip(fig_h * span_x / span_y, 3.5, 9.0)) + 1.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h + 0.6))
    if img is not None:
        ax.imshow(img, origin="upper", zorder=0)
    order = np.argsort(np.abs(scores), kind="mergesort")   # strongest on top
    if np.nanmin(scores) >= 0:
        cmap, norm = SEQUENTIAL_LOCAL, Normalize(0.0, vmax)
    else:
        cmap, norm = DIVERGING, TwoSlopeNorm(0.0, -vmax, vmax)
    zero = scores[order] == 0                                # masked / no co-expression: faint grey
    ax.scatter(x[order][zero], y[order][zero], c="#bdbdbd", alpha=0.35, s=size, linewidths=0, zorder=1,
               rasterized=True)
    keep = order[~zero]
    sc = ax.scatter(x[keep], y[keep], c=scores[keep], cmap=cmap, norm=norm, s=size, linewidths=0, zorder=2,
                    rasterized=True)
    if img is None:
        ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    if pval is None or (isinstance(pval, float) and np.isnan(pval)):
        p = ""
    elif pval == 0 and n_perms > 0:          # no permutation reached the observed value
        p = f", p < {1.0 / n_perms:.2g}"
    else:
        p = f", p = {pval:.3g}"
    ax.set_title(f"{sample_id} · {ligand} → {receptor} local {local_name}\nglobal Moran's R = {morans:.3f}{p}",
                 fontsize=10)
    cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label(f"local {local_name} score", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path

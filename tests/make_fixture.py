"""
A synthetic fixture that follows the workflow contract (conventions.md):

  spatial/<sample_id>.h5ad     st-qc output: log-normalized X (CSR float32), layers["counts"],
                               obs sample_id/x/y/in_tissue/..., obsm["spatial"] in microns,
                               uns["pennsieve"] with platform xenium, kind cell, spot_diameter
  labels/labels.csv            st-celltype output: sample_id, obs_id, cluster, cell_type
  metadata/samples.csv         st-prep output (+ prep.json)

  python -m tests.make_fixture data/fixture [--layout linked|flat|both] [--species human|mouse]

Two samples (s1, s2) of ~2,500 cells each on a jittered 50 x 50 grid (20 um pitch,
1,000 x 1,000 um), 40 genes of which 26 are real ligand/receptor genes of LIANA's
consensus resource. Cell types sit in spatial patches:

  left third    Fibroblast and Tcell interleaved (a mixed stroma/immune patch)
  middle third  Tumor
  right third   Macrophage
  one 3 x 3 block of 9 cells  Rare (fewer than MIN_OBS_PER_LABEL: skipped)

The planted communication: fibroblasts express CXCL12 strongly and T cells CXCR4, and
the two types are mixed in the same patch, so Fibroblast -> Tcell CXCL12^CXCR4 is the
strongest cluster-level interaction and the pair co-localizes spatially (the only
ligand/receptor pair whose expressing cells are neighbours). Tumor cells express
CXCR4 weakly too, but no CXCL12 is near them. Other marker genes give the other
types some weaker interactions (TGFB1 -> TGFBR1/TGFBR2, SPP1 -> CD44, APOE -> LRP1, all
received by Tumor), none of them between mixed cells: each other type sits in its own patch
and meets the next only along a boundary.

--species mouse title-cases the symbols (Cxcl12, Tgfbr1_Tgfbr2 ...), as mouse data looks;
every planted pair exists in LIANA's mouseconsensus resource too.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np

# Real genes of LIANA's consensus resource (ligands, receptors, complex subunits), then filler.
LR_GENES = ["CXCL12", "CXCR4", "TGFB1", "TGFBR1", "TGFBR2", "SPP1", "CD44", "APOE", "LRP1", "COL1A1", "ITGA1",
            "ITGB1", "FN1", "ITGA5", "MIF", "CD74", "CCL5", "CCR5", "CD274", "PDCD1", "CXCL13", "CXCR5", "CCL2",
            "CCR2", "PDGFA", "PDGFRA"]
FILLER = ["ACTB", "GAPDH", "MALAT1", "RPL13", "RPS6", "EEF1A1", "TMSB4X", "FTL", "FTH1", "B2M", "KRT8", "KRT18",
          "VIM", "PTPRC"]
GENES = LR_GENES + FILLER
# Mean counts of the genes that mark a cell type (everything else sits at a low baseline).
PROGRAMS = {
    "Fibroblast": {"CXCL12": 6.0, "COL1A1": 8.0, "FN1": 6.0, "TGFB1": 2.5, "PDGFRA": 2.0, "VIM": 4.0},
    "Tcell": {"CXCR4": 5.0, "CCR5": 2.0, "PDCD1": 1.5, "PTPRC": 4.0, "CCR2": 1.0},
    "Tumor": {"CXCR4": 1.0, "CD44": 2.0, "KRT8": 8.0, "KRT18": 8.0, "CD274": 1.5, "MIF": 3.0, "ITGB1": 2.0, "ITGA1": 1.0,
              "TGFBR1": 1.5, "TGFBR2": 1.5, "LRP1": 1.0},
    "Macrophage": {"SPP1": 6.0, "APOE": 6.0, "CD74": 5.0, "CCL2": 2.0, "PTPRC": 3.0, "CCL5": 1.0},
    "Rare": {"CXCL13": 6.0, "CXCR5": 4.0},
}
BASELINE = {g: 0.08 for g in GENES}
BASELINE.update({g: 1.5 for g in ("ACTB", "GAPDH", "MALAT1", "RPL13", "RPS6", "EEF1A1", "TMSB4X", "FTL", "FTH1", "B2M")})
PITCH_UM, SIDE, CELL_DIAMETER_UM, TARGET_SUM = 20.0, 50, 10.0, 100.0


def mouse_symbol(gene: str) -> str:
    return gene[0].upper() + gene[1:].lower()


def _layout(rng: np.random.Generator, side: int = SIDE, pitch: float = PITCH_UM) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates (um) on a jittered grid, and a label per cell."""
    ii, jj = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
    rows, cols = ii.ravel(), jj.ravel()
    xy = np.column_stack([cols * pitch, rows * pitch]).astype(np.float64)
    xy += rng.uniform(-0.3 * pitch, 0.3 * pitch, xy.shape)
    labels = np.empty(len(rows), dtype=object)
    third = side // 3
    left, middle = cols < third, (cols >= third) & (cols < 2 * third)
    labels[left] = np.where(rng.random(left.sum()) < 0.5, "Fibroblast", "Tcell")
    labels[middle] = "Tumor"
    labels[~(left | middle)] = "Macrophage"
    rare = (rows >= side - 4) & (rows < side - 1) & (cols >= side - 4) & (cols < side - 1)   # 3 x 3 block
    labels[rare] = "Rare"
    return xy, labels


def _counts(rng: np.random.Generator, labels: np.ndarray, genes: list[str]) -> np.ndarray:
    n = len(labels)
    mu = np.tile(np.array([BASELINE[g] for g in genes]), (n, 1))
    for label, program in PROGRAMS.items():
        rows = np.flatnonzero(labels == label)
        for g, m in program.items():
            mu[rows, genes.index(g)] = m
    size = rng.lognormal(0.0, 0.25, n)[:, None]          # per-cell depth
    mu = mu * size
    # Negative binomial via gamma-Poisson (dispersion 0.3).
    lam = rng.gamma(1 / 0.3, mu * 0.3)
    return rng.poisson(lam).astype(np.int32)


def make_sample(sample_id: str, rng: np.random.Generator, species: str, side: int = SIDE, pitch: float = PITCH_UM):
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    genes = GENES if species == "human" else [mouse_symbol(g) for g in GENES]
    xy, labels = _layout(rng, side, pitch)
    counts = _counts(rng, labels, GENES)      # programs use the human symbols; columns renamed below
    n = len(labels)
    obs_names = [f"{sample_id}:cell_{i + 1}" for i in range(n)]
    total = counts.sum(axis=1)
    mt = np.zeros(len(genes), dtype=bool)
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sample_id] * n),
        "x": xy[:, 0], "y": xy[:, 1], "in_tissue": np.ones(n, dtype=np.int64),
        "total_counts": total.astype(np.float32),
        "n_genes_by_counts": (counts > 0).sum(axis=1).astype(np.int64),
        "pct_counts_mt": np.zeros(n, dtype=np.float32),
        "cell_area": rng.normal(80.0, 15.0, n).astype(np.float32),
        "control_counts": np.zeros(n, dtype=np.int64),
        "condition": pd.Categorical(["demo"] * n),
    }, index=obs_names)
    var = pd.DataFrame({"gene_ids": genes, "feature_types": ["Gene Expression"] * len(genes), "mt": mt},
                       index=genes)
    X = sp.csr_matrix(counts.astype(np.float32))
    # st-qc: normalize_total to TARGET_SUM, then log1p; counts kept in layers["counts"].
    scale = TARGET_SUM / np.maximum(total, 1)
    Xn = sp.diags(scale.astype(np.float32)) @ X
    Xn.data = np.log1p(Xn.data).astype(np.float32)
    a = ad.AnnData(X=sp.csr_matrix(Xn, dtype=np.float32), obs=obs, var=var,
                   layers={"counts": sp.csr_matrix(counts.astype(np.int32))})
    a.obsm["spatial"] = xy
    a.uns["pennsieve"] = {"platform": "xenium", "kind": "cell", "coordinate_unit": "micron",
                          "spot_diameter": CELL_DIAMETER_UM, "grid": "none", "source": f"{sample_id}/",
                          "prep_version": "0.1.0", "normalized": True, "target_sum": TARGET_SUM}
    return a, labels


def make_fixture(root: Path, layout: str = "both", species: str = "human", seed: int = 7, side: int = SIDE,
                 pitch: float = PITCH_UM) -> dict[str, Path]:
    """Write the fixture under root/{linked,flat}; returns the layouts' paths. side/pitch scale it up
    (e.g. --side 410 --pitch 12 gives two Xenium-sized samples of ~168 k cells, for memory measurements)."""
    root = Path(root)
    rng = np.random.default_rng(seed)
    samples, label_rows, sample_rows, prep = [], [], [], {"samples": {}}
    cluster_of = {label: f"C{i + 1:02d}" for i, label in enumerate(PROGRAMS)}
    for sid in ("s1", "s2"):
        a, labels = make_sample(sid, rng, species, side, pitch)
        samples.append((sid, a))
        for obs_id, label in zip(a.obs_names, labels):
            label_rows.append({"sample_id": sid, "obs_id": obs_id, "cluster": cluster_of[label], "cell_type": label})
        sample_rows.append({"sample_id": sid, "path": f"{sid}/", "platform": "xenium", "kind": "cell",
                            "n_obs": a.n_obs, "n_genes": a.n_vars, "n_obs_in_tissue": a.n_obs, "condition": "demo"})
        prep["samples"][sid] = {"platform": "xenium", "files": ["cell_feature_matrix.h5", "cells.parquet"],
                                "features_dropped": {}, "coordinate_unit": "micron", "warnings": []}

    layouts = ["linked", "flat"] if layout == "both" else [layout]
    out: dict[str, Path] = {}
    for name in layouts:
        base = root / name
        if base.exists():
            shutil.rmtree(base)
        sub = (lambda port: base / port) if name == "linked" else (lambda port: base)
        for port in ("spatial", "labels", "metadata"):
            sub(port).mkdir(parents=True, exist_ok=True)
        for sid, a in samples:
            a.write_h5ad(sub("spatial") / f"{sid}.h5ad")
        with open(sub("labels") / "labels.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample_id", "obs_id", "cluster", "cell_type"])
            w.writeheader()
            w.writerows(label_rows)
        with open(sub("metadata") / "samples.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sample_rows[0]))
            w.writeheader()
            w.writerows(sample_rows)
        with open(sub("metadata") / "prep.json", "w") as fh:
            json.dump(prep, fh, indent=2)
        out[name] = base
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("root", type=Path)
    p.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    p.add_argument("--species", choices=["human", "mouse"], default="human")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--side", type=int, default=SIDE, help="cells per grid side (default %(default)s)")
    p.add_argument("--pitch", type=float, default=PITCH_UM, help="grid pitch in um (default %(default)s)")
    args = p.parse_args(argv)
    paths = make_fixture(args.root, args.layout, args.species, args.seed, args.side, args.pitch)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

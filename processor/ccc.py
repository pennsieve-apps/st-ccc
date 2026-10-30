"""
The method: ligand-receptor (LR) inference with LIANA+ (Dimitrov et al. 2024).

Per sample, on log-normalized expression and a label column:

1. Cluster level: `liana.mt.rank_aggregate` scores every ligand-receptor pair of the
   resource for every source -> target label pair with several methods (CellPhoneDB,
   Connectome, log2FC, NATMI, SingleCellSignalR) and aggregates their
   ranks (RobustRankAggregate) into `magnitude_rank` (how strongly expressed) and
   `specificity_rank` (how specific to this source/target pair).
2. Spatially informed: `liana.ut.spatial_neighbors` (Gaussian kernel of bandwidth
   BANDWIDTH x spot_diameter) and `liana.mt.bivariate` with the local cosine score
   (ligand and receptor expression co-varying within a neighbourhood, per observation)
   and the global bivariate Moran's R (one number per pair, permutation p-value).

The resource (LIANA's consensus or mouseconsensus, bundled with the package) is
restricted to the pairs whose ligand and receptor genes (every complex subunit) are
in the data, so targeted panels (Xenium, MERSCOPE, CosMx) work.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("st-ccc")

RANK_CUTOFF = 0.05        # magnitude_rank / specificity_rank at or below this count as "top"
MAX_NEIGHBOURS = 100      # liana.ut.spatial_neighbors cap on edges per observation
NZ_PROP = 0.05            # liana.mt.bivariate: ligand and receptor each expressed in >= 5% of observations
LOCAL_NAME = "cosine"     # local score (LIANA's spatial tutorial default), global: bivariate Moran's R
GLOBAL_NAME = "morans"
RESOURCES = {"human": "consensus", "mouse": "mouseconsensus"}
METHOD_COLUMNS = ["lr_means", "cellphone_pvals", "expr_prod", "scaled_weight", "lr_logfc", "spec_weight", "lrscore"]
CLUSTER_COLUMNS = ["sample_id", "source", "target", "ligand_complex", "receptor_complex", "magnitude_rank",
                   "specificity_rank"] + METHOD_COLUMNS
SPATIAL_COLUMNS = ["sample_id", "ligand", "receptor", "global_morans", "global_pval", "n_obs_positive",
                   "n_obs_significant", "mean_local", "ligand_prop", "receptor_prop", "rank"]


# --- resource ---------------------------------------------------------------------------------

def detect_species(var_names) -> str:
    """'mouse' when most symbols have lowercase letters after the first (Cxcl12), else 'human' (CXCL12)."""
    symbols = [str(g) for g in var_names if re.search(r"[A-Za-z]", str(g)[1:])]
    if not symbols:
        return "human"
    mouse_like = sum(1 for g in symbols if re.search(r"[a-z]", g[1:]))
    return "mouse" if mouse_like / len(symbols) > 0.5 else "human"


def resource_name(parameter: str, var_names) -> str:
    """The resource to use: RESOURCE as given, or by species of the gene symbols when 'auto'."""
    if parameter and parameter.lower() != "auto":
        return parameter.lower()
    return RESOURCES[detect_species(var_names)]


def load_resource(name: str) -> pd.DataFrame:
    """LIANA's bundled resource (ligand, receptor; complexes as A_B). Read from the package, no network."""
    import liana as li

    resource = li.rs.select_resource(name)
    return resource.drop_duplicates().reset_index(drop=True)


def covered_resource(resource: pd.DataFrame, var_names) -> pd.DataFrame:
    """The pairs whose ligand and receptor genes (every subunit of a complex) are in the data."""
    genes = set(map(str, var_names))

    def covered(entity: str) -> bool:
        return all(part in genes for part in str(entity).split("_"))

    keep = resource["ligand"].map(covered) & resource["receptor"].map(covered)
    return resource[keep].reset_index(drop=True)


# --- expression and labels --------------------------------------------------------------------

def ensure_log_normalized(adata, target_sum: float | None = None) -> bool:
    """X is log-normalized after st-qc; if it still looks like counts, normalize + log1p. Returns True when done here."""
    pennsieve = adata.uns.get("pennsieve", {}) if isinstance(adata.uns.get("pennsieve", {}), dict) else {}
    if pennsieve.get("normalized"):
        return False
    X = adata.X
    sample = X[: min(adata.n_obs, 2000)]
    data = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
    if data.size and np.allclose(data, np.round(data)) and data.max() > 30:
        import scanpy as sc

        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
        return True
    return False


@dataclass
class LabelStats:
    n_obs: int
    n_labelled: int
    n_unlabelled: int
    kept: dict[str, int] = field(default_factory=dict)
    skipped: dict[str, int] = field(default_factory=dict)


def attach_labels(adata, labels: pd.DataFrame, sample_id: str, column: str, min_obs: int):
    """Join labels.csv onto obs (by obs_id), drop unlabelled observations and labels with fewer than
    min_obs observations. Returns (adata subset, LabelStats)."""
    rows = labels
    if "sample_id" in labels.columns:
        mine = labels["sample_id"].astype(str) == sample_id
        if mine.any():
            rows = labels[mine]
    rows = rows.drop_duplicates("obs_id", keep="first")
    mapping = pd.Series(rows[column].astype(str).values, index=rows["obs_id"].astype(str).values)
    values = adata.obs_names.map(lambda o: mapping.get(o, None))
    labelled = pd.notna(values) & (pd.Series(values).astype(str).str.strip() != "").values
    n_unlabelled = int((~labelled).sum())
    counts = pd.Series(np.asarray(values)[labelled]).value_counts()
    kept = {str(k): int(v) for k, v in counts.items() if v >= min_obs}
    skipped = {str(k): int(v) for k, v in counts.items() if v < min_obs}
    keep = labelled & pd.Series(values).isin(kept).values
    out = adata[keep].copy()
    out.obs[column] = pd.Categorical(np.asarray(values)[keep].astype(str), categories=sorted(kept))
    return out, LabelStats(adata.n_obs, int(labelled.sum()), n_unlabelled, kept, skipped)


# --- cluster-level LR -------------------------------------------------------------------------

def cluster_lr(adata, sample_id: str, column: str, resource: pd.DataFrame, expr_prop: float, n_perms: int,
               min_cells: int, seed: int, n_jobs: int) -> pd.DataFrame:
    """liana.mt.rank_aggregate for one sample; rows sorted by magnitude_rank."""
    import liana as li

    res = li.mt.rank_aggregate(adata, groupby=column, resource=resource, expr_prop=expr_prop, min_cells=min_cells,
                               n_perms=n_perms, seed=seed, n_jobs=n_jobs, use_raw=False, layer=None,
                               inplace=False, verbose=False)
    if res is None or len(res) == 0:
        return pd.DataFrame(columns=CLUSTER_COLUMNS)
    res = res.copy()
    res.insert(0, "sample_id", sample_id)
    for c in CLUSTER_COLUMNS:
        if c not in res.columns:
            res[c] = np.nan
    res = res[CLUSTER_COLUMNS]
    res = res.sort_values(["magnitude_rank", "specificity_rank", "source", "target", "ligand_complex",
                           "receptor_complex"], kind="mergesort").reset_index(drop=True)
    res["rank"] = np.arange(1, len(res) + 1)
    return res


def top_pairs(cluster: pd.DataFrame, n: int) -> list[tuple[str, str]]:
    """Distinct (ligand, receptor) pairs by their best magnitude_rank; all when n <= 0."""
    if cluster.empty:
        return []
    best = (cluster.sort_values("magnitude_rank", kind="mergesort")
            .drop_duplicates(["ligand_complex", "receptor_complex"]))
    if n > 0:
        best = best.head(n)
    return list(zip(best["ligand_complex"].astype(str), best["receptor_complex"].astype(str)))


def source_target_summary(cluster: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Per source -> target: pairs tested, pairs among the sample's top TOP_N, pairs with magnitude_rank
    and specificity_rank <= RANK_CUTOFF, mean magnitude_rank."""
    columns = ["sample_id", "source", "target", "n_interactions", "n_interactions_top", "n_specific",
               "mean_magnitude_rank", "best_magnitude_rank"]
    if cluster.empty:
        return pd.DataFrame(columns=columns)
    c = cluster.copy()
    c["is_top"] = c["rank"] <= top_n
    c["is_specific"] = (c["magnitude_rank"] <= RANK_CUTOFF) & (c["specificity_rank"] <= RANK_CUTOFF)
    g = c.groupby(["sample_id", "source", "target"], observed=True)
    out = g.agg(n_interactions=("rank", "size"), n_interactions_top=("is_top", "sum"),
                n_specific=("is_specific", "sum"), mean_magnitude_rank=("magnitude_rank", "mean"),
                best_magnitude_rank=("magnitude_rank", "min")).reset_index()
    out["n_interactions_top"] = out["n_interactions_top"].astype(int)
    out["n_specific"] = out["n_specific"].astype(int)
    return out.sort_values(["sample_id", "n_interactions_top", "n_specific", "mean_magnitude_rank"],
                           ascending=[True, False, False, True], kind="mergesort")[columns].reset_index(drop=True)


def consensus(frames: list[pd.DataFrame], top_n: int) -> pd.DataFrame:
    """Across samples: for every (source, target, ligand, receptor) the number of samples where it was
    scored (n_samples), where it was among the sample's top TOP_N (n_samples_top), mean and best ranks;
    ranked by n_samples_top, n_samples, mean magnitude_rank."""
    columns = ["rank", "source", "target", "ligand_complex", "receptor_complex", "n_samples", "n_samples_top",
               "mean_magnitude_rank", "mean_specificity_rank", "min_magnitude_rank", "mean_lr_means", "samples"]
    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return pd.DataFrame(columns=columns)
    c = pd.concat(frames, ignore_index=True)
    c["is_top"] = c["rank"] <= top_n
    key = ["source", "target", "ligand_complex", "receptor_complex"]
    g = c.groupby(key, observed=True)
    out = g.agg(n_samples=("sample_id", "nunique"), n_samples_top=("is_top", "sum"),
                mean_magnitude_rank=("magnitude_rank", "mean"), mean_specificity_rank=("specificity_rank", "mean"),
                min_magnitude_rank=("magnitude_rank", "min"), mean_lr_means=("lr_means", "mean"),
                samples=("sample_id", lambda s: ";".join(sorted(set(map(str, s)))))).reset_index()
    out["n_samples_top"] = out["n_samples_top"].astype(int)
    out = out.sort_values(["n_samples_top", "n_samples", "mean_magnitude_rank", "mean_specificity_rank"],
                          ascending=[False, False, True, True], kind="mergesort").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out[columns]


# --- spatially informed LR --------------------------------------------------------------------

def spot_diameter(adata) -> tuple[float, str]:
    """uns['pennsieve']['spot_diameter'], else the median nearest-neighbour distance (with a note)."""
    pennsieve = adata.uns.get("pennsieve", {})
    if isinstance(pennsieve, dict):
        d = pennsieve.get("spot_diameter")
        try:
            if d is not None and float(d) > 0:
                return float(d), "uns.pennsieve.spot_diameter"
        except (TypeError, ValueError):
            pass
    from scipy.spatial import cKDTree

    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    n = min(len(xy), 20000)
    idx = np.random.default_rng(0).choice(len(xy), n, replace=False) if len(xy) > n else np.arange(len(xy))
    dist, _ = cKDTree(xy).query(xy[idx], k=2)
    return float(np.median(dist[:, 1])), "median nearest-neighbour distance (no uns.pennsieve.spot_diameter)"


def spatial_lr(adata, sample_id: str, interactions: list[tuple[str, str]], bandwidth: float, cutoff: float,
               n_perms: int, seed: int):
    """liana.ut.spatial_neighbors + liana.mt.bivariate. Returns (lrdata AnnData with local scores in X and
    local p-values in layers['pvals'], rows for lr_spatial.csv sorted by global Moran's R)."""
    import liana as li

    if not interactions:
        return None, pd.DataFrame(columns=SPATIAL_COLUMNS)
    li.ut.spatial_neighbors(adata, bandwidth=bandwidth, cutoff=cutoff, max_neighbours=MAX_NEIGHBOURS,
                            kernel="gaussian", set_diag=False, spatial_key="spatial", key_added="spatial")
    try:
        lrdata = li.mt.bivariate(adata, interactions=interactions, local_name=LOCAL_NAME, global_name=[GLOBAL_NAME],
                                 n_perms=n_perms, seed=seed, mask_negatives=True, use_raw=False, nz_prop=NZ_PROP,
                                 connectivity_key="spatial_connectivities", verbose=False)
    except ValueError as exc:
        if "non-zero proportions" in str(exc):        # no pair with both genes in >= NZ_PROP of the observations
            return None, pd.DataFrame(columns=SPATIAL_COLUMNS)
        raise
    if lrdata is None or lrdata.n_vars == 0:
        return None, pd.DataFrame(columns=SPATIAL_COLUMNS)
    X = lrdata.X.toarray() if sp.issparse(lrdata.X) else np.asarray(lrdata.X)
    X = np.nan_to_num(X, nan=0.0)
    pv = lrdata.layers.get("pvals") if hasattr(lrdata, "layers") else None
    if pv is not None:
        pv = pv.toarray() if sp.issparse(pv) else np.asarray(pv)
    positive = X > 0
    significant = (positive & (pv <= 0.05)) if pv is not None else positive
    var = lrdata.var
    rows = pd.DataFrame({
        "sample_id": sample_id,
        "ligand": var["ligand"].astype(str).values, "receptor": var["receptor"].astype(str).values,
        "global_morans": var[GLOBAL_NAME].astype(float).values,
        "global_pval": var[f"{GLOBAL_NAME}_pvals"].astype(float).values if f"{GLOBAL_NAME}_pvals" in var else np.nan,
        "n_obs_positive": positive.sum(axis=0).astype(int),
        "n_obs_significant": significant.sum(axis=0).astype(int),
        "mean_local": X.mean(axis=0),
        "ligand_prop": var["ligand_props"].astype(float).values if "ligand_props" in var else np.nan,
        "receptor_prop": var["receptor_props"].astype(float).values if "receptor_props" in var else np.nan,
    }, index=var.index)
    rows = rows.sort_values(["global_morans", "n_obs_positive"], ascending=[False, False], kind="mergesort")
    rows["rank"] = np.arange(1, len(rows) + 1)
    return lrdata, rows.reset_index(drop=True)[SPATIAL_COLUMNS]

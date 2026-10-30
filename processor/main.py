"""
st-ccc: which cell types talk to which, through which ligand-receptor pairs?

Reads the per-sample AnnData files of the `spatial` port (st-qc output: log-normalized
X), the `labels` port (labels.csv from st-celltype, st-domains or cell2location) and
the `metadata` port (samples.csv), and runs LIANA+ per sample for the LABEL_COLUMN
labels. Everything is written under OUTPUT_DIR/<LABEL>/ so several runs (cell types,
domains, niches) can be linked into one folder of the report step:

  lr_cluster.csv              every scored source -> target x ligand -> receptor, per sample
  lr_consensus.csv            the same aggregated across samples
  lr_spatial.csv              global bivariate Moran's R per ligand-receptor pair, per sample
  source_target_summary.csv   interaction counts per source -> target, per sample
  summary.json, provenance.json
  <sample_id>_dotplot.png, <sample_id>_source_target.png, <sample_id>_local_<ligand>_<receptor>.png

See processor/ccc.py for the method.
"""
from __future__ import annotations

import json
import logging
import math
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from processor import threads

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-ccc")
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-ccc"
REPO_ROOT = Path(__file__).resolve().parent.parent
LABEL_RE = re.compile(r"^[a-z0-9-]+$")

DEFAULTS = {
    "LABEL": "",
    "LABEL_COLUMN": "cell_type",
    "RESOURCE": "auto",
    "N_PERMS": "100",
    "EXPR_PROP": "0.1",
    "MIN_OBS_PER_LABEL": "20",
    "TOP_N": "20",
    "SPATIAL": "true",
    "TOP_N_SPATIAL": "50",
    "BANDWIDTH": "3.0",
    "CUTOFF": "0.1",
    "MAX_PLOTS": "6",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def _bool(value: str, name: str) -> bool:
    v = value.strip().lower()
    if v in ("true", "1", "yes", "y", "on"):
        return True
    if v in ("false", "0", "no", "n", "off"):
        return False
    raise ParameterError(f"{name} must be true or false, got {value!r}")


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    p["LABEL_COLUMN"] = raw["LABEL_COLUMN"].strip()
    if not p["LABEL_COLUMN"]:
        raise ParameterError("LABEL_COLUMN must name a column of labels.csv")
    p["LABEL"] = raw["LABEL"].strip() or p["LABEL_COLUMN"].replace("_", "-").lower()
    if not LABEL_RE.match(p["LABEL"]):
        raise ParameterError(f"LABEL must match [a-z0-9-]+, got {p['LABEL']!r}")
    p["RESOURCE"] = raw["RESOURCE"].strip().lower() or "auto"
    p["SPATIAL"] = _bool(raw["SPATIAL"], "SPATIAL")
    try:
        for k in ("N_PERMS", "MIN_OBS_PER_LABEL", "TOP_N", "TOP_N_SPATIAL", "MAX_PLOTS", "SEED", "THREADS"):
            p[k] = int(raw[k])
        for k in ("EXPR_PROP", "BANDWIDTH", "CUTOFF"):
            p[k] = float(raw[k])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["TOP_N_SPATIAL"] < 0 or p["MAX_PLOTS"] < 0 or p["THREADS"] < 0:
        raise ParameterError("TOP_N_SPATIAL, MAX_PLOTS and THREADS must be >= 0")
    if p["N_PERMS"] < 1 or p["MIN_OBS_PER_LABEL"] < 1 or p["TOP_N"] < 1:
        # rank_aggregate needs the CellPhoneDB permutation p-values for its specificity rank
        raise ParameterError("N_PERMS, MIN_OBS_PER_LABEL and TOP_N must be >= 1")
    if not 0 <= p["EXPR_PROP"] <= 1 or not 0 <= p["CUTOFF"] < 1:
        raise ParameterError("EXPR_PROP must be in [0, 1] and CUTOFF in [0, 1)")
    if not p["BANDWIDTH"] > 0 or not math.isfinite(p["BANDWIDTH"]):
        raise ParameterError("BANDWIDTH must be > 0")
    return p


def _rel(path, root: Path):
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _clean_json(obj):
    if isinstance(obj, float):
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):   # numpy scalars
        try:
            return _clean_json(obj.item())
        except (ValueError, TypeError):
            return str(obj)
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def write_table(path: Path, frame, columns=None) -> None:
    frame = frame if columns is None else frame.reindex(columns=columns)
    frame.to_csv(path, index=False, float_format="%.6g")


def tool_versions() -> dict:
    from importlib.metadata import version

    out = {"python": platform.python_version()}
    for name in ("liana", "scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib", "numba"):
        try:
            out[name] = version(name)
        except Exception:   # noqa: BLE001
            out[name] = None
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def _records(frame, n: int) -> list[dict]:
    return [_clean_json(r) for r in frame.head(n).to_dict("records")]


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    n_threads = threads.configure(params["THREADS"])

    # Heavy imports after the thread environment is set.
    import anndata as ad
    import numpy as np
    import pandas as pd

    from processor import ccc as C
    from processor import figures as F
    from processor.inputs import discover

    if params["RESOURCE"] != "auto":
        import liana as li

        if params["RESOURCE"] not in [r.lower() for r in li.rs.show_resources()]:
            log.error("RESOURCE %r is not a LIANA resource (auto, %s)", params["RESOURCE"], ", ".join(li.rs.show_resources()))
            return 2

    found = discover(input_dir)
    if not found.h5ad:
        log.error("no .h5ad files found under %s", input_dir)
        return 1
    if found.labels_csv is None:
        log.error("no labels.csv found under %s", input_dir)
        return 1
    try:
        labels = pd.read_csv(found.labels_csv, dtype=str, keep_default_na=False)
    except (OSError, ValueError) as exc:
        log.error("cannot read %s: %s", found.labels_csv, exc)
        return 1
    if "obs_id" not in labels.columns:
        log.error("%s has no obs_id column", found.labels_csv)
        return 1
    if params["LABEL_COLUMN"] not in labels.columns:
        log.error("LABEL_COLUMN %r is not a column of %s (columns: %s)", params["LABEL_COLUMN"], found.labels_csv,
                  ", ".join(labels.columns))
        return 2
    sample_meta: dict[str, dict] = {}
    if found.samples_csv is not None:
        try:
            table = pd.read_csv(found.samples_csv, dtype=str, keep_default_na=False)
            if "sample_id" in table.columns:
                sample_meta = {r["sample_id"]: dict(r) for _, r in table.iterrows()}
        except (OSError, ValueError) as exc:
            log.warning("cannot read %s: %s", found.samples_csv, exc)
    # samples.csv order first, then any other .h5ad found
    order = [s for s in sample_meta if s in found.h5ad] + [s for s in found.h5ad if s not in sample_meta]
    log.info("%d samples, label column %s, label %s, %d threads", len(order), params["LABEL_COLUMN"], params["LABEL"],
             n_threads)

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)
    column = params["LABEL_COLUMN"]
    resource_used: str | None = None
    species_detected: str | None = None
    resource_full = None
    cluster_frames, st_frames, spatial_frames, statuses, figures = [], [], [], [], []
    sample_summaries: dict[str, dict] = {}
    timings = {"cluster": 0.0, "spatial": 0.0, "figures": 0.0}

    for sid in order:
        path = found.h5ad[sid]
        t = time.monotonic()
        status = {"sample_id": sid, "file": _rel(path, input_dir)}
        try:
            adata = ad.read_h5ad(path)
            if "spatial" not in adata.obsm:
                raise ValueError("no obsm['spatial']")
            adata.var_names_make_unique()
            normalized_here = C.ensure_log_normalized(adata)
            if normalized_here:
                log.warning("%s: X looked like raw counts; normalized and log1p-transformed it here", sid)
            if resource_full is None:
                species_detected = C.detect_species(adata.var_names)
                resource_used = C.resource_name(params["RESOURCE"], adata.var_names)
                resource_full = C.load_resource(resource_used)
                log.info("resource %s (%d pairs; genes look %s)", resource_used, len(resource_full), species_detected)
            resource = C.covered_resource(resource_full, adata.var_names)
            if resource.empty:
                raise ValueError(f"none of the {len(resource_full)} ligand-receptor pairs of {resource_used} "
                                 f"has both genes in the data ({adata.n_vars} genes)")
            sub, lab = C.attach_labels(adata, labels, sid, column, params["MIN_OBS_PER_LABEL"])
            if lab.n_unlabelled:
                log.warning("%s: %d of %d observations have no %s label; dropped", sid, lab.n_unlabelled, lab.n_obs, column)
            if lab.skipped:
                log.info("%s: skipped labels with < %d observations: %s", sid, params["MIN_OBS_PER_LABEL"],
                         ", ".join(f"{k} ({v})" for k, v in lab.skipped.items()))
            if len(lab.kept) < 2:
                raise ValueError(f"only {len(lab.kept)} label(s) with >= {params['MIN_OBS_PER_LABEL']} observations")

            t1 = time.monotonic()
            cluster = C.cluster_lr(sub, sid, column, resource, params["EXPR_PROP"], params["N_PERMS"],
                                   params["MIN_OBS_PER_LABEL"], params["SEED"], n_threads)
            timings["cluster"] += time.monotonic() - t1
            if cluster.empty:
                raise ValueError("no ligand-receptor pair passed the expression filter (EXPR_PROP)")
            st = C.source_target_summary(cluster, params["TOP_N"])
            cluster_frames.append(cluster)
            st_frames.append(st)
            top = cluster.head(params["TOP_N"])
            log.info("%s: %d obs, %d labels, %d LR pairs in the data, %d interactions scored in %.1fs; top: %s",
                     sid, sub.n_obs, len(lab.kept), len(resource), len(cluster), time.monotonic() - t1,
                     "; ".join(f"{r.source}→{r.target} {r.ligand_complex}^{r.receptor_complex} ({r.magnitude_rank:.2g})"
                               for r in top.head(3).itertuples()))
            entry = {
                "n_obs": int(adata.n_obs), "n_obs_used": int(sub.n_obs), "n_genes": int(adata.n_vars),
                "labels": lab.kept, "labels_skipped": lab.skipped, "n_unlabelled": lab.n_unlabelled,
                "n_lr_pairs_in_data": int(len(resource)), "n_interactions": int(len(cluster)),
                "n_interactions_specific": int(((cluster["magnitude_rank"] <= C.RANK_CUTOFF)
                                                & (cluster["specificity_rank"] <= C.RANK_CUTOFF)).sum()),
                "top_interactions": _records(top[["source", "target", "ligand_complex", "receptor_complex",
                                                  "magnitude_rank", "specificity_rank", "lr_means"]], 5),
                "strongest_source_target": _records(st[["source", "target", "n_interactions_top", "n_specific",
                                                        "n_interactions"]], 5),
            }
            for k in ("condition", "batch", "donor", "tissue"):
                if k in sample_meta.get(sid, {}):
                    entry[k] = sample_meta[sid][k]

            # Figures of the cluster-level results.
            t2 = time.monotonic()
            try:
                figures.append(F.dotplot(cluster, sid, params["LABEL"], params["TOP_N"], out / f"{sid}_dotplot.png"))
                figures.append(F.source_target_heatmap(st, sid, params["LABEL"], "n_interactions_top",
                                                       out / f"{sid}_source_target.png"))
            except Exception as exc:   # noqa: BLE001  figures are a convenience; the tables stand without them
                log.exception("%s: figures failed: %s", sid, exc)
            timings["figures"] += time.monotonic() - t2

            # Spatially informed scores.
            if params["SPATIAL"]:
                t3 = time.monotonic()
                diameter, diameter_source = C.spot_diameter(sub)
                bandwidth = params["BANDWIDTH"] * diameter
                pairs = C.top_pairs(cluster, params["TOP_N_SPATIAL"])
                spatial_error = None
                try:
                    lrdata, spatial = C.spatial_lr(sub, sid, pairs, bandwidth, params["CUTOFF"], params["N_PERMS"],
                                                   params["SEED"])
                except Exception as exc:   # noqa: BLE001  keep the cluster-level results of the sample
                    spatial_error = f"{type(exc).__name__}: {exc}"
                    log.error("%s: spatial scores failed: %s", sid, spatial_error)
                    lrdata, spatial = None, pd.DataFrame(columns=C.SPATIAL_COLUMNS)
                timings["spatial"] += time.monotonic() - t3
                spatial_frames.append(spatial)
                entry.update({
                    "spatial": {"bandwidth": bandwidth, "spot_diameter": diameter, "spot_diameter_source": diameter_source,
                                "error": spatial_error,
                                "n_pairs_tested": len(pairs), "n_pairs_scored": int(len(spatial)),
                                "n_pairs_significant": int(((spatial["global_pval"] <= 0.05)
                                                            & (spatial["global_morans"] > 0)).sum()) if len(spatial) else 0,
                                "top_pairs": _records(spatial[["ligand", "receptor", "global_morans", "global_pval",
                                                               "n_obs_positive"]], 5)},
                })
                log.info("%s: spatial: bandwidth %.3g (%s), %d pairs scored in %.1fs; top: %s", sid, bandwidth,
                         diameter_source, len(spatial), time.monotonic() - t3,
                         "; ".join(f"{r.ligand}^{r.receptor} R={r.global_morans:.3f}" for r in spatial.head(3).itertuples()))
                t2 = time.monotonic()
                if lrdata is not None and params["MAX_PLOTS"] > 0:
                    try:
                        X = lrdata.X.toarray() if hasattr(lrdata.X, "toarray") else np.asarray(lrdata.X)
                        names = list(lrdata.var_names)
                        for r in spatial.head(params["MAX_PLOTS"]).itertuples():
                            key = f"{r.ligand}^{r.receptor}"
                            if key not in names:
                                continue
                            scores = np.nan_to_num(X[:, names.index(key)].astype(float), nan=0.0)
                            safe = re.sub(r"[^A-Za-z0-9_.-]", "-", f"{r.ligand}_{r.receptor}")
                            figures.append(F.local_map(sub, scores, sid, r.ligand, r.receptor, float(r.global_morans),
                                                       float(r.global_pval) if r.global_pval == r.global_pval else None,
                                                       C.LOCAL_NAME, out / f"{sid}_local_{safe}.png",
                                                       n_perms=params["N_PERMS"]))
                    except Exception as exc:   # noqa: BLE001
                        log.exception("%s: local score maps failed: %s", sid, exc)
                timings["figures"] += time.monotonic() - t2
                del lrdata
            sample_summaries[sid] = entry
            status.update(status="ok", n_obs=int(adata.n_obs), n_obs_used=int(sub.n_obs), n_labels=len(lab.kept),
                          labels_skipped=lab.skipped, n_interactions=int(len(cluster)),
                          normalized_here=normalized_here, seconds=round(time.monotonic() - t, 2))
            del adata, sub
        except Exception as exc:   # noqa: BLE001  one bad sample is skipped, not fatal
            status.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=round(time.monotonic() - t, 2))
            log.error("%s: failed: %s", sid, exc)
        statuses.append(status)

    if not cluster_frames:
        log.error("every sample failed")
        return 1

    cluster_all = pd.concat(cluster_frames, ignore_index=True)
    st_all = pd.concat(st_frames, ignore_index=True)
    cons = C.consensus(cluster_frames, params["TOP_N"])
    write_table(out / "lr_cluster.csv", cluster_all, C.CLUSTER_COLUMNS + ["rank"])
    write_table(out / "lr_consensus.csv", cons)
    write_table(out / "source_target_summary.csv", st_all)
    if params["SPATIAL"]:
        spatial_all = pd.concat(spatial_frames, ignore_index=True) if spatial_frames else pd.DataFrame(columns=C.SPATIAL_COLUMNS)
        write_table(out / "lr_spatial.csv", spatial_all, C.SPATIAL_COLUMNS)

    # Strongest source -> target pairs over all samples: by interactions in the top sets, then specific ones.
    pairs = (st_all.groupby(["source", "target"], observed=True)
             .agg(n_samples=("sample_id", "nunique"), n_interactions_top=("n_interactions_top", "sum"),
                  n_specific=("n_specific", "sum"), mean_magnitude_rank=("mean_magnitude_rank", "mean")).reset_index()
             .sort_values(["n_interactions_top", "n_specific", "mean_magnitude_rank"], ascending=[False, False, True],
                          kind="mergesort"))
    summary = {
        "label": params["LABEL"], "label_column": column, "resource": resource_used, "species_detected": species_detected,
        "n_samples": len(order), "n_samples_ok": len(cluster_frames), "n_samples_failed": len(order) - len(cluster_frames),
        "n_interactions": int(len(cluster_all)), "n_consensus_interactions": int(len(cons)),
        "n_consensus_in_every_sample": int((cons["n_samples"] == len(cluster_frames)).sum()) if len(cons) else 0,
        "top_interactions": _records(cons, params["TOP_N"]),
        "strongest_source_target": _records(pairs, 10),
        "samples": sample_summaries,
        "spatial": None,
        "parameters": params, "method": {"local_name": C.LOCAL_NAME, "global_name": C.GLOBAL_NAME,
                                         "rank_cutoff": C.RANK_CUTOFF, "max_neighbours": C.MAX_NEIGHBOURS},
        "figures": [f.name for f in figures],
    }
    if params["SPATIAL"] and len(spatial_all):
        sp_cons = (spatial_all.groupby(["ligand", "receptor"], observed=True)
                   .agg(n_samples=("sample_id", "nunique"), mean_global_morans=("global_morans", "mean"),
                        n_samples_significant=("global_pval", lambda s: int((s <= 0.05).sum())),
                        n_obs_positive=("n_obs_positive", "sum")).reset_index()
                   .sort_values(["mean_global_morans"], ascending=False, kind="mergesort"))
        summary["spatial"] = {"n_pairs_scored": int(len(spatial_all)), "top_pairs": _records(sp_cons, params["TOP_N"])}
    write_json(out / "summary.json", summary)
    n_failed = sum(1 for s in statuses if s["status"] != "ok")
    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(),
        "parameters": params, "threads": n_threads,
        "method": {"cluster": "liana.mt.rank_aggregate (CellPhoneDB, Connectome, log2FC, NATMI, SingleCellSignalR; "
                              "RobustRankAggregate) on log-normalized X",
                   "spatial": f"liana.ut.spatial_neighbors (gaussian, bandwidth = BANDWIDTH x spot_diameter, cutoff "
                              f"CUTOFF, set_diag False, max_neighbours {C.MAX_NEIGHBOURS}) + liana.mt.bivariate "
                              f"(local {C.LOCAL_NAME}, global {C.GLOBAL_NAME}, mask_negatives, nz_prop {C.NZ_PROP})",
                   "resource": resource_used, "species_detected": species_detected, "rank_cutoff": C.RANK_CUTOFF},
        "inputs": {"spatial": [_rel(found.h5ad[s], input_dir) for s in order],
                   "labels_csv": _rel(found.labels_csv, input_dir), "samples_csv": _rel(found.samples_csv, input_dir)},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {k: round(v, 2) for k, v in timings.items()} | {"total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(),
        "samples": statuses, "n_ok": len(cluster_frames), "n_failed": n_failed,
        "outputs": sorted(p.name for p in out.iterdir()) + ["provenance.json"],
    })
    head = cons.head(3)
    log.info("done: %d samples (%d failed), %d interactions, consensus top: %s, %.1fs, peak RSS %s MB -> %s",
             len(cluster_frames), n_failed, len(cluster_all),
             "; ".join(f"{r.source}→{r.target} {r.ligand_complex}^{r.receptor_complex} ({r.n_samples_top}/{r.n_samples})"
                       for r in head.itertuples()), time.monotonic() - t0, peak_rss_mb(), out)
    return 0

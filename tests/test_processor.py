"""The synthetic fixture through the processor (both input layouts) and the method's pieces."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from processor import ccc as C
from processor import main as app
from tests.make_fixture import make_fixture

TABLES = {"lr_cluster.csv", "lr_consensus.csv", "lr_spatial.csv", "source_target_summary.csv", "summary.json",
          "provenance.json"}
FAST = {"N_PERMS": "30"}


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


@pytest.fixture(scope="session")
def mouse_root(tmp_path_factory) -> Path:
    return make_fixture(tmp_path_factory.mktemp("mouse"), "linked", species="mouse")["linked"]


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in (FAST | params).items():
        monkeypatch.setenv(name, value)
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_interaction_is_found(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch) == 0
    # Everything under OUTPUT_DIR/<LABEL>/ (LABEL defaults to LABEL_COLUMN with '_' as '-').
    assert [p.name for p in out.iterdir()] == ["cell-type"]
    d = out / "cell-type"
    files = {p.name for p in d.iterdir()}
    assert TABLES <= files
    for sid in ("s1", "s2"):
        assert {f"{sid}_dotplot.png", f"{sid}_source_target.png", f"{sid}_local_CXCL12_CXCR4.png"} <= files

    cluster = rows(d / "lr_cluster.csv")
    assert list(cluster[0])[:7] == ["sample_id", "source", "target", "ligand_complex", "receptor_complex",
                                    "magnitude_rank", "specificity_rank"]
    for sid in ("s1", "s2"):
        first = next(r for r in cluster if r["sample_id"] == sid)          # rank 1 of the sample
        assert (first["source"], first["target"], first["ligand_complex"], first["receptor_complex"]) == \
            ("Fibroblast", "Tcell", "CXCL12", "CXCR4")
        assert first["rank"] == "1"
    # The rare type (9 cells < MIN_OBS_PER_LABEL) is never a source or target.
    assert not any("Rare" in (r["source"], r["target"]) for r in cluster)

    cons = rows(d / "lr_consensus.csv")
    top = cons[0]
    assert (top["source"], top["target"], top["ligand_complex"], top["receptor_complex"]) == \
        ("Fibroblast", "Tcell", "CXCL12", "CXCR4")
    assert top["n_samples"] == "2" and top["n_samples_top"] == "2"

    spatial = rows(d / "lr_spatial.csv")
    for sid in ("s1", "s2"):
        best = next(r for r in spatial if r["sample_id"] == sid)
        assert (best["ligand"], best["receptor"]) == ("CXCL12", "CXCR4")
        assert float(best["global_morans"]) > 0.1 and float(best["global_pval"]) <= 0.05
        # The pair is spatially co-expressed only in the mixed Fibroblast/T-cell patch (left third).
        assert int(best["n_obs_positive"]) > 100

    st = rows(d / "source_target_summary.csv")
    assert {"sample_id", "source", "target", "n_interactions", "n_interactions_top", "n_specific"} <= set(st[0])

    summary = json.loads((d / "summary.json").read_text())
    assert summary["label"] == "cell-type" and summary["resource"] == "consensus"
    assert summary["species_detected"] == "human" and summary["n_samples_ok"] == 2
    assert summary["samples"]["s1"]["labels_skipped"] == {"Rare": 9}
    assert summary["spatial"]["top_pairs"][0]["ligand"] == "CXCL12"

    prov = json.loads((d / "provenance.json").read_text())
    assert prov["app"] == "st-ccc" and prov["version"] == "0.1.0"
    assert prov["tools"]["liana"] and prov["parameters"]["N_PERMS"] == 30
    assert [s["status"] for s in prov["samples"]] == ["ok", "ok"]
    assert prov["inputs"]["labels_csv"].endswith("labels.csv")


def test_mouse_symbols_use_mouseconsensus(mouse_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(mouse_root, out, monkeypatch, SPATIAL="false", LABEL="mouse") == 0
    d = out / "mouse"
    summary = json.loads((d / "summary.json").read_text())
    assert summary["resource"] == "mouseconsensus" and summary["species_detected"] == "mouse"
    cons = rows(d / "lr_consensus.csv")
    assert (cons[0]["ligand_complex"], cons[0]["receptor_complex"]) == ("Cxcl12", "Cxcr4")
    assert not (d / "lr_spatial.csv").exists()
    assert not any(p.name.startswith("s1_local_") for p in d.iterdir())


def test_cluster_column_and_explicit_resource(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, LABEL_COLUMN="cluster", RESOURCE="cellphonedb",
               SPATIAL="false") == 0
    d = out / "cluster"
    assert json.loads((d / "summary.json").read_text())["resource"] == "cellphonedb"
    cluster = rows(d / "lr_cluster.csv")
    assert {r["source"] for r in cluster} <= {"C01", "C02", "C03", "C04"}


def test_raw_counts_are_normalized(fixture_root, tmp_path, monkeypatch):
    """st-prep output (raw counts in X) is accepted and normalized on the fly."""
    import anndata as ad

    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    for f in (src / "spatial").glob("*.h5ad"):
        a = ad.read_h5ad(f)
        a.X = a.layers.pop("counts").astype(np.float32)
        del a.uns["pennsieve"]["normalized"]
        a.write_h5ad(f)
    out = tmp_path / "out"
    assert run(src, out, monkeypatch, SPATIAL="false") == 0
    prov = json.loads((out / "cell-type" / "provenance.json").read_text())
    assert all(s["normalized_here"] for s in prov["samples"])
    cons = rows(out / "cell-type" / "lr_consensus.csv")
    assert cons[0]["ligand_complex"] == "CXCL12"


def test_bad_sample_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    (src / "spatial" / "s2.h5ad").write_bytes(b"not an h5ad")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch, SPATIAL="false") == 0
    prov = json.loads((out / "cell-type" / "provenance.json").read_text())
    assert [s["status"] for s in prov["samples"]] == ["ok", "failed"]
    assert {r["sample_id"] for r in rows(out / "cell-type" / "lr_cluster.csv")} == {"s1"}


@pytest.mark.parametrize("params", [{"LABEL": "Bad_Label"}, {"LABEL_COLUMN": "nope"}, {"N_PERMS": "0"},
                                    {"SPATIAL": "maybe"}, {"RESOURCE": "nope"}, {"BANDWIDTH": "0"}])
def test_bad_parameters_exit_2(fixture_root, tmp_path, monkeypatch, params):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, **params) == 2


def test_no_input_exit_1(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_every_sample_failing_exit_1(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    labels = pd.read_csv(src / "labels" / "labels.csv")
    labels["cell_type"] = "Only"
    labels.to_csv(src / "labels" / "labels.csv", index=False)
    assert run(src, tmp_path / "out", monkeypatch, SPATIAL="false") == 1


# --- pieces ----------------------------------------------------------------------------------

def test_detect_species():
    assert C.detect_species(["CXCL12", "CXCR4", "MT-CO1", "HLA-A"]) == "human"
    assert C.detect_species(["Cxcl12", "Cxcr4", "mt-Co1", "H2-K1"]) == "mouse"
    assert C.resource_name("auto", ["Cxcl12", "Cxcr4"]) == "mouseconsensus"
    assert C.resource_name("CellPhoneDB", ["Cxcl12"]) == "cellphonedb"


def test_covered_resource_requires_every_subunit():
    res = pd.DataFrame({"ligand": ["TGFB1", "TGFB1", "CXCL12"], "receptor": ["TGFBR1_TGFBR2", "TGFBR1_ACVR1", "CXCR4"]})
    kept = C.covered_resource(res, ["TGFB1", "TGFBR1", "TGFBR2", "CXCL12"])
    assert list(kept["receptor"]) == ["TGFBR1_TGFBR2"]


def test_consensus_counts_samples():
    def frame(sid, ranks):
        return pd.DataFrame({"sample_id": sid, "source": "A", "target": "B", "ligand_complex": ["L1", "L2"],
                             "receptor_complex": ["R1", "R2"], "magnitude_rank": ranks, "specificity_rank": [0.1, 0.2],
                             "lr_means": [1.0, 2.0], "rank": [1, 2]})
    cons = C.consensus([frame("a", [0.01, 0.5]), frame("b", [0.02, 0.4])], top_n=1)
    assert list(cons["ligand_complex"]) == ["L1", "L2"]
    assert list(cons["n_samples_top"]) == [2, 0] and list(cons["n_samples"]) == [2, 2]
    assert cons.loc[0, "mean_magnitude_rank"] == pytest.approx(0.015)


def test_spot_diameter_fallback():
    import anndata as ad

    a = ad.AnnData(np.zeros((4, 1), dtype=np.float32))
    a.obsm["spatial"] = np.array([[0, 0], [10, 0], [0, 10], [10, 10]], dtype=float)
    d, source = C.spot_diameter(a)
    assert d == pytest.approx(10.0) and "nearest" in source
    a.uns["pennsieve"] = {"spot_diameter": 12.5}
    assert C.spot_diameter(a) == (12.5, "uns.pennsieve.spot_diameter")

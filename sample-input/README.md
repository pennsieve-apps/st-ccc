# sample-input

This app reads the output of `st-qc` (`spatial/<sample_id>.h5ad`, log-normalized), a
`labels.csv` from `st-celltype` (`cell_type`, `cluster`), `st-domains` (`domain`),
`cell2location` (`dominant_type`) or `st-niches` (`niche`), and `st-prep`'s
`metadata/samples.csv` (+ `prep.json`).

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker (no network) over data/fixture/linked/ -> data/output/cell-type/
```

Two Xenium-like samples (`s1`, `s2`, 2,500 cells each on a jittered 20 µm grid, coordinates in
microns, 40 genes of which 26 are ligand/receptor genes of LIANA's consensus resource). The left
third mixes fibroblasts and T cells, the middle third is tumour, the right third macrophages, and
a 3 × 3 block of a `Rare` type (9 cells) is skipped by `MIN_OBS_PER_LABEL=20`. Planted:
fibroblasts express **CXCL12** and T cells **CXCR4**, in the same patch, so Fibroblast → Tcell
CXCL12 → CXCR4 is rank 1 in both samples and the consensus, and the spatially most co-expressed
pair (global Moran's R ≈ 0.2, p < 0.01). Other marker genes give weaker interactions received
by the tumour (TGFB1 → TGFBR1_TGFBR2, SPP1 → CD44, APOE → LRP1) between cells that only meet
along patch boundaries. `make fixture-mouse` writes the same with mouse symbols (`Cxcl12`), which
selects `mouseconsensus`.

`python -m tests.make_fixture DIR --layout linked --side 410 --pitch 12` writes two
Xenium-sized samples (~168 k cells each) for memory measurements.

## Real data: 10x Xenium FFPE Human Breast Cancer

10x Genomics, *FFPE Human Breast Cancer* (Xenium Onboard Analysis 1.0.1), replicates 1 and 2,
CC BY 4.0 (Janesick et al. 2023, https://doi.org/10.1038/s41467-023-43458-x):
https://www.10xgenomics.com/datasets/ffpe-human-breast-with-custom-add-on-panel-1-standard

1. Run `st-prep` on the two replicate folders (`cell_feature_matrix.h5` + `cells.parquet`) with
   a `samples.csv` (`sample_id, path`), optionally `MAX_OBS_PER_SAMPLE=20000` for a quick run.
2. Run `st-qc` on its output, then `st-celltype` on the qc'd files.
3. Run this app with `spatial/` = st-qc's output, `labels/` = st-celltype's output and
   `metadata/` = st-prep's metadata:

```sh
docker build --platform linux/amd64 -t st-ccc .
docker run --rm --platform linux/amd64 --memory 6g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output -e LABEL_COLUMN=cell_type \
  -v "$PWD/data/xenium:/data/input:ro" -v "$PWD/data/output:/data/output" st-ccc
```

where `data/xenium/spatial/`, `data/xenium/labels/` and `data/xenium/metadata/` hold (or link
to) those outputs. For Visium data, link `st-domains` (`LABEL_COLUMN=domain`) or `cell2location`
(`LABEL_COLUMN=dominant_type`) as `labels/`; the spatial kernel then works in spot diameters
(`uns["pennsieve"]["spot_diameter"]`, full-resolution pixels).

# st-ccc

Which cell types talk to which, and through which ligand–receptor pairs? This app infers
**cell–cell communication** in spatial transcriptomics data with LIANA+ (Dimitrov et al. 2024):
for every pair of labelled groups (a *source* that expresses a ligand and a *target* that
expresses its receptor) it ranks every known ligand–receptor pair, and then asks whether the
strongest pairs are actually **co-expressed by neighbouring observations** in the tissue.

For a biologist: you get, per sample and across samples, a ranked list of interactions such as
"fibroblasts → T cells via CXCL12 → CXCR4", a dot plot of the top interactions, a source × target
heat map of how many top interactions each pair of cell types shares, and, for the pairs with the
strongest spatial signal, maps of the section showing **where** ligand and receptor are
co-expressed. Expression-based inference says that the two cell types *could* talk through a pair;
the spatial score adds that the cells expressing the ligand and the receptor are close enough to.

For a computational analyst: per sample, `liana.mt.rank_aggregate` on log-normalized expression
(CellPhoneDB, Connectome, log2FC, NATMI and SingleCellSignalR scores combined by robust rank
aggregation into `magnitude_rank` and `specificity_rank`) on LIANA's consensus resource; a
consensus table over samples; and `liana.mt.bivariate` on a Gaussian-kernel spatial graph
(`liana.ut.spatial_neighbors`) for the top pairs: global bivariate Moran's R and a local cosine
score per observation, both with permutation p-values.

It is a step of the **C. Single-cell Spatial Profiling** workflow (Xenium, MERSCOPE, CosMx):
`st-prep → st-qc → st-celltype → {st-niches, st-ccc} → st-report`, and works as well on Visium
spots with `st-domains` (`domain`) or `cell2location` (`dominant_type`) labels. It can run more
than once (on cell types, clusters, domains, niches) and several runs can link into one report
input, so **everything is written under `OUTPUT_DIR/<LABEL>/`**.

## Method

1. **Inputs.** Samples are processed one at a time, in `samples.csv` order. `X` is expected to be
   log-normalized (st-qc); a file that still holds integer counts (st-prep output) is normalized
   (`normalize_total` + `log1p`) here, with a warning. Labels are joined from `labels.csv` by
   `obs_id` (rows of the sample); observations without a label are dropped and counted, and labels
   with fewer than `MIN_OBS_PER_LABEL` observations in the sample are skipped for that sample
   (listed in `summary.json`). `Unassigned` (st-celltype) and `Mixed` (cell2location) are kept as
   labels of their own; read their rows with that in mind. A sample needs at least two labels.
2. **Resource.** `RESOURCE=auto` picks LIANA's `consensus` resource (human) or
   `mouseconsensus` when most gene symbols look like mouse symbols (`Cxcl12`, not `CXCL12`).
   It is restricted to the pairs whose ligand and receptor genes (every subunit of a complex such
   as `TGFBR1_TGFBR2`) are in the data, which matters for targeted panels: Xenium's 313-gene breast
   panel covers only 23 of the 4,620 consensus interactions (a whole-transcriptome Visium sample,
   3,544). The resource ships inside the
   `liana` package; nothing is downloaded.
3. **Cluster-level inference** (`liana.mt.rank_aggregate`, `expr_prop=EXPR_PROP`,
   `min_cells=MIN_OBS_PER_LABEL`, `n_perms=N_PERMS`, `seed=SEED`). For every source → target
   label pair, a ligand–receptor pair is scored only if the ligand is expressed in at least
   `EXPR_PROP` of the source's observations and the receptor in `EXPR_PROP` of the target's (every
   subunit). Five methods score it; their ranks are aggregated (RobustRankAggregate, Kolde et al.
   2012) into `magnitude_rank` (how strongly the pair is expressed; from CellPhoneDB's `lr_means`,
   Connectome's and NATMI's `expr_prod` and SingleCellSignalR's `lrscore`) and `specificity_rank`
   (how specific it is to this source/target pair; from CellPhoneDB's permutation p-value,
   Connectome's `scaled_weight`, log2FC's `lr_logfc` and NATMI's `spec_weight`). Lower is better for both; they behave like
   p-values of a rank aggregation, not like calibrated p-values. Rows of a sample are ordered by
   `magnitude_rank` then `specificity_rank`; `rank` is the position.
4. **Consensus across samples.** For every (source, target, ligand, receptor): in how many samples
   it was scored (`n_samples`), in how many it was among that sample's top `TOP_N`
   (`n_samples_top`), mean and best ranks. Ordered by `n_samples_top`, `n_samples`, mean
   magnitude rank.
5. **Spatially informed scores** (`SPATIAL=true`). On the labelled observations of the sample:
   - `liana.ut.spatial_neighbors` with a Gaussian kernel `w = exp(−d² / (2·h²))`, bandwidth
     `h = BANDWIDTH × spot_diameter` (`uns["pennsieve"]["spot_diameter"]`, in coordinate units:
     ~10 µm for cells, the spot diameter in full-resolution pixels for Visium; the median
     nearest-neighbour distance when it is missing), weights below `CUTOFF` dropped (with the
     defaults the neighbourhood reaches 2.15 h ≈ 6.4 cell diameters, ~65 µm for Xenium), at most
     100 neighbours per observation, no self-edge. This kernel graph is LIANA+'s own and replaces
     the squidpy graph of conventions.md "Spatial neighbours" (a weighted kernel is what the
     bivariate scores are defined on).
   - The top `TOP_N_SPATIAL` distinct ligand–receptor pairs of the sample (by best
     `magnitude_rank`; `0` = all scored pairs) go to `liana.mt.bivariate` with
     `mask_negatives=True` and `nz_prop=0.05` (ligand and receptor each expressed in ≥ 5 % of the
     observations). Complexes use the minimum expression of their subunits.
   - **Global**: bivariate Moran's R between ligand and receptor over the graph: > 0 when
     observations with a high ligand have neighbours with a high receptor (and vice versa), 0 when
     unrelated in space, < 0 when they avoid each other. Permutation p-value from `N_PERMS`
     shuffles of the observations (the fraction of permutations reaching the observed value; 0
     means below 1/`N_PERMS`).
   - **Local**: the **cosine similarity** of ligand and receptor expression within each
     observation's kernel-weighted neighbourhood,
     `Σⱼ wᵢⱼ xⱼ yⱼ / √(Σⱼ wᵢⱼ xⱼ² · Σⱼ wᵢⱼ yⱼ²)` (LIANA+'s default local metric and the one used
     in its spatial tutorials; it is bounded in [0, 1] for non-negative expression, cheap, and
     does not need centring, unlike local Moran's or Pearson, which are unstable on sparse
     targeted-panel counts). With `mask_negatives` the score is set to 0 where the
     kernel-weighted z-scores of ligand and receptor are not both positive (only "high–high"
     co-expression counts), and each observation gets a permutation p-value (`N_PERMS` shuffles).
   - Pairs are ordered by global Moran's R. Maps of the local scores are drawn for the
     `MAX_PLOTS` pairs with the highest R.
6. **Source–target summary.** Per sample and source → target: interactions scored, among the top
   `TOP_N`, and "specific" (`magnitude_rank` and `specificity_rank` both ≤ 0.05).

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `spatial` | input (subfolder, required) | `<sample_id>.h5ad` | st-qc output (the AnnData contract: log-normalized `X`, `obsm["spatial"]`, `uns["pennsieve"]["spot_diameter"]`). st-prep output (counts) is normalized here. |
| `labels` | input (subfolder, required) | `labels.csv` | `sample_id, obs_id` + label columns, from st-celltype, st-domains, cell2location or st-niches; `LABEL_COLUMN` picks the column. |
| `metadata` | input (subfolder, required) | `samples.csv`, `prep.json` | st-prep's metadata; sets the sample order and copies `condition`/`batch`/`donor`/`tissue` into `summary.json`. Without it, every `.h5ad` found is a sample. |
| `ccc` | output | `<LABEL>/…` (below) | The only output, so it is the whole of `OUTPUT_DIR`; the app writes into `OUTPUT_DIR/<LABEL>/`. |

Without links the same files are found by name anywhere under `INPUT_DIR` (symlinks followed).

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `LABEL` | string | (empty) | Name of the run and of the output folder, `[a-z0-9-]+`. Empty: `LABEL_COLUMN` with `_` as `-` (`cell_type` → `cell-type`). Anything else fails the run (exit 2). |
| `LABEL_COLUMN` | string | `cell_type` | Column of `labels.csv` whose values are the communicating groups (`cell_type`, `cluster`, `domain`, `dominant_type`, `niche`). |
| `RESOURCE` | string | `auto` | `auto` = `consensus` (human) or `mouseconsensus` (mouse symbols); or any resource LIANA ships (`cellphonedb`, `cellchatdb`, `connectomedb2020`, ...). |
| `N_PERMS` | integer | `100` | Permutations for the CellPhoneDB p-values (label shuffles; feed the specificity rank) and for the global and local spatial p-values (observation shuffles). ≥ 1. |
| `EXPR_PROP` | number | `0.1` | Minimum fraction of a label's observations expressing the ligand (source) / receptor (target), every subunit. |
| `MIN_OBS_PER_LABEL` | integer | `20` | Labels with fewer observations in a sample are skipped there (also LIANA's `min_cells`). |
| `TOP_N` | integer | `20` | Interactions per sample counted as "top" (dot plot, `n_interactions_top`, `n_samples_top`). |
| `SPATIAL` | boolean | `true` | Compute the spatially informed scores, `lr_spatial.csv` and the local maps. |
| `TOP_N_SPATIAL` | integer | `50` | Distinct ligand–receptor pairs per sample (by best magnitude rank) that get spatial scores; `0` = all scored pairs (slower). |
| `BANDWIDTH` | number | `3.0` | Gaussian kernel bandwidth in spot / cell diameters. |
| `CUTOFF` | number | `0.1` | Kernel weights below this are dropped (sets the neighbourhood radius, 2.15 × bandwidth at 0.1). |
| `MAX_PLOTS` | integer | `6` | Local score maps per sample. |
| `SEED` | integer | `42` | Seed of the permutations. |
| `THREADS` | integer | `0` | Threads for BLAS / numba / LIANA (`0` = the container's CPU quota, at most 8). |

## Outputs (`OUTPUT_DIR/<LABEL>/`)

| File | Contents |
|---|---|
| `lr_cluster.csv` | Every scored interaction, per sample: `sample_id, source, target, ligand_complex, receptor_complex, magnitude_rank, specificity_rank, lr_means, cellphone_pvals, expr_prod, scaled_weight, lr_logfc, spec_weight, lrscore, rank`. |
| `lr_consensus.csv` | Across samples: `rank, source, target, ligand_complex, receptor_complex, n_samples, n_samples_top, mean_magnitude_rank, mean_specificity_rank, min_magnitude_rank, mean_lr_means, samples` (`;`-separated). |
| `lr_spatial.csv` | (`SPATIAL=true`) per sample and pair: `sample_id, ligand, receptor, global_morans, global_pval, n_obs_positive` (local score > 0), `n_obs_significant` (local score > 0 and local p ≤ 0.05), `mean_local, ligand_prop, receptor_prop, rank`. |
| `source_target_summary.csv` | `sample_id, source, target, n_interactions, n_interactions_top, n_specific, mean_magnitude_rank, best_magnitude_rank`. |
| `summary.json` | `label`, `label_column`, `resource`, `species_detected`, sample counts, `top_interactions` (consensus head), `strongest_source_target`, per sample (`labels`, `labels_skipped`, `n_unlabelled`, `n_lr_pairs_in_data`, `n_interactions`, `n_interactions_specific`, top interactions, spatial bandwidth and top pairs), `spatial.top_pairs` (mean R over samples), `parameters`, `figures`. |
| `<sample_id>_dotplot.png` | Top `TOP_N` interactions: rows ligand → receptor, columns source → target, colour = mean expression (`lr_means`), size = specificity. |
| `<sample_id>_source_target.png` | Source × target heat map of the number of top-`TOP_N` interactions. |
| `<sample_id>_local_<ligand>_<receptor>.png` | Each observation at its position (y down), coloured by the local cosine score, with global R and p in the title; on the lowres H&E for Visium. |
| `provenance.json` | App and tool versions, parameters, threads, method, inputs seen, per-sample status (observations used, labels skipped, seconds, error), timings, peak RSS, outputs. |

A sample that fails is logged and skipped (a failure of the spatial part alone keeps the sample's
cluster-level results); the run fails with exit 1 when no input is found or every sample failed,
and with exit 2 on bad parameters or a `LABEL_COLUMN` missing from `labels.csv`.

## Results

### Synthetic fixture

Two samples of 2,491 labelled cells (`tests/make_fixture.py`, see `sample-input/README.md`), a
planted Fibroblast → T cell CXCL12 → CXCR4 interaction in a mixed patch: rank 1 in both samples
(`magnitude_rank` 0.03 / 0.02) and in the consensus (top in 2/2 samples), and the spatially
strongest pair (global R 0.20 / 0.22, p < 0.01; R 0.21 at 168 k cells per sample). The `Rare`
label (9 cells) is skipped; mouse symbols select `mouseconsensus` and find `Cxcl12 → Cxcr4`.

### Xenium FFPE human breast cancer

RESULTS_XENIUM

## Runtime

RUNTIME_SECTION

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python     # pytest on the synthetic fixture, both input layouts, no Docker
make run PYTHON=.venv/bin/python      # Docker (no network, 3 GB) over data/fixture/linked -> data/output/cell-type/
```

By hand:

```sh
docker build --platform linux/amd64 -t st-ccc .
docker run --rm --platform linux/amd64 --memory 6g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output -e LABEL_COLUMN=cell_type \
  -v "$PWD/data/fixture/linked:/data/input:ro" -v "$PWD/data/output:/data/output" st-ccc
```

## Limitations

- **Inference, not measurement.** Ligand and receptor mRNA in the right cells is necessary, not
  sufficient, for signalling: protein levels, secretion, processing, co-receptors and the
  downstream response are not observed. Treat the output as ranked hypotheses.
- **Resource and panel coverage.** Only interactions in the resource can be found, and on a
  targeted panel only those whose every gene is on the panel; many pathways are invisible to a
  300-gene panel. Interactions that are not protein ligand → receptor (metabolites, ECM
  mechanics, contact via adhesion complexes not in the resource) are not covered.
- **Labels set the resolution.** Communication is between the label groups you provide. Mis-typed
  cells, segmentation spill-over between touching cells (common in imaging-based data: a T cell
  next to a fibroblast picks up fibroblast transcripts) and `Unassigned` groups all produce
  interactions; spill-over in particular inflates the spatial scores of pairs whose genes mark
  neighbouring cell types.
- **Ranks are not p-values.** `magnitude_rank` and `specificity_rank` are aggregate ranks; there is
  no multiple-testing correction across thousands of source/target/pair combinations. The global
  and local spatial p-values are per pair and per observation, uncorrected, and limited to
  1/`N_PERMS`.
- **The spatial scores are label-free.** Moran's R and the local score say that ligand and
  receptor are co-expressed nearby, not which cell types express them; read them next to the
  cluster-level table. They are computed on labelled observations only.
- **Bandwidth matters.** Paracrine ranges differ by molecule; one kernel (default ~65 µm reach for
  cells, ~6 spots for Visium) is a compromise. Spot-level data (Visium) mixes several cells per
  spot, so autocrine and paracrine signals cannot be told apart.
- **Across samples** the consensus counts samples; it does not test conditions against each other
  (no differential communication).

## Licences

Apache-2.0 for this app (see `LICENSE`). LIANA+ is BSD-3-Clause; its ligand–receptor resources
ship inside the package. Every bundled component and its licence: `THIRD_PARTY.md`. The image has
no GPL code.

## Citations

LIANA+, the framework this app runs (rank aggregation, spatial bivariate scores, resources):

> Dimitrov D, Schäfer PSL, Farr E, Rodriguez-Mier P, Lobentanzer S, Badia-i-Mompel P, Dugourd A,
> Tanevski J, Ramirez Flores RO, Saez-Rodriguez J. **LIANA+ provides an all-in-one framework for
> cell–cell communication inference.** Nat Cell Biol 26(9):1613–1622 (2024).
> https://doi.org/10.1038/s41556-024-01469-w, PMID 39223377.

LIANA, the method comparison that introduced the consensus resource and the rank aggregate:

> Dimitrov D, Türei D, Garrido-Rodriguez M, Burmedi PL, Nagai JS, Boys C, Ramirez Flores RO,
> Kim H, Szalai B, Costa IG, Valdeolivas A, Dugourd A, Saez-Rodriguez J. **Comparison of methods
> and resources for cell-cell communication inference from single-cell RNA-Seq data.** Nat Commun
> 13(1):3224 (2022). https://doi.org/10.1038/s41467-022-30755-0, PMID 35680885.

Robust rank aggregation:

> Kolde R, Laur S, Adler P, Vilo J. **Robust rank aggregation for gene list integration and
> meta-analysis.** Bioinformatics 28(4):573–580 (2012). https://doi.org/10.1093/bioinformatics/btr709,
> PMID 22247279.

The methods aggregated:

> Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R. **CellPhoneDB: inferring cell–cell
> communication from combined expression of multi-subunit ligand–receptor complexes.** Nat Protoc
> 15(4):1484–1506 (2020). https://doi.org/10.1038/s41596-020-0292-x, PMID 32103204.

> Raredon MSB, Yang J, Garritano J, Wang M, Kushnir D, Schupp JC, Adams TS, Greaney AM, Leiby KL,
> Kaminski N, Kluger Y, Levchenko A, Niklason LE. **Computation and visualization of cell–cell
> signaling topologies in single-cell systems data using Connectome.** Sci Rep 12(1):4187 (2022).
> https://doi.org/10.1038/s41598-022-07959-x, PMID 35264704.

> Hou R, Denisenko E, Ong HT, Ramilowski JA, Forrest ARR. **Predicting cell-to-cell communication
> networks using NATMI.** Nat Commun 11(1):5011 (2020). https://doi.org/10.1038/s41467-020-18873-z,
> PMID 33024107.

> Cabello-Aguilar S, Alame M, Kon-Sun-Tack F, Fau C, Lacroix M, Colinge J. **SingleCellSignalR:
> inference of intercellular networks from single-cell transcriptomics.** Nucleic Acids Res
> 48(10):e55 (2020). https://doi.org/10.1093/nar/gkaa183, PMID 32196115.

Moran's statistic (the global bivariate R generalises it):

> Moran PAP. **Notes on continuous stochastic phenomena.** Biometrika 37(1-2):17–23 (1950).
> https://doi.org/10.1093/biomet/37.1-2.17.

Scanpy and anndata (normalization, data structures):

> Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data analysis.**
> Genome Biol 19(1):15 (2018). https://doi.org/10.1186/s13059-017-1382-0, PMID 29409532.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** J Open Source Softw 9(101):4371 (2024). https://doi.org/10.21105/joss.04371.

Demo data:

> Janesick A, Shelansky R, Gottscho AD, Wagner F, Williams SR, Rouault M, Beliakoff G,
> Morrison CA, Oliveira MF, Sicherman JT, Kohlway A, Abousoud J, Drennon TY, Mohabbat SH,
> 10x Development Teams, Taylor SEB. **High resolution mapping of the tumor microenvironment using
> integrated single-cell, spatial and in situ analysis.** Nat Commun 14(1):8353 (2023).
> https://doi.org/10.1038/s41467-023-43458-x, PMID 38114474. Data: 10x Genomics, *FFPE Human
> Breast Cancer* (Xenium), CC BY 4.0.

Every DOI was checked on Crossref (title, authors, journal, year, volume, pages) and every PMID on
PubMed (esearch by DOI, then esummary: title, author list, year); the Moran and anndata papers have
no PMID. `CITATION.cff` lists all of them with their full author lists.

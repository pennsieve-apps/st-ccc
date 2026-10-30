# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the versions pinned in `requirements.txt` (the complete
`pip freeze` of the image); their licences apply to the image. Nothing in the image is GPL or
LGPL (this image has no igraph / leidenalg), so it can be redistributed under these permissive
licences. Licences were read from each installed package's metadata and, for LIANA+, from its
`LICENSE` file.

| Component | What | Licence | Source |
|---|---|---|---|
| LIANA+ (`liana`) 1.10.0 | `rank_aggregate` (CellPhoneDB, Connectome, log2FC, NATMI, SingleCellSignalR scores + robust rank aggregation), `spatial_neighbors`, `bivariate` (global Moran's R, local cosine); ships the ligand-receptor resources (below). | BSD-3-Clause (Copyright (c) 2025, Daniel Dimitrov) | https://github.com/saezlab/liana-py |
| Scanpy 1.12.5 | Normalization when X still holds counts; used by LIANA+. | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata 0.13.4, mudata 0.4.1 | Reading `.h5ad`; LIANA+'s data structures. | BSD-3-Clause | https://github.com/scverse/anndata |
| NumPy 2.5.4, SciPy 1.18.1, pandas 2.3.3, h5py 3.16.0 | Arrays, sparse matrices, KD-tree, tables, HDF5. | BSD-3-Clause (NumPy also bundles 0BSD / MIT / Zlib / CC0-1.0 parts) | PyPI |
| Matplotlib 3.11.2 | PNG figures. | Matplotlib licence (PSF-based, BSD-compatible) | https://matplotlib.org |
| numba 0.68.1, llvmlite 0.50.0 | JIT kernels used by LIANA+ and Scanpy. | BSD-2-Clause; llvmlite bundles LLVM (Apache-2.0 WITH LLVM-exception) | https://numba.pydata.org |
| scikit-learn 1.9.1, statsmodels 0.15.0, umap-learn 0.5.12, pynndescent 0.7.0, networkx 3.7, seaborn 0.13.2, joblib, threadpoolctl, patsy, cloudpickle | Dependencies of Scanpy / LIANA+. | BSD-3-Clause / BSD-2-Clause | PyPI |
| plotnine 0.15.8, mizani 0.14.6 | Dependencies of LIANA+ (plotting; not used for our figures). | MIT / BSD-3-Clause | PyPI |
| docrep 0.3.2 | Docstring templating, a LIANA+ dependency (installed from its sdist; pure Python). | Apache-2.0 | https://github.com/Chilipp/docrep |
| tqdm, session-info2, legacy-api-wrap, fast-array-utils, certifi | Dependencies. | MPL-2.0 (tqdm: MPL-2.0 AND MIT), unmodified | PyPI |
| zarr, numcodecs, formulaic, narwhals, natsort, pydantic(-core, -settings), annotated-types, typing-inspection, donfig, interface_meta, array-api-compat, fonttools, pyparsing, PyYAML, six, pytz, pillow, msgspec, python-dotenv, wrapt, contourpy, cycler, kiwisolver, packaging, python-dateutil, typing_extensions, tzdata, google-crc32c, scverse-misc | Dependencies of the above. | MIT / BSD / MIT-CMU (Pillow) / Apache-2.0 or BSD (packaging, python-dateutil) / PSF-2.0 (typing_extensions) / Apache-2.0 (tzdata, google-crc32c) | PyPI |
| python:3.12.11-slim-bookworm | Base image (Debian 12 + CPython 3.12), pinned by digest. | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Ligand-receptor resources (data, inside the `liana` package)

`RESOURCE=auto` uses LIANA's **consensus** resource (4,620 human interactions, complexes as
`A_B`) or **mouseconsensus** (3,989 mouse interactions) from `liana/resource/omni_resource.csv`,
which is part of the BSD-3-Clause `liana` wheel. Nothing is downloaded at run time (the image is
built and tested with `--network none`). The consensus resource was assembled by the LIANA
authors through OmniPath from the curated interactions of CellPhoneDB, CellChatDB, ICELLNET,
connectomeDB2020 and CellTalkDB (Dimitrov et al. 2022, cited in the README). The other resources
selectable with `RESOURCE` come from the same file. If you redistribute results for commercial
use, check the terms of the original databases as well.

## Not bundled

- **Demo data** (not in the image or the repository): 10x Genomics, *FFPE Human Breast Cancer*
  (Xenium, Rep1 and Rep2), CC BY 4.0 (Janesick et al. 2023, cited in the README), processed by
  st-prep, st-qc and st-celltype.

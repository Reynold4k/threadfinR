<p align="center">
  <img src="https://raw.githubusercontent.com/Reynold4k/Threadfin/main/docs/_static/threadfin_overview.png" width="720" alt="Threadfin turns paired single-cell RNA and receptor sequencing into a clone embedding, one point per clone">
</p>

<h1 align="center">Threadfin for R</h1>

<p align="center"><b>Which B-cell and T-cell clones behave alike, and what explains it: paired single-cell RNA and receptor data in, clone-level answers out</b></p>

<p align="center">
  <a href="https://github.com/Reynold4k/threadfinR/actions/workflows/R-CMD-check.yaml"><img src="https://img.shields.io/github/actions/workflow/status/Reynold4k/threadfinR/R-CMD-check.yaml?branch=main&style=flat-square&label=R-CMD-check" alt="R-CMD-check"></a>
  <a href="https://pypi.org/project/threadfin/"><img src="https://img.shields.io/pypi/v/threadfin?style=flat-square&color=2a78d6&label=threadfin%20(Python)" alt="threadfin on PyPI"></a>
  <a href="https://github.com/Reynold4k/threadfinR/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/licence-MIT-2a78d6?style=flat-square" alt="MIT licence"></a>
</p>

<p align="center"><sub>The R interface to the Python package <a href="https://github.com/Reynold4k/Threadfin">Threadfin</a>. Every analysis runs in Threadfin itself through reticulate, so the numbers are those of the Python package.</sub></p>

## Get started

```r
install.packages("remotes")
remotes::install_github("Reynold4k/threadfinR")
threadfin::install_threadfin()     # once: a Python environment "r-threadfin" with Threadfin and Dandelion from PyPI
```

Restart R. You need, for the same cells, the raw gene counts, the receptor contigs and the donor and sample of
each cell. One call gives the answers:

```r
library(threadfin)

res <- run_threadfin(
  "b_cells.h5ad",                           # raw counts, one row per cell
  bcr = "filtered_contig_annotations.csv",  # Cell Ranger contigs, or an AIRR .tsv
  donor_key = "donor",                      # the individual each cell came from
  sample_key = "sample",                    # the sample or library each cell came from
  n_jobs = 8
)
summary(res)                                # the plain-language report
plot_overview(res)                          # clone map, coherence, stability, programme composition
res$clones                                  # one row per clone with two or more cells
```

A Seurat object or a SingleCellExperiment works the same way, with the contigs as a data frame, and the results go
back to it:

```r
res <- run_threadfin(seurat_obj, bcr = read.csv("filtered_contig_annotations.csv"),
                     donor_key = "donor", sample_key = "sample", test = c("isotype", "antigen_binding"))
seurat_obj$clone_programme <- res$cells$clone_programme[match(colnames(seurat_obj), res$cells$cell)]
```

No data at hand: `simulate_repertoire()` makes a repertoire with a known answer, and `vignette("threadfin")` runs
everything on it.

## What it answers

The cells of one clone are spread over states: germinal-centre B cells, ASCs (antibody-secreting cells) and memory
B cells. Two clones can have the same average state and still differ in how
their cells are spread. Threadfin describes every clone by that distribution, compares clones rather than cells,
and tests with clones as the unit, so large clones cannot dominate and P values are not inflated.

| Question | In the result | Step by step |
|---|---|---|
| Are cells of one clone more alike than cells of different clones? | `res$coherence` | `clonal_coherence()` |
| How similar are clones to one another? | `res$clones` (`x`, `y`), `plot_clone_map()` | `embed_clones()` |
| Which clones behave alike? | `res$programmes`, `res$composition` | `find_programmes()` |
| Does a label (antigen binding, isotype, sort gate) explain how clones differ? | `res$profile_tests`, `res$tests` | `profile_association()`, `association_test()` |
| Do clones keep their state over time or across tissues? | `res$memory` (with `time_key`) | `clonal_memory()` |
| How much less mutated is a clone's ASC output than its own germinal-centre cells? | | `export_lag()` |
| Which genes are shared within clones, or mark the clone landscape? | | `gene_heritability()`, `programme_markers()`, `landscape_genes()` |

## What you get

| Function | Runs | Returns |
|---|---|---|
| `run_threadfin()` | `threadfin.run` | `summary`, `coherence`, `clones`, `programmes`, `composition`, `profile_tests`, `tests`, `memory`, `cells`, `methods`, as data frames and text |
| `as_anndata()` | | the AnnData Threadfin works on, from an `.h5ad` file, a Seurat object or a SingleCellExperiment |
| `read_contigs()`, `define_clones_dandelion()` | `threadfin.read_contigs`, `define_clones_dandelion` | contigs with AIRR names; one clone id per cell |
| `prepare_embedding()`, `clone_profiles()`, `clonal_coherence()`, `embed_clones()`, `find_programmes()`, ... | the steps of `threadfin.pp` and `threadfin.tl` | data frames and lists |
| `plot_overview()`, `plot_clone_map()`, `plot_threadfin()`, ... | `ThreadfinResult.plot`, `threadfin.pl` | the figure on the R device, in a file (PDF, PNG, SVG) or as a raster |
| `save_threadfin()`, `read_threadfin()` | `to_csv`, `write_h5ad`, `ThreadfinResult.from_adata` | the clone table, report, Methods paragraph and AnnData; the result again |

T cells work the same way (pass the TCR contigs as `bcr`); clones defined elsewhere go in a `clone_id` column of
the cell metadata, without `bcr`. `py_object(res)` is the Python result, for anything not shown in R.

## The whole Python package

`tf` is Threadfin itself, under the name its documentation uses (`import threadfin as tf`), so every function of
the [Python API](https://github.com/Reynold4k/Threadfin/blob/main/docs/API.md) is one `$` away, and whole numbers
reach integer arguments as integers:

```r
adata <- tf$sim$simulate_repertoire(scenario = "bifurcation", random_state = 0)   # a repertoire with a known answer
tf$tl$map_quality(features, xy)             # score any clone map: clone profiles and their coordinates
tf$tl$lineage_heritability(y, forest)       # do close relatives on a lineage tree share their state?
```

## Python environment

`install_threadfin()` makes a virtual environment (conda when no Python 3.10 or newer is found) and installs
`threadfin[dandelion]` (4.6 or later) from PyPI with pip. To use a Python you manage yourself, set
`RETICULATE_PYTHON` to it before the first threadfin call. `threadfin_available()` says whether Threadfin loads.

Clone ids are the same in every session. reticulate starts Python with R's environment, so unless `PYTHONHASHSEED=0`
was set before Python started (for example in `.Renviron`), Threadfin calls the clones in a separate Python process
started with seed 0 (about 2 s more), as it does in Python; the clone ids are identical either way.

## Learn more

* `vignette("threadfin")`: the analysis from R on a simulated repertoire, no download needed.
* [Threadfin documentation](https://github.com/Reynold4k/Threadfin/blob/main/docs/README.md) and
  [tutorials](https://github.com/Reynold4k/Threadfin/tree/main/docs/tutorials): the methods, every function and
  the analyses in depth.

## How to cite

Please cite Threadfin as in its [CITATION.cff](https://github.com/Reynold4k/Threadfin/blob/main/CITATION.cff), and
the tools it builds on: Dandelion (Suo *et al.*, *Nat. Biotechnol.* 2024), Scanpy (Wolf *et al.*, *Genome Biol.*
2018) and UMAP (McInnes *et al.*, *J. Open Source Softw.* 2018).

## Licence

MIT, see [LICENSE.md](LICENSE.md).

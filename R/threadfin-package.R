#' threadfin: clone-level analysis of paired single-cell RNA and receptor data from R
#'
#' An R interface to the Python package [Threadfin](https://github.com/Reynold4k/Threadfin): each B-cell or T-cell
#' clone is described by the distribution of its cells' transcriptional states, so that clones, not cells, can be
#' compared, mapped and tested. The analyses run in Threadfin itself through reticulate, so their results are those
#' of the Python package.
#'
#' * Set up once: [install_threadfin()]; check with [threadfin_available()].
#' * One call: [run_threadfin()] (an `.h5ad` file, a Seurat object or a SingleCellExperiment, with Cell Ranger or
#'   AIRR contigs), then [plot_overview()], [save_threadfin()] and [read_threadfin()].
#' * Step by step: [as_anndata()], [define_clones_dandelion()], [prepare_embedding()], [clone_profiles()],
#'   [clonal_coherence()], [embed_clones()], [find_programmes()], [profile_association()], [association_test()],
#'   [clonal_memory()], [export_lag()], [gene_heritability()], [programme_markers()].
#' * Figures: [plot_threadfin()] and its [shortcuts][threadfin_plots].
#' * Everything else: [tf], the Python package itself.
#'
#' @keywords internal
#' @importFrom Matrix sparseMatrix
"_PACKAGE"

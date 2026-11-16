#' Gene expression for Threadfin: an AnnData from a file, a Seurat object or a SingleCellExperiment
#'
#' Threadfin works on an AnnData of raw counts, one row per cell. [run_threadfin()] converts R objects on its own;
#' call `as_anndata()` to work step by step ([clone_profiles()], [clonal_coherence()], ...), because the steps
#' write their results into the AnnData, as in Python. R objects are converted directly, without SeuratDisk or
#' zellkonverter: the raw counts become a sparse `X` (cells x genes), the cell metadata `obs`, and each reduction
#' an `obsm` embedding.
#'
#' @param x The cells' gene expression:
#'   * an `.h5ad` file, read with `anndata.read_h5ad`;
#'   * a Seurat object (v4 or v5): the raw counts of `assay` (layer `layer`; split v5 layers are joined), the
#'     metadata and the reductions;
#'   * a SingleCellExperiment (or another SummarizedExperiment): the `assay` counts, `colData()` and
#'     `reducedDims()`;
#'   * a matrix of raw counts, genes x cells (a `dgCMatrix` or a base matrix), with the cell metadata in `obs`;
#'   * a Python AnnData, returned as it is.
#' @param assay Seurat assay (default: the default assay) or SummarizedExperiment assay (default `"counts"`) with
#'   the raw counts.
#' @param layer Seurat layer (v5) or slot (v4) with the raw counts.
#' @param reductions Which reductions (Seurat) or `reducedDims` (SingleCellExperiment) to carry over: `TRUE` (all),
#'   `FALSE` (none) or their names. Each becomes `obsm["X_<name>"]` (names that start with `X_` are kept), so a
#'   Seurat `"pca"` is `basis = "X_pca"`.
#' @param obs For a matrix: the cell metadata, one row per cell (row names are cell names, or rows are in the
#'   order of the columns).
#' @return A Python AnnData. Counts that are whole numbers below 2^24 are stored as float32 (exactly, as scanpy
#'   reads 10x files), others as float64.
#' @seealso [run_threadfin()], [simulate_repertoire()]
#' @examplesIf threadfin_available() && requireNamespace("SeuratObject", quietly = TRUE)
#' counts <- Matrix::rsparsematrix(50, 20, 0.3, rand.x = function(n) rpois(n, 3) + 1)
#' dimnames(counts) <- list(paste0("Gene", 1:50), paste0("cell", 1:20))
#' so <- SeuratObject::CreateSeuratObject(counts)
#' so$donor <- rep(c("d1", "d2"), each = 10)
#' adata <- as_anndata(so)
#' adata
#' @export
as_anndata <- function(x, assay = NULL, layer = "counts", reductions = TRUE, obs = NULL) {
  .require_threadfin()
  if (is.character(x)) {
    if (length(x) != 1L) stop("x must be one .h5ad file", call. = FALSE)
    path <- path.expand(x)
    if (!file.exists(path)) stop("no such file: ", x, call. = FALSE)
    if (!grepl("\\.h5ad$", path, ignore.case = TRUE)) {
      stop("a file must be an .h5ad file (read other formats into a Seurat object or a SingleCellExperiment)",
        call. = FALSE
      )
    }
    return(.glue(TRUE)$read_h5ad(normalizePath(path)))
  }
  if (.is_py(x)) {
    return(x)
  }
  parts <- if (inherits(x, "Seurat")) {
    .seurat_parts(x, assay = assay, layer = layer, reductions = reductions)
  } else if (methods::is(x, "SummarizedExperiment")) {
    .se_parts(x, assay = assay, reductions = reductions)
  } else if (is.matrix(x) || methods::is(x, "Matrix")) {
    list(counts = x, obs = obs, obsm = NULL)
  } else {
    stop("x must be an .h5ad file, a Seurat object, a SingleCellExperiment, a counts matrix or a Python AnnData ",
      "(got ", class(x)[1], ")",
      call. = FALSE
    )
  }
  .anndata(parts$counts, obs = parts$obs, obsm = parts$obsm)
}

#' Receptor contigs as Threadfin reads them
#'
#' Runs `threadfin.read_contigs()`: the contig-level receptor table (one row per BCR or TCR chain) from Cell
#' Ranger's `filtered_contig_annotations.csv`, an AIRR rearrangement `.tsv` or a data frame, with AIRR column names
#' (`cell_id`, `v_call`, `j_call`, `junction`, ...) and only the contigs Cell Ranger calls cells and high-confidence.
#' Use it to prepare contigs for [define_clones_dandelion()] (Tutorial 2, section 6).
#'
#' @param bcr A file path, a data frame (read from one of these files; it goes to Python as pandas reads the
#'   file), or a Python object (a pandas table, a scirpy AnnData or MuData).
#' @param obs_names Cell names of the expression data: when a table has both a 10x `barcode` and an AIRR `cell_id`
#'   column, the one that matches more of them names the cells.
#' @return A data frame, one row per contig.
#' @seealso [define_clones_dandelion()], [run_threadfin()]
#' @examplesIf threadfin_available()
#' # contigs <- read_contigs("filtered_contig_annotations.csv")
#' # contigs <- contigs[contigs$cell_id %in% colnames(seurat_obj), ]
#' @export
read_contigs <- function(bcr, obs_names = NULL) {
  .frame(.call_py("read_contigs", list(.bcr_arg(bcr)), .kwargs(list(obs_names = .names_arg(obs_names)))))
}

#' A simulated repertoire with a known answer
#'
#' Runs `threadfin.sim.simulate_repertoire()`: an AnnData with ground-truth clones, programmes and cell states, for
#' trying Threadfin without data (Tutorial 1). It already has clone ids (`obs$clone_id`) and an embedding
#' (`obsm["X_pca"]`), so run it with `run_threadfin(adata, donor_key = "donor", sample_key = "context",
#' basis = "X_pca")`.
#'
#' @param ... Arguments of `simulate_repertoire()`, e.g. `n_clones` (3000), `n_donors` (4), `n_programmes` (4),
#'   `scenario` (`"default"`, `"null"` or `"bifurcation"`), `n_timepoints` and `memory` (a time course),
#'   `random_state` (0).
#' @return A Python AnnData: `X` = `obsm["X_pca"]` = the embedding; `obs` columns `clone_id`, `donor`, `context`,
#'   `true_programme`, `true_state`, `clone_size` (and `timepoint` with `n_timepoints > 1`).
#' @seealso [run_threadfin()]
#' @examplesIf threadfin_available()
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' adata
#' @export
simulate_repertoire <- function(...) {
  .call_py(c("sim", "simulate_repertoire"), list(), .kwargs(list(...)), convert = TRUE)
}

# ------------------------------------------------------------------------------------------------ arguments

# The expression argument of run_threadfin(): a path, an R object or a Python AnnData, as an AnnData.
.adata_arg <- function(data, assay = NULL, layer = "counts", reductions = TRUE) {
  if (inherits(data, "threadfin_result")) {
    stop("data is already a Threadfin result; pass the AnnData (as_anndata()) or the original object", call. = FALSE)
  }
  if (is.matrix(data) || methods::is(data, "Matrix")) {
    stop("a counts matrix needs its cell metadata (donor, sample): pass as_anndata(counts, obs = metadata)",
      call. = FALSE
    )
  }
  as_anndata(data, assay = assay, layer = layer, reductions = reductions)
}

# The receptor argument: NULL, one file path, a data frame (typed as pandas reads the file) or a Python object.
.bcr_arg <- function(bcr) {
  if (is.null(bcr) || .is_py(bcr)) {
    return(bcr)
  }
  if (is.character(bcr)) {
    if (length(bcr) != 1L) {
      stop("bcr must be one file; read several libraries into one data frame (prefix their barcodes) and pass it",
        call. = FALSE
      )
    }
    path <- path.expand(bcr)
    if (!file.exists(path)) stop("no such file: ", bcr, call. = FALSE)
    return(normalizePath(path))
  }
  if (is.data.frame(bcr)) {
    return(.contigs_to_py(bcr))
  }
  stop("bcr must be a file path (filtered_contig_annotations.csv or an AIRR .tsv), a data frame of contigs or a ",
    "Python object",
    call. = FALSE
  )
}

# A vector of names as a Python list (also for one name), NULL stays None.
.names_arg <- function(x) if (is.null(x)) NULL else as.list(as.character(x))

# The AnnData of an analysis step: a Python AnnData or the AnnData of a run_threadfin() result.
.adata_of <- function(x) {
  if (inherits(x, "threadfin_result")) {
    return(.py_get(py_object(x), "adata"))
  }
  if (.is_py(x)) {
    return(x)
  }
  stop("x must be a Python AnnData (as_anndata(), simulate_repertoire()) or a run_threadfin() result: the steps ",
    "write their results into the AnnData",
    call. = FALSE
  )
}

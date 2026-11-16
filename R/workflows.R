#' Clone-level analysis in one call
#'
#' Runs `threadfin.run()`, the complete Threadfin analysis: clones are called from the receptor contigs within each
#' donor (Dandelion, paired heavy and light chains; for T cells beta and alpha), an expression embedding is built
#' without the receptor genes, one profile is estimated per clone with two or more cells, clonal coherence is
#' tested against clone labels shuffled within sampling units, the clone map is laid out, clonal programmes are
#' looked for, and the labels in `test` are tested with clones as the unit.
#'
#' The AnnData that Threadfin works on keeps every result (`py_object(res)$adata`). When `data` is a Python
#' AnnData, it is changed in place, as in Python; files and R objects are converted to a new AnnData first.
#'
#' @param data The cells' gene expression with raw counts: an `.h5ad` file, a Seurat object, a
#'   SingleCellExperiment, or a Python AnnData (for a counts matrix, `as_anndata(counts, obs = metadata)`). Its cell
#'   metadata needs the donor and sample of each cell.
#' @param bcr The receptor contigs (BCR or TCR), one row per chain: a Cell Ranger `filtered_contig_annotations.csv`
#'   or an AIRR `.tsv` file, a data frame read from one of them (with the same cell names as `data`), or a Python
#'   object (a pandas table, a scirpy AnnData or MuData). `NULL` when the metadata already holds clone ids
#'   (`clone_key`, default `"clone_id"`).
#' @param donor_key Metadata column naming the individual each cell came from: clones are called within donors and
#'   labels are shuffled within donors.
#' @param sample_key Metadata column naming the sample or library of each cell: every clone is described relative
#'   to the cells of its sample, and the coherence test shuffles clone labels within it. Never a tissue, time point
#'   or sort gate you want to compare (those go in `test`); then give the library as `strata_key`.
#' @param ... Further arguments of `threadfin.run()`, by name: `strata_key`, `batch_key`, `state_key`,
#'   `composition_weight`, `time_key`, `test` (labels to test, e.g. `c("isotype", "antigen_binding")`),
#'   `clone_key`, `basis` (an existing embedding; for R objects the reduction name, e.g. `"pca"`, or its `obsm`
#'   key), `counts_layer`, `representation`, `map_n_neighbors`, `map_min_dist`, `n_perm`, `n_boot`,
#'   `random_state`, `n_jobs` (processes; the result is the same for any number), `clone_options` (e.g.
#'   `list(identity = 0.9)`), `verbose`. See the
#'   [function reference](https://github.com/Reynold4k/Threadfin/blob/main/docs/API.md#tfrun).
#' @inheritParams as_anndata
#' @return A `threadfin_result`, a list of
#'   * `summary`: the plain-language report (`summary(res)` prints it);
#'   * `coherence`: one row: the clonal ICC (`icc`), its mean under shuffled clone labels (`null_mean`), `excess`,
#'     `p_value`, the strata and numbers of clones and cells;
#'   * `clones`: one row per profiled clone (two or more cells): size, reliability, donor, map position `x`, `y`,
#'     programme and bootstrap confidence;
#'   * `programmes`: one row per programme: clones, cells, mean reliability and confidence, `stability`;
#'   * `composition`: share of each cell state in each programme (with `state_key`), or `NULL`;
#'   * `profile_tests`: one row per label of `test`: the share of clone-profile variance it explains (`r2`), the
#'     same for shuffled labels (`null_mean`), `excess_r2`, `p_value`, number of clones, sample-level flag;
#'   * `tests`: per label, the programme effects (one row per programme and level, with odds ratios and FDR);
#'   * `memory`: the clonal memory index with its interval, p-value and per-pair table (with `time_key`), or `NULL`;
#'   * `cells`: one row per cell: its clone, receptor columns and programme (to add to a Seurat object or
#'     SingleCellExperiment);
#'   * `methods`: a Methods paragraph with the settings and package versions used;
#'   * `settings` and `skipped` (analyses that could not be run on these data, with the reason).
#'
#'   The Python `ThreadfinResult` is kept for everything else ([py_object()]); results pass to Python functions as
#'   it.
#' @seealso [plot_overview()], [save_threadfin()], [profile_association()] and the other steps, [tf]
#' @examplesIf threadfin_available()
#' \donttest{
#' # a simulated repertoire with a known answer; it has clone ids and an embedding already
#' adata <- simulate_repertoire(n_clones = 600, random_state = 0)
#' res <- run_threadfin(adata, donor_key = "donor", sample_key = "context", basis = "X_pca",
#'                      state_key = "true_state", test = "true_programme")
#' summary(res)
#' res$coherence
#' head(res$clones)
#' res$profile_tests[, c("label", "r2", "null_mean", "p_value")]
#' }
#' # your own data:
#' # res <- run_threadfin("b_cells.h5ad", bcr = "filtered_contig_annotations.csv",
#' #                      donor_key = "donor", sample_key = "sample")
#' # res <- run_threadfin(seurat_obj, bcr = read.csv("filtered_contig_annotations.csv"),
#' #                      donor_key = "donor", sample_key = "sample", n_jobs = 8)
#' @export
run_threadfin <- function(data, bcr = NULL, donor_key = NULL, sample_key = NULL, ..., assay = NULL,
                          layer = "counts", reductions = TRUE) {
  .require_threadfin()
  dots <- list(...)
  if (!is.null(dots[["basis"]]) && !.is_py(data) && !is.character(data)) dots[["basis"]] <- .obsm_key(dots[["basis"]])
  adata <- .adata_arg(data, assay = assay, layer = layer, reductions = reductions)
  kw <- .kwargs(c(list(donor_key = donor_key, sample_key = sample_key), dots))
  py <- .call_py("run", list(adata, .bcr_arg(bcr)), kw)
  .threadfin_result(py)
}

# The R view of a ThreadfinResult: its tables as data frames, its texts, and the Python object.
.threadfin_result <- function(py) {
  coherence <- .result_value(.py_get(py, "coherence"))
  profile_tests <- .result_value(.py_get(py, "profile_tests"))
  skipped <- unlist(.result_value(.py_get(py, "skipped")))
  memory <- .py_get(py, "memory")
  structure(
    list(
      summary = .value(.py_method(py, "summary")),
      coherence = .scalar_table(coherence),
      clones = .frame(.py_get(py, "clones")),
      programmes = .frame(.py_get(py, "programmes")),
      composition = .frame(.py_get(py, "composition")),
      profile_tests = .rbind_fill(lapply(profile_tests, .scalar_table)),
      tests = .result_value(.py_get(py, "tests")),
      memory = if (.is_none(memory)) NULL else .dict_result(memory, "threadfin_memory"),
      cells = .frame(.glue()$run_columns(py)),
      methods = .value(.py_method(py, "methods")),
      settings = .result_value(.py_get(py, "settings")),
      skipped = if (is.null(skipped)) character() else skipped
    ),
    class = "threadfin_result", py = py
  )
}

# The entries of a result that are single values, as a one-row data frame (NULL for NULL).
.scalar_table <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- unclass(x)
  keep <- vapply(x, function(v) is.atomic(v) && length(v) == 1L && is.null(dim(v)), logical(1))
  out <- as.data.frame(x[keep], stringsAsFactors = FALSE, optional = TRUE)
  rownames(out) <- NULL
  out
}

# Rows of several one-row tables, with the union of their columns (NA where a table has no such column).
.rbind_fill <- function(tables) {
  tables <- Filter(Negate(is.null), tables)
  if (!length(tables)) {
    return(NULL)
  }
  cols <- unique(unlist(lapply(tables, names)))
  rows <- lapply(tables, function(t) {
    for (c in setdiff(cols, names(t))) t[[c]] <- NA
    t[cols]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# A dict result as an R list with a class (print and as.data.frame methods), e.g. clonal_coherence().
.dict_result <- function(py, class) {
  structure(.result_value(py), class = c(class, "threadfin_dict"))
}

#' The Python object behind a Threadfin result
#'
#' A [run_threadfin()] or [read_threadfin()] result keeps its Python `ThreadfinResult`, so that anything the R
#' result does not show is one `$` away: `py_object(res)$adata` is the AnnData with every result in
#' `uns["threadfin"]`, `py_object(res)$coherence$null` the shuffled ICCs, `py_object(res)$profile_tests` the full
#' label tests. Results also pass to Python functions as this object (`tf$pl$overview(py_object(res)$adata)`).
#'
#' @param x A Threadfin result.
#' @return The Python object (its tables convert to data frames when accessed).
#' @examplesIf threadfin_available()
#' \donttest{
#' res <- run_threadfin(simulate_repertoire(n_clones = 300, random_state = 0),
#'                      donor_key = "donor", sample_key = "context", basis = "X_pca",
#'                      map_n_neighbors = 15, map_min_dist = 0.5, verbose = FALSE)
#' py_object(res)
#' py_object(res)$adata
#' }
#' @export
py_object <- function(x) {
  py <- attr(x, "py", exact = TRUE)
  if (is.null(py)) stop("not a Threadfin result (no Python object attached)", call. = FALSE)
  if (reticulate::py_is_null_xptr(py)) {
    stop("the Python object of this result is gone (a new R session?); run the analysis again, or read it back ",
      "with read_threadfin()",
      call. = FALSE
    )
  }
  py
}

#' @exportS3Method reticulate::r_to_py
r_to_py.threadfin_result <- function(x, convert = FALSE) py_object(x)

#' A Methods paragraph for a Threadfin analysis
#'
#' Runs `ThreadfinResult.methods()`: what was run, with the settings and package versions of this analysis and the
#' references it relies on, as a starting point for the Methods section of a paper (`res$methods` holds the same
#' text). When the clones were called before the run with [define_clones_dandelion()], pass those clone ids: their
#' recorded settings are then described too.
#'
#' @param x A [run_threadfin()] or [read_threadfin()] result.
#' @param clone_ids The clone ids returned by [define_clones_dandelion()], or `NULL`.
#' @return The paragraph, a character string.
#' @examplesIf threadfin_available()
#' \donttest{
#' res <- run_threadfin(simulate_repertoire(n_clones = 300, random_state = 0),
#'                      donor_key = "donor", sample_key = "context", basis = "X_pca",
#'                      map_n_neighbors = 15, map_min_dist = 0.5, verbose = FALSE)
#' cat(threadfin_methods(res))
#' }
#' @export
threadfin_methods <- function(x, clone_ids = NULL) {
  if (!inherits(x, "threadfin_result")) stop("x must be a Threadfin result", call. = FALSE)
  ids <- if (is.null(clone_ids)) NULL else .series_to_py(clone_ids)
  .value(.py_method(py_object(x), "methods", ids))
}

# A named vector (e.g. clone ids) as a pandas Series, with its attrs (the clone-calling record) restored.
.series_to_py <- function(x) {
  if (.is_py(x)) {
    return(x)
  }
  if (is.null(names(x))) stop("clone_ids must be a named vector (cell -> clone id)", call. = FALSE)
  vals <- as.list(as.character(unname(x)))
  vals[is.na(x)] <- list(NULL)
  attrs <- attr(x, "attrs", exact = TRUE)
  .glue()$series_from_r(vals, names(x), attr(x, "name", exact = TRUE), if (length(attrs)) .to_py(attrs) else NULL)
}

#' Save a Threadfin result, and read it back
#'
#' `save_threadfin()` writes what Tutorial 2 saves, each file by Threadfin itself: the clone table
#' (`<prefix>_clones.csv`, `result.clones.to_csv`), the report (`<prefix>_summary.txt`), the Methods paragraph
#' (`<prefix>_methods.txt`) and the AnnData with every result and the record of the run (`<prefix>_result.h5ad`,
#' `adata.write_h5ad`). `read_threadfin()` rebuilds the result from such an `.h5ad` file
#' (`ThreadfinResult.from_adata`), so that its summary, tables and figures are available again.
#'
#' @param x A [run_threadfin()] result.
#' @param dir Folder to write to (created when missing).
#' @param prefix Start of the file names.
#' @param what Which files to write: any of `"clones"`, `"summary"`, `"methods"`, `"h5ad"`.
#' @param file An `.h5ad` file written by `save_threadfin()` (or by `adata.write_h5ad()` after `threadfin.run()`).
#' @return `save_threadfin()`: the paths written, invisibly. `read_threadfin()`: a `threadfin_result` (see
#'   [run_threadfin()]).
#' @examplesIf threadfin_available()
#' \donttest{
#' res <- run_threadfin(simulate_repertoire(n_clones = 300, random_state = 0),
#'                      donor_key = "donor", sample_key = "context", basis = "X_pca",
#'                      map_n_neighbors = 15, map_min_dist = 0.5, verbose = FALSE)
#' files <- save_threadfin(res, tempdir())
#' basename(files)
#' again <- read_threadfin(files[["h5ad"]])
#' identical(again$summary, res$summary)
#' }
#' @export
save_threadfin <- function(x, dir = ".", prefix = "threadfin", what = c("clones", "summary", "methods", "h5ad")) {
  if (!inherits(x, "threadfin_result")) stop("x must be a Threadfin result", call. = FALSE)
  what <- match.arg(what, several.ok = TRUE)
  dir <- path.expand(dir)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  py <- py_object(x)
  paths <- c(
    clones = paste0(prefix, "_clones.csv"), summary = paste0(prefix, "_summary.txt"),
    methods = paste0(prefix, "_methods.txt"), h5ad = paste0(prefix, "_result.h5ad")
  )[what]
  paths[] <- file.path(normalizePath(dir), paths)
  glue <- .glue()
  for (w in what) {
    switch(w,
      clones = glue$method(.py_get(py, "clones"), "to_csv", list(paths[["clones"]])),
      summary = glue$write_text(paths[["summary"]], .py_method(py, "summary")),
      methods = glue$write_text(paths[["methods"]], .py_method(py, "methods")),
      h5ad = glue$method(.py_get(py, "adata"), "write_h5ad", list(paths[["h5ad"]]))
    )
  }
  invisible(paths)
}

#' @rdname save_threadfin
#' @export
read_threadfin <- function(file) {
  .require_threadfin()
  adata <- as_anndata(file)
  py <- .glue()$result_from_adata(adata)
  .py_warnings()
  .threadfin_result(py)
}

#' @export
print.threadfin_result <- function(x, ...) {
  cat(.py_repr(x), "\n")
  n <- if (is.null(x$clones)) 0L else nrow(x$clones)
  cat(sprintf("%s profiled clones, %s programmes; summary(x) prints the report.\n", format(n, big.mark = ","),
              if (is.null(x$programmes)) 0L else nrow(x$programmes)))
  if (length(x$skipped)) cat("Not run:", paste(names(x$skipped), collapse = "; "), "\n")
  cat("Elements:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.threadfin_result <- function(object, ...) structure(object$summary, class = "threadfin_summary")

#' @export
print.threadfin_summary <- function(x, ...) {
  cat(unclass(x), sep = "\n")
  invisible(x)
}

#' @export
print.threadfin_dict <- function(x, ...) {
  title <- switch(class(x)[1],
    threadfin_coherence = "Clonal coherence",
    threadfin_profile_association = "Profile association",
    threadfin_memory = "Clonal memory",
    threadfin_export_lag = "Export lag",
    "Threadfin result"
  )
  cat(title, "\n")
  tab <- .scalar_table(x)
  if (!is.null(tab) && ncol(tab)) print(tab, row.names = FALSE, digits = 4)
  other <- setdiff(names(x), names(tab))
  if (length(other)) cat("Also:", paste(other, collapse = ", "), "\n")
  invisible(x)
}

#' @export
as.data.frame.threadfin_dict <- function(x, row.names = NULL, optional = FALSE, ...) .scalar_table(x)

.py_repr <- function(x) {
  py <- attr(x, "py", exact = TRUE)
  if (is.null(py) || reticulate::py_is_null_xptr(py)) {
    "<the Python object of this result is gone; run the analysis again>"
  } else {
    reticulate::py_str(py)
  }
}

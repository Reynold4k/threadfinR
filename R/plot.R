#' Draw any Threadfin figure
#'
#' Calls a `threadfin.pl` function (matplotlib) and shows the figure on the R graphics device, writes it to a file,
#' or returns it as a raster. With `file`, the file is the one Threadfin itself writes (the function's own `save`
#' argument), in the format of its extension (`.pdf`, `.png`, `.svg`, ...).
#'
#' @param fun Name of a `threadfin.pl` function (`"clone_map"`, `"coherence"`, `"stability"`,
#'   `"programme_composition"`, `"association"`, `"memory"`, `"label_effects"`, `"heritability"`, `"overview"`,
#'   `"cells"`, `"signature_heatmap"`; see `names(tf$pl)`) or a Python callable that returns a matplotlib Figure or
#'   Axes.
#' @param ... Arguments of the plotting function: unnamed ones are passed by position, named ones by name. A
#'   [run_threadfin()] result stands for its AnnData; data frames and lists of results (e.g. from
#'   [association_test()] or [clonal_memory()]) are converted back to pandas tables and dicts.
#' @param file File to write; `NULL`: none (with `output = "file"`, a temporary PNG).
#' @param dpi Resolution. `NULL`: in `file`, Threadfin's own (300 dpi for raster formats); on the device, the
#'   resolution at which the figure fills it (enlarge the device for more detail); in a raster, the figure's own.
#' @param output `"device"` draws the figure on the current R graphics device and returns `file` invisibly;
#'   `"file"` only writes the file and returns its path; `"raster"` returns the image as a raster for further R
#'   graphics.
#' @return The path of `file` (`"device"`: invisibly, `NULL` without a file), or a raster (`"raster"`).
#' @seealso [plot_overview()], [plot_clone_map()] and the other shortcuts.
#' @examplesIf threadfin_available()
#' \donttest{
#' res <- run_threadfin(simulate_repertoire(n_clones = 300, random_state = 0),
#'                      donor_key = "donor", sample_key = "context", basis = "X_pca",
#'                      map_n_neighbors = 15, map_min_dist = 0.5)
#' plot_threadfin("coherence", res)
#' f <- plot_threadfin("stability", res, file = tempfile(fileext = ".pdf"), output = "file")
#' file.exists(f)
#' }
#' @export
plot_threadfin <- function(fun, ..., file = NULL, dpi = NULL, output = c("device", "file", "raster")) {
  output <- match.arg(output)
  .require_threadfin()
  if (is.character(fun)) {
    if (length(fun) != 1L) stop("fun must be one function name", call. = FALSE)
    fn <- .glue()$attr(.glue()$attr(tf, "pl"), fun)
    if (.is_none(fn)) stop("threadfin.pl has no function '", fun, "'", call. = FALSE)
  } else if (.is_py(fun)) {
    fn <- fun
  } else {
    stop("fun must be the name of a threadfin.pl function or a Python callable", call. = FALSE)
  }
  dots <- list(...)
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  args <- unname(lapply(dots[!nzchar(nm)], .plot_arg))
  kwargs <- .kwargs(lapply(dots[nzchar(nm)], .plot_arg))
  .render(fn, args, kwargs, file = file, dpi = dpi, output = output)
}

#' Threadfin's figures
#'
#' Shortcuts of [plot_threadfin()] for the figures of a Threadfin analysis:
#' * `plot_overview()`: the four-panel overview of a [run_threadfin()] result (`ThreadfinResult.plot()`): the
#'   clone map with one panel per programme, the coherence test, programme stability and programme composition;
#' * `plot_clone_map()`: the clone map, one point per clone (size = number of cells), coloured by a clone-table
#'   or label column (`threadfin.pl.clone_map`; categories as one panel each);
#' * `plot_label_effects()`: how much of the difference between clones each label explains
#'   (`threadfin.pl.label_effects`);
#' * `plot_association()`: programme effects of one label, as a forest plot (`threadfin.pl.association`);
#' * `plot_memory()`: clonal memory, the change of a clone between snapshots against other clones
#'   (`threadfin.pl.memory`).
#'
#' @param x A [run_threadfin()] result; for `plot_clone_map()` also an AnnData after [embed_clones()]; for
#'   `plot_label_effects()` also a named list of [profile_association()] results; for `plot_association()` also an
#'   [association_test()] table; for `plot_memory()` also a [clonal_memory()] result.
#' @param color Column to colour the clones by (`NULL`: the programme).
#' @param label The label of `x$tests` to draw (default: the first).
#' @param level The level of a categorical label to draw (default: the rarest).
#' @param ... Further arguments of the Python plotting function (e.g. `facet`, `ncols`, `title`).
#' @inheritParams plot_threadfin
#' @return As [plot_threadfin()].
#' @examplesIf threadfin_available()
#' \donttest{
#' res <- run_threadfin(simulate_repertoire(n_clones = 300, random_state = 0), donor_key = "donor",
#'                      sample_key = "context", basis = "X_pca", state_key = "true_state",
#'                      test = "true_programme", map_n_neighbors = 15, map_min_dist = 0.5)
#' plot_overview(res)
#' plot_clone_map(res, color = "n_cells")
#' plot_label_effects(res)
#' pdf_file <- plot_overview(res, file = tempfile(fileext = ".pdf"), output = "file")
#' }
#' @name threadfin_plots
NULL

#' @rdname threadfin_plots
#' @export
plot_overview <- function(x, file = NULL, dpi = NULL, output = c("device", "file", "raster")) {
  output <- match.arg(output)
  if (inherits(x, "threadfin_result")) {
    return(.render(.py_get(py_object(x), "plot"), list(), .kwargs(list()), file = file, dpi = dpi, output = output))
  }
  plot_threadfin("overview", x, file = file, dpi = dpi, output = output)
}

#' @rdname threadfin_plots
#' @export
plot_clone_map <- function(x, color = NULL, ..., file = NULL, dpi = NULL, output = c("device", "file", "raster")) {
  plot_threadfin("clone_map", x, color = color, ..., file = file, dpi = dpi, output = output)
}

#' @rdname threadfin_plots
#' @export
plot_label_effects <- function(x, ..., file = NULL, dpi = NULL, output = c("device", "file", "raster")) {
  results <- if (inherits(x, "threadfin_result")) .py_get(py_object(x), "profile_tests") else x
  plot_threadfin("label_effects", results, ..., file = file, dpi = dpi, output = output)
}

#' @rdname threadfin_plots
#' @export
plot_association <- function(x, label = NULL, level = NULL, ..., file = NULL, dpi = NULL,
                             output = c("device", "file", "raster")) {
  if (inherits(x, "threadfin_result")) {
    tests <- .py_get(py_object(x), "tests")
    labels <- as.character(unlist(.value(.glue()$keys(tests))))
    if (!length(labels)) stop("the result has no programme tests (test = ..., with two or more programmes)", call. = FALSE)
    label <- .null_or(label, labels[[1]])
    if (!label %in% labels) stop("no programme test of '", label, "' (tests: ", paste(labels, collapse = ", "), ")",
                                 call. = FALSE)
    x <- .glue()$item(tests, label)
  }
  plot_threadfin("association", x, level = level, ..., file = file, dpi = dpi, output = output)
}

#' @rdname threadfin_plots
#' @export
plot_memory <- function(x, ..., file = NULL, dpi = NULL, output = c("device", "file", "raster")) {
  if (inherits(x, "threadfin_result")) {
    x <- .py_get(py_object(x), "memory")
    if (.is_none(x)) stop("the result has no clonal memory (run_threadfin(..., time_key = ))", call. = FALSE)
  }
  plot_threadfin("memory", x, ..., file = file, dpi = dpi, output = output)
}

# An argument of a plotting function: a run result stands for its AnnData.
.plot_arg <- function(a) {
  if (inherits(a, "threadfin_result")) {
    return(.py_get(py_object(a), "adata"))
  }
  if (is.null(a)) NULL else .to_py(a)
}

.render <- function(fn, args, kwargs, file = NULL, dpi = NULL, output = "device") {
  if (!is.null(file)) file <- normalizePath(path.expand(file), mustWork = FALSE)
  if (output == "file" && is.null(file)) file <- tempfile("threadfin-", fileext = ".png")
  png <- if (output == "file") NULL else tempfile("threadfin-", fileext = ".png")
  target <- NULL
  if (output == "device" && is.null(dpi)) {
    # render at the device's own resolution, so that the image is drawn without resampling
    px <- tryCatch(grDevices::dev.size("px"), error = function(e) NULL)
    if (length(px) == 2L && all(is.finite(px) & px > 0)) target <- as.list(as.numeric(px))
  }
  info <- .glue(TRUE)$render(fn, args, kwargs, save = file, png = png, dpi = dpi, target = target)
  .py_warnings()
  switch(output,
    file = file,
    raster = .read_raster(png),
    device = {
      .draw_png(png, info$background)
      invisible(file)
    }
  )
}

.read_raster <- function(path) {
  if (requireNamespace("png", quietly = TRUE)) {
    return(png::readPNG(path, native = TRUE))
  }
  grDevices::as.raster(reticulate::py_to_r(.glue()$read_png(path)), max = 255)
}

.draw_png <- function(path, background = "white") {
  img <- .read_raster(path)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = background, col = NA))
  grid::grid.raster(img, interpolate = TRUE)
}

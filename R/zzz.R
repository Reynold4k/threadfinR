#' The Python package Threadfin
#'
#' `tf` is the Python package threadfin, loaded the first time it is used, under the name the Threadfin
#' documentation gives it (`import threadfin as tf`). Every Python function is reachable with `$`: `tf$run()`,
#' `tf$tl$clone_profiles()`, `tf$tl$tune_clone_map()`, `tf$pl$clone_map()`, `tf$sim$simulate_repertoire()`, ...
#' with the arguments and documentation of the
#' [Python function reference](https://github.com/Reynold4k/Threadfin/blob/main/docs/API.md). reticulate converts
#' the arguments; in addition, a whole number given where a function takes an integer (`random_state = 0`,
#' `n_jobs = 8`) arrives as an integer, and a single string given where it takes a list of strings arrives as a
#' one-element list. Results are Python objects; their tables convert to data frames when you access them.
#'
#' The wrapper functions ([run_threadfin()], [as_anndata()], the analysis steps such as [clone_profiles()] and
#' [profile_association()], [plot_threadfin()], [save_threadfin()]) cover what the tutorials use, with R objects in
#' and out.
#'
#' Functions that take a list of names where R has a single string (`export_lag(source_states = )`,
#' `gene_heritability(genes = )`) need `list("GC")` when called through `tf`; the wrappers do that for you.
#' The R package tensorflow also exports an object called `tf`; when both are attached, `threadfin::tf` is
#' Threadfin.
#'
#' @format A Python module (a thin proxy of `threadfin`, see `system.file("python", package = "threadfin")`),
#'   imported by reticulate on first use.
#' @seealso [install_threadfin()], [threadfin_available()], `reticulate::py_help()` for a function's Python help.
#' @examplesIf threadfin_available()
#' tf$`__version__`
#' # a repertoire with a known answer
#' adata <- tf$sim$simulate_repertoire(n_clones = 300, random_state = 0)
#' adata
#' head(names(tf$tl))
#' @export
tf <- NULL

.state <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .state$python_path <- normalizePath(system.file("python", package = pkgname, lib.loc = libname), mustWork = FALSE)
  # Registered first and with priority, so that reticulate looks for threadfin when it chooses a Python: it then
  # finds the environment made by install_threadfin() ("r-threadfin") by its own rule for "r-<module>" environments.
  .state$threadfin <- reticulate::import(
    "threadfin",
    delay_load = list(priority = 10L, environment = "r-threadfin")
  )
  tf <<- reticulate::import(
    "threadfin_r",
    delay_load = list(
      before_load = .add_python_path,
      on_load = .remove_python_path,
      on_error = function(e) stop(.no_threadfin_message(conditionMessage(e)), call. = FALSE)
    )
  )
  invisible()
}

# The folder of the proxy module goes on sys.path while it is imported (once; a failed import leaves it there, so
# that a second attempt reports the real cause again).
.add_python_path <- function() {
  sys <- reticulate::import("sys", convert = FALSE)
  if (!.state$python_path %in% reticulate::py_to_r(sys$path)) sys$path$insert(0L, .state$python_path)
  invisible()
}

.remove_python_path <- function() {
  sys <- reticulate::import("sys", convert = FALSE)
  if (.state$python_path %in% reticulate::py_to_r(sys$path)) sys$path$remove(.state$python_path)
  invisible()
}

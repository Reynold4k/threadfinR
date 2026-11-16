#' Install the Python package Threadfin
#'
#' Creates a Python environment (by default the virtual environment `"r-threadfin"`, which this package finds on its
#' own) and installs Threadfin from PyPI into it with pip, with Dandelion for clone calling. Run it once, then
#' restart R.
#'
#' @param envname Name (or path) of the environment. `"r-threadfin"` is found automatically; for another name, call
#'   `reticulate::use_virtualenv(envname)` (or `reticulate::use_condaenv()`) before using the package.
#' @param method `"auto"` (a virtual environment when a Python >= 3.10 that can create one is found, else conda),
#'   `"virtualenv"` or `"conda"`.
#' @param extras Extras of the Python package: `"dandelion"` (the default) installs Dandelion 0.5.x, which calls
#'   clones from receptor contigs; `NULL` installs Threadfin alone (then pass clone ids yourself).
#' @param version Version requirement of threadfin, e.g. `">=4.6"` or `"4.6.1"` (exactly this version); `NULL`
#'   installs the latest release.
#' @param python_version Python version for a new environment: a constraint for virtual environments (default
#'   `">=3.10"`), a version for conda (default `"3.12"`).
#' @param ... Passed to [reticulate::virtualenv_install()] or [reticulate::conda_install()], e.g.
#'   `pip_options = "--upgrade"`.
#' @return The pip requirement that was installed, invisibly.
#' @seealso [threadfin_available()], [tf]
#' @examples
#' \dontrun{
#' install_threadfin()                    # threadfin with Dandelion, from PyPI
#' install_threadfin(version = "4.6.1")   # exactly this version
#' }
#' @export
install_threadfin <- function(envname = "r-threadfin", method = c("auto", "virtualenv", "conda"),
                              extras = "dandelion", version = ">=4.6", python_version = NULL, ...) {
  method <- match.arg(method)
  spec <- .pip_spec(extras, version)
  if (method == "auto") method <- .install_method(envname, python_version)
  if (method == "virtualenv") {
    if (!reticulate::virtualenv_exists(envname)) {
      python <- .venv_starter(.null_or(python_version, ">=3.10"))
      reticulate::virtualenv_create(envname, python = python, packages = NULL)
    }
    reticulate::virtualenv_install(envname, spec, ...)
  } else {
    reticulate::conda_install(envname, spec, python_version = .null_or(python_version, "3.12"), pip = TRUE, ...)
  }
  hint <- if (identical(envname, "r-threadfin")) {
    "library(threadfin) finds it on its own."
  } else if (method == "conda") {
    sprintf("Use it with reticulate::use_condaenv(\"%s\") before calling threadfin functions.", envname)
  } else {
    sprintf("Use it with reticulate::use_virtualenv(\"%s\") before calling threadfin functions.", envname)
  }
  if (reticulate::py_available(initialize = FALSE)) hint <- paste(hint, "Restart R to use it.")
  message(sprintf("Installed %s into the %s environment \"%s\". %s", spec, method, envname, hint))
  invisible(spec)
}

#' Is the Python package Threadfin available?
#'
#' Loads Threadfin in the Python that reticulate uses (see [install_threadfin()]) and reports whether that worked.
#' It never downloads or builds a Python of its own.
#'
#' @return `TRUE` or `FALSE`.
#' @seealso [install_threadfin()]
#' @examples
#' threadfin_available()
#' @export
threadfin_available <- function() {
  if (isTRUE(.state$ready)) {
    return(TRUE)
  }
  if (!reticulate::py_available(initialize = FALSE) && is.na(Sys.getenv("RETICULATE_USE_MANAGED_VENV", NA))) {
    # Without a configured Python, reticulate would download uv and build a managed environment to answer this.
    Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")
    on.exit(Sys.unsetenv("RETICULATE_USE_MANAGED_VENV"), add = TRUE)
  }
  isTRUE(tryCatch(.require_threadfin(), error = function(e) FALSE))
}

.pip_spec <- function(extras = "dandelion", version = ">=4.6") {
  if (!is.null(extras)) {
    extras <- as.character(extras)
    if (anyNA(extras) || !all(grepl("^[A-Za-z0-9_.-]+$", extras))) {
      stop("extras must be names of threadfin extras, e.g. \"dandelion\"", call. = FALSE)
    }
  }
  spec <- paste0("threadfin", if (length(extras)) paste0("[", paste(extras, collapse = ","), "]"))
  if (is.null(version)) {
    return(spec)
  }
  if (length(version) != 1L || is.na(version)) stop("version must be one requirement string", call. = FALSE)
  if (!nzchar(version)) {
    return(spec)
  }
  version <- gsub("\\s+", "", version)
  if (grepl("^[0-9]", version)) version <- paste0("==", version)
  paste0(spec, version)
}

.install_method <- function(envname, python_version = NULL) {
  if (reticulate::virtualenv_exists(envname)) {
    return("virtualenv")
  }
  if (isTRUE(tryCatch(reticulate::condaenv_exists(envname), error = function(e) FALSE))) {
    return("conda")
  }
  if (length(.venv_starters(.null_or(python_version, ">=3.10")))) {
    return("virtualenv")
  }
  if (!is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))) {
    return("conda")
  }
  "virtualenv" # .venv_starter() explains what is missing
}

# Pythons of the requested version that can create virtual environments (venv and ensurepip present).
.venv_starters <- function(version) {
  found <- tryCatch(reticulate::virtualenv_starter(version, all = TRUE), error = function(e) NULL)
  paths <- if (is.data.frame(found)) as.character(found$path) else character()
  paths[vapply(paths, .can_make_venv, logical(1))]
}

.venv_starter <- function(version) {
  ok <- .venv_starters(version)
  if (length(ok)) {
    return(ok[[1]])
  }
  stop(
    "install_threadfin() needs Python ", version, " with its venv and ensurepip modules, and found none. ",
    "Install one with reticulate::install_python(\"3.12\") and run install_threadfin() again, ",
    "or use install_threadfin(method = \"conda\").",
    call. = FALSE
  )
}

.can_make_venv <- function(python) {
  status <- tryCatch(
    suppressWarnings(system2(python, c("-c", shQuote("import venv, ensurepip")), stdout = FALSE, stderr = FALSE)),
    error = function(e) 1L
  )
  identical(as.integer(status), 0L)
}

.no_threadfin_message <- function(detail = NULL) {
  if (length(detail) && grepl("install_threadfin()", detail, fixed = TRUE)) {
    return(detail)
  }
  python <- tryCatch(
    if (reticulate::py_available(initialize = FALSE)) reticulate::py_config()$python else NULL,
    error = function(e) NULL
  )
  paste0(
    "The Python package threadfin is not available",
    if (length(python)) paste0(" in the Python that reticulate uses (", python, ")"),
    ".\nInstall it once with threadfin::install_threadfin() and restart R, or point reticulate at a Python that ",
    "has it: Sys.setenv(RETICULATE_PYTHON = \"/path/to/python\") before the first threadfin call.",
    if (length(detail) && nzchar(detail)) paste0("\nDetails: ", detail)
  )
}

# Import threadfin (through the proxy `tf`) and the glue module once; clear error otherwise.
.require_threadfin <- function() {
  if (isTRUE(.state$ready)) {
    return(invisible(TRUE))
  }
  res <- tryCatch(reticulate::py_get_attr(tf, "__name__"), error = function(e) e)
  if (inherits(res, "error")) stop(.no_threadfin_message(conditionMessage(res)), call. = FALSE)
  .state$glue <- reticulate::import("threadfin_r.glue", convert = FALSE)
  .state$glue_r <- reticulate::import("threadfin_r.glue", convert = TRUE)
  .state$ready <- TRUE
  invisible(TRUE)
}

# The glue module (data movement and figures only; see inst/python/threadfin_r/glue.py).
.glue <- function(convert = FALSE) {
  .require_threadfin()
  if (convert) .state$glue_r else .state$glue
}

# Call a Threadfin function by its path, e.g. .call_py(c("tl", "clone_profiles"), list(adata), .kwargs(...)), with
# R numbers coerced as the Python signature asks and Threadfin's warnings re-issued as R warnings. The result stays
# a Python object; convert = TRUE gives it reticulate's conversion on access (for an AnnData handed to the user).
.call_py <- function(path, args = list(), kwargs = NULL, convert = FALSE) {
  if (is.null(kwargs)) kwargs <- .kwargs(list())
  out <- .glue(convert)$call_path(as.list(path), unname(args), kwargs)
  .py_warnings()
  out
}

# Threadfin's warnings of the last call, as R warnings (other packages' warnings are printed by Python).
.py_warnings <- function() {
  for (msg in as.character(unlist(reticulate::py_to_r(.glue()$last_warnings())))) warning(msg, call. = FALSE)
  invisible()
}

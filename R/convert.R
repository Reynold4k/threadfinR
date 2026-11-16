# Moving data between R and Python. Nothing here computes: values are carried over unchanged, with R's missing
# values kept missing (reticulate turns NA in integer / logical columns into -2147483648 / FALSE, so those columns go
# over as doubles and are restored on the Python side; see inst/python/threadfin_r/glue.py).

.null_or <- function(x, y) if (is.null(x)) y else x

.is_py <- function(x) inherits(x, "python.builtin.object")

.is_none <- function(x) is.null(x) || inherits(x, "python.builtin.NoneType")

# ------------------------------------------------------------------------------------------------- R -> Python

# Any argument for a Python call: the run result gives its Python object, data frames go over with missing values
# and their pandas index intact, lists become dicts (named) or lists, named vectors become dicts; everything else
# is left to reticulate.
.to_py <- function(x) {
  if (inherits(x, "threadfin_result")) {
    return(py_object(x))
  }
  if (.is_py(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    return(.df_to_py(x))
  }
  if (is.list(x) && !is.object(x)) {
    items <- lapply(unclass(x), .to_py)
    nm <- names(x)
    if (!is.null(nm) && all(nzchar(nm))) {
      return(reticulate::py_dict(nm, unname(items), convert = FALSE))
    }
    return(reticulate::r_to_py(unname(items)))
  }
  if (is.list(x)) {
    return(.to_py(unclass(x)))
  }
  if (is.atomic(x) && !is.null(names(x)) && length(x) && !is.factor(x)) {
    return(as.list(x))
  }
  x
}

# Named arguments for a Python call, as a dict (an empty one when there are none); NULL becomes None.
.kwargs <- function(args) {
  if (!length(args)) {
    return(reticulate::dict())
  }
  nm <- names(args)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("additional arguments must be named: they are passed by name to the Python function", call. = FALSE)
  }
  reticulate::py_dict(nm, unname(lapply(args, function(a) if (is.null(a)) NULL else .to_py(a))), convert = FALSE)
}

.has_rownames <- function(df) .row_names_info(df, type = 1L) > 0L

# A data frame as a pandas DataFrame. Row names, or `index`, become the index; a table that came from pandas gets
# its index back from the columns named in attr(df, "index").
.df_to_py <- function(df, index = NULL) {
  set_index <- attr(df, "index", exact = TRUE)
  df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
  if (is.null(index) && .has_rownames(df)) index <- rownames(df)
  if (!is.null(set_index) && !all(set_index %in% names(df))) set_index <- NULL
  logical <- character()
  for (j in seq_along(df)) {
    v <- df[[j]]
    if (is.logical(v) && anyNA(v)) {
      df[[j]] <- as.numeric(v)
      logical <- c(logical, names(df)[j])
    } else if (is.integer(v) && !is.factor(v) && anyNA(v)) {
      df[[j]] <- as.numeric(v)
    } else if (is.list(v)) {
      df[[j]] <- vapply(v, function(e) if (is.null(e)) NA_character_ else paste(as.character(e), collapse = ";"), "")
    } else if (inherits(v, c("Date", "POSIXt"))) {
      df[[j]] <- ifelse(is.na(v), NA_character_, as.character(v))
    }
  }
  rownames(df) <- NULL
  .glue()$frame_from_r(
    reticulate::r_to_py(df),
    index = if (is.null(index)) NULL else as.character(index),
    logical = if (length(logical)) as.list(logical) else NULL,
    set_index = if (length(set_index)) as.list(set_index) else NULL
  )
}

# A receptor contig table (Cell Ranger filtered_contig_annotations.csv or AIRR rearrangements, read into R) as
# pandas reads the same table from a file: every column goes over as the text it would be in a file (logical values
# in the AIRR spelling T / F, numbers so that they read back to the same value) and is read by pandas.read_csv, the
# reader Threadfin uses for these files.
.contigs_to_py <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
  if (!ncol(df)) stop("the receptor table has no columns", call. = FALSE)
  if (anyDuplicated(names(df))) stop("the receptor table has duplicated column names", call. = FALSE)
  cols <- lapply(df, .as_file_text)
  .glue()$contigs_from_r(unname(cols), names(df))
}

.as_file_text <- function(v) {
  out <- if (is.factor(v)) {
    as.character(v)
  } else if (is.logical(v)) {
    ifelse(v, "T", "F")
  } else if (is.integer(v)) {
    as.character(v)
  } else if (is.double(v)) {
    ok <- is.finite(v)
    if (all(v[ok] == round(v[ok])) && all(abs(v[ok]) < 1e15)) sprintf("%.0f", v) else sprintf("%.17g", v)
  } else if (is.list(v)) {
    vapply(v, function(e) if (is.null(e)) NA_character_ else paste(as.character(e), collapse = ";"), "")
  } else {
    as.character(v)
  }
  out[is.na(v)] <- NA_character_
  if (is.double(v)) out[is.nan(v) | is.na(v)] <- NA_character_
  as.list(out)
}

# A sparse dgCMatrix from any matrix-like expression object (base, Matrix, DelayedArray, ...).
.as_dgc <- function(m) {
  if (inherits(m, "dgCMatrix")) {
    return(m)
  }
  if (is.data.frame(m)) m <- as.matrix(m)
  if (is.matrix(m)) storage.mode(m) <- "double"
  out <- tryCatch(
    methods::as(methods::as(methods::as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix"),
    error = function(e) tryCatch(methods::as(m, "dgCMatrix"), error = function(e) NULL)
  )
  if (!inherits(out, "dgCMatrix")) {
    stop("cannot convert the expression matrix (class ", class(m)[1], ") to a sparse dgCMatrix", call. = FALSE)
  }
  out
}

# The obsm key of an R embedding: Seurat reduction "pca" -> "X_pca"; names that start with "X_" are kept.
.obsm_key <- function(name) ifelse(startsWith(name, "X_"), name, paste0("X_", name))

# Embeddings (cells x dimensions) in the order of `cells`; embeddings that miss a cell are left out with a message.
.embeddings_for <- function(emb, cells) {
  out <- list()
  for (nm in names(emb)) {
    e <- as.matrix(emb[[nm]])
    if (!is.numeric(e) || is.null(rownames(e)) || !all(cells %in% rownames(e))) {
      message("embedding '", nm, "' does not cover every cell; it is left out")
      next
    }
    e <- e[cells, , drop = FALSE]
    storage.mode(e) <- "double"
    out[[.obsm_key(nm)]] <- unname(e)
  }
  out
}

# The AnnData Threadfin takes: raw counts (cells x genes), obs metadata, embeddings in obsm.
.anndata <- function(counts, obs = NULL, obsm = NULL) {
  m <- .as_dgc(counts)
  cells <- colnames(m)
  genes <- rownames(m)
  if (is.null(cells) || anyDuplicated(cells)) {
    stop("the expression matrix needs unique cell names (its column names)", call. = FALSE)
  }
  if (is.null(genes)) stop("the expression matrix needs gene names (its row names)", call. = FALSE)
  if (!is.null(obs)) {
    obs <- as.data.frame(obs, optional = TRUE)
    if (!is.null(rownames(obs)) && .has_rownames(obs)) {
      i <- match(cells, rownames(obs))
      if (anyNA(i)) stop(sum(is.na(i)), " cells have no row in the cell metadata", call. = FALSE)
      obs <- obs[i, , drop = FALSE]
    } else if (nrow(obs) != length(cells)) {
      stop("the cell metadata has ", nrow(obs), " rows for ", length(cells), " cells", call. = FALSE)
    }
    obs <- obs[vapply(obs, function(v) is.atomic(v) && is.null(dim(v)), logical(1))]
  }
  obs_py <- if (!is.null(obs) && ncol(obs)) .df_to_py(obs, index = cells) else NULL
  obsm_py <- if (length(obsm)) reticulate::py_dict(names(obsm), unname(obsm), convert = FALSE) else NULL
  .glue(TRUE)$anndata_from_r(
    .np(m@x), .np(m@i), .np(m@p), ncol(m), nrow(m), cells, genes,
    obs = obs_py, obsm = obsm_py
  )
}

# An R vector as a 1-D NumPy array (not a Python list).
.np <- function(v) {
  dim(v) <- length(v)
  reticulate::r_to_py(v)
}

# ---------------------------------------------------------------------------------------------- Seurat objects

.so <- function(name) getExportedValue("SeuratObject", name)

.seurat_counts <- function(x, assay, layer) {
  if (exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    a <- x[[assay]]
    layers <- .so("Layers")(a, search = layer)
    if (!length(layers)) {
      stop("assay '", assay, "' has no layer '", layer, "' (layers: ", paste(.so("Layers")(a), collapse = ", "), ")",
        call. = FALSE
      )
    }
    if (length(layers) > 1L) a <- .so("JoinLayers")(a, layers = layer)
    return(.so("LayerData")(a, layer = layer))
  }
  .so("GetAssayData")(x, assay = assay, slot = layer)
}

.seurat_parts <- function(x, assay = NULL, layer = "counts", reductions = TRUE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("reading a Seurat object needs the R package SeuratObject", call. = FALSE)
  }
  assay <- .null_or(assay, .so("DefaultAssay")(x))
  m <- .seurat_counts(x, assay, layer)
  order <- intersect(colnames(x), colnames(m)) # joined v5 layers can come back grouped: keep the object's order
  if (!identical(order, colnames(m))) m <- m[, order, drop = FALSE]
  cells <- colnames(m)
  meta <- x[[]]
  emb <- list()
  if (!isFALSE(reductions)) {
    red <- .so("Reductions")(x)
    if (is.character(reductions)) red <- intersect(reductions, red)
    for (r in red) emb[[r]] <- .so("Embeddings")(x, reduction = r)
  }
  list(counts = m, obs = meta[cells, , drop = FALSE], obsm = .embeddings_for(emb, cells))
}

# ------------------------------------------------------------------------ SingleCellExperiment / SummarizedExperiment

.se_parts <- function(x, assay = NULL, reductions = TRUE) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("reading a SingleCellExperiment needs the R package SummarizedExperiment", call. = FALSE)
  }
  if (is.null(assay)) {
    assay <- if ("counts" %in% SummarizedExperiment::assayNames(x)) "counts" else 1L
  }
  m <- SummarizedExperiment::assay(x, assay, withDimnames = TRUE)
  cells <- colnames(x)
  if (is.null(cells) || anyDuplicated(cells)) {
    stop("the SingleCellExperiment needs unique column names (cell barcodes)", call. = FALSE)
  }
  obs <- as.data.frame(SummarizedExperiment::colData(x), optional = TRUE)
  obs <- obs[vapply(obs, function(v) is.atomic(v) && is.null(dim(v)), logical(1))]
  rownames(obs) <- cells
  emb <- list()
  if (!isFALSE(reductions) && methods::is(x, "SingleCellExperiment")) {
    rd <- getExportedValue("SingleCellExperiment", "reducedDims")(x)
    red <- names(rd)
    if (is.character(reductions)) red <- intersect(reductions, red)
    for (r in red) {
      e <- as.matrix(rd[[r]])
      if (is.null(rownames(e))) rownames(e) <- cells
      emb[[r]] <- e
    }
  }
  list(counts = m, obs = obs, obsm = .embeddings_for(emb, cells))
}

# ------------------------------------------------------------------------------------------------- Python -> R

# A pandas table as a data frame. A labelled index becomes the first column(s) and its name is kept in
# attr(df, "index"), so that the table goes back to pandas with its index; plain DataFrame.attrs are kept in
# attr(df, "attrs").
.frame <- function(x) {
  if (.is_none(x)) {
    return(NULL)
  }
  df <- .r_value(.glue(TRUE)$frame_value(x))
  if (!is.data.frame(df)) stop("could not convert a pandas table to a data frame", call. = FALSE)
  df
}

.finish_frame <- function(df, logical = NULL, index = NULL, attrs = NULL) {
  if (is.null(df)) {
    return(NULL)
  }
  if (!is.data.frame(df)) stop("could not convert a pandas table to a data frame", call. = FALSE)
  attr(df, "pandas.index") <- NULL
  for (col in unlist(logical)) df[[col]] <- as.logical(df[[col]])
  for (j in seq_along(df)) if (is.double(df[[j]])) df[[j]][is.nan(df[[j]])] <- NA_real_
  rownames(df) <- NULL
  index <- as.character(unlist(index))
  if (length(index)) attr(df, "index") <- index
  if (length(attrs)) attr(df, "attrs") <- attrs
  df
}

# A Python value (str, number, dict, list) as R; NULL for None.
.value <- function(x) {
  if (.is_none(x)) {
    return(NULL)
  }
  reticulate::py_to_r(x)
}

# Any Threadfin result (dict, table, Series, arrays) as plain R values: dicts become named lists, tables data
# frames (as .frame), Series named vectors, arrays vectors / matrices (glue.value_for_r).
.result_value <- function(x) .r_value(.glue(TRUE)$value_for_r(x))

.r_value <- function(v) {
  if (is.null(v)) {
    return(NULL)
  }
  if (is.list(v) && !is.data.frame(v) && !.is_py(v)) {
    nm <- names(v)
    if (!is.null(nm) && "__frame__" %in% nm) {
      return(.finish_frame(v[["__frame__"]], v[["__logical__"]], v[["__index__"]], .r_value(v[["__attrs__"]])))
    }
    if (!is.null(nm) && "__series__" %in% nm) {
      out <- .r_vector(v[["__series__"]])
      names(out) <- as.character(unlist(v[["__names__"]]))
      attrs <- .r_value(v[["__attrs__"]])
      if (length(attrs)) attr(out, "attrs") <- attrs
      return(out)
    }
    out <- lapply(v, .r_value)
    names(out) <- nm
    return(out)
  }
  if (is.array(v) && length(dim(v)) == 1L) v <- as.vector(v)
  if (is.double(v)) v[is.nan(v)] <- NA_real_
  v
}

.r_vector <- function(v) {
  if (is.list(v)) v <- unlist(lapply(v, function(e) if (is.null(e)) NA else e))
  if (is.array(v)) v <- as.vector(v)
  if (is.double(v)) v[is.nan(v)] <- NA_real_
  v
}

# Attribute / method of a Python object, unconverted.
.py_get <- function(obj, name) .glue()$attr(obj, name)

.py_method <- function(obj, name, ...) {
  args <- list(...)
  .glue()$method(obj, name, unname(lapply(args, .to_py)))
}

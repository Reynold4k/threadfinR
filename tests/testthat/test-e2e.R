# End to end on real data (opt-in): the files of Tutorial 5 (GSE287123, made by its steps 1-6) and the result of
# threadfin.run() on them in a separate Python process. Set THREADFIN_E2E_DIR to a folder with
#   gse287123_rbd.h5ad          raw counts with donor, sample and the sort gates of each cell
#   gse287123_rbd_contigs.tsv   the AIRR contigs of those cells
#   python_reference.pkl        the Python result: every table and text (pandas.to_pickle of a dict)
# and THREADFIN_E2E_NJOBS to the n_jobs of the Python run (default 1; the result does not depend on it).

e2e_dir <- function() {
  dir <- Sys.getenv("THREADFIN_E2E_DIR")
  testthat::skip_if(!nzchar(dir), "THREADFIN_E2E_DIR is not set (real data, opt-in)")
  dir
}

e2e_args <- function() {
  list(donor_key = "donor", sample_key = "donor", strata_key = "sample",
       test = c("division_gate", "rbd_bait", "zone_gate"),
       n_jobs = as.integer(Sys.getenv("THREADFIN_E2E_NJOBS", "1")), verbose = FALSE)
}

# The Python result in the shape of a threadfin_result (the same conversions as the R wrapper).
python_reference <- function(dir) {
  if (!exists("e2e_reference", envir = .fixtures, inherits = FALSE)) {
    pd <- reticulate::import("pandas", convert = FALSE)
    ref <- pd$read_pickle(file.path(dir, "python_reference.pkl"))
    rv <- function(key) threadfin:::.result_value(threadfin:::.glue()$item(ref, key))
    out <- list(
      summary = rv("summary"),
      coherence = threadfin:::.scalar_table(rv("coherence")),
      clones = rv("clones"),
      programmes = rv("programmes"),
      composition = rv("composition"),
      profile_tests = threadfin:::.rbind_fill(lapply(rv("profile_tests"), threadfin:::.scalar_table)),
      tests = rv("tests"),
      memory = NULL,
      cells = rv("cells"),
      methods = rv("methods"),
      settings = rv("settings")
    )
    assign("e2e_reference", out, envir = .fixtures)
  }
  get("e2e_reference", envir = .fixtures, inherits = FALSE)
}

e2e_r_objects <- function(dir) {
  a <- as_anndata(file.path(dir, "gse287123_rbd.h5ad"))
  glue <- threadfin:::.glue()
  mx <- reticulate::py_to_r(glue$matrix_for_r(glue$attr(a, "X")))
  obs <- threadfin:::.frame(glue$attr(a, "obs"))
  cells <- obs[[1]]
  obs <- obs[-1]
  rownames(obs) <- cells
  genes <- as.character(reticulate::py_to_r(glue$method(glue$attr(a, "var_names"), "tolist")))
  counts <- methods::new("dgCMatrix",
    x = as.numeric(mx$x), i = as.integer(mx$i), p = as.integer(mx$p),
    Dim = c(as.integer(mx$n_vars), as.integer(mx$n_obs)), Dimnames = list(genes, cells)
  )
  list(counts = counts, obs = obs)
}

test_that("GSE287123 from the files: the R tables equal those of threadfin.run in Python", {
  dir <- e2e_dir()
  skip_if_no_dandelion()
  ref <- python_reference(dir)
  res <- quiet(do.call(run_threadfin, c(list(file.path(dir, "gse287123_rbd.h5ad"),
                                             bcr = file.path(dir, "gse287123_rbd_contigs.tsv")), e2e_args())))
  expect_same_result(res, ref, c(result_tables, "cells", "summary", "methods"))
  # and the numbers of Tutorial 5
  expect_identical(nrow(res$clones), 1881L)
  expect_identical(res$settings$clone_calling$n_clones, 8736L)
  expect_identical(round(res$coherence$icc, 3), 0.169)
  expect_identical(round(res$coherence$null_mean, 3), 0.05)
  expect_identical(nrow(res$programmes), 3L)
  expect_identical(res$programmes$n_clones, c(1076L, 743L, 62L))
  expect_identical(round(res$profile_tests$r2, 3), c(0.085, 0.023, 0.037))
  expect_identical(res$settings$clone_calling$dandelion_hash_seed, 0L)
  assign("e2e_files", res, envir = .fixtures)
})

test_that("GSE287123 as a Seurat object and a SingleCellExperiment with an AIRR data frame: the same tables", {
  dir <- e2e_dir()
  skip_if_no_dandelion()
  skip_if_not_installed("SeuratObject")
  ref <- python_reference(dir)
  parts <- e2e_r_objects(dir)
  airr <- utils::read.delim(file.path(dir, "gse287123_rbd_contigs.tsv"))
  so <- SeuratObject::CreateSeuratObject(parts$counts, meta.data = parts$obs)
  res <- quiet(do.call(run_threadfin, c(list(so, bcr = airr), e2e_args())))
  expect_same_result(res, ref, c(result_tables, "summary"))
  expect_identical(res$cells$clone_id, ref$cells$clone_id)
  expect_identical(res$cells$clone_programme, ref$cells$clone_programme)
  rm(so)
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = parts$counts),
                                                    colData = S4Vectors::DataFrame(parts$obs, check.names = FALSE))
  res2 <- quiet(do.call(run_threadfin, c(list(sce, bcr = airr), e2e_args())))
  expect_same_result(res2, ref, c(result_tables, "summary"))
  expect_identical(res2$cells$clone_id, ref$cells$clone_id)
})

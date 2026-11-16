test_that("run_threadfin returns Threadfin's tables as data frames and keeps the Python result", {
  skip_if_no_threadfin()
  res <- fixture("run")
  expect_s3_class(res, "threadfin_result")
  expect_type(res$summary, "character")
  expect_match(res$summary, "Threadfin summary")
  expect_s3_class(res$coherence, "data.frame")
  expect_identical(nrow(res$coherence), 1L)
  expect_true(all(c("icc", "null_mean", "excess", "p_value", "n_perm", "strata_key") %in% names(res$coherence)))
  expect_s3_class(res$clones, "data.frame")
  expect_identical(names(res$clones)[1], "clone_id")
  expect_true(all(c("n_cells", "reliability", "x", "y", "clone_programme") %in% names(res$clones)))
  expect_true(all(res$clones$n_cells >= 2))
  expect_s3_class(res$programmes, "data.frame")
  expect_true(all(c("programme", "n_clones", "stability") %in% names(res$programmes)))
  expect_s3_class(res$composition, "data.frame")
  expect_equal(unname(rowSums(res$composition[-1])), rep(1, nrow(res$composition)), tolerance = 1e-8)
  expect_identical(res$profile_tests$label, "true_programme")
  expect_true(all(c("r2", "null_mean", "excess_r2", "p_value", "n_clones", "sample_level") %in% names(res$profile_tests)))
  expect_named(res$tests, "true_programme")
  expect_s3_class(res$tests$true_programme, "data.frame")
  expect_null(res$memory)
  expect_identical(names(res$cells), c("cell", "clone_id", "clone_programme", "clone_programme_assigned"))
  expect_identical(nrow(res$cells), as.integer(reticulate::py_to_r(fixture("sim")$n_obs)))
  expect_match(res$methods, "Clonal analyses used Threadfin")
  expect_identical(res$settings$basis, "X_pca")
  expect_type(res$skipped, "character")
  # the numbers are those of the Python result
  py <- py_object(res)
  expect_equal(res$coherence$icc, reticulate::py_to_r(py$coherence)$icc)
  expect_identical(res$clones, threadfin:::.frame(threadfin:::.py_get(py, "clones")))
  expect_identical(res$summary, reticulate::py_to_r(py$summary()))
  # results pass to Python as their Python object
  expect_identical(reticulate::py_to_r(reticulate::r_to_py(res)$summary()), res$summary)
  expect_output(print(res), "ThreadfinResult")
  expect_output(print(summary(res)), "clonal coherence")
})

test_that("the one call gives what Threadfin's own run gives on the same data", {
  skip_if_no_threadfin()
  py <- quiet(tf$run(sim_copy(), donor_key = "donor", sample_key = "context", basis = "X_pca", state_key = "true_state",
                     test = "true_programme", verbose = FALSE))
  direct <- threadfin:::.threadfin_result(py)
  expect_same_result(fixture("run"), direct, c(result_tables, "cells", "summary", "settings"))
})

test_that("Seurat objects, SingleCellExperiments and h5ad files give the same result as the AnnData", {
  skip_if_no_threadfin()
  parts <- fixture("parts")
  ref <- fixture("run")
  args <- list(donor_key = "donor", sample_key = "context", state_key = "true_state", test = "true_programme",
               verbose = FALSE)
  h5 <- tempfile(fileext = ".h5ad")
  sim_copy()$write_h5ad(h5)
  expect_same_result(quiet(do.call(run_threadfin, c(list(h5, basis = "X_pca"), args))), ref)
  skip_if_not_installed("SeuratObject")
  so <- SeuratObject::CreateSeuratObject(parts$counts, meta.data = parts$obs)
  so[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = parts$emb, key = "PC_", assay = "RNA")
  expect_same_result(quiet(do.call(run_threadfin, c(list(so, basis = "pca"), args))), ref)
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = parts$counts), colData = S4Vectors::DataFrame(parts$obs, check.names = FALSE),
    reducedDims = list(PCA = parts$emb)
  )
  expect_same_result(quiet(do.call(run_threadfin, c(list(sce, basis = "PCA"), args))), ref)
})

test_that("a run with time points reports clonal memory", {
  skip_if_no_threadfin()
  a <- simulate_repertoire(n_clones = 400L, n_timepoints = 2L, memory = 0.8, random_state = 0L)
  res <- quiet(run_threadfin(a, donor_key = "donor", sample_key = "context", basis = "X_pca", time_key = "timepoint",
                             map_n_neighbors = 15L, map_min_dist = 0.5, n_boot = 10L, verbose = FALSE))
  expect_s3_class(res$memory, "threadfin_memory")
  expect_true(is.numeric(res$memory$memory_index))
  expect_length(res$memory$memory_index_ci, 2L)
  expect_s3_class(res$memory$pairs, "data.frame")
  expect_equal(res$memory$memory_index, reticulate::py_to_r(py_object(res)$memory)$memory_index)
  expect_output(print(res$memory), "Clonal memory")
  expect_s3_class(as.data.frame(res$memory), "data.frame")
})

test_that("Threadfin's warnings are R warnings, and its errors R errors", {
  skip_if_no_threadfin()
  a <- simulate_repertoire(n_clones = 300L, random_state = 0L)
  msgs <- character()
  withCallingHandlers(
    run_threadfin(a, donor_key = "donor", basis = "X_pca", map_n_neighbors = 15L, map_min_dist = 0.5,
                  n_boot = 5L, verbose = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No sample_key was given", msgs, fixed = TRUE)))
  expect_error(run_threadfin(sim_copy(), donor_key = "no_such_column", basis = "X_pca", verbose = FALSE),
               "no_such_column")
})

test_that("a result saves to files and reads back", {
  skip_if_no_threadfin()
  res <- fixture("run")
  dir <- tempfile("saved-")
  files <- save_threadfin(res, dir)
  expect_named(files, c("clones", "summary", "methods", "h5ad"))
  expect_true(all(file.exists(files)))
  expect_identical(readLines(files[["summary"]], warn = FALSE), strsplit(res$summary, "\n")[[1]])
  clones <- utils::read.csv(files[["clones"]])
  expect_identical(clones$clone_id, res$clones$clone_id)
  again <- quiet(read_threadfin(files[["h5ad"]]))
  expect_identical(again$summary, res$summary)
  expect_identical(again$methods, res$methods)
  # the same values; an h5ad file keeps dicts in the order of their keys and text columns as categories
  for (k in c("coherence", "profile_tests")) expect_equal(again[[k]][names(res[[k]])], res[[k]], label = k)
  for (k in c("clones", "programmes", "composition")) expect_equal(again[[k]], res[[k]], label = k)
  expect_identical(as.character(again$cells$clone_id), res$cells$clone_id)
  expect_identical(as.character(again$cells$clone_programme), as.character(res$cells$clone_programme))
  only <- save_threadfin(res, dir, prefix = "x", what = "methods")
  expect_identical(basename(only), "x_methods.txt")
})

test_that("threadfin_methods() is the Methods paragraph", {
  skip_if_no_threadfin()
  res <- fixture("run")
  expect_identical(threadfin_methods(res), res$methods)
  expect_error(threadfin_methods(res, clone_ids = "x"), "named vector")
})

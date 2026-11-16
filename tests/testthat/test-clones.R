# Clone calling with Dandelion and the whole one-call path from files and R objects, on the synthetic files of
# helper-threadfin.R (raw counts and receptor contigs).

run_args <- list(donor_key = "donor", sample_key = "sample", test = "state", n_perm = 50L, n_boot = 5L,
                 verbose = FALSE)

test_that("define_clones_dandelion calls the same clones from a data frame and from the file", {
  skip_if_no_dandelion()
  files <- fixture("files")
  airr <- utils::read.delim(files$paths[["airr"]])
  airr$donor <- files$obs[airr$cell_id, "donor"]
  ids <- quiet(define_clones_dandelion(airr, donor_key = "donor", verbose = FALSE))
  expect_type(ids, "character")
  expect_setequal(names(ids), rownames(files$obs))
  expect_identical(length(unique(ids)), 80L)                       # 40 clones of three cells, 40 of one
  expect_true(all(table(ids) %in% c(1L, 3L)))
  info <- attr(ids, "attrs")
  expect_identical(info$dandelion_hash_seed, 0L)
  expect_identical(info$n_clones, 80L)
  expect_identical(info$dandelion_subprocess, Sys.getenv("PYTHONHASHSEED") != "0")
  # the same contigs from the file (donor added), with Threadfin's reader
  tsv <- tempfile(fileext = ".tsv")
  utils::write.table(airr, tsv, sep = "\t", quote = FALSE, na = "", row.names = FALSE)
  from_file <- quiet(define_clones_dandelion(tsv, donor_key = "donor", verbose = FALSE))
  expect_identical(as.vector(from_file[names(ids)]), as.vector(ids))
  # the clone ids of a Python call
  py <- quiet(tfpy()$define_clones_dandelion(tfpy()$read_contigs(tsv), donor_key = "donor", verbose = FALSE))
  expect_identical(as.vector(threadfin:::.result_value(py)[names(ids)]), as.vector(ids))
  # read_contigs gives AIRR names, also for Cell Ranger's csv
  rc <- read_contigs(files$paths[["tenx"]])
  expect_true(all(c("cell_id", "v_call", "j_call", "junction", "junction_aa", "locus") %in% names(rc)))
  expect_identical(nrow(rc), nrow(files$tenx))
})

test_that("the one call from files equals Threadfin's own run on the same files", {
  skip_if_no_dandelion()
  files <- fixture("files")
  res <- quiet(do.call(run_threadfin, c(list(files$paths[["h5ad"]], bcr = files$paths[["airr"]]), run_args)))
  expect_s3_class(res, "threadfin_result")
  expect_identical(nrow(res$clones), 40L)
  expect_true("bcr_v_call" %in% names(res$cells))
  expect_identical(res$settings$receptor, "BCR")
  expect_identical(res$settings$clone_calling$dandelion_hash_seed, 0L)
  py <- quiet(do.call(tf$run, c(list(as_anndata(files$paths[["h5ad"]]), bcr = files$paths[["airr"]]), run_args)))
  expect_same_result(res, threadfin:::.threadfin_result(py), c(result_tables, "cells", "summary"))
  expect_gt(res$coherence$icc, res$coherence$null_mean)
})

test_that("Seurat and SingleCellExperiment objects with contig data frames give the result of the files", {
  skip_if_no_dandelion()
  skip_if_not_installed("SeuratObject")
  files <- fixture("files")
  ref <- quiet(do.call(run_threadfin, c(list(files$paths[["h5ad"]], bcr = files$paths[["airr"]]), run_args)))
  so <- SeuratObject::CreateSeuratObject(Matrix::Matrix(files$counts, sparse = TRUE), meta.data = files$obs)
  airr <- utils::read.delim(files$paths[["airr"]])
  res <- quiet(do.call(run_threadfin, c(list(so, bcr = airr), run_args)))
  expect_same_result(res, ref)
  # a data frame of contigs gives clone ids only (Threadfin summarises receptors per cell from a file)
  expect_identical(res$cells$clone_id, ref$cells$clone_id)
  expect_identical(res$cells$clone_programme, ref$cells$clone_programme)
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = files$counts),
                                                    colData = S4Vectors::DataFrame(files$obs))
  expect_same_result(quiet(do.call(run_threadfin, c(list(sce, bcr = airr), run_args))), ref)
  # Cell Ranger's csv as a file and as a data frame
  tenx_file <- quiet(do.call(run_threadfin, c(list(files$paths[["h5ad"]], bcr = files$paths[["tenx"]]), run_args)))
  tenx_df <- quiet(do.call(run_threadfin, c(list(sce, bcr = utils::read.csv(files$paths[["tenx"]])), run_args)))
  expect_same_result(tenx_df, tenx_file)
  expect_identical(tenx_file$clones$n_cells, ref$clones$n_cells)
})

test_that("the result does not depend on n_jobs (worker processes under reticulate)", {
  skip_if_no_dandelion()
  skip_on_cran()
  files <- fixture("files")
  one <- quiet(do.call(run_threadfin, c(list(files$paths[["h5ad"]], bcr = files$paths[["airr"]]), run_args)))
  two <- quiet(do.call(run_threadfin, c(list(files$paths[["h5ad"]], bcr = files$paths[["airr"]], n_jobs = 2L),
                                        run_args[names(run_args) != "n_jobs"])))
  expect_same_result(two, one)
})

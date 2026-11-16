test_that("data frames cross to pandas and back with missing values and the index intact", {
  skip_if_no_threadfin()
  df <- data.frame(
    s = c("x", NA, "z"), i = c(1L, NA, 3L), b = c(TRUE, NA, FALSE), d = c(1.5, NA, 2),
    f = factor(c("u", NA, "v")), j = 1:3, row.names = c("r1", "r2", "r3"), stringsAsFactors = FALSE
  )
  p <- threadfin:::.df_to_py(df)
  expect_match(reticulate::py_str(p$dtypes), "boolean")
  back <- threadfin:::.frame(p)
  expect_identical(attr(back, "index"), "index")
  expect_identical(back$index, c("r1", "r2", "r3"))
  expect_identical(back$s, df$s)
  expect_equal(back$i, c(1, NA, 3))
  expect_identical(back$b, df$b)
  expect_equal(back$d, df$d)
  expect_identical(as.character(back$f), as.character(df$f))
  expect_equal(back$j, 1:3)
  # a table from pandas goes back with its index (also after keeping some rows)
  again <- threadfin:::.df_to_py(back[back$j > 1, ])
  expect_identical(unlist(reticulate::py_to_r(again$index$tolist())), c("r2", "r3"))
  expect_false("index" %in% unlist(reticulate::py_to_r(again$columns$tolist())))
})

test_that("nested Threadfin results become lists, vectors and data frames", {
  skip_if_no_threadfin()
  reticulate::py_run_string(paste(
    "import numpy as np, pandas as pd",
    "_tfr_obj = {'a': 1, 'b': np.float64(2.5), 'ci': (0.1, 0.9), 'none': None, 'arr': np.arange(3),",
    "            'zero_d': np.array(4.0), 'names': ['p', 'q'], 'nested': {'k': [1.0, float('nan')]},",
    "            'tab': pd.DataFrame({'x': [1, 2], 'y': ['u', None]}, index=pd.Index(['r1', 'r2'], name='clone_id')),",
    "            'records': [{'level': 'a', 'n': 1}, {'level': 'b', 'n': 2}],",
    "            'series': pd.Series([0.5, 0.25], index=['c1', 'c2'])}",
    "_tfr_obj['tab'].attrs['sample_level'] = True",
    sep = "\n"
  ))
  v <- threadfin:::.result_value(reticulate::py_eval("_tfr_obj", convert = FALSE))
  expect_identical(v$a, 1L)
  expect_identical(v$b, 2.5)
  expect_identical(v$ci, c(0.1, 0.9))
  expect_null(v$none)
  expect_equal(v$arr, 0:2)
  expect_identical(v$zero_d, 4)
  expect_identical(v$names, c("p", "q"))
  expect_identical(v$nested$k, c(1, NA))
  expect_s3_class(v$tab, "data.frame")
  expect_identical(v$tab$clone_id, c("r1", "r2"))
  expect_identical(v$tab$y, c("u", NA))
  expect_identical(attr(v$tab, "index"), "clone_id")
  expect_identical(attr(v$tab, "attrs"), list(sample_level = TRUE))
  expect_s3_class(v$records, "data.frame")
  expect_identical(v$records$level, c("a", "b"))
  expect_identical(v$series, c(c1 = 0.5, c2 = 0.25))
})

test_that("receptor contigs reach pandas as pandas reads the file", {
  skip_if_no_threadfin()
  files <- fixture("files")
  reticulate::py_run_string(paste(
    "def _tfr_same(a, b):",
    "    out = {}",
    "    for c in a.columns:",
    "        x, y = a[c], b[c]",
    "        if x.dtype == bool or y.dtype == bool:   # logical values arrive in the AIRR spelling T / F",
    "            out[c] = bool((x.astype(str).str.upper().str[0] == y.astype(str).str.upper().str[0]).all())",
    "        else:",
    "            out[c] = bool(x.equals(y)) and str(x.dtype) == str(y.dtype)",
    "    return out",
    sep = "\n"
  ))
  pd <- reticulate::import("pandas", convert = FALSE)
  main <- reticulate::import_main(convert = TRUE)
  airr <- utils::read.delim(files$paths[["airr"]])
  same <- main$`_tfr_same`(pd$read_csv(files$paths[["airr"]], sep = "\t", low_memory = FALSE),
                           threadfin:::.contigs_to_py(airr))
  expect_true(all(unlist(same)), label = paste(names(same)[!unlist(same)], collapse = ", "))
  tenx <- utils::read.csv(files$paths[["tenx"]])
  same <- main$`_tfr_same`(pd$read_csv(files$paths[["tenx"]], low_memory = FALSE), threadfin:::.contigs_to_py(tenx))
  expect_true(all(unlist(same)), label = paste(names(same)[!unlist(same)], collapse = ", "))
  # factors, doubles and missing values are written as a file would hold them
  txt <- threadfin:::.as_file_text
  expect_identical(unlist(txt(factor(c("a", NA)))), c("a", NA))
  expect_identical(unlist(txt(c(1, 20, NA))), c("1", "20", NA))
  expect_identical(as.numeric(unlist(txt(c(0.1, 1 / 3)))), c(0.1, 1 / 3))
  expect_identical(unlist(txt(c(TRUE, NA, FALSE))), c("T", NA, "F"))
})

test_that("matrices of any kind become dgCMatrix", {
  m <- matrix(c(0, 1, 2, 0), 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  s <- threadfin:::.as_dgc(m)
  expect_s4_class(s, "dgCMatrix")
  expect_equal(as.matrix(s), m)
  expect_s4_class(threadfin:::.as_dgc(methods::as(s, "TsparseMatrix")), "dgCMatrix")
  expect_identical(threadfin:::.obsm_key(c("pca", "X_umap", "PCA")), c("X_pca", "X_umap", "X_PCA"))
})

same_anndata <- function() {
  reticulate::py_eval(paste0(
    "lambda a, b, cols: bool((a.X != b.X).nnz == 0 and str(a.X.dtype) == str(b.X.dtype) and ",
    "list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names) and ",
    "all(a.obs[c].astype(str).tolist() == b.obs[c].astype(str).tolist() for c in cols))"
  ))
}

test_that("a counts matrix, a Seurat object and a SingleCellExperiment give the AnnData of the h5ad file", {
  skip_if_no_threadfin()
  files <- fixture("files")
  ref <- as_anndata(files$paths[["h5ad"]])
  same <- same_anndata()
  cols <- list("donor", "sample", "state")
  a <- as_anndata(files$counts, obs = files$obs)
  expect_true(same(a, ref, cols))
  expect_identical(reticulate::py_eval("lambda a: str(a.X.dtype)")(a), "float32")
  skip_if_not_installed("SeuratObject")
  so <- SeuratObject::CreateSeuratObject(Matrix::Matrix(files$counts, sparse = TRUE), meta.data = files$obs)
  expect_true(same(as_anndata(so), ref, cols))
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = files$counts),
                                                    colData = S4Vectors::DataFrame(files$obs))
  expect_true(same(as_anndata(sce), ref, cols))
})

test_that("Seurat reductions and reducedDims become obsm embeddings", {
  skip_if_no_threadfin()
  skip_if_not_installed("SeuratObject")
  parts <- fixture("parts")
  so <- SeuratObject::CreateSeuratObject(parts$counts, meta.data = parts$obs)
  so[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = parts$emb, key = "PC_", assay = "RNA")
  a <- as_anndata(so)
  obsm_keys <- function(a) unlist(threadfin:::.value(threadfin:::.glue()$keys(a$obsm)))
  expect_true("X_pca" %in% obsm_keys(a))
  same_emb <- reticulate::py_eval("lambda a, b: bool((a.obsm['X_pca'] == b.obsm['X_pca']).all())")
  expect_true(same_emb(a, fixture("sim")))
  expect_false("X_pca" %in% obsm_keys(as_anndata(so, reductions = FALSE)))
  # Seurat v5 split layers are joined first
  skip_if_not_installed("SeuratObject", "5.0.0")
  so2 <- so
  so2[["RNA"]] <- split(so2[["RNA"]], f = so2$donor)
  expect_gt(length(SeuratObject::Layers(so2[["RNA"]], search = "counts")), 1L)
  m <- threadfin:::.seurat_counts(so2, "RNA", "counts")
  expect_equal(m[, parts$cells], parts$counts)
})

test_that("unclear input stops with a clear message", {
  skip_if_no_threadfin()
  expect_error(as_anndata(list(1)), "h5ad file, a Seurat object")
  expect_error(as_anndata("no_such_file.h5ad"), "no such file")
  expect_error(as_anndata(c("a.h5ad", "b.h5ad")), "one .h5ad")
  expect_error(threadfin:::.bcr_arg(1), "bcr must be")
  expect_error(threadfin:::.bcr_arg(c("a.csv", "b.csv")), "one file")
  expect_error(run_threadfin(fixture("run"), donor_key = "donor"), "already a Threadfin result")
  expect_error(run_threadfin(fixture("files")$counts, donor_key = "donor"), "needs its cell metadata")
})

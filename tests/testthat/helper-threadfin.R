# Tests that need the Python package skip cleanly without it.
skip_if_no_threadfin <- function() {
  testthat::skip_if_not(threadfin_available(), "Python package threadfin not available")
}

skip_if_no_dandelion <- function() {
  skip_if_no_threadfin()
  testthat::skip_if_not(reticulate::py_module_available("dandelion"), "Python package dandelion not available")
}

# The Python package without reticulate's conversion: direct Python calls whose results the tests convert with the
# package's own converters, to compare with what the wrappers return.
tfpy <- function() reticulate::import("threadfin", convert = FALSE)

# Threadfin's warnings arrive as R warnings; tests that are not about them silence them.
quiet <- function(expr) suppressWarnings(expr)

# Simulated repertoires of the Python package and results on them, computed once per test run.
.fixtures <- new.env(parent = emptyenv())

fixture <- function(name) {
  if (!exists(name, envir = .fixtures, inherits = FALSE)) {
    value <- switch(name,
      sim = simulate_repertoire(n_clones = 600, random_state = 0L),
      run = quiet(run_threadfin(sim_copy(), donor_key = "donor", sample_key = "context", basis = "X_pca",
                                state_key = "true_state", test = "true_programme", verbose = FALSE)),
      parts = sim_parts(),
      time = time_course(),
      files = user_files(),
      stop("unknown fixture ", name)
    )
    assign(name, value, envir = .fixtures)
  }
  get(name, envir = .fixtures, inherits = FALSE)
}

# A fresh copy of the simulated repertoire (the analyses write into their AnnData).
sim_copy <- function() fixture("sim")$copy()

# The simulated repertoire as R objects: obs (data frame), the embedding (cells x 20, row names = cells) and a
# counts matrix (genes x cells) that the analysis does not use when basis is given.
sim_parts <- function() {
  a <- fixture("sim")
  glue <- threadfin:::.glue()
  obs <- threadfin:::.frame(glue$attr(a, "obs"))
  cells <- obs[[1]]
  obs <- obs[-1]
  rownames(obs) <- cells
  emb <- reticulate::py_to_r(glue$item(glue$attr(a, "obsm"), "X_pca"))
  dimnames(emb) <- list(cells, paste0("PC_", seq_len(ncol(emb))))
  counts <- Matrix::Matrix(round(abs(t(emb)) * 10), sparse = TRUE)
  dimnames(counts) <- list(paste0("g", seq_len(nrow(counts))), cells)
  list(obs = obs, emb = emb, counts = counts, cells = cells)
}

# A simulated time course with profiles, a clone map and programmes, for the step-by-step functions.
time_course <- function() {
  a <- simulate_repertoire(n_clones = 400L, n_timepoints = 2L, memory = 0.8, random_state = 0L)
  quiet({
    clone_profiles(a, basis = "X_pca", context_key = "context", donor_key = "donor", representation = "kernel",
                   verbose = FALSE)
    embed_clones(a, n_neighbors = 15L, min_dist = 0.5, verbose = FALSE)
    find_programmes(a, verbose = FALSE)
  })
  a
}

# A small synthetic experiment as a user has it on disk: raw counts with donor and sample (.h5ad), and the
# receptor contigs as an AIRR table (.tsv) and as a Cell Ranger filtered_contig_annotations.csv. 40 clones of
# three cells and 40 cells of their own, in two donors with two samples each; cells of a clone share a state.
user_files <- function() {
  withr_seed <- function(seed, code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    code
  }
  withr_seed(1L, {
    codon <- c(A = "GCT", C = "TGT", D = "GAT", E = "GAA", F = "TTT", G = "GGT", H = "CAT", I = "ATT", K = "AAA",
               L = "CTG", M = "ATG", N = "AAT", P = "CCT", Q = "CAG", R = "CGT", S = "TCT", T = "ACT", V = "GTT",
               W = "TGG", Y = "TAT")
    aa <- setdiff(names(codon), c("C", "W"))
    junction <- function(n) {
      x <- paste0("C", paste(sample(aa, n, replace = TRUE), collapse = ""), "W")
      c(aa = x, nt = paste(codon[strsplit(x, "")[[1]]], collapse = ""))
    }
    n_clones <- 40L
    sizes <- c(rep(3L, n_clones), rep(1L, 40L))
    clone <- rep(seq_along(sizes), sizes)
    n_cells <- length(clone)
    donor_of_clone <- rep(c("D1", "D2"), length.out = length(sizes))
    cells <- sprintf("cell%03d", seq_len(n_cells))
    donor <- donor_of_clone[clone]
    sample <- paste0(donor, "_S", ifelse(stats::runif(n_cells) < 0.5, 1, 2))
    state <- sample(1:3, length(sizes), replace = TRUE)[clone]
    genes <- c(sprintf("Gene%03d", 1:300), "Ighm", "Ighg1", "Igkc", "Iglc1", "Ighv1-72", "Igkv3-2")
    base <- stats::rgamma(length(genes), 2, 1)
    lift <- matrix(1, length(genes), 3)
    for (s in 1:3) lift[sample(300, 40), s] <- 4
    mu <- base * lift[, state] * rep(stats::runif(n_cells, 0.7, 1.3), each = length(genes))
    counts <- matrix(stats::rpois(length(mu), mu), length(genes), n_cells, dimnames = list(genes, cells))
    obs <- data.frame(donor = donor, sample = sample, state = paste0("S", state), row.names = cells)
    vh <- c("IGHV1-1*01", "IGHV1-2*02", "IGHV3-23*01", "IGHV4-34*01")
    jh <- c("IGHJ1*01", "IGHJ4*02", "IGHJ6*02")
    vk <- c("IGKV1-1*01", "IGKV3-20*01", "IGKV1-5*03")
    jk <- c("IGKJ1*01", "IGKJ2*01")
    rows <- list()
    for (k in seq_along(sizes)) {
      h <- junction(sample(10:14, 1))
      l <- junction(sample(8:10, 1))
      v1 <- sample(vh, 1)
      j1 <- sample(jh, 1)
      v2 <- sample(vk, 1)
      j2 <- sample(jk, 1)
      for (cell in cells[clone == k]) {
        rows[[length(rows) + 1L]] <- data.frame(
          cell_id = cell, sequence_id = paste0(cell, "_contig_", 1:2), locus = c("IGH", "IGK"),
          v_call = c(v1, v2), d_call = c("IGHD1-1*01", NA), j_call = c(j1, j2), c_call = c("IGHM", "IGKC"),
          junction = c(h[["nt"]], l[["nt"]]), junction_aa = c(h[["aa"]], l[["aa"]]), productive = TRUE,
          umi_count = c(12L, 7L), consensus_count = c(240L, 150L)
        )
      }
    }
    airr <- do.call(rbind, rows)
    tenx <- data.frame(
      barcode = airr$cell_id, is_cell = "true", contig_id = airr$sequence_id, high_confidence = "true",
      length = nchar(airr$junction) + 300L, chain = airr$locus, v_gene = sub("\\*.*", "", airr$v_call),
      d_gene = ifelse(is.na(airr$d_call), "None", sub("\\*.*", "", airr$d_call)),
      j_gene = sub("\\*.*", "", airr$j_call), c_gene = airr$c_call, full_length = "true", productive = "true",
      cdr3 = airr$junction_aa, cdr3_nt = airr$junction, reads = airr$consensus_count, umis = airr$umi_count
    )
    dir <- tempfile("threadfin-files-")
    dir.create(dir)
    paths <- c(h5ad = file.path(dir, "cells.h5ad"), airr = file.path(dir, "airr_rearrangement.tsv"),
               tenx = file.path(dir, "filtered_contig_annotations.csv"))
    a <- as_anndata(Matrix::Matrix(counts, sparse = TRUE), obs = obs)
    a$write_h5ad(paths[["h5ad"]])
    utils::write.table(airr, paths[["airr"]], sep = "\t", quote = FALSE, na = "", row.names = FALSE)
    utils::write.csv(tenx, paths[["tenx"]], quote = FALSE, row.names = FALSE)
    list(paths = paths, counts = counts, obs = obs, airr = airr, tenx = tenx)
  })
}

# The tables of two Threadfin results that must agree when the same analysis ran on the same data.
result_tables <- c("coherence", "clones", "programmes", "composition", "profile_tests", "tests", "memory")

expect_same_result <- function(object, expected, tables = result_tables) {
  for (k in tables) testthat::expect_identical(object[[k]], expected[[k]], label = paste0("$", k))
}

png_magic <- function(path) {
  identical(readBin(path, "raw", 8L), as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))
}

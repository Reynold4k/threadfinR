# Each step gives what the Python function gives on a copy of the same AnnData.

tc <- function() fixture("time")$copy()

test_that("clone_profiles, clonal_coherence, embed_clones and find_programmes equal their Python functions", {
  skip_if_no_threadfin()
  frame <- threadfin:::.frame
  a <- simulate_repertoire(n_clones = 400L, n_timepoints = 2L, memory = 0.8, random_state = 0L)
  b <- a$copy()
  p_r <- quiet(clone_profiles(a, basis = "X_pca", context_key = "context", donor_key = "donor",
                              representation = "kernel", verbose = FALSE))
  p_py <- quiet(frame(tfpy()$tl$clone_profiles(b, basis = "X_pca", context_key = "context", donor_key = "donor",
                                           representation = "kernel", verbose = FALSE)))
  expect_identical(p_r, p_py)
  expect_identical(names(p_r)[1:3], c("clone_id", "n_cells", "reliability"))
  coh <- clonal_coherence(a, strata_key = "context", n_perm = 50L, verbose = FALSE)
  coh_py <- threadfin:::.result_value(tfpy()$tl$clonal_coherence(b, strata_key = "context", n_perm = 50L,
                                                                  verbose = FALSE))
  expect_s3_class(coh, "threadfin_coherence")
  expect_identical(coh$icc, coh_py$icc)
  expect_identical(coh$p_value, coh_py$p_value)
  expect_identical(unclass(coh), coh_py)
  expect_identical(nrow(as.data.frame(coh)), 1L)
  expect_output(print(coh), "Clonal coherence")
  m_r <- quiet(embed_clones(a, n_neighbors = 15L, min_dist = 0.5, verbose = FALSE))
  m_py <- quiet(frame(tfpy()$tl$embed_clones(b, n_neighbors = 15L, min_dist = 0.5, verbose = FALSE)))
  expect_identical(m_r, m_py)
  expect_true(all(c("x", "y") %in% names(m_r)))
  pr_r <- quiet(find_programmes(a, n_boot = 10L, verbose = FALSE))
  pr_py <- quiet(frame(tfpy()$tl$find_programmes(b, n_boot = 10L, verbose = FALSE)))
  expect_identical(pr_r, pr_py)
  expect_identical(names(pr_r)[1], "programme")
})

test_that("label tests and clone labels equal their Python functions", {
  skip_if_no_threadfin()
  a <- tc()
  b <- tc()
  pa <- quiet(profile_association(a, "true_programme", n_perm = 200L, verbose = FALSE))
  pa_py <- quiet(threadfin:::.result_value(tfpy()$tl$profile_association(b, "true_programme", n_perm = 200L,
                                                                     verbose = FALSE)))
  expect_s3_class(pa, "threadfin_profile_association")
  expect_identical(unclass(pa), pa_py)
  expect_s3_class(pa$clone_size, "data.frame")
  at <- quiet(association_test(a, "true_programme", n_perm = 200L, verbose = FALSE))
  at_py <- quiet(threadfin:::.frame(tfpy()$tl$association_test(b, "true_programme", n_perm = 200L, verbose = FALSE)))
  expect_identical(at, at_py)
  expect_true(all(c("programme", "level", "odds_ratio", "pvalue", "fdr") %in% names(at)))
  expect_identical(attr(at, "attrs")$strata_key, "donor")
  cl <- clone_labels(a, "true_state", how = "fraction:S1")
  expect_type(cl, "double")
  expect_true(all(startsWith(names(cl), "D")))
  expect_identical(cl, threadfin:::.result_value(tfpy()$tl$clone_labels(b, "true_state", how = "fraction:S1")))
})

test_that("clonal memory, programme composition and the gene-level summaries equal their Python functions", {
  skip_if_no_threadfin()
  a <- tc()
  b <- tc()
  rv <- threadfin:::.result_value
  mem <- quiet(clonal_memory(a, "timepoint", order = c("t0", "t1"), n_null = 50L, n_boot = 50L, verbose = FALSE))
  mem_py <- quiet(rv(tfpy()$tl$clonal_memory(b, "timepoint", order = list("t0", "t1"), n_null = 50L, n_boot = 50L,
                                         verbose = FALSE)))
  expect_identical(unclass(mem), mem_py)
  expect_identical(mem$levels, c("t0", "t1"))
  comp <- programme_composition(a, "true_state", weight = "clones")
  expect_identical(comp, threadfin:::.frame(tfpy()$tl$programme_composition(b, "true_state", weight = "clones")))
  # the simulated repertoire has no genes: its 20 embedding dimensions stand in for them (layer = NULL: X)
  pm <- quiet(programme_markers(a, layer = NULL, n_top = 3L, verbose = FALSE))
  expect_identical(pm, quiet(threadfin:::.frame(tfpy()$tl$programme_markers(b, layer = NULL, n_top = 3L,
                                                                        verbose = FALSE))))
  lg <- quiet(landscape_genes(a, layer = NULL, n_top = 5L, verbose = FALSE))
  expect_identical(lg, quiet(threadfin:::.frame(tfpy()$tl$landscape_genes(b, layer = NULL, n_top = 5L, verbose = FALSE))))
  gh <- quiet(gene_heritability(a, layer = NULL, n_perm = 20L, n_boot = 20L, verbose = FALSE))
  gh_py <- quiet(tfpy()$tl$gene_heritability(b, layer = NULL, n_perm = 20L, n_boot = 20L, verbose = FALSE))
  expect_identical(gh, threadfin:::.frame(gh_py))
  expect_identical(attr(gh, "index"), "gene")
  genes <- gh$gene
  sets <- list(first = genes[1:5], second = genes[6:12])
  gs <- quiet(geneset_heritability(a, sets))
  gs_py <- quiet(threadfin:::.frame(tfpy()$tl$geneset_heritability(gh_py, lapply(sets, as.list))))
  expect_identical(gs, gs_py)
  expect_error(geneset_heritability(tc(), sets), "gene_heritability")
  expect_error(geneset_heritability(a, list("x")), "named list")
})

test_that("export_lag passes single state names as lists", {
  skip_if_no_threadfin()
  a <- simulate_repertoire(n_clones = 400L, random_state = 0L)
  b <- a$copy()
  # stand-in columns on the simulated cells: a germinal-centre state and its output, mutation and isotype
  reticulate::py_run_string(paste(
    "import numpy as np",
    "def _tfr_states(a):",
    "    rng = np.random.default_rng(0)",
    "    a.obs['compartment'] = np.where(a.obs['true_state'].isin(['S0', 'S1', 'S2']), 'GC', 'ASC')",
    "    a.obs['mutation_frequency'] = rng.uniform(0, 0.1, a.n_obs)",
    "    a.obs['isotype'] = 'IGHG'",
    sep = "\n"
  ))
  main <- reticulate::import_main()
  main$`_tfr_states`(a)
  main$`_tfr_states`(b)
  lag <- export_lag(a, state_key = "compartment", source_states = "GC", output_states = "ASC")
  lag_py <- threadfin:::.result_value(tfpy()$tl$export_lag(b, state_key = "compartment", source_states = list("GC"),
                                                       output_states = list("ASC")))
  expect_s3_class(lag, "threadfin_export_lag")
  expect_identical(unclass(lag), lag_py)
  expect_identical(lag$settings$source_states, "GC")
  expect_length(lag$median_lag_ci, 2L)
  expect_s3_class(lag$pairs, "data.frame")
})

test_that("prepare_embedding builds the embedding without receptor genes", {
  skip_if_no_threadfin()
  files <- fixture("files")
  a <- as_anndata(files$paths[["h5ad"]])
  b <- as_anndata(files$paths[["h5ad"]])
  expect_identical(quiet(prepare_embedding(a, n_top_genes = 200L, n_comps = 10L, integrate = NULL, verbose = FALSE)), a)
  quiet(tfpy()$pp$prepare_embedding(b, n_top_genes = 200L, n_comps = 10L, integrate = NULL, verbose = FALSE))
  same <- reticulate::py_eval(paste0(
    "lambda a, b: bool((a.obsm['X_threadfin'] == b.obsm['X_threadfin']).all() and ",
    "a.var['threadfin_receptor_gene'].equals(b.var['threadfin_receptor_gene']))"
  ))
  expect_true(same(a, b))
  var <- a$var
  ig <- stats::setNames(var$threadfin_receptor_gene, rownames(var))
  expect_true(all(ig[c("Ighm", "Igkc", "Ighv1-72")]))
  expect_false(any(ig[c("Gene001", "Gene002")]))
})

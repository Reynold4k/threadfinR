# The analysis steps of threadfin.tl and threadfin.pp that the tutorials use, one R function each. Every step works
# on a Python AnnData (as_anndata(), simulate_repertoire()) or on the AnnData of a run_threadfin() result, and
# writes its results into that AnnData, as in Python: later steps read them from there (Tutorial 4).

#' Clone profiles
#'
#' Runs `threadfin.tl.clone_profiles()`, step 1 of Tutorial 4: one profile per clone with two or more cells, from
#' the cells' coordinates in an expression embedding, centred on their sampling context, summarised by the clone's
#' distribution of states (`representation = "kernel"`, as [run_threadfin()] uses) or its mean, and shrunk towards
#' the average profile the fewer cells the clone has. Clones are never filtered on reliability.
#'
#' @param x A Python AnnData ([as_anndata()], [simulate_repertoire()]) or a [run_threadfin()] result. The results
#'   are written into its AnnData (`uns["threadfin"]`), where the later steps find them.
#' @param ... Arguments of the Python function, by name; for `clone_profiles()`: `clone_key` (`"clone_id"`),
#'   `basis` (`"X_pca"`), `context_key` (the sampling context, e.g. the sample), `donor_key`, `representation`
#'   (`"mean"` or `"kernel"`), `n_features`, `n_components`, `bandwidth`, `smooth`, `min_cells` (2),
#'   `random_state`, `verbose`.
#' @return A data frame, one row per profiled clone: `clone_id`, `n_cells`, `reliability` and, when available,
#'   `donor` and `n_contexts`.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' profiles <- clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                            representation = "kernel")
#' head(profiles)
#' clonal_coherence(adata, strata_key = "context")
#' @export
clone_profiles <- function(x, ...) {
  .frame(.call_py(c("tl", "clone_profiles"), list(.adata_of(x)), .kwargs(list(...))))
}

#' Clonal coherence
#'
#' Runs `threadfin.tl.clonal_coherence()`, step 2 of Tutorial 4: are cells of one clone more alike than cells of
#' different clones? The share of the variance of cell states explained by clone identity (the clonal ICC), against
#' the same statistic with clone labels shuffled among the cells of each sampling unit (`strata_key`).
#'
#' @inheritParams clone_profiles
#' @param ... Arguments of `clonal_coherence()`: `strata_key` (the finest sampling unit: library, sort gate or
#'   tissue sample; default the context of the profiles), `n_perm` (200), `random_state`, `verbose`.
#' @return A list: `icc`, `null_mean`, `null_sd`, `null_q95`, `excess` (observed minus shuffled), `p_value`,
#'   `n_perm`, `strata_key`, `n_clones`, `n_cells`, the ICC of each embedding dimension (`icc_per_feature`) and the
#'   shuffled ICCs (`null`). `as.data.frame()` gives the single values as one row.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor")
#' coh <- clonal_coherence(adata, strata_key = "context", n_perm = 100)
#' coh$icc
#' as.data.frame(coh)
#' @export
clonal_coherence <- function(x, ...) {
  .dict_result(.call_py(c("tl", "clonal_coherence"), list(.adata_of(x)), .kwargs(list(...))), "threadfin_coherence")
}

#' The clone map
#'
#' Runs `threadfin.tl.embed_clones()`, step 3 of Tutorial 4: a two-dimensional layout of every profiled clone, one
#' point per clone, with UMAP settings chosen for the dataset so that clones that are connected in profile space
#' are not broken into islands (Tutorial 3). The map is a display: no statistic depends on it.
#'
#' @inheritParams clone_profiles
#' @param ... Arguments of `embed_clones()`: `n_neighbors` and `min_dist` (`"auto"`, or numbers to set them by
#'   hand), `spread`, `coherence_target`, `tune_kwargs`, `random_state`, `n_jobs` (processes for the automatic
#'   search; the map is the same for any number), `verbose`.
#' @return The clone table (as [clone_profiles()]) with the map coordinates `x` and `y`.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' \donttest{
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                representation = "kernel")
#' map <- embed_clones(adata, n_neighbors = 15, min_dist = 0.5)
#' head(map[, c("clone_id", "n_cells", "x", "y")])
#' }
#' @export
embed_clones <- function(x, ...) {
  .frame(.call_py(c("tl", "embed_clones"), list(.adata_of(x)), .kwargs(list(...))))
}

#' Clonal programmes
#'
#' Runs `threadfin.tl.find_programmes()`, step 4 of Tutorial 4: groups of clones with similar profiles (Leiden
#' communities of the clone graph), each split tested against a single Gaussian group, with the stability of every
#' programme from resampling the cells of each clone. Programme labels go to `obs$clone_programme`.
#'
#' @inheritParams clone_profiles
#' @param ... Arguments of `find_programmes()`: `resolution` (`"auto"`), `resolutions`, `n_neighbors`,
#'   `weight_by_reliability`, `n_boot` (30), `stability_threshold`, `min_programme_size`, `min_clones`,
#'   `test_splits`, `alpha`, `n_null`, `assign_remaining`, `min_posterior`, `random_state`, `n_jobs`, `key_added`,
#'   `verbose`.
#' @return A data frame, one row per programme: `programme`, `n_clones`, `n_cells`, `mean_reliability`,
#'   `mean_confidence`, `stability` (bootstrap Jaccard: >= 0.75 stable, < 0.6 do not interpret).
#' @family analysis steps
#' @examplesIf threadfin_available()
#' \donttest{
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                representation = "kernel")
#' find_programmes(adata)
#' }
#' @export
find_programmes <- function(x, ...) {
  .frame(.call_py(c("tl", "find_programmes"), list(.adata_of(x)), .kwargs(list(...))))
}

#' Does a label explain how clones differ?
#'
#' `profile_association()` runs `threadfin.tl.profile_association()`: the share of the variation between clone
#' profiles explained by a cell label collapsed to one value per clone (as in PERMANOVA for categories, a linear
#' trend for numbers), against the label shuffled among clones of the same donor and of similar size. It needs no
#' programmes. `association_test()` runs `threadfin.tl.association_test()`: which programmes differ in the label,
#' with clones as replicates (odds ratios within donor x clone-size strata for categories, rank effects for
#' numbers). `clone_labels()` runs `threadfin.tl.clone_labels()`: the label collapsed to one value per clone.
#'
#' A label that is the same for every cell of each sample (a tissue, sort gate or visit sequenced in its own
#' library) compares samples, not independent clones; the result says so (`sample_level`).
#'
#' @inheritParams clone_profiles
#' @param label_key The `obs` column (cell metadata column) to test.
#' @param ... Arguments of the Python function. `profile_association()`: `how` (`"auto"`, `"majority"`,
#'   `"mean"`, `"median"` or `"fraction:<level>"`, e.g. `"fraction:IGHA"`), `min_frac`, `strata_key`,
#'   `size_strata`, `sample_key`, `n_perm` (2000), `random_state`, `verbose`. `association_test()`: the same and
#'   `programme_key`. `clone_labels()`: `how`, `min_frac`, `clone_key`.
#' @return `profile_association()`: a list with `r2`, `null_mean`, `excess_r2`, `p_value`, `n_clones`, `kind`,
#'   the correlation of each level with the first profile axes (`axis_correlation`), the clone sizes behind the
#'   label (`clone_size`, a data frame) and the sample-level check (`sample_level`, `n_samples`);
#'   `as.data.frame()` gives the single values as one row. `association_test()`: a data frame, one row per
#'   programme and level (categories) or per programme (numbers), with FDR; its attribute `"attrs"` holds the
#'   strata and the sample-level check. `clone_labels()`: a vector named by clone id.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' \donttest{
#' adata <- simulate_repertoire(n_clones = 600, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                representation = "kernel")
#' pa <- profile_association(adata, "true_programme", n_perm = 500)
#' pa
#' head(clone_labels(adata, "true_state", how = "fraction:S1"))
#' find_programmes(adata)
#' head(association_test(adata, "true_programme", n_perm = 500))
#' }
#' @export
profile_association <- function(x, label_key, ...) {
  .dict_result(
    .call_py(c("tl", "profile_association"), list(.adata_of(x), label_key), .kwargs(list(...))),
    "threadfin_profile_association"
  )
}

#' @rdname profile_association
#' @export
association_test <- function(x, label_key, ...) {
  .frame(.call_py(c("tl", "association_test"), list(.adata_of(x), label_key), .kwargs(list(...))))
}

#' @rdname profile_association
#' @export
clone_labels <- function(x, label_key, ...) {
  .result_value(.call_py(c("tl", "clone_labels"), list(.adata_of(x), label_key), .kwargs(list(...))))
}

#' Do clones keep their state over time or across tissues?
#'
#' Runs `threadfin.tl.clonal_memory()`, step 6 of Tutorial 4: the memory index compares how much a clone changes
#' between its consecutive levels of `time_key` (time points, or tissues) with how much it differs from other clones
#' of its donor seen at the same two levels, with each snapshot's own sampling noise removed. 1: clones keep
#' everything that distinguishes them; 0: a clone is no more similar to its earlier self than to another clone.
#'
#' @inheritParams clone_profiles
#' @param time_key The `obs` column with the time point (or tissue) of each cell.
#' @param order The order of its levels; `NULL`: the order of a factor's levels (an ordered factor in R), else
#'   numbers in the labels (`d7` before `d28`), else alphabetical, with a warning.
#' @param ... Arguments of `clonal_memory()`: `reference`, `noise`, `programme_key`, `min_cells`, `min_pairs`,
#'   `n_null`, `n_boot`, `random_state`, `verbose`.
#' @return A list: `memory_index` with its bootstrap interval `memory_index_ci`, `p_value`, `n_clones`, `n_pairs`,
#'   `mean_change` and `mean_change_random`, the per-pair table `pairs` and, with programmes, `persistence`,
#'   `persistence_null` and the `transitions` tables. `as.data.frame()` gives the single values as one row.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' \donttest{
#' adata <- simulate_repertoire(n_clones = 400, n_timepoints = 2, memory = 0.8, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                representation = "kernel")
#' mem <- clonal_memory(adata, "timepoint", n_null = 100, n_boot = 100)
#' mem$memory_index
#' mem$memory_index_ci
#' }
#' @export
clonal_memory <- function(x, time_key, order = NULL, ...) {
  kw <- c(list(...), if (!is.null(order)) list(order = .names_arg(order)))
  .dict_result(.call_py(c("tl", "clonal_memory"), list(.adata_of(x), time_key), .kwargs(kw)), "threadfin_memory")
}

#' How much less mutated is a clone's output than its own germinal-centre cells?
#'
#' Runs `threadfin.tl.export_lag()`: for every clone captured both in the germinal centre (`source_states`) and in
#' its output (`output_states`, e.g. ASCs, antibody-secreting cells), the lag is the mean mutation of its
#' output cells minus that of its germinal-centre cells; the result is the median lag across clones with a
#' distribution-free 95 % interval and a sign test, and with `time_key` the lag in days at the germinal centre's
#' rate of maturation.
#'
#' @inheritParams clone_profiles
#' @param state_key The `obs` column with each cell's state.
#' @param source_states,output_states The levels of `state_key` that are the germinal centre and its output.
#' @param within `obs` columns that must match between the two compartments of a clone (the visit, the tissue).
#' @param ... Arguments of `export_lag()`: `clone_key`, `mutation_key` (`"mutation_frequency"`), `time_key`,
#'   `isotype_key`, `switched_only`, `min_cells`, `donor_key`, `clock`.
#' @return A list: `median_lag`, `median_lag_ci`, `sign_test_pvalue`, `n_clones`, the per-clone table `pairs`, and
#'   with `time_key` the maturation slope and the lag in days; with `donor_key`, `per_donor`.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' # lag <- export_lag(adata, state_key = "cell_state", source_states = "GC",
#' #                   output_states = "ASC", within = "timepoint", time_key = "timepoint")
#' # lag$median_lag; lag$median_lag_ci
#' @export
export_lag <- function(x, state_key, source_states, output_states, within = character(), ...) {
  kw <- c(list(state_key = state_key, source_states = .names_arg(source_states),
               output_states = .names_arg(output_states), within = .names_arg(within)), list(...))
  .dict_result(.call_py(c("tl", "export_lag"), list(.adata_of(x)), .kwargs(kw)), "threadfin_export_lag")
}

#' Cell states in each programme
#'
#' Runs `threadfin.tl.programme_composition()`: the share of each cell state in each programme (rows sum to 1),
#' counting every cell once (`weight = "cells"`) or every clone once (`weight = "clones"`).
#'
#' @inheritParams clone_profiles
#' @param state_key The `obs` column with the cell states.
#' @param ... Arguments of `programme_composition()`: `key` (the programme column), `weight`, `clone_key`.
#' @return A data frame, one row per programme and one column per state.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' \donttest{
#' adata <- simulate_repertoire(n_clones = 300, random_state = 0)
#' clone_profiles(adata, basis = "X_pca", context_key = "context", donor_key = "donor",
#'                representation = "kernel")
#' find_programmes(adata)
#' programme_composition(adata, "true_state", weight = "clones")
#' }
#' @export
programme_composition <- function(x, state_key, ...) {
  .frame(.call_py(c("tl", "programme_composition"), list(.adata_of(x), state_key), .kwargs(list(...))))
}

#' Genes behind the clone landscape
#'
#' `programme_markers()` runs `threadfin.tl.programme_markers()`: genes that distinguish each programme, tested
#' with clones as replicates (with `state_key` and `state`, within one cell state). `landscape_genes()` runs
#' `threadfin.tl.landscape_genes()`: genes whose expression changes smoothly across the clones (Moran's I on the
#' clone graph), with no programmes needed. Both describe the landscape the profiles were built from; they are
#' not independent tests of it.
#'
#' @inheritParams clone_profiles
#' @param ... Arguments of the Python function: `layer` (`"log_norm"`, written when Threadfin builds the
#'   embedding), `context_key`, `min_frac_expressed`, `n_top`, `state_key`, `state`, `min_cells_in_state`,
#'   `verbose`; `programme_markers()` also `programme_key`, `landscape_genes()` also `n_neighbors`.
#' @return `programme_markers()`: a data frame `programme, gene, mean_diff, auc, pvalue, fdr, rank` (and `state`).
#'   `landscape_genes()`: a data frame `gene, morans_i, z, pvalue, fdr, corr_axis1, rank`.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' # after run_threadfin() on counts (which writes the "log_norm" layer):
#' # programme_markers(res)
#' # landscape_genes(res, n_top = 30)
#' @export
programme_markers <- function(x, ...) {
  .frame(.call_py(c("tl", "programme_markers"), list(.adata_of(x)), .kwargs(list(...))))
}

#' @rdname programme_markers
#' @export
landscape_genes <- function(x, ...) {
  .frame(.call_py(c("tl", "landscape_genes"), list(.adata_of(x)), .kwargs(list(...))))
}

#' Which genes are shared within clones?
#'
#' `gene_heritability()` runs `threadfin.tl.gene_heritability()`: the clonal ICC of every expressed gene, with a
#' within-stratum permutation p-value and a clone bootstrap. `geneset_heritability()` runs
#' `threadfin.tl.geneset_heritability()` on that table: are named gene sets more clonal than genes of similar
#' expression? The table and its clone bootstrap are taken from the AnnData, where `gene_heritability()` stored
#' them.
#'
#' @inheritParams clone_profiles
#' @param genes Genes to test (default: those expressed in at least `min_frac_expressed` of the cells).
#' @param gene_sets A named list of gene vectors, e.g. `list(GC = c("Bcl6", "Aicda"), ASC = c("Prdm1", "Xbp1"))`.
#' @param ... Arguments of the Python function. `gene_heritability()`: `clone_key`, `context_key`, `strata_key`,
#'   `layer` (`"log_norm"`), `min_frac_expressed`, `exclude_receptor_genes`, `n_perm`, `n_boot`, `chunk_size`,
#'   `random_state`, `verbose`. `geneset_heritability()`: `n_bins`, `seed`.
#' @return `gene_heritability()`: a data frame, one row per gene (`gene`, `icc`, `pvalue`, `fdr`, `excess_icc`,
#'   ...), sorted by `excess_icc`. `geneset_heritability()`: a data frame, one row per gene set, with the
#'   difference from the expression-matched background, its standard error, interval, p-value and FDR.
#' @family analysis steps
#' @examplesIf threadfin_available()
#' # after run_threadfin() on counts:
#' # genes <- gene_heritability(res)
#' # geneset_heritability(res, list(GC = c("Bcl6", "Aicda", "Mki67"),
#' #                                ASC = c("Prdm1", "Xbp1", "Jchain")))
#' @export
gene_heritability <- function(x, genes = NULL, ...) {
  kw <- c(list(...), if (!is.null(genes)) list(genes = .names_arg(genes)))
  .frame(.call_py(c("tl", "gene_heritability"), list(.adata_of(x)), .kwargs(kw)))
}

#' @rdname gene_heritability
#' @export
geneset_heritability <- function(x, gene_sets, ...) {
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || !all(nzchar(names(gene_sets)))) {
    stop("gene_sets must be a named list of gene vectors", call. = FALSE)
  }
  adata <- .adata_of(x)
  uns <- .glue()$item(.py_get(adata, "uns"), "threadfin")
  table <- tryCatch(.glue()$item(uns, "gene_heritability"), error = function(e) NULL)
  if (is.null(table)) stop("run gene_heritability() on this AnnData first", call. = FALSE)
  boot <- tryCatch(.glue()$item(uns, "gene_heritability_bootstrap"), error = function(e) NULL)
  sets <- reticulate::py_dict(names(gene_sets), lapply(gene_sets, .names_arg), convert = FALSE)
  .frame(.call_py(c("tl", "geneset_heritability"), list(table, sets), .kwargs(c(list(bootstrap = boot), list(...)))))
}

#' An expression embedding without the receptor genes
#'
#' Runs `threadfin.pp.prepare_embedding()` (Tutorial 5, step 8): counts normalised to 10,000 per cell and
#' log-transformed (kept as `layers["log_norm"]`), the most variable genes without immunoglobulin and T-cell
#' receptor genes, scaling, PCA and optionally Harmony over technical batches, stored as `obsm["X_threadfin"]`.
#' [run_threadfin()] does this itself when no `basis` is given.
#'
#' @inheritParams clone_profiles
#' @param ... Arguments of `prepare_embedding()`: `counts_layer`, `batch_key`, `n_top_genes` (3000), `n_comps` (30),
#'   `exclude_receptor_genes`, `integrate` (`"harmony"` or `NULL`), `key_added`, `random_state`, `verbose`.
#' @return The AnnData, invisibly (changed in place).
#' @family analysis steps
#' @examplesIf threadfin_available()
#' # adata <- as_anndata(seurat_obj)
#' # prepare_embedding(adata, n_top_genes = 3000, n_comps = 30)
#' @export
prepare_embedding <- function(x, ...) {
  adata <- .adata_of(x)
  .call_py(c("pp", "prepare_embedding"), list(adata), .kwargs(list(...)))
  invisible(if (inherits(x, "threadfin_result")) adata else x)
}

#' Call clones with Dandelion
#'
#' Runs `threadfin.define_clones_dandelion()`, Threadfin's clone definition (Tutorial 2, section 6): cells belong
#' to one clone when they use the same V and J genes on both chains (heavy and light; beta and alpha), with
#' junctions of the same length and at least 85 % identical in amino acids for B cells, identical for T cells.
#' Called within each donor (`donor_key`), so no clone spans two individuals. Clones are called under Python's
#' string-hash seed 0 (in a separate Python process when this one was started without `PYTHONHASHSEED=0`), so the
#' clone ids are the same in every session.
#'
#' @param contigs The contig-level receptor table (one row per chain): a data frame with AIRR or Cell Ranger
#'   columns (e.g. from [read_contigs()]), a file path (read with [read_contigs()]), or a pandas table.
#' @param ... Arguments of `define_clones_dandelion()`: `donor_key` (a contig column naming the donor),
#'   `use_light_chain`, `identity`, `cell_key`, `productive_only`, `require_heavy_chain`, `check_contigs`, `out_col`,
#'   `verbose`, `hash_seed`.
#' @return A character vector of clone ids named by cell (missing for cells without a heavy chain); its attribute
#'   `"attrs"` records the settings, the clone counts and the hash seed used. Pass it to [threadfin_methods()] to
#'   describe the clone calling.
#' @seealso [read_contigs()], [run_threadfin()] (which calls clones itself when given `bcr`)
#' @examplesIf threadfin_available()
#' # contigs <- read_contigs("filtered_contig_annotations.csv")
#' # contigs$donor <- donor_of_cell[contigs$cell_id]
#' # ids <- define_clones_dandelion(contigs, donor_key = "donor")
#' # seurat_obj$clone_id <- unname(ids[colnames(seurat_obj)])
#' @export
define_clones_dandelion <- function(contigs, ...) {
  if (is.character(contigs)) {
    contigs <- .call_py("read_contigs", list(.bcr_arg(contigs)))
  } else if (is.data.frame(contigs)) {
    contigs <- .contigs_to_py(contigs)
  }
  .result_value(.call_py("define_clones_dandelion", list(contigs), .kwargs(list(...))))
}

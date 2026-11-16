# threadfin 0.1.0

* First version: the R interface to the Python package Threadfin (>= 4.6).
* `install_threadfin()` and `threadfin_available()` set up and check the Python side; `tf` is the Python package
  itself.
* `run_threadfin()` runs `threadfin.run()` on an h5ad file, a Seurat object (v4 and v5), a SingleCellExperiment or
  an AnnData, with Cell Ranger or AIRR contigs as files or data frames, and returns its tables as data frames.
* The steps of Threadfin's tutorials: `as_anndata()`, `read_contigs()`, `define_clones_dandelion()`,
  `prepare_embedding()`, `clone_profiles()`, `clonal_coherence()`, `embed_clones()`, `find_programmes()`,
  `profile_association()`, `association_test()`, `clone_labels()`, `clonal_memory()`, `export_lag()`,
  `programme_composition()`, `programme_markers()`, `landscape_genes()`, `gene_heritability()`,
  `geneset_heritability()` and `simulate_repertoire()`.
* Figures: `plot_overview()`, `plot_clone_map()`, `plot_label_effects()`, `plot_association()`, `plot_memory()`
  and `plot_threadfin()` for any `threadfin.pl` function; `save_threadfin()` and `read_threadfin()`.

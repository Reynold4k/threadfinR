test_that("tf reaches every Python function, with R numbers coerced", {
  skip_if_no_threadfin()
  expect_true(all(c("run", "tl", "pp", "pl", "sim", "read_contigs", "ThreadfinResult") %in% names(tf)))
  expect_true(all(c("clone_profiles", "embed_clones", "tune_clone_map", "map_quality") %in% names(tf$tl)))
  expect_match(tf$`__version__`, "^4\\.")
  a <- tf$sim$simulate_repertoire(n_clones = 600, random_state = 0)     # R doubles where Python takes int
  b <- fixture("sim")
  same <- reticulate::py_eval("lambda a, b: bool((a.obsm['X_pca'] == b.obsm['X_pca']).all() and a.obs.equals(b.obs))")
  expect_true(same(a, b))
})

test_that("arguments take the types of the Python signature", {
  skip_if_no_threadfin()
  reticulate::py_run_string(paste(
    "from typing import Literal, Optional, Sequence",
    "def _tfr_f(n: int, x: float, seed: 'int | None' = None, k=3, names: Sequence[str] = ('a',),",
    "           which: 'str | Sequence[str]' = 'a', cols: 'Sequence[str] | Literal[\"all\"] | None' = None,",
    "           m: Optional[int] = None, z: 'float | None' = None, bins: 'bool | tuple[int, ...] | None' = None,",
    "           within=(), sets: 'dict[str, list[str]]' = None, auto: 'int | str' = 'auto'):",
    "    return [type(v).__name__ for v in (n, x, seed, k, names, which, cols, m, z, bins, within)] + \\",
    "           [type(sets['s']).__name__, type(auto).__name__, type(bins[0]).__name__]",
    "def _tfr_apply(f, out):",
    "    return f(*out[0], **out[1])",
    sep = "\n"
  ))
  main <- reticulate::import_main(convert = FALSE)
  out <- threadfin:::.glue()$coerce_call(
    main$`_tfr_f`, list(),
    reticulate::dict(n = 5, x = 2, seed = 0, k = 4, names = "b", which = "c", cols = "all", m = 7, z = 1,
                     bins = list(2, 5, 20), within = "visit", sets = reticulate::dict(s = "Bcl6"), auto = 30)
  )
  types <- reticulate::py_to_r(main$`_tfr_apply`(main$`_tfr_f`, out))
  expect_identical(unlist(types), c("int", "float", "int", "int", "list", "str", "str", "int", "float", "list",
                                    "list", "list", "int", "int"))
})

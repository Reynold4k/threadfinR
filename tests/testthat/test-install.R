test_that("pip requirements are built from extras and versions", {
  spec <- threadfin:::.pip_spec
  expect_equal(spec(), "threadfin[dandelion]>=4.6")
  expect_equal(spec(NULL), "threadfin>=4.6")
  expect_equal(spec("dandelion", "4.6.1"), "threadfin[dandelion]==4.6.1")
  expect_equal(spec(c("dandelion", "data"), NULL), "threadfin[dandelion,data]")
  expect_equal(spec(NULL, ">= 4.6, < 5"), "threadfin>=4.6,<5")
  expect_error(spec("dandelion; rm -rf"), "extras")
  expect_error(spec("dandelion", c(">=4.6", "<5")), "one requirement")
})

test_that("install_threadfin validates its method", {
  expect_error(install_threadfin(method = "pip"), "should be one of")
})

test_that("threadfin_available answers TRUE or FALSE", {
  ok <- threadfin_available()
  expect_type(ok, "logical")
  expect_length(ok, 1L)
})

test_that("the missing-Python message names install_threadfin()", {
  msg <- threadfin:::.no_threadfin_message("No module named 'threadfin'")
  expect_match(msg, "install_threadfin()", fixed = TRUE)
  expect_match(msg, "No module named 'threadfin'", fixed = TRUE)
  expect_identical(threadfin:::.no_threadfin_message(msg), msg)
})

test_that("results without their Python object are refused clearly", {
  expect_error(py_object(list(1)), "not a Threadfin result")
  expect_error(threadfin_methods(list()), "Threadfin result")
  expect_error(save_threadfin(data.frame()), "Threadfin result")
  expect_error(threadfin:::.adata_of(data.frame()), "Python AnnData")
})

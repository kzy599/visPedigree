library(testthat)
library(visPedigree)

test_that("pedmat default behavior is method='A'", {
  tped <- tidyped(small_ped)
  
  # Default call
  res <- pedmat(tped)
  
  # Check class and attributes
  expect_true(inherits(res, "Matrix") || inherits(res, "pedmat"))
  expect_equal(attr(res, "method"), "A")
  
  # Check dimensions
  n <- nrow(tped)
  expect_equal(nrow(res), n)
  expect_equal(ncol(res), n)
  
  # Check content (diagonal should be >= 1)
  d <- Matrix::diag(res)
  expect_true(all(d >= 1))
})

test_that("pedmat handles method='f' correctly", {
  tped <- tidyped(small_ped, inbreed = TRUE)
  
  # Calculate f using pedmat
  f_vec <- pedmat(tped, method = "f")
  
  # Check type
  expect_true(is.numeric(f_vec))
  expect_false(is.matrix(f_vec))
  expect_equal(length(f_vec), nrow(tped))
  
  # Compare with tidyped inbreed column
  # Sort both to ensure matching order
  expected_f <- tped$f
  names(expected_f) <- tped$Ind
  
  # Align 
  f_vec_sorted <- f_vec[tped$Ind]
  expect_equal(f_vec_sorted, expected_f, tolerance = 1e-8)
})

test_that("pedmat supports compact mode", {
  # Use small_ped but ensure it has siblings to merge
  tped <- tidyped(small_ped)
  
  # Run compact mode
  res_compact <- pedmat(tped, method = "A", compact = TRUE)
  
  # Check attributes
  ci <- attr(res_compact, "call_info")
  expect_true(ci$compact)
  
  # If there are siblings, n_compact < n_original
  # small_ped has siblings E and P? Wait, lets check small_ped structure or structure of result
  if (ci$n_compact < ci$n_original) {
    expect_lt(nrow(res_compact), nrow(tped))
  }
  
  # Check attributes existence
  expect_false(is.null(attr(res_compact, "compact_map")))
  expect_false(is.null(attr(res_compact, "family_summary")))
  
  # Expand back
  res_expanded <- expand_pedmat(res_compact)
  expect_equal(nrow(res_expanded), nrow(tped))
})

test_that("pedmat argument threads is accepted", {
  tped <- tidyped(small_ped)
  # Just check it runs without error
  expect_error(pedmat(tped, threads = 1), NA)
})

test_that("vismat runs with defaults", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped)
  
  # Should run without error and return a trellis object (lattice)
  # We can't easily check the plot content in unit tests but can check object class
  p <- vismat(A, showgraph = FALSE) # Assuming vismat returns the plot object invisibly or visibly
  
  # Lattice plots are class "trellis"
  expect_s3_class(p, "trellis")
})

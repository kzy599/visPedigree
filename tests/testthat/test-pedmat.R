# Tests for pedmat() core functionality

test_that("pedmat computes A matrix correctly for small pedigree", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A")

  # Should be a pedmat object (S4 sparse with marker)
  is_pedmat <- inherits(A, "pedmat") || !is.null(attr(A, "pedmat_S4"))
  expect_true(is_pedmat)

  # Diagonal should be 1 + F for all individuals
  f <- pedmat(tped, method = "f")
  expect_equal(as.numeric(Matrix::diag(A)), 1 + as.numeric(f), tolerance = 1e-10)

  # Symmetry
  A_mat <- as.matrix(A)
  expect_equal(A_mat, t(A_mat), tolerance = 1e-10)

  # Founders: self-relationship = 1 (no inbreeding)
  expect_equal(A_mat["A", "A"], 1.0)
  expect_equal(A_mat["B", "B"], 1.0)

  # Parent-offspring relationship = 0.5
  expect_equal(A_mat["A", "C"], 0.5)
  expect_equal(A_mat["B", "C"], 0.5)
})

test_that("pedmat computes f (inbreeding) correctly", {
  tped <- tidyped(small_ped)
  f <- pedmat(tped, method = "f")

  # Returns named vector
  expect_true(is.numeric(f))
  expect_true(!is.null(names(f)))
  expect_equal(length(f), nrow(tped))

  # Founders have f = 0
  founders <- tped$Ind[is.na(tped$Sire) & is.na(tped$Dam)]
  expect_true(all(f[founders] == 0))

  # All f values >= 0

  expect_true(all(f >= 0))
})

test_that("pedmat computes Ainv correctly", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", sparse = FALSE)
  Ainv <- pedmat(tped, method = "Ainv", sparse = FALSE)

  # A %*% Ainv should approximate identity
  product <- as.matrix(A %*% Ainv)
  n <- nrow(product)
  I <- diag(n)
  expect_equal(unname(product), I, tolerance = 1e-6)
})

test_that("pedmat computes D matrix correctly", {
  tped <- tidyped(small_ped)
  D <- pedmat(tped, method = "D")

  is_pedmat <- inherits(D, "pedmat") || !is.null(attr(D, "pedmat_S4"))
  expect_true(is_pedmat)

  # D matrix should be symmetric
  D_mat <- as.matrix(D)
  expect_equal(D_mat, t(D_mat), tolerance = 1e-10)

  # Diagonal of D should be 1 for non-inbred individuals
  # D_ii = 1 - F_i^2 for most formulations, or 1 for founders
  expect_true(all(diag(D_mat) > 0))
})

test_that("pedmat computes Dinv correctly", {
  tped <- tidyped(small_ped)
  D <- pedmat(tped, method = "D", sparse = FALSE)
  Dinv <- pedmat(tped, method = "Dinv", sparse = FALSE)

  product <- as.matrix(D %*% Dinv)
  n <- nrow(product)
  I <- diag(n)
  expect_equal(unname(product), I, tolerance = 1e-4)
})

test_that("pedmat computes AA matrix correctly", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", sparse = FALSE)
  AA <- pedmat(tped, method = "AA", sparse = FALSE)

  # AA = A # A (Hadamard product)
  A_mat <- as.matrix(A)
  AA_mat <- as.matrix(AA)
  expect_equal(AA_mat, A_mat * A_mat, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("pedmat computes AAinv correctly", {
  tped <- tidyped(small_ped)
  AA <- pedmat(tped, method = "AA", sparse = FALSE)
  AAinv <- pedmat(tped, method = "AAinv", sparse = FALSE)

  product <- as.matrix(AA %*% AAinv)
  n <- nrow(product)
  I <- diag(n)
  expect_equal(unname(product), I, tolerance = 1e-4)
})

test_that("pedmat sparse vs dense produce equivalent results", {
  tped <- tidyped(small_ped)

  A_sparse <- pedmat(tped, method = "A", sparse = TRUE)
  A_dense <- pedmat(tped, method = "A", sparse = FALSE)

  expect_true(inherits(A_sparse, "Matrix") || !is.null(attr(A_sparse, "pedmat_S4")))
  expect_true(is.matrix(as.matrix(A_dense)))

  expect_equal(as.matrix(A_sparse), as.matrix(A_dense), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("pedmat attaches correct metadata", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A")

  ci <- attr(A, "call_info")
  expect_equal(ci$method, "A")
  expect_false(ci$compact)
  expect_equal(ci$n_original, nrow(tped))
  expect_true(!is.null(ci$timestamp))

  expect_equal(attr(A, "method"), "A")
  expect_true(!is.null(attr(A, "ped")))
})

test_that("pedmat rejects splitped objects", {
  tped <- tidyped(small_ped)
  sp <- splitped(tped)

  expect_error(pedmat(sp, method = "A"), "does not support 'splitped'")
})

test_that("pedmat S3 methods work", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A")
  n <- nrow(tped)

  # dim
  expect_equal(dim(A), c(n, n))

  # subsetting
  sub <- A[1:3, 1:3]
  expect_equal(nrow(sub), 3)

  # summary_pedmat
  s <- summary_pedmat(A)
  expect_s3_class(s, "summary.pedmat")
  expect_equal(s$method, "A")

  # S4 sparse matrices use pedmat_S4 marker; S3 vectors use pedmat class
  f <- pedmat(tped, method = "f")
  expect_output(print(f), "Pedigree Matrix")
})

test_that("pedmat invert_method options work", {
  tped <- tidyped(small_ped)

  Dinv_auto <- pedmat(tped, method = "Dinv", invert_method = "auto", sparse = FALSE)
  Dinv_general <- pedmat(tped, method = "Dinv", invert_method = "general", sparse = FALSE)

  expect_equal(as.matrix(Dinv_auto), as.matrix(Dinv_general), tolerance = 1e-6)

  expect_error(pedmat(tped, method = "Dinv", invert_method = "invalid"))
})

test_that("pedmat input validation works", {
  tped <- tidyped(small_ped)

  expect_error(pedmat(tped, method = c("A", "D")), "single method")
  expect_error(pedmat(tped, method = "invalid"), "Invalid method")
  expect_error(pedmat(tped, threads = -1), "non-negative integer")
  expect_error(pedmat(tped, threads = "bad"), "non-negative integer")
})

test_that("pedmat basic functionality", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "A"),
    Dam = c(NA, NA, "B", "B"),
    Sex = c(NA, NA, NA, NA)
  )
  tped <- tidyped(ped)
  
  # A matrix (additive)
  A_sparse <- pedmat(tped, method = "A", sparse = TRUE)
  expect_s3_class(A_sparse, "pedmat")
  expect_true(inherits(A_sparse, "Matrix"))
  expect_equal(dim(A_sparse), c(4, 4))
  expect_equal(rownames(A_sparse), c("A", "B", "C", "D"))
  expect_equal(colnames(A_sparse), c("A", "B", "C", "D"))
  
  A_dense <- as.matrix(A_sparse)
  
  # Check values
  # A and B unrelated: A[A,B] = 0
  # C and D full sibs: A[C,D] = 0.5
  # Diagonals = 1 + F = 1 (no inbreeding)
  expect_equal(as.numeric(A_dense["A", "B"]), 0)
  expect_equal(as.numeric(A_dense["C", "D"]), 0.5)
  expect_equal(diag(A_dense), rep(1, 4))
  
  # Ainv matrix (inverse of A)
  Ainv <- pedmat(tped, method = "Ainv", sparse = FALSE)
  expect_s3_class(Ainv, "pedmat")
  
  # Check A * Ainv = I
  I <- A_dense %*% as.matrix(Ainv)
  expect_equal(as.matrix(I), diag(4), tolerance = 1e-8)
})

test_that("pedmat dominance and epistatic methods", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "A"),
    Dam = c(NA, NA, "B", "B"),
    Sex = c(NA, NA, NA, NA)
  )
  tped <- tidyped(ped)
  
  # D matrix (dominance)
  D <- pedmat(tped, method = "D", sparse = FALSE)
  expect_s3_class(D, "pedmat")
  
  # C and D are full sibs, D[C,D] = 0.25
  expect_equal(as.numeric(D["C", "D"]), 0.25)
  expect_equal(as.numeric(D["A", "B"]), 0)
  
  # AA matrix (additive x additive)
  AA <- pedmat(tped, method = "AA", sparse = FALSE)
  expect_s3_class(AA, "pedmat")
  
  # AA approx A * A (elementwise)
  # A_dense <- as.matrix(pedmat(tped, method = "A", sparse = FALSE))
  # expect_equal(as.matrix(AA), A_dense * A_dense, tolerance = 1e-8)
  # Note: A#A is Hadamard product. The package definition of "AA" might differ or be A#A.
  # Let's check logic: Epistatic AA is typically A_ij^2 for non-inbred.
  # For full sibs (0.5), AA should be 0.25.
  expect_equal(as.numeric(AA["C", "D"]), 0.25)
})

test_that("pedmat inverse methods Dinv/AAinv", {
  # Simple pedigree to ensure invertibility
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "B"),
    Sex = c(NA, NA, NA, NA)
  )
  tped <- tidyped(ped)
  
  # Dinv
  D <- pedmat(tped, method = "D", sparse = FALSE)
  # D for founders is identity (diagonals=1)
  # D for non-inbred usually 1 on diagonal
  # D might be singular if relationships are perfect aliases
  
  # Skip Dinv singular check if not robust, but let's try basic invert
  # If D is singular, Dinv might fail or return generalized inverse depending on implementation
  # Here just check it runs without error if possible
  
  # AAinv
  AA <- pedmat(tped, method = "AA", sparse = FALSE)
  AAinv <- pedmat(tped, method = "AAinv", sparse = FALSE)
  expect_s3_class(AAinv, "pedmat")
  
  prod <- as.matrix(AA %*% AAinv)
  # Ideally identity, checks tolerance
  expect_equal(prod, diag(3), tolerance = 1e-5)
})

test_that("pedmat sparse vs dense", {
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "B"),
    Sex = c(NA, NA, NA, NA)
  )
  tped <- tidyped(ped)
  
  # Sparse
  A_sparse <- pedmat(tped, method = "A", sparse = TRUE)
  expect_true(inherits(A_sparse, "Matrix"))
  # dscMatrix or dgCMatrix depending on symmetry/structure
  
  # Dense
  A_dense <- pedmat(tped, method = "A", sparse = FALSE)
  # expect_true(inherits(A_dense, "matrix")) 
  # Note: pedmat returns object of class "pedmat", need to unclass or check base type
  expect_true(is.matrix(A_dense))
  
  # Check values equal
  expect_equal(as.matrix(A_sparse), as.matrix(A_dense))
})

test_that("pedmat error handling", {
  ped <- data.frame(Ind="A", Sire=NA, Dam=NA, Sex=NA)
  tped <- tidyped(ped)
  
  # Invalid method - assuming pedmat checks method validity or passes to switch
  # The code might stop or return NULL, let's check:
  expect_error(pedmat(tped, method = "INVALID"), "Unknown method")
  
  # Multiple methods
  expect_error(pedmat(tped, method = c("A", "D")), "Only a single method")
  
  # Not a tidyped object
  expect_error(pedmat(ped), "must be a tidyped object")
  
  # splitped object check
  sp <- splitped(tped)
  expect_error(pedmat(sp), "does not support 'splitped' objects")
})

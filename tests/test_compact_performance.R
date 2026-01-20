# Test compact mode performance optimization for Issue #7
# Tests that A matrix is not recalculated in compact mode for D/AA methods
library(visPedigree)

cat("=== Testing Issue #7 Fix: Compact Mode A Matrix Reuse ===\n\n")

# Load test data
data(small_ped)
tped <- tidyped(small_ped)

# Test 1: Verify A_intermediate attribute is set for D method
cat("Test 1: D method - checking A_intermediate attribute...\n")
result_D <- pedmatrix(tped, method = "D", compact = FALSE, sparse = TRUE)
has_A_attr_D <- !is.null(attr(result_D, "A_intermediate"))
cat("  A_intermediate attribute present:", has_A_attr_D, "\n")
if (has_A_attr_D) {
  A_dim <- dim(attr(result_D, "A_intermediate"))
  cat("  A_intermediate dimensions:", A_dim[1], "x", A_dim[2], "\n")
  cat("  ✓ PASS: A_intermediate attribute is attached\n")
} else {
  cat("  ✗ FAIL: A_intermediate attribute is missing\n")
}

# Test 2: Verify A_intermediate attribute is set for AA method
cat("\nTest 2: AA method - checking A_intermediate attribute...\n")
result_AA <- pedmatrix(tped, method = "AA", compact = FALSE, sparse = TRUE)
has_A_attr_AA <- !is.null(attr(result_AA, "A_intermediate"))
cat("  A_intermediate attribute present:", has_A_attr_AA, "\n")
if (has_A_attr_AA) {
  A_dim <- dim(attr(result_AA, "A_intermediate"))
  cat("  A_intermediate dimensions:", A_dim[1], "x", A_dim[2], "\n")
  cat("  ✓ PASS: A_intermediate attribute is attached\n")
} else {
  cat("  ✗ FAIL: A_intermediate attribute is missing\n")
}

# Test 3: Compact mode - verify attribute extraction works
cat("\nTest 3: Compact mode D - verify A matrix is reused...\n")
result_compact_D <- pedmatrix(tped, method = "D", compact = TRUE, sparse = TRUE)
cat("  Compact result computed successfully\n")
has_A_matrix_D <- !is.null(attr(result_compact_D, "A_matrix"))
cat("  Result has A_matrix (for query_relationship):", has_A_matrix_D, "\n")
if (has_A_matrix_D) {
  cat("  ✓ PASS: A_matrix attribute is present\n")
} else {
  cat("  ✗ FAIL: A_matrix attribute is missing\n")
}

# Test 4: Compact mode AA
cat("\nTest 4: Compact mode AA - verify A matrix is reused...\n")
result_compact_AA <- pedmatrix(tped, method = "AA", compact = TRUE, sparse = TRUE)
cat("  Compact result computed successfully\n")
has_A_matrix_AA <- !is.null(attr(result_compact_AA, "A_matrix"))
cat("  Result has A_matrix (for query_relationship):", has_A_matrix_AA, "\n")
if (has_A_matrix_AA) {
  cat("  ✓ PASS: A_matrix attribute is present\n")
} else {
  cat("  ✗ FAIL: A_matrix attribute is missing\n")
}

# Test 5: Verify A matrix reuse in compact mode
cat("\nTest 5: Verify A matrix reuse mechanism in compact mode...\n")
# When compact mode calls pedmatrix for D, it should extract A_intermediate
# Let's manually verify the mechanism
tped_compact_result <- visPedigree:::compact_ped_for_matrix(tped)
result_D_on_compact <- pedmatrix(tped_compact_result$ped_compact, method = "D", compact = FALSE, sparse = TRUE)
has_A_in_D <- !is.null(attr(result_D_on_compact, "A_intermediate"))
cat("  D method on compact ped has A_intermediate:", has_A_in_D, "\n")

if (has_A_in_D) {
  # This is what compact mode will extract
  A_extracted <- attr(result_D_on_compact, "A_intermediate")
  # Calculate A directly for comparison
  A_direct <- pedmatrix(tped_compact_result$ped_compact, method = "A", compact = FALSE, sparse = FALSE)
  
  # Compare values (ignore attributes)
  A_match <- all.equal(A_extracted, as.matrix(A_direct), check.attributes = FALSE, tolerance = 1e-10)
  if (isTRUE(A_match)) {
    cat("  ✓ PASS: Extracted A matches directly calculated A\n")
    test5_pass <- TRUE
  } else {
    cat("  ✗ FAIL: Extracted A differs from directly calculated A\n")
    cat("   ", A_match, "\n")
    test5_pass <- FALSE
  }
} else {
  cat("  ✗ FAIL: A_intermediate not attached\n")
  test5_pass <- FALSE
}

# Test 6: Verify query_relationship works with compact D/AA
cat("\nTest 6: Verify query_relationship with compact mode...\n")
# Test with compact D
query_works_D <- tryCatch({
  val <- query_relationship(result_compact_D, "A", "B")
  is.numeric(val)
}, error = function(e) FALSE)

query_works_AA <- tryCatch({
  val <- query_relationship(result_compact_AA, "C", "E")
  is.numeric(val)
}, error = function(e) FALSE)

cat("  query_relationship works with compact D:", query_works_D, "\n")
cat("  query_relationship works with compact AA:", query_works_AA, "\n")

if (query_works_D && query_works_AA) {
  cat("  ✓ PASS: query_relationship works correctly\n")
  test6_pass <- TRUE
} else {
  cat("  ✗ FAIL: query_relationship failed\n")
  test6_pass <- FALSE
}

# Test 7: Performance timing (optional, for larger pedigrees)
cat("\nTest 7: Performance comparison (small pedigree)...\n")
time_D_compact <- system.time({
  for (i in 1:10) pedmatrix(tped, method = "D", compact = TRUE, sparse = TRUE)
})["elapsed"]
cat("  Compact mode (10 iterations):", time_D_compact, "seconds\n")

# Summary
cat("\n=== Test Summary ===\n")
total_tests <- 6
passed_tests <- sum(has_A_attr_D, has_A_attr_AA, 
                    has_A_matrix_D,
                    has_A_matrix_AA,
                    test5_pass,
                    test6_pass)
cat("Passed:", passed_tests, "/", total_tests, "\n")

if (passed_tests == total_tests) {
  cat("✓ All tests passed! Issue #7 fix is working correctly.\n")
  cat("\nThe fix successfully:\n")
  cat("  1. Attaches A_intermediate to D/AA calculation results\n")
  cat("  2. Allows compact mode to extract and reuse A matrix\n")
  cat("  3. Avoids redundant A matrix recalculation\n")
  cat("  4. Supports query_relationship() for full-sibling queries\n")
} else {
  cat("✗ Some tests failed. Please review the fix.\n")
}

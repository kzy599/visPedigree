test_that("pedmat threads argument works correctly", {
  # Load small data
  tped <- tidyped(small_ped)
  
  # 1. Invalid input handling
  expect_error(pedmat(tped, threads = "invalid"), "'threads' must be a single non-negative integer")
  expect_error(pedmat(tped, threads = -1), "'threads' must be a single non-negative integer")
  
  # 2. OpenMP availability check
  # We can't easily force OpenMP on/off, but we can check if it warns when threads > 0 
  # and OpenMP is reported as unavailable
  
  # Check if OpenMP is available (internal helper)
  has_openmp <- visPedigree:::cpp_openmp_available()
  
  if (!has_openmp) {
    # Expect warning if threads requested but not available
    expect_warning(
      pedmat(tped, method = "A", threads = 2),
      "OpenMP is not available"
    )
    
    # Expect NO warning if threads=0 or 1 (default/serial)
    # Note: threads=1 might still trigger code paths but usually doesn't warn if implementation is smart,
    # but the current implementation warns if threads > 0 and no OpenMP.
    # Actually code says: if (threads > 0) { if (!available) warning... }
    # So threads=1 SHOULD warn if OpenMP is missing according to current logic? 
    # Let's verify code:
    # if (threads > 0) {
    #   if (!cpp_openmp_available()) {
    #     warning(...)
    #   }
    # }
    expect_warning(
        pedmat(tped, method="A", threads = 1),
        "OpenMP is not available"
    )
    
    expect_silent(pedmat(tped, method="A", threads = 0))
    
  } else {
    # If OpenMP IS available
    expect_silent(pedmat(tped, method = "A", threads = 2))
    expect_silent(pedmat(tped, method = "A", threads = 0))
    
    # We can't easily test speedup on small_ped, but we can verify it runs
    res <- pedmat(tped, method = "A", threads = 2)
    # Check if it is a valid pedmat object (S3 class or S4 with marker)
    is_pedmat <- inherits(res, "pedmat") || !is.null(attr(res, "pedmat_S4"))
    expect_true(is_pedmat)
  }
})

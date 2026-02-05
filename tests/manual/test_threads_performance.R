
pkgload::load_all(".")
library(data.table)

# Load the large pedigree
load("data/big_family_size_ped.rda")

# Prepare tidyped
tped <- tidyped(big_family_size_ped)
print(paste("Pedigree size:", nrow(tped)))

# Function to run performance test
run_test <- function(n_threads) {
  cat(sprintf("\nRunning pedmat(method='Ainv') with threads = %d...\n", n_threads))
  
  # Measure time
  time_taken <- system.time({
    # We only test method="Ainv" as that's the one with the OpenMP implementation optimized
    # Note: method='A' uses Armadillo which might or might not use OpenMP internally depending on BLAS,
    # but our custom OpenMP code is in cpp_build_ainv_triplets called by 'Ainv'.
    res <- pedmat(tped, method = "Ainv", threads = n_threads)
  })
  
  print(time_taken)
  return(time_taken["elapsed"])
}

# Run with 1 thread
time_1 <- run_test(1)

# Run with 4 threads (or max available if less)
# Check available cores
n_cores <- parallel::detectCores()
n_use <- min(4, n_cores)

if (n_use > 1) {
    time_multi <- run_test(n_use)
    
    # Report speedup
    speedup <- time_1 / time_multi
    cat(sprintf("\nSpeedup: %.2fx (1 vs %d threads)\n", speedup, n_use))
    
    if (speedup > 1.1) {
      cat("SUCCESS: Multi-threading provided speedup.\n")
    } else {
      cat("WARNING: No significant speedup observed.\n")
    }
} else {
    cat("Skipping multi-thread test: Only 1 core detected.\n")
}

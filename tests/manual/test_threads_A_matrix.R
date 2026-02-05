
pkgload::load_all(".")
library(data.table)

# Load the deep pedigree (larger dense matrices amplify parallel benefits)
# Note: huge pedigrees (like big_family_size_ped) might be too slow for dense A matrix
# So we use deep_ped (4396 individuals) which creates a ~4400x4400 matrix (approx 19M elements)
load("data/deep_ped.rda")
tped <- tidyped(deep_ped)
print(paste("Pedigree size:", nrow(tped)))

# Function to run performance test
run_test_A <- function(n_threads) {
  cat(sprintf("\nRunning pedmat(method='A') with threads = %d...\n", n_threads))
  
  # Measure time
  time_taken <- system.time({
    # method='A' uses cpp_calculate_A which returns an Armadillo matrix.
    # Parallelization here depends on whether Armadillo/BLAS uses OpenMP.
    # The pure C++ implementation of A calculation in cpp_calculate_A seems serial in the provided code,
    # but let's verify if OpenMP threads setting affects it (e.g. via BLAS downstream or if I missed OMP implementation)
    res <- pedmat(tped, method = "A", threads = n_threads, sparse = FALSE) # Dense for max computation
  })
  
  print(time_taken)
  return(time_taken["elapsed"])
}

# Run with 1 thread
time_1 <- run_test_A(1)

# Run with 4 threads
n_cores <- parallel::detectCores()
n_use <- min(4, n_cores)

if (n_use > 1) {
    time_multi <- run_test_A(n_use)
    
    # Report speedup
    speedup <- time_1 / time_multi
    cat(sprintf("\nSpeedup: %.2fx (1 vs %d threads)\n", speedup, n_use))
    
    if (speedup > 1.1) {
      cat("SUCCESS: Multi-threading provided speedup.\n")
    } else {
      cat("WARNING: No significant speedup observed for A matrix.\n")
      cat("Possible reason: The calculation of A (cpp_calculate_A) might be serial code.\n") 
    }
}

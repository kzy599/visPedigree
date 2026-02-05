# Use installed package instead of pkgload to avoid dylib issues on macOS
library(visPedigree)
library(data.table)

# Load the deep pedigree (4399 individuals)
load("data/deep_ped.rda")
tped <- tidyped(deep_ped)
print(paste("Pedigree size:", nrow(tped)))

# Function to run performance test
run_test_D <- function(n_threads) {
  cat(sprintf("\nRunning pedmat(method='D') with threads = %d...\n", n_threads))
  
  # Measure time
  time_taken <- system.time({
    res <- pedmat(tped, method = "D", threads = n_threads, sparse = FALSE)
  })
  
  print(time_taken)
  return(time_taken["elapsed"])
}

# Run with 1 thread
time_1 <- run_test_D(1)

# Run with 4 threads
n_cores <- parallel::detectCores()
n_use <- min(4, n_cores)

if (n_use > 1) {
    time_multi <- run_test_D(n_use)
    
    # Report speedup
    speedup <- time_1 / time_multi
    cat(sprintf("\nSpeedup: %.2fx (1 vs %d threads)\n", speedup, n_use))
    
    if (speedup > 1.5) {
      cat("SUCCESS: Significant speedup achieved with multi-threading!\n")
    } else if (speedup > 1.1) {
      cat("MODERATE: Some speedup observed.\n")
    } else {
      cat("WARNING: No significant speedup observed.\n")
    }
} else {
    cat("Skipping multi-thread test: Only 1 core detected.\n")
}

library(visPedigree)
library(nadiv)
library(data.table)

# Test with different pedigrees
test_pedigrees <- list(
  small = small_ped,
  simple = simple_ped,
  deep = deep_ped
)

cat("=================================================================\n")
cat("Comparing visPedigree (parallel) vs nadiv for D matrix calculation\n")
cat("=================================================================\n\n")

for (ped_name in names(test_pedigrees)) {
  cat(sprintf("\n--- Testing with %s pedigree ---\n", ped_name))
  
  ped_data <- test_pedigrees[[ped_name]]
  
  # Prepare for visPedigree
  tped <- tidyped(ped_data)
  n <- nrow(tped)
  cat(sprintf("Pedigree size: %d individuals\n", n))
  
  # visPedigree: D matrix (1 thread)
  cat("\n[visPedigree] Calculating D matrix (1 thread)...\n")
  time_vis_1 <- system.time({
    D_vis_1 <- pedmat(tped, method = "D", threads = 1, sparse = FALSE)
  })
  cat(sprintf("  Time: %.3f seconds\n", time_vis_1["elapsed"]))
  
  # visPedigree: D matrix (4 threads)
  cat("\n[visPedigree] Calculating D matrix (4 threads)...\n")
  time_vis_4 <- system.time({
    D_vis_4 <- pedmat(tped, method = "D", threads = 4, sparse = FALSE)
  })
  cat(sprintf("  Time: %.3f seconds\n", time_vis_4["elapsed"]))
  cat(sprintf("  Speedup: %.2fx\n", time_vis_1["elapsed"] / time_vis_4["elapsed"]))
  
  # nadiv: D matrix
  cat("\n[nadiv] Calculating D matrix...\n")
  # Prepare pedigree for nadiv (needs numeric ID, Sire, Dam)
  ped_nadiv <- data.frame(
    ID = seq_len(n),
    Sire = match(tped$Sire, tped$Ind, nomatch = 0),
    Dam = match(tped$Dam, tped$Ind, nomatch = 0)
  )
  
  time_nadiv <- system.time({
    # nadiv uses makeD() function
    D_nadiv <- as.matrix(makeD(ped_nadiv)$D)
  })
  cat(sprintf("  Time: %.3f seconds\n", time_nadiv["elapsed"]))
  
  # Compare results
  cat("\n[Comparison]\n")
  
  # Remove dimnames for comparison (order might differ)
  D_vis_clean <- as.matrix(D_vis_4)
  dimnames(D_vis_clean) <- NULL
  dimnames(D_nadiv) <- NULL
  
  # Check dimensions
  if (all(dim(D_vis_clean) == dim(D_nadiv))) {
    cat(sprintf("  Dimensions match: %d x %d\n", nrow(D_vis_clean), ncol(D_vis_clean)))
  } else {
    cat("  WARNING: Dimensions do NOT match!\n")
    cat(sprintf("    visPedigree: %d x %d\n", nrow(D_vis_clean), ncol(D_vis_clean)))
    cat(sprintf("    nadiv: %d x %d\n", nrow(D_nadiv), ncol(D_nadiv)))
  }
  
  # Check if values match
  max_diff <- max(abs(D_vis_clean - D_nadiv))
  cat(sprintf("  Max absolute difference: %.2e\n", max_diff))
  
  if (max_diff < 1e-10) {
    cat("  ✓ Results match perfectly!\n")
  } else if (max_diff < 1e-6) {
    cat("  ✓ Results match within numerical tolerance\n")
  } else {
    cat("  ✗ Results differ significantly!\n")
  }
  
  # Performance summary
  cat("\n[Performance Summary]\n")
  cat(sprintf("  visPedigree (1 thread): %.3f s\n", time_vis_1["elapsed"]))
  cat(sprintf("  visPedigree (4 threads): %.3f s (%.2fx speedup)\n", 
              time_vis_4["elapsed"], time_vis_1["elapsed"] / time_vis_4["elapsed"]))
  cat(sprintf("  nadiv: %.3f s\n", time_nadiv["elapsed"]))
  cat(sprintf("  visPedigree vs nadiv: %.2fx %s\n", 
              time_nadiv["elapsed"] / time_vis_4["elapsed"],
              ifelse(time_vis_4["elapsed"] < time_nadiv["elapsed"], "faster", "slower")))
  
  cat("\n")
}

cat("=================================================================\n")
cat("All tests completed!\n")
cat("=================================================================\n")

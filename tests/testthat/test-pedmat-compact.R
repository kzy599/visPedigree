# Test script for compact matrix functionality
# Source this file to test the compact pedmat feature

library(visPedigree)
library(data.table)

cat("Testing Compact Pedigree Matrix\n")
cat(strrep("=", 60), "\n\n")

# Create a test pedigree with large full-sibling families
create_test_ped <- function() {
  n_sires <- 5
  n_dams <- 10
  offspring_per_family <- 50
  
  ped_list <- list()
  ind_counter <- 1
  
  # Create sires
  sires <- paste0("S", sprintf("%02d", 1:n_sires))
  ped_list[[length(ped_list) + 1]] <- data.table(
    Ind = sires,
    Sire = NA,
    Dam = NA,
    Sex = "male"
  )
  
  # Create dams
  dams <- paste0("D", sprintf("%02d", 1:n_dams))
  ped_list[[length(ped_list) + 1]] <- data.table(
    Ind = dams,
    Sire = NA,
    Dam = NA,
    Sex = "female"
  )
  
  # Create offspring (full siblings)
  for (s in 1:n_sires) {
    for (d in 1:n_dams) {
      offspring_ids <- paste0("O", sprintf("%04d", 
                              (ind_counter):(ind_counter + offspring_per_family - 1)))
      ind_counter <- ind_counter + offspring_per_family
      
      ped_list[[length(ped_list) + 1]] <- data.table(
        Ind = offspring_ids,
        Sire = sires[s],
        Dam = dams[d],
        Sex = sample(c("male", "female"), offspring_per_family, replace = TRUE)
      )
    }
  }
  
  ped <- rbindlist(ped_list)
  return(ped)
}

# Test 1: Create and tidy pedigree
cat("Test 1: Creating test pedigree...\n")
ped <- create_test_ped()
ped_tidy <- tidyped(ped, addnum = TRUE)
cat("  Original pedigree size:", nrow(ped_tidy), "individuals\n")
cat("  Founders:", sum(is.na(ped_tidy$Sire) & is.na(ped_tidy$Dam)), "\n")
cat("  Offspring:", sum(!is.na(ped_tidy$Sire) & !is.na(ped_tidy$Dam)), "\n\n")

# Test 2: Calculate inbreeding without compact
cat("Test 2: Calculate inbreeding WITHOUT compact...\n")
time_normal <- system.time({
  f_normal <- pedmat(ped_tidy, method = "f", compact = FALSE)
})
cat("  Time:", round(time_normal["elapsed"], 3), "seconds\n")
cat("  Result length:", length(f_normal), "\n")
cat("  Mean f:", round(mean(f_normal), 6), "\n\n")

# Test 3: Calculate inbreeding WITH compact
cat("Test 3: Calculate inbreeding WITH compact...\n")
time_compact <- system.time({
  f_compact <- pedmat(ped_tidy, method = "f", compact = TRUE)
})
cat("  Time:", round(time_compact["elapsed"], 3), "seconds\n\n")

print(f_compact)
cat("\n")

# Test 4: Verify correctness
cat("Test 4: Verify compact results...\n")
compact_map <- attr(f_compact, "compact_map")
family_summary <- attr(f_compact, "family_summary")
cat("  Compact map rows:", nrow(compact_map), "\n")
cat("  Family summary rows:", nrow(family_summary), "\n")
cat("  Result length:", length(f_compact), "\n")

# Check if full siblings have same f value
test_family <- family_summary[1]
family_inds <- compact_map[FamilyID == test_family$FamilyID, Ind]
if (length(family_inds) > 1) {
  rep_ind <- compact_map[FamilyID == test_family$FamilyID & IsRepresentative == TRUE, RepIndNum]
  f_values <- f_normal[family_inds]
  cat("  Testing family", test_family$FamilyID, "with", length(family_inds), "siblings:\n")
  cat("    All f values identical:", length(unique(f_values)) == 1, "\n")
  cat("    f value:", unique(f_values)[1], "\n")
  cat("    Compact result f:", f_compact[rep_ind], "\n")
  cat("    Values match:", all.equal(unique(f_values)[1], f_compact[rep_ind]), "\n")
}
cat("\n")

# Test 5: Query functions
cat("Test 5: Test query functions...\n")
ind1 <- compact_map[IsCompacted == TRUE][1, Ind]
ind2 <- compact_map[IsCompacted == TRUE][2, Ind]
cat("  Query individual '", ind1, "':\n", sep = "")
f_query <- query_relationship(f_compact, ind1)
cat("    Result:", f_query, "\n")
cat("    Matches original:", all.equal(f_normal[ind1], f_query), "\n\n")

# Test 6: Expand matrix
cat("Test 6: Test expand function...\n")
f_expanded <- expand_pedmat(f_compact)
cat("  Expanded length:", length(f_expanded), "\n")
cat("  Matches original length:", length(f_expanded) == length(f_normal), "\n")
cat("  Values match:", all.equal(as.numeric(f_normal), as.numeric(f_expanded)), "\n\n")

# Test 7: Summary
cat("Test 7: Summary method...\n")
summary(f_compact)
cat("\n")

# Test 8: A matrix (if pedigree not too large)
if (nrow(ped_tidy) < 1000) {
  cat("Test 8: Test A matrix with compact...\n")
  time_A <- system.time({
    A_compact <- pedmat(ped_tidy, method = "A", compact = TRUE, sparse = TRUE)
  })
  A_mat <- A_compact
  class(A_mat) <- setdiff(class(A_mat), "pedmat")
  A_compact_map <- attr(A_compact, "compact_map")
  cat("  Time:", round(time_A["elapsed"], 3), "seconds\n")
  cat("  Matrix dimensions:", nrow(A_mat), "x", ncol(A_mat), "\n")
  cat("  Matrix class:", class(A_mat)[1], "\n")
  
  # Verify diagonal (should be 1 + f)
  diag_values <- Matrix::diag(A_mat)
  rep_inds <- A_compact_map[IsRepresentative == TRUE, IndNum]
  rep_inds <- rep_inds[!is.na(rep_inds)]  # Remove NA entries
  expected_diag <- 1 + f_normal[A_compact_map[IsRepresentative == TRUE & !is.na(IndNum), Ind]]
  cat("  Diagonal check (1+f):", all.equal(as.numeric(diag_values), as.numeric(expected_diag)), "\n\n")
}

cat(strrep("=", 60), "\n")
cat("All tests completed!\n")

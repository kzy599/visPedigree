#!/usr/bin/env Rscript
# Integration test for pedanalysis functions on real pedigree data
library(data.table)
devtools::load_all(".")

# ============================================================
# Setup: Read CSV and build tidyped
# ============================================================
csv_path <- "inst/extdata/g2016G00to2025G09Ped.csv"
ped_raw <- fread(csv_path, colClasses = "character")

extract_year <- function(id) {
  m <- regmatches(id, regexpr("(\\d{2})[BG]", id))
  ifelse(length(m) > 0 & nchar(m) > 0,
         as.integer(substr(m, 1, 2)) + 2000L,
         NA_integer_)
}
ped_raw[, Year := sapply(AnimalID, extract_year)]
ped_raw[SireID == "0", SireID := NA_character_]
ped_raw[DamID == "0", DamID := NA_character_]
ped_raw <- unique(ped_raw, by = c("AnimalID", "SireID", "DamID"))

tped <- tidyped(ped_raw, addnum = TRUE, inbreed = TRUE, genmethod = "bottom")
tped[, Year := as.integer(Year)]

# Backfill Year from parents
year_lookup <- setNames(tped$Year, tped$Ind)
repeat {
  na_idx <- which(is.na(tped$Year))
  if (length(na_idx) == 0) break
  filled <- 0L
  for (i in na_idx) {
    sire_year <- if (!is.na(tped$Sire[i])) year_lookup[tped$Sire[i]] else NA_integer_
    dam_year <- if (!is.na(tped$Dam[i])) year_lookup[tped$Dam[i]] else NA_integer_
    parent_years <- c(sire_year, dam_year)
    parent_years <- parent_years[!is.na(parent_years)]
    if (length(parent_years) > 0) {
      new_year <- max(parent_years) + 1L
      tped$Year[i] <- new_year
      year_lookup[tped$Ind[i]] <- new_year
      filled <- filled + 1L
    }
  }
  if (filled == 0L) break
}

cat("Setup complete:", nrow(tped), "individuals,", sum(is.na(tped$Year)), "Year NAs\n\n")

# ============================================================
# Test 1: pedgenint(timevar = "Year")
# ============================================================
cat("=== Test 1: pedgenint(timevar='Year') ===\n")
tryCatch({
  gen_int <- pedgenint(tped, timevar = "Year")
  cat("Result:\n")
  print(gen_int)
  cat("\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ============================================================
# Test 2: pedne(by = "Year")
# ============================================================
cat("=== Test 2: pedne(by='Year') ===\n")
tryCatch({
  ne_result <- pedne(tped, by = "Year")
  cat("Result:\n")
  print(ne_result)
  cat("\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ============================================================
# Test 3: pedstats(timevar = "Year")
# ============================================================
cat("=== Test 3: pedstats(timevar='Year') ===\n")
tryCatch({
  stats_result <- pedstats(tped, timevar = "Year")
  print(stats_result)
  cat("\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ============================================================
# Test 4: pedcontrib(cohort = last gen)
# ============================================================
cat("=== Test 4: pedcontrib(cohort=last gen) ===\n")
tryCatch({
  max_gen <- max(tped$Gen)
  cohort_ids <- tped[Gen == max_gen, Ind]
  cat("Cohort size (Gen", max_gen, "):", length(cohort_ids), "\n")
  
  # Use a smaller cohort for speed (sample 500)
  if (length(cohort_ids) > 500) {
    set.seed(42)
    cohort_ids <- sample(cohort_ids, 500)
    cat("Subsampled to 500 for speed\n")
  }
  
  t0 <- proc.time()
  contrib_result <- pedcontrib(tped, cohort = cohort_ids, mode = "both", top = 20)
  elapsed <- (proc.time() - t0)["elapsed"]
  cat("Elapsed:", elapsed, "sec\n")
  print(contrib_result)
  cat("\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ============================================================
# Test 5: vispstat()
# ============================================================
cat("=== Test 5: vispstat() ===\n")
tryCatch({
  stats_result2 <- pedstats(tped, timevar = "Year")
  
  # Test genint plot
  p1 <- vispstat(stats_result2, type = "genint")
  cat("genint plot created successfully\n")
  
  # Test ecg plot
  p2 <- vispstat(stats_result2, type = "ecg")
  cat("ecg plot created successfully\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ============================================================
# Test 6: pedrel() with family sampling
# ============================================================
cat("=== Test 6: pedrel(by='Gen') with family sampling ===\n")
tryCatch({
  # Sample 1 individual per family
  set.seed(42)
  sampled <- tped[, .SD[sample(.N, 1)], by = Family]
  cat("Sampled", nrow(sampled), "from", uniqueN(tped$Family), "families\n")
  
  # Rebuild tidyped on sampled subset
  sampled_tped <- tidyped(sampled[, .(Ind, Sire, Dam)], addnum = TRUE, genmethod = "bottom")
  cat("Sampled tidyped:", nrow(sampled_tped), "individuals\n")
  
  # Check max group size
  group_sizes <- sampled_tped[, .N, by = Gen]
  cat("Max group size:", max(group_sizes$N), "\n")
  
  if (max(group_sizes$N) <= 25000) {
    rel_result <- pedrel(sampled_tped, by = "Gen")
    cat("Result:\n")
    print(rel_result)
  } else {
    cat("Groups still too large for dense A matrix (>25K). Using sample param.\n")
    rel_result <- pedrel(sampled_tped, by = "Gen", sample = 5000)
    cat("Result:\n")
    print(rel_result)
  }
  cat("\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

cat("=== All tests completed ===\n")

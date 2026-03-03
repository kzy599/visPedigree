# tests/test_G09_integration.R
# G09_prunedPed.csv 集成测试: 27,635 条记录，含完整 Sex/Year/Breed/FamilyID
# 运行方式: devtools::load_all("."); source("tests/test_G09_integration.R")

suppressMessages({
  devtools::load_all(".")
  library(data.table)
})

csv_path <- "inst/extdata/G09_prunedPed.csv"
cat("========================================\n")
cat("G09 Integration Test\n")
cat("========================================\n\n")

# ── Step 0: Read and prepare data ────────────────────────────────
cat("=== Step 0: Data preparation ===\n")
raw <- fread(csv_path)
cat("Raw rows:", nrow(raw), "\n")

# Rename SexID -> Sex (tidyped expects "Sex")
setnames(raw, "SexID", "Sex")

# Convert Sire/Dam "0" -> NA
raw[SireID == "0", SireID := NA_character_]
raw[DamID == "0", DamID := NA_character_]

# Rename to standard names
setnames(raw, c("AnimalID", "SireID", "DamID"), c("Ind", "Sire", "Dam"))

cat("Founders (both parents NA):", sum(is.na(raw$Sire) & is.na(raw$Dam)), "\n")
cat("Sex table:\n"); print(table(raw$Sex))
cat("Year range:", range(raw$Year), "\n")
cat("Breed table:\n"); print(table(raw$Breed))
cat("\n")

# ── Step 1: tidyped ──────────────────────────────────────────────
cat("=== Step 1: tidyped() ===\n")
tped <- tryCatch({
  suppressMessages(tidyped(raw))
}, error = function(e) { cat("ERROR:", e$message, "\n"); NULL })

if (is.null(tped)) stop("tidyped failed. Cannot continue.")
cat("tidyped rows:", nrow(tped), "\n")
cat("Gen range:", range(tped$Gen), "\n")
cat("Columns:", paste(names(tped), collapse = ", "), "\n")

# Check extra columns preserved
stopifnot("Year" %in% names(tped))
stopifnot("Breed" %in% names(tped))
cat("PASS: tidyped created with extra columns\n\n")

# ── Test 1: pedgenint (7 pathways) ───────────────────────────────
cat("=== Test 1: pedgenint(timevar='Year') ===\n")
tryCatch({
  gi <- suppressMessages(pedgenint(tped, timevar = "Year"))
  print(gi)
  stopifnot(nrow(gi) == 7)
  stopifnot(all(c("SS","SD","DS","DD","SO","DO","Average") %in% gi$Pathway))
  stopifnot(all(gi$N > 0))
  stopifnot(all(gi$Mean > 0))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 2: pedgenint by Breed ───────────────────────────────────
cat("=== Test 2: pedgenint(timevar='Year', by='Breed') ===\n")
tryCatch({
  gi_breed <- suppressMessages(pedgenint(tped, timevar = "Year", by = "Breed"))
  cat("Groups:", length(unique(gi_breed$Group)), "\n")
  cat("Rows:", nrow(gi_breed), "\n")
  print(head(gi_breed, 14))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 3: pedne(by='Year') ────────────────────────────────────
cat("=== Test 3: pedne(by='Year') ===\n")
tryCatch({
  ne <- suppressMessages(pedne(tped, by = "Year"))
  print(ne)
  stopifnot(all(c("Cohort","N","MeanF","DeltaF","Ne") %in% names(ne)))
  # Ne can be NA when DeltaF=0 (no inbreeding increase in that cohort)
  ne_valid <- ne[!is.na(Ne)]
  stopifnot(all(ne_valid$Ne > 0 & is.finite(ne_valid$Ne)))
  n_na <- sum(is.na(ne$Ne))
  if (n_na > 0) {
    cat("Note:", n_na, "cohort(s) with DeltaF=0, Ne set to NA\n")
  }
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 4: pedecg ──────────────────────────────────────────────
cat("=== Test 4: pedecg() ===\n")
tryCatch({
  ecg <- pedecg(tped)
  cat("ECG range:", range(ecg$ECG), "\n")
  cat("FullGen range:", range(ecg$FullGen), "\n")
  cat("MaxGen range:", range(ecg$MaxGen), "\n")
  stopifnot(nrow(ecg) == nrow(tped))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 5: pedsubpop by Year ────────────────────────────────────
cat("=== Test 5: pedsubpop(by='Year') ===\n")
tryCatch({
  sp_year <- pedsubpop(tped, by = "Year")
  print(sp_year)
  stopifnot(nrow(sp_year) == length(unique(tped$Year)))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 6: pedsubpop by Breed ───────────────────────────────────
cat("=== Test 6: pedsubpop(by='Breed') ===\n")
tryCatch({
  sp_breed <- pedsubpop(tped, by = "Breed")
  print(sp_breed)
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 7: pedsubpop connected components ───────────────────────
cat("=== Test 7: pedsubpop() [connected components] ===\n")
tryCatch({
  sp_comp <- pedsubpop(tped)
  print(sp_comp)
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 8: pedinbreed_class ─────────────────────────────────────
cat("=== Test 8: pedinbreed_class() ===\n")
tryCatch({
  inb <- pedinbreed_class(tped)
  print(inb[])
  stopifnot(abs(sum(inb$Percentage) - 100) < 0.01)
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 9: pedcontrib ───────────────────────────────────────────
cat("=== Test 9: pedcontrib(mode='both') ===\n")
tryCatch({
  # Use 2025 cohort (largest year group)
  cand_ids <- tped[Year == 2025, Ind]
  cat("Candidate cohort (Year=2025):", length(cand_ids), "\n")
  # Sample if too many for speed
  if (length(cand_ids) > 500) {
    set.seed(42)
    cand_ids <- sample(cand_ids, 500)
  }
  cont <- suppressMessages(pedcontrib(tped, cohort = cand_ids, mode = "both", top = 10))
  cat("Ne_f:", cont$summary$Ne_f, "\n")
  cat("Ne_a:", cont$summary$Ne_a, "\n")
  cat("Top 5 founders:\n")
  print(head(cont$founders, 5))
  cat("Top 5 ancestors:\n")
  print(head(cont$ancestors, 5))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 10: pedancestry(labelvar='Breed') ───────────────────────
cat("=== Test 10: pedancestry(labelvar='Breed') ===\n")
tryCatch({
  anc <- suppressMessages(pedancestry(tped, labelvar = "Breed"))
  cat("Columns:", paste(names(anc), collapse = ", "), "\n")
  cat("Rows:", nrow(anc), "\n")
  # All rows should sum to 1.0 (100% ancestry accounted for)
  row_sums <- rowSums(anc[, -"Ind"], na.rm = TRUE)
  cat("Row sum range:", range(row_sums), "\n")
  stopifnot(all(abs(row_sums - 1.0) < 1e-6))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 11: pedrel(by='Year') ───────────────────────────────────
cat("=== Test 11: pedrel(by='Year', compact=TRUE) ===\n")
tryCatch({
  rel_year <- suppressMessages(suppressWarnings(pedrel(tped, by = "Year", compact = TRUE)))
  print(rel_year)
  stopifnot(all(c("Year","NTotal","NUsed","MeanRel") %in% names(rel_year)))
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 12: pedstats + vispstat ─────────────────────────────────
cat("=== Test 12: pedstats(timevar='Year') ===\n")
tryCatch({
  st <- suppressMessages(pedstats(tped, timevar = "Year"))
  cat("Summary:\n"); print(st$summary)
  cat("\nGeneration intervals:\n"); print(st$gen_intervals)
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

# ── Test 13: pedpartial ─────────────────────────────────────────
cat("=== Test 13: pedpartial (top 3 ancestors of a random inbred) ===\n")
tryCatch({
  # Find an inbred individual
  inb_ped <- inbreed(tped)
  inbred_inds <- inb_ped[f > 0.05, Ind]
  cat("Individuals with F > 0.05:", length(inbred_inds), "\n")
  if (length(inbred_inds) > 0) {
    target <- inbred_inds[1]
    cat("Target individual:", target, "F =", inb_ped[Ind == target, f], "\n")
    # Get actual ancestors of the target individual via pedcontrib
    target_cont <- suppressMessages(pedcontrib(tped, cohort = target, mode = "ancestor", top = 3))
    top_anc <- target_cont$ancestors$Ind[1:min(3, nrow(target_cont$ancestors))]
    cat("Actual top ancestors:", paste(top_anc, collapse = ", "), "\n")

    part <- pedpartial(tped, ancestors = top_anc)
    cat("Partial inbreeding decomposition columns:", paste(names(part), collapse = ", "), "\n")
    cat("Target row:\n")
    print(part[Ind == target])

    # The partial inbreeding values should be >= 0
    anc_cols <- setdiff(names(part), c("Ind", "f"))
    stopifnot(all(part[, lapply(.SD, function(x) all(x >= 0)), .SDcols = anc_cols]))
  } else {
    cat("No highly inbred individuals found. Skipping.\n")
  }
  cat("PASS\n\n")
}, error = function(e) cat("ERROR:", e$message, "\n\n"))

cat("========================================\n")
cat("G09 Integration Test Complete\n")
cat("========================================\n")

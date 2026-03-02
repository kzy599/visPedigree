# Tests for compact pedigree matrix functionality

# Helper to create a pedigree with large full-sibling families
create_fullsib_ped <- function(n_sires = 3, n_dams = 4, offspring_per_family = 10) {
  ped_list <- list()
  ind_counter <- 1

  sires <- paste0("S", sprintf("%02d", seq_len(n_sires)))
  ped_list[[length(ped_list) + 1]] <- data.table::data.table(
    Ind = sires, Sire = NA, Dam = NA, Sex = "male"
  )

  dams <- paste0("D", sprintf("%02d", seq_len(n_dams)))
  ped_list[[length(ped_list) + 1]] <- data.table::data.table(
    Ind = dams, Sire = NA, Dam = NA, Sex = "female"
  )

  for (s in seq_len(n_sires)) {
    for (d in seq_len(n_dams)) {
      ids <- paste0("O", sprintf("%04d", ind_counter:(ind_counter + offspring_per_family - 1)))
      ind_counter <- ind_counter + offspring_per_family
      sexes <- rep(c("male", "female"), length.out = offspring_per_family)
      ped_list[[length(ped_list) + 1]] <- data.table::data.table(
        Ind = ids, Sire = sires[s], Dam = dams[d], Sex = sexes
      )
    }
  }
  data.table::rbindlist(ped_list)
}

test_that("compact pedigree reduces dimensions", {
  ped <- create_fullsib_ped(n_sires = 3, n_dams = 4, offspring_per_family = 20)
  tped <- tidyped(ped, addnum = TRUE)
  n_orig <- nrow(tped)

  f_compact <- pedmat(tped, method = "f", compact = TRUE)

  ci <- attr(f_compact, "call_info")
  expect_true(ci$compact)
  expect_equal(ci$n_original, n_orig)
  expect_true(ci$n_compact < n_orig)
})

test_that("compact f matches non-compact f for all individuals", {
  ped <- create_fullsib_ped(n_sires = 2, n_dams = 3, offspring_per_family = 5)
  tped <- tidyped(ped, addnum = TRUE)

  f_full <- pedmat(tped, method = "f", compact = FALSE)
  f_compact <- pedmat(tped, method = "f", compact = TRUE)

  # Expand compact and compare
  f_expanded <- expand_pedmat(f_compact)
  expect_equal(length(f_expanded), length(f_full))
  expect_equal(as.numeric(f_expanded), as.numeric(f_full), tolerance = 1e-10)
})

test_that("compact A matrix expand matches non-compact A", {
  ped <- create_fullsib_ped(n_sires = 2, n_dams = 2, offspring_per_family = 5)
  tped <- tidyped(ped, addnum = TRUE)

  A_full <- pedmat(tped, method = "A", compact = FALSE, sparse = FALSE)
  A_compact <- pedmat(tped, method = "A", compact = TRUE, sparse = FALSE)

  A_expanded <- expand_pedmat(A_compact)

  expect_equal(dim(A_expanded), dim(A_full))
  expect_equal(as.matrix(A_expanded), as.matrix(A_full), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("compact_map and family_summary are well-formed", {
  ped <- create_fullsib_ped(n_sires = 2, n_dams = 3, offspring_per_family = 8)
  tped <- tidyped(ped, addnum = TRUE)

  f_compact <- pedmat(tped, method = "f", compact = TRUE)

  compact_map <- attr(f_compact, "compact_map")
  family_summary <- attr(f_compact, "family_summary")
  compact_stats <- attr(f_compact, "compact_stats")

  # compact_map should cover all original individuals
  expect_true(nrow(compact_map) >= nrow(tped))

  # family_summary should have one row per full-sib family
  expect_true(nrow(family_summary) > 0)
  expect_true("FamilyID" %in% names(family_summary))
  expect_true("RepInd" %in% names(family_summary))

  # compact_stats should report compression
  expect_equal(compact_stats$n_original, nrow(tped))
  expect_true(compact_stats$compression_ratio < 1)
  expect_true(compact_stats$memory_saved_pct > 0)
})

test_that("query_relationship works for full-siblings in compact mode", {
  tped <- tidyped(small_ped)

  A_compact <- pedmat(tped, method = "A", compact = TRUE)
  A_full <- pedmat(tped, method = "A", compact = FALSE)

  # C, D, E are full-sibs (AxB)
  val_compact <- query_relationship(A_compact, "C", "D")
  val_full <- query_relationship(A_full, "C", "D")
  expect_equal(val_compact, val_full, tolerance = 1e-10)
})

test_that("query_relationship returns NA for unknown individuals", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", compact = TRUE)
  expect_true(is.na(query_relationship(A, "NONEXISTENT", "A")))
})

test_that("expand_pedmat returns unchanged for non-compact", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", compact = FALSE)
  A2 <- expand_pedmat(A)
  expect_equal(dim(A), dim(A2))
})

test_that("compact mode works with D matrix", {
  ped <- create_fullsib_ped(n_sires = 2, n_dams = 2, offspring_per_family = 5)
  tped <- tidyped(ped, addnum = TRUE)

  D_compact <- pedmat(tped, method = "D", compact = TRUE)
  D_full <- pedmat(tped, method = "D", compact = FALSE)

  # Query a specific pair
  val_compact <- query_relationship(D_compact, "S01", "D01")
  val_full <- query_relationship(D_full, "S01", "D01")
  expect_equal(val_compact, val_full, tolerance = 1e-6)
})

test_that("summary_pedmat works for compact matrices", {
  ped <- create_fullsib_ped(n_sires = 2, n_dams = 3, offspring_per_family = 10)
  tped <- tidyped(ped, addnum = TRUE)

  A_compact <- pedmat(tped, method = "A", compact = TRUE)
  s <- summary_pedmat(A_compact)

  expect_s3_class(s, "summary.pedmat")
  expect_true(s$compact)
  expect_true(!is.null(s$compact_stats))
  expect_output(print(s), "Compaction Results")
})

test_that("pedigree with no full-sibs returns uncompacted", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "A"),
    Dam = c(NA, NA, "B", NA),
    Sex = c("male", "female", "male", "female"),
    stringsAsFactors = FALSE
  )
  tped <- tidyped(ped, addnum = TRUE)

  f_compact <- pedmat(tped, method = "f", compact = TRUE)
  ci <- attr(f_compact, "call_info")

  # No full-sibs to compact, so n_compact should equal n_original
  expect_true(ci$n_compact <= ci$n_original)
})

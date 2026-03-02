library(testthat)
library(visPedigree)
library(data.table)

# Create a minimal testing dataset
test_ped_dt <- data.table::data.table(
  Ind = c("A", "B", "C", "D", "E", "F"),
  Sire = c(NA, NA, "A", "A", "C", "C"),
  Dam = c(NA, NA, "B", "B", "D", "D"),
  Sex = c("male", "female", "male", "female", "male", "female"),
  Gen = c(1, 1, 2, 2, 3, 3),
  Group = c("G1", "G1", "G2", "G2", "G3", "G3")
)
# C and D are full siblings (F _ offspring = 0.25). 
# E and F are offspring of full siblings.

test_ped <- suppressMessages(tidyped(test_ped_dt, cand = c("E", "F")))

test_that("pedrel calculates relation correctly and handles boundaries", {
  # N < 2 case for G1 (after subsetting by Gen, if one was dropped, or just test N=1)
  dt_n1 <- data.table::data.table(Ind = c("A", "B", "C"), Sire = c(NA, NA, "A"), Dam = c(NA, NA, "B"), Sex = c("male", "female", "male"), Gen = c(1, 1, 2))
  ped_n1 <- suppressMessages(tidyped(dt_n1))
  
  res_n1 <- suppressWarnings(pedrel(ped_n1, by = "Gen"))
  expect_equal(res_n1[Gen == 2, NUsed], 1)
  expect_true(is.na(res_n1[Gen == 2, MeanRel]))
  
  # For test_ped, group by Gen
  rel <- pedrel(test_ped, by = "Gen")
  
  # Gen 1: A, B (Unrelated) -> 0
  expect_equal(rel[Gen == 1, MeanRel], 0)
  
  # Gen 2: C, D (full sibs). a_ij = 0.5. Mean of off-diagonals = 0.5
  expect_equal(rel[Gen == 2, MeanRel], 0.5)
  
  # Test sample size
  rel_sample <- suppressWarnings(pedrel(test_ped, by = "Gen", sample = 1))
  expect_true(all(is.na(rel_sample$MeanRel)))
})

test_that("pedgenint keeps unweighted 4-pathway logic for Average", {
  test_ped$BirthYear <- c(2000, 2001, 2005, 2006, 2010, 2012)
  # SS: C-A (2005-2000=5), E-C (2010-2005=5)
  # SD: D-A (2006-2000=6), F-C (2012-2005=7)
  # DS: C-B (2005-2001=4), E-D (2010-2006=4)
  # DD: D-B (2006-2001=5), F-D (2012-2006=6)
  
  suppressMessages(
    genint_res <- pedgenint(test_ped, timevar = "BirthYear", unit = "year")
  )
  
  # Pathways:
  # SS mean = 5
  # SD mean = 6.5
  # DS mean = 4
  # DD mean = 5.5
  # Average mean = (5 + 6.5 + 4 + 5.5) / 4 = 21 / 4 = 5.25
  
  avg_res <- genint_res[Pathway == "Average"]
  expect_equal(avg_res$Mean, 5.25)
  expect_true(!is.na(avg_res$SD))
})

test_that("pedcontrib Ne_f and Ne_a compute on full sets despite top cutoff", {
  cont_all <- suppressMessages(pedcontrib(test_ped, mode = "both", top = 100))
  cont_top1 <- suppressMessages(pedcontrib(test_ped, mode = "both", top = 1))
  
  # Ne_f and Ne_a should be identical because they're based on full arrays before truncation
  expect_equal(cont_all$summary$Ne_f, cont_top1$summary$Ne_f)
  expect_equal(cont_all$summary$Ne_a, cont_top1$summary$Ne_a)
  
  # Ensure reported reflects top
  expect_equal(cont_top1$summary$n_founders_reported, 1)
  expect_equal(cont_all$summary$n_founders_reported, length(cont_all$founders$Ind))
  expect_true(cont_top1$summary$n_founders_total > 1) 
})

test_that("pedpartial and pedancestry run without addnum=TRUE initially", {
  # Create tidyped without IndNum (addnum = FALSE)
  ped_nonum <- suppressMessages(tidyped(test_ped_dt, addnum = FALSE))
  
  expect_false("IndNum" %in% names(ped_nonum))
  
  # Should not crash and automatically handle missing nums
  part_res <- pedpartial(ped_nonum, ancestors = c("A", "B"), top = 2)
  expect_true(all(c("A", "B") %in% names(part_res)))
  expect_equal(nrow(part_res), nrow(ped_nonum))
  
  # pedancestry with labels
  ped_nonum$Label <- ifelse(ped_nonum$Ind %in% c("A", "B"), "Base", NA)
  anc_res <- pedancestry(ped_nonum, labelvar = "Label")
  expect_true("Base" %in% names(anc_res))
  # All non-founders descend completely from A and B, so Base should be 1.0 for all
  expect_true(all(anc_res$Base == 1.0))
})

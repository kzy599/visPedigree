library(visPedigree)
library(data.table)

# 1. Create a dummy pedigree with date strings
ped_data <- data.table(
  Ind = c("A", "B", "C", "D", "E"),
  Sire = c(NA, NA, "A", "A", "C"),
  Dam = c(NA, NA, "B", "B", "D"),
  Sex = c("male", "female", "male", "female", "male"),
  BirthDate = c("2020-01-01", "2020-03-01", "2020-12-01", "2021-02-01", "2022-06-01")
)

tped <- tidyped(ped_data)

# Test 1: Calculate generation interval in months
cat("\n--- Test 1: Monthly Intervals ---\n")
gi_month <- pedgenint(tped, timevar = "BirthDate", unit = "month")
print(gi_month)

# Test 2: Calculate with cycle_length (e.g., 6 months per generation)
cat("\n--- Test 2: Generation Equivalents (6 months/cycle) ---\n")
gi_gen <- pedgenint(tped, timevar = "BirthDate", unit = "month", cycle_length = 6)
print(gi_gen)

# Test 3: Standard pedstats integration
cat("\n--- Test 3: pedstats with Years ---\n")
stats <- pedstats(tped, timevar = "BirthDate", unit = "year")
print(stats)

# Test 4: Plot check (visual inspection in real R, but here we check the object)
cat("\n--- Test 4: Plot Metadata Check ---\n")
p <- vispstat(stats, type = "genint")
cat("Plot ylab:", p$y.limits[[1]], "\n") # Note: lattice structure varies, but we check if it runs

# Test 5: pedne integration
cat("\n--- Test 5: pedne ---\n")
ne_res <- pedne(tped, by = "BirthDate")
print(ne_res)

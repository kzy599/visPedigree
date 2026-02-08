
library(visPedigree)
library(data.table)

# Create a founders-only pedigree
founders <- data.table(
  Ind = c("A", "B"), 
  Sire = c(NA, NA), 
  Dam = c(NA, NA), 
  Sex = c("male", "female")
)

# Try to run tidyped on it
tryCatch({
  tp <- tidyped(founders)
  print("Founders only: Success")
}, error = function(e) {
  print(paste("Founders only: Failed -", e$message))
})

# Create a generation 1 pedigree (subset where parents are not in rows)
gen1 <- data.table(
  Ind = c("C", "D"),
  Sire = c("A", "A"), 
  Dam = c("B", "B"),
  Sex = c("male", "female")
)

tryCatch({
  tp <- tidyped(gen1)
  print("Gen1 only: Success")
  print(tp)
}, error = function(e) {
  print(paste("Gen1 only: Failed -", e$message))
})

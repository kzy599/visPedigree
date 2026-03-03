devtools::load_all(".")
library(data.table)
raw <- fread("inst/extdata/G09_prunedPed.csv")
setnames(raw, "SexID", "Sex")
raw[SireID == "0", SireID := NA_character_]
raw[DamID == "0", DamID := NA_character_]
setnames(raw, c("AnimalID","SireID","DamID"), c("Ind","Sire","Dam"))
tped <- suppressMessages(tidyped(raw))

target <- "MRYGY18G020A306"
cat("Target:", target, "\n")
cat("Sire:", tped[Ind == target, Sire], "\n")
cat("Dam:", tped[Ind == target, Dam], "\n")

sire_id <- tped[Ind == target, Sire]
dam_id <- tped[Ind == target, Dam]
cat("Sire parents:", tped[Ind == sire_id, Sire], tped[Ind == sire_id, Dam], "\n")
cat("Dam parents:", tped[Ind == dam_id, Sire], tped[Ind == dam_id, Dam], "\n")

# Check if this is a simple case where parents share a parent
gs_s <- tped[Ind == sire_id, Sire]
gs_d <- tped[Ind == sire_id, Dam]
gd_s <- tped[Ind == dam_id, Sire]
gd_d <- tped[Ind == dam_id, Dam]
cat("Paternal GS:", gs_s, "GD:", gs_d, "\n")
cat("Maternal GS:", gd_s, "GD:", gd_d, "\n")

# Find common ancestor
common <- intersect(c(gs_s, gs_d), c(gd_s, gd_d))
cat("Common grandparents:", paste(common, collapse=", "), "\n")

# pedcontrib to find top ancestors
cont <- suppressMessages(pedcontrib(tped, cohort = target, mode = "ancestor", top = 5))
cat("\nTop ancestors (pedcontrib):\n")
print(cont$ancestors)

# Find common ancestors between sire and dam paths
# An individual is inbred when its parents share ancestors
sire_anc <- suppressMessages(pedcontrib(tped, cohort = sire_id, mode = "founder", top = 50))
dam_anc <- suppressMessages(pedcontrib(tped, cohort = dam_id, mode = "founder", top = 50))
cat("\nSire top founders:\n")
print(head(sire_anc$founders, 10))
cat("\nDam top founders:\n")
print(head(dam_anc$founders, 10))

# Common founders
common_founders <- intersect(sire_anc$founders_full$Ind, dam_anc$founders_full$Ind)
cat("\nCommon founders between sire and dam:", length(common_founders), "\n")
if (length(common_founders) > 0) {
  cat(paste(head(common_founders, 10), collapse=", "), "\n")
  # pedpartial with common founders
  cat("\npedpartial with common founders:\n")
  top_cf <- head(common_founders, 5)
  part <- pedpartial(tped, ancestors = top_cf)
  cat("Target row:\n")
# Also trace back great-grandparents manually
cat("\nGreat-grandparents (sire side):\n")
for (gp in c(gs_s, gs_d)) {
  cat(" ", gp, "-> Sire:", tped[Ind == gp, Sire], "Dam:", tped[Ind == gp, Dam], "\n")
}
cat("Great-grandparents (dam side):\n")
for (gp in c(gd_s, gd_d)) {
  cat(" ", gp, "-> Sire:", tped[Ind == gp, Sire], "Dam:", tped[Ind == gp, Dam], "\n")
}

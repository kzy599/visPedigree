
library(visPedigree)
library(data.table)

# 1. 准备数据
ped_data <- data.table(
   Ind = c("A", "B", "C", "D", "E", "F", "G"),
   Sire = c(NA, NA, "A", "A", "C", "E", "E"),
   Dam = c(NA, NA, "B", "B", "D", "B", "F"),
   Sex = c("male", "female", "male", "female", "male", "female", "male")
 )
tped <- tidyped(ped_data)

# Case 1: showf=TRUE (Internal calculation). This usually worked? 
# Wait, if internal calculation happens, tidyped with inbreed is run inside prepare.
# Let's see: ped2igraph calls `finalize_graph`. `finalize_graph` modifies label.
# So even showf=TRUE would be broken without my fix, because layout used the modified label.
# Maybe the user was wrong about Case 1 being optimized? Or maybe Case 1 didn't have large enough F to trigger modification?
# Or maybe the labels matched coincidentally?
# Ah, if f=0, label is NOT modified!
# "Note: Inbreeding coefficients of 0 are not shown in the graph."
# So if the user's example has f=0 for all, labels are clean.
# But for inbred individuals, f > 0, labels modified => logic breakage.

# Let's make sure we have some Non-Zero inbreeding to test.
# Pedigree:
# A, B Founders
# C = A x B
# D = A x B (Full sib to C)
# E = C x D (Full sib mating! f(E) should be 0.25)
# F = E x B (Backcross)
# G = E x F (Inbred)

# E has f = 0.25.
# So E's label will be "E\n[0.25]".
# E is parent of F (Sire=E, Dam=B).
# F is Generation 3 (A/B=1, C/D=2, E=3, F/G=4?). No, generations logic:
# Gen 1: A, B
# Gen 2: C, D
# Gen 3: E
# Gen 4: F
# Gen 5: G

# When laying out F (Gen 4), we look for Sire E (Gen 3).
# Sire label is "E".
# Parent node E label is "E\n[0.25]".
# Match failed -> F not centered under E.

# My fix ensures Parent node has Ind="E". Match succeeds.

cat("Running visped with showf=TRUE...\n")
tryCatch({
    # We can't easily check visual layout in script, but we can ensure it runs without error.
    # And we can check internal structure if we could access invisible output.
    res <- visped(inbreed(tped), compact = FALSE, showf = TRUE, showgraph = FALSE, file = tempfile(fileext = ".pdf"))
    
    # We can verify that 'Ind' column exists in the nodes and matches regex of clean IDs
    nodes <- res$g |> igraph::as_data_frame(what="vertices")
    
    if (!"Ind" %in% names(nodes)) {
        stop("Fix failed: Ind column not preserved in graph nodes.")
    }
    
    cat("Success: Ind column preserved. Graph generation passed.\n")
}, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
})

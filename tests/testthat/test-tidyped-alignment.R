library(testthat)
library(data.table)

# Load package functions if not available
# if (!exists("tidyped")) {
#   devtools::load_all()
# }

test_that("Mate Alignment: Push Down Logic (pmax)", {
  # Case 1: Founder (F) mates with Descendant (D) at Gen 3.
  # A -> B -> D
  # F -> D (F is valid parent of D? No, F and D are mates? No wait)
  # Correction: Mate alignment is about optimizing shared children.
  # If A and B are parents of C.
  # A is Gen 1 (Founder).
  # B is Gen 3 (Child of ...).
  # C is Gen 4.
  # Old Logic (pmin): Pull B up to Gen 1. 
  #   -> If B has parents at Gen 2, B(1) < Parents(2) -> ERROR.
  # New Logic (pmax): Push A down to Gen 3.
  #   -> A(3), B(3), C(4). OK.
  
  ped <- data.table(
    Ind  = c("G1", "G2", "G3", "F_Mate", "Kind"),
    Sire = c(NA,   "G1", "G2", NA,       "G3"),     # G1 -> G2 -> G3. F_Mate x G3 -> Kind
    Dam  = c(NA,   NA,   NA,   NA,       "F_Mate")
  )
  
  # Expected generations (bottom-up or top-down, let's use bottom)
  # Bottom-up tends to push things deep.
  # Initial calculation:
  # G1 (Gen 1 or similar distance from bottom)
  # Kind is leaf.
  
  # Let's inspect generations directly.
  # G1 (Gen 1)
  # G2 (Gen 2)
  # G3 (Gen 3)
  # F_Mate (Gen 1 initially if treated as founder)
  # Kind (Gen 4)
  
  res <- tidyped(ped, genmethod = "bottom")
  setkey(res, Ind)
  
  # Check Baseline Generations
  expect_gt(res["G3", Gen], res["G2", Gen])
  expect_gt(res["G2", Gen], res["G1", Gen])
  
  # Check Mate Alignment
  # F_Mate should be aligned with G3
  expect_equal(res["F_Mate", Gen], res["G3", Gen], 
               info = "Founder Mate should be pushed down to Spouse's generation")
  
  # Check Topology
  # Kind must be > Parents
  expect_gt(res["Kind", Gen], res["F_Mate", Gen])
  expect_gt(res["Kind", Gen], res["G3", Gen])
})


test_that("Mate Alignment: Topological Constraint (Block Push)", {
  # Case 2: Founder (F) mates with G3 (Deep).
  # But F ALREADY has a child C at Gen 2.
  # If we push F to Gen 3 (to match G3), F(3) > C(2). VIOLATION.
  # Logic should block the push. F stays at 1.
  
  ped <- data.table(
    Ind  = c("G1", "G2", "G3", "F", "C", "Kind"),
    Sire = c(NA,   "G1", "G2", NA,  NA,  "G3"),
    Dam  = c(NA,   NA,   NA,   NA,  "F", "F") 
  )
  # Relationships:
  # G1 -> G2 -> G3
  # F -> C (F is parent of C)
  # G3 x F -> Kind (F is parent of Kind)
  
  # Initial Gen estimate:
  # G1: 1
  # G2: 2
  # G3: 3
  # F: 1 (Founder)
  # C: 2 (Child of F)
  # Kind: 4 (Child of G3(3) and F(1/3?))
  
  # Conflict:
  # Mate Alignment wants F and G3 to be same Gen.
  # G3 is Gen 3.
  # F ideal is Gen 3.
  # But F has child C at Gen 2.
  # F cannot be >= C. F must be < C.
  # F max valid Gen is 1 (if C is 2).
  # Result: F stays at 1. G3 stays at 3.
  
  res <- tidyped(ped, genmethod = "top")
  setkey(res, Ind)
  
  # F should NOT be pushed to 3
  expect_lt(res["F", Gen], res["C", Gen])
  expect_equal(res["F", Gen], 1, info = "F constrained to Gen 1 by child C")
  expect_equal(res["G3", Gen], 3, info = "G3 remains at Gen 3")
  
  # Mates are misaligned
  expect_true(res["F", Gen] != res["G3", Gen], info = "Mates remain misaligned due to constraints")
})

test_that("Reproduce Bug with Parents X, U, V", {
  # The small_ped bug.
  # X (Child) was Gen 1. U,V (Parents) were Gen 4.
  # Fixed logic should ensure X > U, V.
  
  # Load real simplified data if possible, or mock it closely.
  # The bug was triggered by X mating with a founder or similar.
  # Mock:
  # U(4) x V(4) -> X(5)
  # X(5) x N(1) -> Offspring
  # Old Bug: N(1) pulls X to 1.
  # New Fix: X(5) pushes N to 5.
  
  ped <- data.table(
    Ind = c("U", "V", "X", "N", "O"),
    Sire = c(NA, NA, "U", "X", "X"), # U->X, X->O
    Dam = c(NA, NA, "V", "N", "N"),  # V->X, N->O
    GenInit = c(4, 4, 5, 1, 6) # Simulated initial state or structural implication
  )
  # We provide structure, tidyped calculates Gen.
  # We need depth for U,V to be 4.
  # Chain: G1->G2->G3->U
  
  ped_chain <- data.table(
    Ind = c("A","B","C","U", "V", "X", "N", "O"),
    Sire = c(NA,"A","B","C", NA,  "U", NA,  "X"),
    Dam =  c(NA,NA,NA, NA,  NA,  "V", NA,  "N")
  )
  # A(1)->B(2)->C(3)->U(4).
  # V(1/4?) - lets assume V is deep or aligned.
  # X(5) from U(4).
  # N(1).
  # O(6) from X(5)xN(1).
  
  res <- tidyped(ped_chain, genmethod = "bottom")
  setkey(res, Ind)
  
  # Check X vs U
  expect_gt(res["X", Gen], res["U", Gen])
  
  # Check Mate Alignment (X and N)
  # X is 5 (approx). N should be pushed to 5.
  expect_equal(res["N", Gen], res["X", Gen], info = "Founder N aligned to Mate X")
})

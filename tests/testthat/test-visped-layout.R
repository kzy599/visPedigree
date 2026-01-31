
test_that("visped handles isolated nodes in layout correctly", {
  skip_if_not_installed("visPedigree")
  
  # This specific test case (small_ped with bottom sorting + specific highlight)
  # created a graph where some nodes existed in the node table but not in the edge list
  # causing a mismatch when assigning x coordinates back to the node table.
  # The fix was to ensure graph_from_data_frame includes 'vertices' argument.
  
  tped_bottom <- tidyped(small_ped, genmethod = "bottom")
  expect_error(visped(tped_bottom, highlight = "Z1"), NA)
  
  # Test with another method to be sure
  tped_top <- tidyped(small_ped, genmethod = "top")
  expect_error(visped(tped_top, highlight = "Z1"), NA)
  
  # Expanded test cases with other datasets
  
  # 1. simple_ped
  tped_simple_bottom <- tidyped(simple_ped, genmethod = "bottom")
  expect_error(visped(tped_simple_bottom), NA)
  
  # 2. deep_ped (larger dataset)
  tped_deep <- tidyped(deep_ped, genmethod = "bottom")
  # Taking a subset via cand to test layout logic without overwhelming plot
  tped_deep_sub <- tidyped(deep_ped, cand="K110550H", trace="up", tracegen=3, genmethod="bottom")
  expect_error(visped(tped_deep_sub), NA)
  
  # 3. loop_ped (with loops - although tidyped cleans/warns, visped should layout cleanly whatever tidyped produces)
  # loop_ped causes error in tidyped by default, skip full fail test here, assume clean input
  
  # 4. Edge case: very small layout - requires at least some parentage to be built otherwise tidyped stops
  # Create a minimal valid pedigree
  mini_ped <- data.frame(Ind=c("A", "B"), Sire=c(NA, "A"), Dam=c(NA, NA))
  tmini <- tidyped(mini_ped)
  expect_error(visped(tmini), NA)
})

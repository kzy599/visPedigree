# visPedigree Optimization Roadmap

This document outlines the planned architectural and functional optimizations for future versions of the `visPedigree` package.

## Phase 1: Core Architecture Upgrade (S3 Class System & Validation)
- [ ] **Define S3 Class**: Add `class(res) <- c("tidyped", class(res))` to the return value of `tidyped()`.
- [ ] **Implement Constructor**: Create an internal `new_tidyped()` function to ensure consistent object creation.
- [ ] **Implement Validator**: Write `validate_tidyped()` to centrally check for required pedigree columns (Ind, Sire, Dam) and data types.
- [ ] **Define Print Method**: Implement `print.tidyped()` to provide users with a concise summary of the pedigree (e.g., total individuals, generations, sex ratio).

## Phase 2: Modular Refactoring (Decoupling Logic & Presentation)
- [ ] **Extract Calculation Logic**: Encapsulate coordinate calculation and layout generation from `visped()` into a standalone function `prepare_ped_graph()`.
- [ ] **Unify Plotting Engine**: Refactor `visped()` to act as a wrapper for `prepare_ped_graph()`, handling output devices (PDF/Screen) and argument dispatching.
- [ ] **Introduce Argument Passing (`...`)**: Add `...` to all exported functions and use `utils::modifyList()` and `do.call()` to safely pass arguments to underlying functions (e.g., `igraph::plot` or `pdf()`).

## Phase 3: User Experience Optimization (Pipes & Methods)
- [ ] **Pipe-Friendly Design**: Ensure the first argument of all processing functions (`sortped`, `numped`, `tidyped`) is `ped`, and they return the modified object.
- [ ] **Implement Plot Method**: Define `plot.tidyped()` so users can visualize directly via `plot(ped_tidy)`.
- [ ] **Enhance Error Messaging**: Use `validate_tidyped()` in all entry points to provide human-readable error messages (e.g., identifying missing columns or invalid characters).

## Phase 4: Performance & Code Quality
- [ ] **Optimize Memory Usage**: Review internal code to use `data.table` reference modifications (`set()` family) where possible without breaking original data, improving speed for large pedigrees.
- [ ] **Comprehensive Unit Testing**: Write test cases using `testthat` specifically for the new S3 methods and validators.
- [ ] **Update Documentation**: Clearly document the use of `...` in `roxygen2` comments and update examples to demonstrate pipe operations.

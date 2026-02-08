# Pedigree Analysis Module Design (v2.0)

## 1. Overview
The goal is to implement a comprehensive pedigree analysis pipeline comparable to RelaX2 and optiSel, but integrated into the `visPedigree` ecosystem using `tidyped` objects.

## 2. Core API

### 2.1 Main Entry Point: `pedstats`
*   **Function**: `pedstats(ped, tasks = NULL, ...)`
*   **Description**: A wrapper function that runs selected or all analysis modules.
*   **Returns**: An S3 object of class `pedstats`.

### 2.2 Functional Modules

#### A. Depth & Completeness
*   **`pedecg(ped, by = NULL)`**
    *   Calculates Equivalent Complete Generations (ECG).
    *   Returns: `data.table(Ind, ECG, FullGen, MaxGen)`.
*   **`pedcompl(ped, maxgen = 5, by = NULL)`**
    *   Calculates Pedigree Completeness per generation.
    *   Returns: `data.table(Group, Gen, Completeness)`.

#### B. Population & Diversity
*   **`pedne(ped, timevar = "BirthYear", method = "Gutierrez")`**
    *   Calculates effective population size ($N_e$) based on individual rate of inbreeding ($\Delta F$).
    *   Returns: `data.table(Period, Ne, DeltaF, MeanF, N)`.
*   **`pedgenint(ped, timevar = "BirthYear")`**
    *   Calculates generation intervals for 4 pathways (SS, SD, DS, DD).
    *   Returns: `data.table(Pathway, Interval, N, SE)`.
*   **`pedsubpop(ped)`**
    *   Analyzes subpopulation structure (wraps `splitped`).
    *   Returns: `data.table(Subpop, Size, MeanF, Ne, ...)`.

#### C. Contributions
*   **`pedcontrib(ped, cohort = NULL, mode = c("founder", "ancestor"))`**
    *   Calculates genetic contributions from founders ($f_e$) or influential ancestors ($f_a$).
    *   **Algorithm**: 
        *   `founder`: Gene dropping / probability flow.
        *   `ancestor`: Boichard's iterative algorithm (requires C++).
    *   Returns: 
        *   `summary`: $f_e, f_a, f_e/f_a$
        *   `contributors`: `data.table(Ind, Contrib, MarginalContrib, Rank)`

## 3. Data Structure: `pedstats` Object

The `pedstats` object is a list containing the following components (keys in snake_case):

```r
list(
  # 1. Summary Table (Descriptive Statistics)
  # Columns: Group, N, MeanF, MeanECG, Ne, Fe, Fa
  summary = data.table(...), 

  # 2. Trends (for plotting)
  # Columns: Year, MeanF, DeltaF, Ne
  trend = data.table(...),

  # 3. Contributors (Founders/Ancestors)
  # Columns: Ind, Contrib, CumContrib, Rank
  founders = data.table(...),
  ancestors = data.table(...),

  # 4. Intervals
  # Columns: Pathway, Interval, N
  interval = data.table(...),

  # 5. Metadata
  params = list(time_var = ..., cohort = ...)
)
```

## 4. Implementation Strategy

1.  **C++ Extensions (`src/`)**: 
    *   Implement `cpp_calculate_ecg` (Top-down DP).
    *   Implement `cpp_boichard_algo` for $f_a$.
2.  **R Interface**:
    *   Ensure all functions utilize `data.table` in-place modification (`:=`) where possible for performance.
    *   Functions should detect grouping columns in `tidyped` (e.g., `Breed`) specified by `by` argument.

## 5. Visualization (`visstats`)
*   `visstats(x, type = "trend")`: Plot F and Ne trends.
*   `visstats(x, type = "contrib")`: Plot top ancestors.
*   `visstats(x, type = "interval")`: Plot generation intervals.

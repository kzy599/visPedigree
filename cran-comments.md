## Test environments
* local macOS Sequoia 15.4, R 4.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The only remaining note is the standard "New submission" note. All other non-standard files and directories (.github, todo, etc.) have been added to .Rbuildignore.

## Resubmission

Responded to CRAN review comments:

- Added missing \value tags for plot.tidyped, print.tidyped, and print.summary.tidyped.
- Removed commented-out example line in tidyped by replacing it with an executable try() call.

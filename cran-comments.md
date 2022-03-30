## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

There is no ERROR nor WARNING from R CMD check. 

The only NOTE is because of multiple local function definitions are given for 'solveForbeta0' with different formal arguments
to accommodate different kinds of covariate structure, and this is not significant issue in our view.

This is a first submission.

This document is obtained using `usethis:: use_cran_comments(open = rlang::is_interactive())`.

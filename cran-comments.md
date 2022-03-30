## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 1 warnings | 2 note

There is no ERROR from R CMD check. The only WARNING on Windows is because:

```
"'qpdf' is needed for checks on size reduction of PDFs".
```

The first NOTE is because of multiple local function definitions are given for 'solveForbeta0' with different formal arguments
to accommodate different kinds of covariate structure.

The second NOTE is because this is a first submission.

This document is obtained using `usethis:: use_cran_comments(open = rlang::is_interactive())`.

## Submission

This is an update of RNentropy from version 1.2.3 (the current CRAN release) to
1.3.3.

The update adds an analysis of changes in relative isoform usage (isoform
switching) across samples: the new functions `RN_iso_calc()`,
`RNentropy_iso_switch()` and `RN_iso_select()`, and the example dataset
`RN_IsoSwitch_Example_S7`. The existing functions are unchanged. The citation
file and several help pages were also corrected.

## Test environments

* Local: Ubuntu 22.04.5 LTS, R 4.6.1
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release and oldrel-1)
* win-builder: <!-- fill in after running devtools::check_win_devel() -->

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.

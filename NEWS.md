# RNentropy 1.3.2

- `RN_iso_select` now combines the sample p-values of each TESTED gene into a
  gene-level p-value with the Simes method and compares it with the threshold,
  without correction across genes. The threshold argument `lpv_t` is renamed
  `gene_pv`, the `method` argument and the `lpv_adj` component are removed, and
  the gene-level p-values are returned as `gene_pv`. `selected` reports the
  gene-level p-value (`ISO_GENE_PV`) next to the uncorrected per-sample
  p-values.
- `RN_iso_calc` and `RNentropy_iso_switch` return the per-sample p-values as
  `pv` (formerly `lpv`), with columns named `ISO_PV_<sample>` (formerly
  `ISO_LPV_<sample>`).
- Renamed the sample columns of `RN_IsoSwitch_Example_S7` and of the compact S7
  example file to plain tissue names (`brain`, `heart`, `kidney`, `liver`,
  `lung`, `muscle`).

# RNentropy 1.3.1

- Added `RN_iso_select` to apply multiple-testing correction to isoform-switch
  local p-values and retain genes significant in at least one sample.
- Added runnable file-based examples for the standard and isoform-switch
  workflows using a shared synthetic regression dataset.

# RNentropy 1.3.0

- Initialized the public GitHub development repository from the verified CRAN
  1.2.3 release.
- Added continuous integration and repository provenance metadata.
- Added `RN_iso_calc` and `RNentropy_iso_switch` for the analysis of isoform
  switching, with explicit testability information and results validated against
  the independent C++ implementation.
- Added the complete historical S7 isoform-switch example as the compressed
  `RN_IsoSwitch_Example_S7` package dataset.

# RNentropy 1.2.3

- Fixed small bug in .RN_design_check preventing compatibility with R >=4.2.0
- Added reference to new publication
- Removed unneeded initialization of .Random.seed
- Fixed a few typos in the documentation
- Added README.md, NEWS.md, and cran-comments.md files

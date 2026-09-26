# Repository provenance

This repository was initialized from the original RNentropy 1.2.3 package
development tree. Its three historical commits from April 2022 are preserved,
and commit `66d7151` is marked by the `v1.2.3` tag.

The package payload at that tag was compared with the official
[CRAN RNentropy 1.2.3 package](https://cran.r-project.org/package=RNentropy).
All substantive package files matched; the differences in `DESCRIPTION` were
CRAN-generated repository metadata and line wrapping, and `data/datalist` was a
build-generated file.

The baseline also passed `R CMD check --no-manual --no-vignettes` with status
OK under R 4.6.1 on Linux on 26 August 2026.

The `main` branch begins post-release development at version 1.2.3.9000. The
independent historical C++ implementation is preserved separately at
[Federico77z/RNentropy_cpp](https://github.com/Federico77z/RNentropy_cpp).

The isoform-switch implementation added on the `isoswitch` branch was recovered
from the latest unpublished R development trees and adapted to the current
package. Its eligibility rules, numerical results, and untestable-sample states
were validated against `RNentropy_iso_switch` in the archived C++ repository.
The package includes a compact regression fixture derived from the historical
S7 example and the complete S7 input as the compressed
`RN_IsoSwitch_Example_S7` package dataset. Both drop the `_s7_1` suffix from
the historical sample column names (e.g. `brain_s7_1` becomes `brain`); their
values are unchanged. The C++ output used as the expected result keeps the
historical names. The complete S12 and S13 examples
remain outside the package repository because of their size.

The complete comparison was repeated on 27 August 2026 with the default
`min_expr = 1` and `pseudocount = 0.01`:

| Example | Genes | Tested | NOTEST(ISO) | NOTEST(EXP) | Numeric scores | Maximum absolute difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| S7 | 28,426 | 4,058 | 15,306 | 9,062 | 24,133 | 4.982e-5 |
| S12 | 28,426 | 5,325 | 15,306 | 7,795 | 31,733 | 4.546e-5 |
| S13 | 28,426 | 4,766 | 15,306 | 8,354 | 28,359 | 4.605e-5 |

All gene classifications and all `-`/`*` marker positions matched. Numerical
scores were compared with the C++ `main.res` files; the differences shown above
are consistent with their six-digit output precision. The historical
`summary.res` files were not used because rows containing `-` or `*` omit those
fields and shift the remaining values.

The synthetic file used for the package's file-based examples was copied
unchanged from `tests/data/synthetic_expression.tsv` in the archived C++
repository. It is a regression and documentation fixture rather than a
biological dataset. Its SHA-256 checksum is
`daf83d525ac6b43bdaf011d940a524d0c8e2e188d9cff942fa5a4f9732f04594`.

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

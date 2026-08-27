
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNentropy

<!-- badges: start -->

[![R-CMD-check](https://github.com/Federico77z/RNentropy/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Federico77z/RNentropy/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is the implementation of a method based on information theory
devised for the identification of genes showing a significant variation
of expression across multiple conditions. Given expression estimates
from any number of RNA-Seq samples and conditions it identifies genes or
transcripts with a significant variation of expression across all the
conditions studied, together with the samples in which they are over- or
under-expressed. [Zambelli F. et
al. (2018)](https://doi.org/10.1093/nar/gky055).

A detailed walk-through on how to use RNentropy is available at
[Zambelli F., Pavesi G.
(2021)](https://doi.org/10.1007/978-1-0716-1307-8_6).

The `main` branch is the development line based on the CRAN 1.2.3
release. The exact release sources are preserved by the
[`v1.2.3`](https://github.com/Federico77z/RNentropy/tree/v1.2.3) tag.

## Installation

Install the stable release from CRAN:

``` r
install.packages("RNentropy")
```

Install the development version from GitHub:

``` r
pak::pak("Federico77z/RNentropy")
```

## Example

This is a basic example showing how to use RNentropy. Please see
[Zambelli F., Pavesi G.
(2021)](https://doi.org/10.1007/978-1-0716-1307-8_6) for more info.

``` r
library(RNentropy)
# basic example code
##load expression values and experiment design
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
#Run RNentropy 
Results <- RN_calc(RN_Brain_Example_tpm, RN_Brain_Example_design)
#select only genes with significant changes of expression
Results <- RN_select(Results)
#Compute the Point Mutual information Matrix
Results <- RN_pmi(Results)
```

## Isoform-switch analysis

The package includes the complete historical S7 data used to validate
the isoform-switch implementation. It contains expression values for
81,314 transcripts from 28,426 genes across six tissues. Load the data
frame and run the analysis with `RN_iso_calc()`:

``` r
data("RN_IsoSwitch_Example_S7")
S7_Results <- RN_iso_calc(
  RN_IsoSwitch_Example_S7,
  gene.col = "GENE_ID"
)
S7_Results$gene_status
S7_Results$sample_status
S7_Results$lpv
```

For transcript expression stored in a text file,
`RNentropy_iso_switch()` provides the corresponding file-based
interface; see `?RNentropy_iso_switch` for its input options.

## Independent C++ implementation

The independently developed C++ reference implementation, including the
historical isoform-switch executable, is archived at
[Federico77z/RNentropy\_cpp](https://github.com/Federico77z/RNentropy_cpp).

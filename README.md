
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNentropy

<!-- badges: start -->
<!-- badges: end -->

This is sn implementation of a method based on information theory
devised for the identification of genes showing a significant variation
of expression across multiple conditions. Given expression estimates
from any number of RNA-Seq samples and conditions it identifies genes or
transcripts with a significant variation of expression across all the
conditions studied, together with the samples in which they are over- or
under-expressed. Zambelli et al. (2018) <doi:10.1093/nar/gky055>.

A detailed walk-through on how to use RNentropy is available at
<doi:10.1007/978-1-0716-1307-8_6>

## Installation

You can install the development version of RNentropy like so:

``` r
install.packages("RNentropy")
```

## Example

This is a basic example showing how to use RNentropy. Please see
<doi:10.1007/978-1-0716-1307-8_6> for more info.

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

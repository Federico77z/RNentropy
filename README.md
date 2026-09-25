
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

The `main` branch contains the 1.3.1 release. The exact CRAN 1.2.3
release sources are preserved by the
[`v1.2.3`](https://github.com/Federico77z/RNentropy/tree/v1.2.3) tag.

## Installation

Install the stable release from CRAN:

``` r
install.packages("RNentropy")
```

Install version 1.3.1 from GitHub:

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

### Starting from an expression file

`RNentropy()` reads normalized expression values directly from a text
file. The package includes a small synthetic example with three
conditions and two replicates per condition. The first three columns
contain transcript IDs, gene IDs, and comments; the remaining columns
contain expression values.

``` r
Example_file <- system.file(
  "extdata", "RNentropy_example.tsv",
  package = "RNentropy"
)
Sample_names <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
Example_design <- matrix(
  c(1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1),
  nrow = 6, byrow = TRUE,
  dimnames = list(Sample_names, c("A", "B", "C"))
)

File_Results <- RNentropy(
  Example_file,
  tr.col = 1,
  design = Example_design,
  header = FALSE,
  skip.col = c(2, 3),
  col.names = Sample_names
)
File_Results <- RN_select(File_Results)
File_Results$selected
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
S7_Results <- RN_iso_select(S7_Results)
S7_Results$gene_status
S7_Results$sample_status
S7_Results$lpv
S7_Results$gpv_adj
S7_Results$selected
```

`RN_iso_select()` combines the finite sample p-values of each gene with
status `TESTED` into one gene-level p-value with the Simes method,
applies multiple-testing correction across genes, and retains genes
whose adjusted gene-level p-value is at or below the threshold. Its
`gpv_t` argument is expressed as a raw p-value and defaults to 0.01; the
correction method defaults to Benjamini-Hochberg (`BH`). The selected
table reports the gene-level p-value (`ISO_GPV`), its adjusted value
(`CORR_ISO_GPV`), and the uncorrected per-sample isoform-switch scores
(`ISO_LPV_<sample>`), which show in which samples the change occurs.
This selection step is specific to the R package and has no counterpart
in the historical C++ isoform-switch executable.

### Starting from a transcript expression file

The same synthetic file can be analyzed for changes in relative isoform
usage. Here the gene-ID column is retained, while the comments column is
skipped.

``` r
Example_file <- system.file(
  "extdata", "RNentropy_example.tsv",
  package = "RNentropy"
)
Sample_names <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
Example_design <- matrix(
  c(1, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1),
  nrow = 6, byrow = TRUE,
  dimnames = list(Sample_names, c("A", "B", "C"))
)

Iso_File_Results <- RNentropy_iso_switch(
  Example_file,
  tr.col = 1,
  gene.col = 2,
  design = Example_design,
  header = FALSE,
  skip.col = 3,
  col.names = c("GENE_ID", Sample_names)
)
Iso_File_Results <- RN_iso_select(Iso_File_Results)
Iso_File_Results$gene_status
Iso_File_Results$selected
```

This compact dataset is synthetic and intended to demonstrate input
handling and expected analysis behavior. The complete S7 dataset above
provides the historical biological example. See `?RNentropy_iso_switch`
for all file-input options.

## Independent C++ implementation

The independently developed C++ reference implementation, including the
historical isoform-switch executable, is archived at
[Federico77z/RNentropy\_cpp](https://github.com/Federico77z/RNentropy_cpp).

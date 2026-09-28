\name{RNentropy-package}
\alias{RNentropy-package}
\docType{package}
\title{
\packageTitle{RNentropy}
}
\description{
\packageDescription{RNentropy}
}
\author{
\packageAuthor{RNentropy}

Maintainer: \packageMaintainer{RNentropy}
}
\details{
The package provides two analyses.

Variation of expression: \code{\link{RNentropy}} (from a file) or
\code{\link{RN_calc}} (from a data.frame) compute a global p-value for each gene
or transcript, testing whether its expression varies across all samples, and a
local p-value for each sample, testing whether it is over- or under-expressed.
\code{\link{RN_select}} selects genes with a significant global p-value and
summarises the local p-values of each condition, and \code{\link{RN_pmi}}
computes the point mutual information between conditions.

Isoform switching: \code{\link{RNentropy_iso_switch}} (from a file) or
\code{\link{RN_iso_calc}} (from a data.frame) compute, for each gene and each
sample, a p-value for a change in relative isoform usage with respect to the
samples of the other conditions. \code{\link{RN_iso_select}} combines them
into a gene-level p-value and selects the genes with a significant change.
}

\references{
Zambelli F, Mastropasqua F, Picardi E, D'Erchia AM, Pesole G, Pavesi G (2018).
RNentropy: an entropy-based tool for the detection of significant variation of
gene expression across multiple RNA-Seq experiments. \emph{Nucleic Acids
Research}, 46(8), e46. \doi{10.1093/nar/gky055}

Zambelli F, Pavesi G (2021). Using RNentropy to Detect Significant Variation in
Gene Expression Across Multiple RNA-Seq or Single-Cell RNA-Seq Samples. In
\emph{RNA Bioinformatics}, Methods in Molecular Biology, 2284, 77--96.
\doi{10.1007/978-1-0716-1307-8_6}
}
\keyword{ package }

\examples{
#load expression values and experiment design
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
#compute statistics and p-values (considering only a subset of genes due to
#examples running time limit of CRAN). 
Results <- RN_calc(RN_Brain_Example_tpm[1:10000,], RN_Brain_Example_design)
#select only genes with significant changes of expression
Results <- RN_select(Results)
#Compute the Point Mutual information Matrix
Results <- RN_pmi(Results)

#load expression values and experiment design
data("RN_BarresLab_FPKM", "RN_BarresLab_design")
#compute statistics and p-values (considering only a subset of genes due to
#examples running time limit of CRAN)
Results_B <- RN_calc(RN_BarresLab_FPKM[1:10000,], RN_BarresLab_design)
#select only genes with significant changes of expression
Results_B <- RN_select(Results_B)
#Compute the Point Mutual information matrix
Results_B <- RN_pmi(Results_B)

#isoform-switch analysis on a subset of the S7 example
Iso_file <- system.file("extdata", "isoswitch_s7.tsv", package = "RNentropy")
Iso_results <- RNentropy_iso_switch(Iso_file, tr.col = "TR_ID",
  gene.col = "GENE_ID")
Iso_results <- RN_iso_select(Iso_results)
Iso_results$selected
}

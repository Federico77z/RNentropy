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
\references{doi = {10.1093/nar/gky055}
doi = {10.1007/978-1-0716-1307-8_6}
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
}

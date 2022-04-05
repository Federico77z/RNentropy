\name{RN_BarresLab_FPKM}
\alias{RN_BarresLab_FPKM}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
RN_BarresLab_FPKM
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Data expression values (mean FPKM) for different cerebral cortex cell from Zhang et al: "An RNA-Seq transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex", J Neurosci 2014 Sep. 
RN_BarresLab_design is the corresponding design matrix.
}
\usage{data("RN_BarresLab_FPKM")}
\format{
  A data frame with 12978 observations on the following 7 variables.
  \describe{
    \item{\code{Astrocyte}}{a numeric vector}
    \item{\code{Neuron}}{a numeric vector}
    \item{\code{OPC}}{a numeric vector}
    \item{\code{NFO }}{a numeric vector}
    \item{\code{MO}}{a numeric vector}
    \item{\code{Microglia}}{a numeric vector}
    \item{\code{Endothelial}}{a numeric vector}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Zhang Y, Chen K, Sloan SA, Bennett ML, Scholze AR, O'Keeffe S, Phatnani HP,
Guarnieri P, Caneda C, Ruderisch N, Deng S, Liddelow SA, Zhang C, Daneman R,
Maniatis T, Barres BA, Wu JQ. An RNA-sequencing transcriptome and splicing
database of glia, neurons, and vascular cells of the cerebral cortex. J Neurosci. 2014 Sep 3;34(36):11929-47. doi: 10.1523/JNEUROSCI.1860-14.2014. Erratum in: J
Neurosci. 2015 Jan 14;35(2):846-6. PubMed PMID: 25186741; PubMed Central PMCID:
PMC4152602.

Data were downloaded from https://web.stanford.edu/group/barres_lab/brain_rnaseq.html
}
\examples{
data(RN_BarresLab_FPKM)
}
\keyword{datasets}

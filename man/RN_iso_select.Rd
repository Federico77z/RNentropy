\name{RN_iso_select}
\alias{RN_iso_select}
\title{
Select genes with significant changes in relative isoform usage.
}
\description{
Select genes with an adjusted isoform-switch p-value lower than or equal to a
user-defined threshold in at least one sample.
}
\usage{
RN_iso_select(Results, lpv_t = 0.01, method = "BH")
}
\arguments{
  \item{Results}{
The output of \code{RN_iso_calc} or \code{RNentropy_iso_switch}.
}
  \item{lpv_t}{
Threshold for adjusted local p-values, expressed as a raw p-value between 0 and
1 (Default: 0.01).
}
  \item{method}{
Multiple-testing correction method. Available methods are those accepted by
\code{p.adjust}; type \code{p.adjust.methods} to see the list. The default is
\code{"BH"} (Benjamini-Hochberg).
}
}

\details{
Correction is applied as one family across all finite sample p-values belonging
to genes with \code{gene_status} equal to \code{"TESTED"}. Genes that were not
tested and sample p-values that could not be calculated are excluded from the
correction family.

This selection step is specific to the R package. The historical C++
isoform-switch implementation produces unfiltered results and has no
corresponding selector.
}

\value{
The original input list containing two additional components:
\item{lpv_adj}{A gene by sample matrix containing adjusted p-values as
\code{-log10(adjusted p-value)}. Ineligible values are reported as \code{NA}.}
\item{selected}{A data.frame containing TESTED genes with an adjusted p-value
less than or equal to \code{lpv_t} in at least one sample. Gene identifiers are
stored as row names. The data.frame contains gene status, original
isoform-switch LPVs, and corrected LPVs named
\code{CORR_ISO_LPV_<sample>}. Rows are ordered by the strongest corrected
signal.}
}

\author{
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}

\examples{
Isoforms <- data.frame(
  gene = c("gene_1", "gene_1", "gene_2", "gene_2"),
  sample_1 = c(10, 1, 5, 5),
  sample_2 = c(1, 10, 5, 5)
)
Results <- RN_iso_calc(Isoforms, gene.col = "gene")
Results <- RN_iso_select(Results)
Results$selected
}
\keyword{ RNentropy }
\keyword{ Select }

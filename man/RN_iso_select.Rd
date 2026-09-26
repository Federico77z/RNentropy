\name{RN_iso_select}
\alias{RN_iso_select}
\title{
Select genes with significant changes in relative isoform usage
}
\description{
Combine the isoform-switch p-values of each gene into one gene-level p-value
and select genes with a gene-level p-value lower than or equal to a
user-defined threshold.
}
\usage{
RN_iso_select(Results, gpv_t = 0.01)
}
\arguments{
  \item{Results}{
The output of \code{RN_iso_calc} or \code{RNentropy_iso_switch}.
}
  \item{gpv_t}{
Threshold for gene-level p-values, expressed as a raw p-value between 0 and 1
(Default: 0.01).
}
}

\details{
For each gene with \code{gene_status} equal to \code{"TESTED"}, the finite
sample p-values \eqn{p_{(1)} \le \dots \le p_{(k)}}{p(1) <= ... <= p(k)} are
combined with the Simes method into the gene-level p-value
\eqn{\min_i k p_{(i)} / i}{min(k * p(i) / i)}. It tests the null hypothesis
that relative isoform usage does not change in any sample, accounting for the
multiple sample tests of the gene, and it is valid under independence and under
the positive dependence expected between samples of the same gene.

Gene-level p-values are compared directly with \code{gpv_t}; no correction is
applied across genes. Genes that were not tested, and TESTED genes without any
finite sample p-value, receive no gene-level p-value and are never selected.
Sample p-values are not corrected: the original isoform-switch LPVs are
reported next to the gene-level value to show in which samples the change
occurs.

This selection step is specific to the R package. The historical C++
isoform-switch implementation produces unfiltered results and has no
corresponding selector.
}

\value{
The original input list containing two additional components:
\item{gpv}{A named vector containing, for each gene, the Simes gene-level
p-value as \code{-log10(p-value)}. Genes without a gene-level p-value are
reported as \code{NA}.}
\item{selected}{A data.frame containing TESTED genes with a gene-level p-value
less than or equal to \code{gpv_t}. Gene identifiers are stored as row names.
The data.frame contains gene status, the gene-level p-value (\code{ISO_GPV}),
and the original isoform-switch LPVs (\code{ISO_LPV_<sample>}), all as
\code{-log10(p-value)}. Rows are ordered by decreasing \code{ISO_GPV}.}
}

\references{
Simes R. J. (1986) An improved Bonferroni procedure for multiple tests of
significance. \emph{Biometrika}, 73(3), 751--754.
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
Results$gpv
Results$selected
}
\keyword{ RNentropy }
\keyword{ Select }

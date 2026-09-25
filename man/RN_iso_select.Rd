\name{RN_iso_select}
\alias{RN_iso_select}
\title{
Select genes with significant changes in relative isoform usage
}
\description{
Combine the isoform-switch p-values of each gene into one gene-level p-value,
correct these across genes, and select genes with an adjusted gene-level
p-value lower than or equal to a user-defined threshold.
}
\usage{
RN_iso_select(Results, gpv_t = 0.01, method = "BH")
}
\arguments{
  \item{Results}{
The output of \code{RN_iso_calc} or \code{RNentropy_iso_switch}.
}
  \item{gpv_t}{
Threshold for adjusted gene-level p-values, expressed as a raw p-value between
0 and 1 (Default: 0.01).
}
  \item{method}{
Multiple-testing correction method applied to the gene-level p-values.
Available methods are those accepted by \code{p.adjust}; type
\code{p.adjust.methods} to see the list. The default is \code{"BH"}
(Benjamini-Hochberg).
}
}

\details{
For each gene with \code{gene_status} equal to \code{"TESTED"}, the finite
sample p-values \eqn{p_{(1)} \le \dots \le p_{(k)}}{p(1) <= ... <= p(k)} are
combined with the Simes method into the gene-level p-value
\eqn{\min_i k p_{(i)} / i}{min(k * p(i) / i)}. It tests the null hypothesis
that relative isoform usage does not change in any sample, and it is valid
under independence and under the positive dependence expected between samples
of the same gene.

Correction is then applied as one family across all gene-level p-values.
Genes that were not tested, and TESTED genes
without any finite sample p-value, are excluded from the correction family.
Sample p-values are not corrected: the original isoform-switch LPVs are
reported next to the adjusted gene-level value to show in which samples the
change occurs.

}

\value{
The original input list containing three additional components:
\item{gpv}{A named vector containing, for each gene, the Simes gene-level
p-value as \code{-log10(p-value)}. Genes excluded from the correction family
are reported as \code{NA}.}
\item{gpv_adj}{A named vector containing the adjusted gene-level p-values as
\code{-log10(adjusted p-value)}, with \code{NA} for excluded genes.}
\item{selected}{A data.frame containing TESTED genes with an adjusted
gene-level p-value less than or equal to \code{gpv_t}. Gene identifiers are
stored as row names. The data.frame contains gene status, the gene-level
p-value (\code{ISO_GPV}), the adjusted gene-level p-value
(\code{CORR_ISO_GPV}), and the original isoform-switch LPVs
(\code{ISO_LPV_<sample>}), all as \code{-log10(p-value)}. Rows are ordered by
decreasing \code{CORR_ISO_GPV}, then by decreasing \code{ISO_GPV}.}
}

\references{
Simes R. J. (1986) An improved Bonferroni procedure for multiple tests of
significance. \emph{Biometrika}, 73(3), 751--754.

Benjamini Y., Hochberg Y. (1995) Controlling the false discovery rate: a
practical and powerful approach to multiple testing. \emph{Journal of the Royal
Statistical Society, Series B}, 57(1), 289--300.
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
Results$gpv_adj
Results$selected
}
\keyword{ RNentropy }
\keyword{ Select }

\name{RN_iso_calc}
\alias{RN_iso_calc}
\title{
RN_iso_calc
}
\description{
This function calculates local p-values for changes in relative isoform usage,
output as -log10(p-value). At least two isoforms and two samples with expression
strictly greater than \code{min_expr} are required for a gene to be tested.
}
\usage{
RN_iso_calc(X, gene.col, design = NULL, min_expr = 1,
  pseudocount = 0.01)
}
\arguments{
  \item{X}{
data.frame with one row for each transcript, a column containing the
corresponding gene identifier, and numeric expression columns. It may contain
additional non numeric columns.
}
  \item{gene.col}{
The column name or number in \code{X} containing the gene identifiers.
}
  \item{design}{
The RxC design matrix where R (rows) corresponds to the number of numeric
columns (samples) in \code{X} and C (columns) to the number of conditions. It
must be a binary matrix with one and only one '1' for every row. Replicates of
the sample being tested are excluded from the reference distribution. design
defaults to a square matrix of independent samples (diagonal = 1, everything
else = 0).
}
  \item{min_expr}{
Minimum expression value used to determine whether a gene can be tested.
Expression must be strictly greater than this value (Default = 1).
}
  \item{pseudocount}{
Positive pseudocount added to the isoform frequencies before re-normalization
(Default = 0.01).
}
}

\value{
A list containing:
\item{expr}{The original expression data.frame.}
\item{design}{The experimental design matrix.}
\item{pv}{A gene by sample matrix containing -log10(p-values), with columns
named \code{ISO_PV_<sample>}. Values that cannot be calculated are reported as
\code{NA}.}
\item{gene_status}{A named vector reporting \code{TESTED},
\code{NOTEST(ISO)}, or \code{NOTEST(EXP)} for each gene.}
\item{sample_status}{A gene by sample matrix reporting \code{TESTED},
\code{NO_CURRENT_EXPRESSION}, or \code{NO_REFERENCE_EXPRESSION}. The last two
values correspond to '-' and '*' in the C++ implementation, respectively.}
}

\author{
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}

\examples{
Isoforms <- data.frame(
  gene = c("gene_1", "gene_1", "gene_2"),
  sample_1 = c(10, 1, 4),
  sample_2 = c(1, 10, 5)
)
Results <- RN_iso_calc(Isoforms, gene.col = "gene")
Results$gene_status
Results$pv
}
\keyword{ RNentropy }
\keyword{ Run }

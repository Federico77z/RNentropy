\name{RNentropy_iso_switch}
\alias{RNentropy_iso_switch}
\title{
RNentropy_iso_switch
}
\description{
This function runs the RNentropy isoform-switch analysis on a file containing
transcript expression values and the corresponding gene identifiers.
}
\usage{
RNentropy_iso_switch(file, tr.col, gene.col, design = NULL,
  header = TRUE, skip = 0, skip.col = NULL, col.names = NULL,
  min_expr = 1, pseudocount = 0.01)
}
\arguments{
  \item{file}{
Expression data file. It must contain a column with univocal transcript
identifiers and a column with the corresponding gene identifiers.
}
  \item{tr.col}{
The column name or number (1 based) in the file containing the univocal
transcript identifiers.
}
  \item{gene.col}{
The column name or number (1 based) in the original file containing the gene
identifiers.
}
  \item{design}{
The design matrix. Refer to the help for RN_iso_calc for further info.
}
  \item{header}{
Set this to FALSE if \code{file} does not have an header row. If header is TRUE
the read.table function used by RNentropy will try to assign names to each
column using the header (Default = TRUE).
}
  \item{skip}{
Number of rows to skip at the beginning of \code{file}. Rows beginning with a
'#' token will be skipped even if skip is set to 0.
}
  \item{skip.col}{
Columns that will not be imported from \code{file} (1 based) and not included
in the subsequent analysis.
}
  \item{col.names}{
Assign names to the imported columns. The number of names must correspond to
the columns effectively imported (thus, no name for the skipped ones). Also,
tr.col is not considered an imported column.
}
  \item{min_expr}{
Minimum expression value used to determine whether a gene can be tested
(Default = 1).
}
  \item{pseudocount}{
Positive pseudocount added to the isoform frequencies before re-normalization
(Default = 0.01).
}
}

\value{
Refer to the help for RN_iso_calc.
}

\author{
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}

\examples{
File <- system.file("extdata", "isoswitch_s7.tsv", package = "RNentropy")
Results <- RNentropy_iso_switch(File, tr.col = "TR_ID",
  gene.col = "GENE_ID")
Results$gene_status
}
\keyword{ RNentropy }
\keyword{ Run }

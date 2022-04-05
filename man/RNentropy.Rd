\name{RNentropy}
\alias{RNentropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RNentropy 
}
\description{
This function runs RNentropy on a file containing normalized expression values.
}
\usage{
RNentropy(file, tr.col, design, header = TRUE, skip = 0, skip.col = NULL, 
col.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
Expression data file. It must contain a column with univocal identifiers for each row (typically transcript or gene ids).
}
  \item{tr.col}{
%%     ~~Describe \code{tr.col} here~~
The column # (1 based) in the file containing the univocal identifiers (usually transcript or gene ids).
}
  \item{design}{
%%     ~~Describe \code{design} here~~
The RxC design matrix where R (rows) corresponds to the number of numeric columns (samples) in 'file' and C (columns) to the number of conditions. It must be a binary matrix with one and only one '1' for every row, corresponding to the condition (column) for which the sample corresponding to the row has to be considered a biological ot technical replicate. See the example 'RN_Brain_Example_design' for the design matrix of 'RN_Brain_Example_tpm' which has three replicates for three conditions (three rows) for a total of nine samples (nine rows).
}
  \item{header}{
%%     ~~Describe \code{header} here~~
Set this to FALSE if 'file' does not have an header row. If header is TRUE the read.table function used by RNentropy will try to assign names to each column using the header. Default = TRUE
}
  \item{skip}{
%%     ~~Describe \code{skip} here~~
Number of rows to skip at the beginning of 'file'. Rows beginning with a '#' token will be skipped even if skip is set to 0. 
}
  \item{skip.col}{
%%     ~~Describe \code{skip.col} here~~
Columns that will not be imported from 'file' (1 based) and not included in the subsequent analysis. Useful if you have more than annotation column, e.g. gene name in the first, transcript ID in the second columns.
}
  \item{col.names}{
%%     ~~Describe \code{col.names} here~~
Assign names to the imported columns. The number of names must correspond to the columns effectively imported (thus, no name for the skipped ones). Also, tr.col is not considered an imported column. Useful to assign sample names to expression columns.
}
}

\value{
Refer to the help for RN_calc.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
%%  ~~who you are~~
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}


\examples{
#load expression values and experiment design
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
#compute statistics and p-values (considering only a subset of genes due to
#examples running time limit of CRAN) 
Results <- RN_calc(RN_Brain_Example_tpm[1:10000,], RN_Brain_Example_design)
#select only genes with significant changes of expression
Results <- RN_select(Results)


## The function is currently defined as
function (file, tr.col, design, header = TRUE, skip = 0, skip.col = NULL, 
    col.names = NULL) 
{
    TABLE <- read.table(file, row.names = tr.col, header = header, 
        skip = skip, blank.lines.skip = TRUE, comment.char = "#")
    if (!is.null(skip.col)) {
        TABLE <- .RN_delete_col(TABLE, tr.col, skip.col)
    }
    if (!is.null(col.names)) {
        colnames(TABLE) <- col.names
    }
    return(RN_calc(TABLE, design))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")


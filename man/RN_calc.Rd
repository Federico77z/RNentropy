\name{RN_calc}
\alias{RN_calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
RN_calc
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes both global and local p-values, and returns the results in a list containing for each gene the original expression values and the associated global and local p-values (as  -log10(p-value)).
}
\usage{
RN_calc(X, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
data.frame with expression values. It may contain additional non numeric columns (eg. a column with gene names).
}
  \item{design}{
%%     ~~Describe \code{design} here~~
The RxC design matrix where R (rows) corresponds to the number of numeric columns (samples) in 'file' and C (columns) to the number of conditions. It must be a binary matrix with one and only one '1' for every row, corresponding to the condition (column) for which the sample corresponding to the row has to be considered a biological ot technical replicate. See the example 'RN_Brain_Example_design' for the design matrix of 'RN_Brain_Example_tpm' which has three replicates for three conditions (three rows) for a total of nine samples (nine rows). design defaults to a square matrix of independent samples (diagonal = 1, everything else = 0)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{gpv}{-log10 of the global p-values}
\item{lpv}{-log10 of the local p-values}
\item{c_like}{results formatted as in the output of the C++ implementation of RNentropy.}
\item{res}{The results data.frame with the original expression values and the associated -log10 of global and local p-values.}
\item{design}{the experimental design matrix}
}

\author{
%%  ~~who you are~~
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
#compute statistics and p-values (considering only a subset of genes due to
#examples running time limit of CRAN)
Results <- RN_calc(RN_Brain_Example_tpm[1:10000,], RN_Brain_Example_design)

## The function is currently defined as

function(X, design = NULL)
{
  if(is.null(design))
  {
    design <- .RN_default_design(sum(sapply(X, is.numeric)))
  }
  
	Results <- list(expr = X, design = design)

	GPV <- RN_calc_GPV(X, bind = FALSE)
	LPV <- RN_calc_LPV(X, design = design, bind = FALSE)	

	TABLE = cbind(X,'---',GPV,'---',LPV)	

	Results$gpv <- GPV
	Results$lpv <- LPV
	Results$c_like <- TABLE
	Results$res <- cbind(X, GPV, LPV)
	
	return(Results)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")
\keyword{ Run }% __ONLY ONE__ keyword per line

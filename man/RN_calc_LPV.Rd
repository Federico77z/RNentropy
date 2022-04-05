\name{RN_calc_LPV}
\alias{RN_calc_LPV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
RN_calc_LPV
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates local p-values from expression data, output as -log10(p-value).
}
\usage{
RN_calc_LPV(X, design, bind = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
data.frame with expression values. It may contain additional non numeric columns (eg. a column with gene names)
}
  \item{design}{
%%     ~~Describe \code{design} here~~
The design matrix. Refer to the help for RN_calc for further info.
}
  \item{bind}{
%%     ~~Describe \code{bind} here~~
See Value (Default = TRUE)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If bind is TRUE the function returns a data.frame with the original expression values from data.frame X and attached columns with the -log10(p-values) computed for each sample, otherwise only the matrix of -log10(p-values) is returned.
}

\author{
%%  ~~who you are~~
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}

\examples{
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
LPV <- RN_calc_LPV(RN_Brain_Example_tpm, RN_Brain_Example_design)
## The function is currently defined as

function(X, design = NULL, bind = TRUE)
{
  if(is.null(design))
  {
    design <- .RN_default_design(sum(sapply(X, is.numeric)))
  }
  
	rnums <- sapply(X, is.numeric)

	.RN_design_check(X, design, rnums)
	RL <- .RN_get_replicate_list(design)
	
	l1 <- rep("LOC_LPV", length(X[rnums]))
	l2 <- colnames(X[rnums])
	PV <- apply(X[rnums],  1, '.RN_calc_LPV_row', RL = RL)
	PV <- t(PV)
	colnames(PV) <- paste(l1, l2, sep="_")

	if(bind)
	{
		return (cbind(X, PV))
	}
	else
		return (PV)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")


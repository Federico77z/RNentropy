\name{RN_calc_GPV}
\alias{RN_calc_GPV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
RN_calc_GPV 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates global p-values from expression data, represented as -log10(p-values). 
}
\usage{
RN_calc_GPV(X, bind = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
data.frame with expression values. It may contain additional non numeric columns (eg. a column with gene names)
}
  \item{bind}{
%%     ~~Describe \code{bind} here~~
See Value (Default: TRUE)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If bind is TRUE the function returns a data.frame with the original expression values from 'X' and an attached column with the -log10() of the global p-value, otherwise only the numeric vector of -log10(p-values) is returned.
}

\author{
%%  ~~who you are~~
Giulio Pavesi - Dep. of Biosciences, University of Milan

Federico Zambelli - Dep. of Biosciences, University of Milan
}


\examples{
data("RN_Brain_Example_tpm")
GPV <- RN_calc_GPV(RN_Brain_Example_tpm)

## The function is currently defined as
function (X, bind = TRUE) 
{
    rnums <- sapply(X, is.numeric)
    GL_LPV <- apply(X[rnums], 1, ".RN_calc_GPV_row")
    if (bind) {
        GPV <- cbind(X, GL_LPV)
        return(GPV)
    }
    else {
        return(GL_LPV)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")


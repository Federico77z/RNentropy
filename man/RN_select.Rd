\name{RN_select}
\alias{RN_select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select transcripts/genes with significant p-values.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Select transcripts with global p-value lower than an user defined threshold and provide a summary of over- or under-expression according to local p-values. 
}
\usage{
RN_select(Results, gpv_t = 0.01, lpv_t = 0.01, method = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Results}{
%%     ~~Describe \code{Results} here~~
The output of RNentropy or RN_calc.
}
  \item{gpv_t}{
%%     ~~Describe \code{gpv_t} here~~
Threshold for global p-value. (Default: 0.01)
}
  \item{lpv_t}{
%%     ~~Describe \code{lpv_t} here~~
Threshold for local p-value. (Default: 0.01)
}
  \item{method}{
%%     ~~Describe \code{method} here~~
Multiple test correction method. Available methods are the ones of p.adjust. Type p.adjust.methods to see the list. Default: BH (Benjamini & Hochberg)
}
}

\value{
The original input containing

\item{gpv}{-log10 of the global p-values}
\item{lpv}{-log10 of the local p-values}
\item{c_like}{results formatted as in the output of the C++ implementation of RNentropy.}
\item{res}{The results data.frame containing the original expression values together with the -log10 of global and local p-values.}
\item{design}{The experimental design matrix.}

and a new dataframe

\item{selected}{Transcripts/genes with a corrected global p-value lower than gpv_t. For each condition it will contain a column where values can be -1,0,1 or NA. 1 means that all the replicates of this condition have expression value higher than the average and local p-value <= lpv_t (thus the corresponding gene will be over-expressed in this condition).  -1 means that all the replicates of this condition  have expression value lower than the average and local p-value <= lpv_t (thus the corresponding gene will be under-expressed in this condition). 0 means that at least one of the replicates has a local p-value > lpv_t. NA means that the local p-values of the replicates are not consistent for this condition, that is, at least one replicate results to be over-expressed and at least one results to be under-expressed.}
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
data("RN_Brain_Example_tpm", "RN_Brain_Example_design")
#compute statistics and p-values (considering only a subset of genes due to
#examples running time limit of CRAN)
Results <- RN_calc(RN_Brain_Example_tpm[1:10000,], RN_Brain_Example_design)
Results <- RN_select(Results)

## The function is currently defined as
function (Results, gpv_t = 0.01, lpv_t = 0.01, method = "BH") 
{
    lpv_t <- -log10(lpv_t)
    gpv_t <- -log10(gpv_t)
    Results$gpv_bh <- -log10(p.adjust(10^-Results$gpv, method = method))
    true_rows <- (Results$gpv_bh >= gpv_t)
    design_b <- t(Results$design > 0)
    Results$lpv_sel <- data.frame(row.names = rownames(Results$lpv)[true_rows])
    for (d in seq_along(design_b[, 1])) {
        col <- apply(Results$lpv[true_rows, ], 1, ".RN_select_lpv_row", 
            design_b[d, ], lpv_t)
        Results$lpv_sel <- cbind(Results$lpv_sel, col)
        colnames(Results$lpv_sel)[length(Results$lpv_sel)] <- paste("condition", 
            d, sep = "_")
    }
    lbl <- Results$res[, !sapply(Results$res, is.numeric)]
    Results$selected <- cbind(lbl[true_rows], Results$gpv[true_rows], 
	Results$gpv_bh[true_rows], Results$lpv_sel)
    colnames(Results$selected) <- c(names(which(!sapply(Results$res, 
	is.numeric))), "GL_LPV", "Corr. GL_LPV", colnames(Results$lpv_sel))
    Results$selected <- Results$selected[order(Results$selected[,3], decreasing=TRUE),]
    Results$lpv_sel <- NULL
    return(Results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")
\keyword{ Select }% __ONLY ONE__ keyword per line

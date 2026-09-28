\name{RN_pmi}
\alias{RN_pmi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute point mutual information matrix for the experimental conditions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute point mutual information for experimental conditions from the over-expressed genes identified by RN_select.
}
\usage{
RN_pmi(Results)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Results}{
%%     ~~Describe \code{Results} here~~
The output of RNentropy, RN_calc or RN_select. If RN_select has not already run on the results it will be invoked by RN_pmi using default arguments.
}
}

\value{
The original input containing

\item{gpv}{-log10 of the Global p-values}
\item{lpv}{-log10 of the Local p-values}
\item{c_like}{a table similar to the one you obtain running the C++ implementation of RNentropy.}
\item{res}{The results data.frame containing the original expression values together with the -log10 of Global and Local p-values.}
\item{design}{The experimental design matrix.}

\item{selected}{Transcripts/genes with a corrected Global p-value lower than gpv_t. Each condition N gets a condition_N column which values can be -1,0,1 or NA. 1 means that all the replicates of this condition seems to be consistently over-expressed w.r.t the overall expression of the transcript in all the conditions (that is, all the replicates of condition N have a positive Local p-value <= lpv_t).  -1 means that all the replicates of this condition seems to be consistently under-expressed w.r.t the overall expression of the transcript in all the conditions (that is, all the replicates of condition N have a negative Local p-value and abs(Local p-values) <= lpv_t). 0 means that one or more replicates have an abs(Local p-value) > lpv_t. NA means that the Local p-values of the replicates are not consistent for this condition.}

And two new matrices:
\item{pmi}{Point mutual information matrix of the conditions.}
\item{npmi}{Normalized point mutual information matrix of the conditions.}

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
Results <- RN_pmi(Results)
}

\keyword{ RNentropy }% use one of  RShowDoc("KEYWORDS")
\keyword{ PMI }% __ONLY ONE__ keyword per line

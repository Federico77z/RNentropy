\name{RN_IsoSwitch_Example_S7}
\alias{RN_IsoSwitch_Example_S7}
\docType{data}
\title{
RN_IsoSwitch_Example_S7
}
\description{
Transcript expression values from six tissues in the complete historical S7
example used to validate the RNentropy isoform-switch implementation. RefSeq
transcript identifiers are used as row names.
}
\usage{data("RN_IsoSwitch_Example_S7")}
\format{
  A data frame with 81314 observations on the following 7 variables.
  \describe{
    \item{\code{GENE_ID}}{a character vector containing gene identifiers}
    \item{\code{brain}}{a numeric vector}
    \item{\code{heart}}{a numeric vector}
    \item{\code{kidney}}{a numeric vector}
    \item{\code{liver}}{a numeric vector}
    \item{\code{lung}}{a numeric vector}
    \item{\code{muscle}}{a numeric vector}
  }
}
\source{
Historical S7 benchmark preserved with the independent C++ implementation of
RNentropy at
\url{https://github.com/Federico77z/RNentropy_cpp}.
}
\examples{
data(RN_IsoSwitch_Example_S7)
}
\keyword{datasets}

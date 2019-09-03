\name{untableSet}
\alias{untableSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
data.frame conversion
}
\description{
Converts a frequency table into a raw data.frame.
}
\usage{untableSet(data, freqvar = "Freq")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The \code{data.frame} including a frequency variable "Freq".
}

\item{freqvar}{
    The name of the frequency/weights variable which is used to expand \code{V}.
}
}

\value{
A \code{data.frame}.
}
\references{
  Alexander Pilhoefer \emph{ New approaches in visualization of categorical data:
       R-package extracat}\cr
  Journal of Statistical Software, submitted March 2010
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}



\examples{
require(MASS)

hs2 = untableSet(housing)
summary(hs2)
}

\keyword{ frequency table }


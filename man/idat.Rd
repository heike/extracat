\name{idat}
\alias{idat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
indicator dataframe
}
\description{
Converts all factor variables in a dataframe to a set of binary variables.}
\usage{
idat(x, allcat = FALSE, keep = "Freq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
dataframe
}
\item{allcat}{
Whether or not to keep all categories or leave the last one out.
}
\item{keep}{
Variables which are kept unchanged such as a frequency variable.
}
}

\value{
dataframe
}


\seealso{
\link[extracat]{Burt}, \link[extracat]{imat}
}
\examples{
require(MASS)
idat(housing)
}


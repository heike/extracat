\name{imat}
\alias{imat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
indicator variables
}
\description{
converts a single categorical variable into indicator variables
}
\usage{
imat(x, allcat = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A factor variable.
}
\item{allcat}{
Whether or not to keep all categories or leave the last one out.
}
}

\value{
matrix
}

\examples{
require(MASS)
imat(housing$Type)
}


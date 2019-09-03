\name{cohen}
\alias{cohen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cohens Kappa for rectangular matrices
}
\description{
Cohen's Kappa for quadratic and non-quadratic matrices using L1-weights.
}
\usage{
cohen(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix with non-negative entries.
}
}

\value{
Cohen's Kappa
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{kendalls}, \link[extracat]{BCI}, \link[extracat]{WBCI}
}
\examples{
a <- arsim(2000,c(12,12),6)
cohen(a)
cohen(optile(a))
}

\name{WBCI}
\alias{WBCI}
\title{
The Weighted Bertin Classification Index
}
\description{
The weighted Bertin Classification Criterion using weights according to the Hamming distance is normalized by means of the independence case.
}
\usage{
WBCI(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix.
}
}
\value{
The criterion value.
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{kendalls}
}
\examples{
M <-arsim(1000, c(12,12), 3)
BCI(M)
WBCI(M)


M2 <- optile(M, iter = 100)
BCI(M2)
WBCI(M2)

M3 <- optile(M, fun = "WBCC", iter = 100)
BCI(M3)
WBCI(M3)
}


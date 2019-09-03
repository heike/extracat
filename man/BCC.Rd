\name{BCC}
\alias{BCC}


\title{
The Bertin Classification Criterion
}
\description{
Computes the Bertin Classification Criterion for a contingency table of any dimensions.
}
\usage{
BCC(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data matrix, table or array.
}
}
\details{
The BCC counts the number of observation pairs which differ in all variables but are not fully concordant, (i.e. neither of the two observations of each pair is larger than the other in all variables).
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
BCC(M)

M2 <- optile(M, iter = 100)
BCC(M2)



}


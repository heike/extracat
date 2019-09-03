\name{kendalls}
\alias{kendalls}

\title{
Kendalls Tau for a matrix
}
\description{
Computes Kendalls Tau for a two-way table or matrix.
}
\usage{
kendalls(x)
}
\arguments{
  \item{x}{
A two-way table or matrix.
}
}
\details{
Kendalls tau is a rank-correlation coefficient.
}
\value{
numeric between -1 and +1.
}

\author{
Alexander Pilhoefer
}

\examples{
M <- arsim(300,c(8,8),3)
kendalls(M)
kendalls(optile(M))
}

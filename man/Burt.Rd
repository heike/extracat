\name{Burt}
\alias{Burt}
\title{
Burt atrix
}
\description{
The Burt matrix is a quadratic matrix where each row and column corresponds to a category in one of the variables. The entries of the matrix are the frequencies of the corresponding combination of categories.
}
\usage{
Burt(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A dataframe with factor variables or a contingency table.
}
}

\value{
A matrix.
}

\author{
Alexander Pilhoefer
}



\seealso{
\link[extracat]{idat}, \link[extracat]{imat}
}
\examples{
require(MASS)
Burt(housing)
th <- xtabs(Freq~Sat+Infl+Type, data = housing)
Burt(th)

}


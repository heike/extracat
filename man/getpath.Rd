\name{getpath}
\alias{getpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Path extracation from quickfechner objects
}
\description{
Uses the path.matrix to obtain the shortest paths of the quickfechner object.
}
\usage{
getpath(fm, pm = NULL, from = 1, to = nrow(fm))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fm}{
A Fechner matrix computed via \link[extracat]{quickfechner}
}
  \item{pm}{
The path.matrix. Only necessary if the corresponding attribute in fm is missing.
}
  \item{from}{
Start index.
}
  \item{to}{
Final index.
}
}
\details{
The path.matrix is defined as follows: The entry of the i-th row and j-th column is the index of the last node visited before j on the shortest path from i to j.
}
\value{
A vector of indices defining the shortest path from i to j in the original matrix passed to the quickfechner function.
}

\author{
Alexander Pilhoefer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#not a distance matrix, but a similarity matrix in some sense
cx <- 1-abs(cor(olives[-c(1,2,11)]))

cx2 <- quickfechner(cx)
getpath(cx2,from=1,to=5)

}


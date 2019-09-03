\name{getIs}
\alias{getIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ME reordering for biclust objects.
}
\description{
Basically this is an auxiliary function used by \link[extracat]{heattile}:
It extracts the cluster indices from a biclust object and optimizes the order of the rows and columns in the data matrix with respect to these clusters. Uses the Measure of Effectiveness as an optimization criterion. See \link[extracat]{ME} and \link[extracat]{optME}.
}
\usage{
getIs(biclust, dim, nstart = 20, solver = "nn", adjust.dist = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biclust}{
The biclust object.
}
  \item{dim}{
The dimension of the matrix.
}
  \item{nstart}{
Number of starting points for the TSP solver in \link[extracat]{optME}.
}
  \item{solver}{
The TSP solver to use with \link[extracat]{optME}: See \link[TSP]{solve_TSP}.
}
  \item{adjust.dist}{
If \code{TRUE} the ME values used as a distance matrix for the TSP are slightly adjusted by adding a the hamming distance divided by a constant.
This keeps identical cases together and also preserves the orders within such groups.
}
}
\details{
The algorithm first computes an indicator matrix for each cluster and then combines these matrices to a 3D table.
Then for the first and the second dimension the category orders are optimized with respect to \link[extracat]{ME}.
The optimization is done via \link[extracat]{optME} which uses a TSP solver.
}
\value{
The cluster indices with respect to the optimized row and column orders in form of a list.
The optimized orders for the data matrix are attached as an attribute \code{attr(x, "orders")}.
}

\author{
Alexander Pilhoefer
}



\seealso{
 \link[extracat]{heattile},  \link[extracat]{getIs2}
}


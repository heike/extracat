\name{optME}
\alias{optME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizing ME
}
\description{
Computes optimal category orders for each dimension separately. Uses a TSP solver to achieve the best ME value.}
\usage{
optME(x, dims = NULL, nstart = 1,
 solver = "nearest_insertion",
  return.table = TRUE, adjust.dist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix, table or array.
}
  \item{dims}{
Which dimensions to reorder. The dimensions are reordered independently.
}
  \item{nstart}{
The number of different starting points for the TSP solver. If \code{nstart} is greater or equal to the number of cities in a dimension, the solver uses each city once.
}
  \item{solver}{
Should be one of \code{"nn", "nearest_insertion", "cheapest_insertion","farthest_insertion"}.  See \link[TSP]{solve_TSP}.
}
  \item{return.table}{
Whether or not to return the optimized table. If \code{FALSE} only the new category orders are returned. If \code{TRUE} the new orders are attached to the table as an attribute \code{"orders"}.
}
  \item{adjust.dist}{
If \code{TRUE} a small proportion of the euclidean distances between the category profiles (e.g. rows) is added to the ME distance value.
The idea is to keep identical profiles together which is otherwise not guaranteed, since the ME values can be identical even if the profiles aren't.
}
}
\details{
Each dimension is optimized separately via a TSP solver.
}
\value{
The passed object as a table with optimized category orders.
}

\seealso{
\link[extracat]{ME}
}
\examples{
a <- arsim(2000,c(8,9,10),3,0.2)
ME(a)
a2<-optME(a)
ME(a2)
}



\name{sortandcut}
\alias{sortandcut}

\title{
Sort-and-Cut Reordering
}
\description{
An implementation of the sort-and-cut algorithm which is a mixture of the top-down-partitioning algorithm used by \link[extracat]{cfluctile} and the reordering techniques available through \link[extracat]{optile}. See details.
}
\usage{
sortandcut(x,iter=20, tau0 = NULL, fun = "BCC", method = "WBCI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix or 2D table.
}
  \item{iter}{
The number of random initial orderings for \link[extracat]{optile}.
}
  \item{tau0}{
The minimum criterion value for a new cut. See \link[extracat]{cfluctile}.
}
  \item{fun}{
The reordering function used by \link[extracat]{optile}. Currently \code{"BCC"},  \code{"barysort"} and \code{"preclass"} are available.
}
  \item{method}{
The method argument for cfluctile which defines the criterion used to find an optimal partition.
}
}
\details{
The algorithm sorts a matrix using \link[extracat]{optile} and cuts the reordered matrix once using \link[extracat]{cfluctile} with \code{nsplit = 1}.
Then the same procedure is applied to the resulting submatrices at the top left and the bottom right. The partitioning stops when the best cut leads to a criterion below \code{tau0}.
}
\value{
The reordered matrix. The row and column order vectors are attached as an attributes \code{attr(x ,"orders")}.
}



\seealso{
\link[extracat]{cfluctile}
}
\examples{
M <- arsim(12000,c(30,40),7,noise=0.3)
c1 <- cfluctile(M1<-optile(M, iter = 20))
c2 <- cfluctile(M2<-sortandcut(M))

}

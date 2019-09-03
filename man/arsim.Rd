\name{arsim}
\alias{arsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
block-structured arrays}
\description{
Generates an array or matrix that includes k fully separated block-clusters.
}
\usage{
arsim(n, dim, k, noise = 0, shuffle = TRUE, v = 0.1, minc = 1,
 exp.prop = NULL, min.prop = 1/dim/4, noise.type = "s",
 dimnames=list(LETTERS,1:max(dim)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of observations in the array.
}
  \item{dim}{
The dimension of the array.
}
  \item{k}{
The number of clusters. 1 for no clusters.
}
  \item{noise}{
The proportion of noise among the observations. There are two choices for \code{noise.type}.
}
  \item{shuffle}{
Whether or not to shuffle the original category orders randomly.
}
  \item{v}{
A variability parameter for the assignment of the observations to the block clusters. Small values lea
}
  \item{minc}{
The minimum number of categories each cluster must have in each variable. E.g. \code{minc = 2} means, that each block cluster covers at least 2 categories in each dimension.
}
  \item{exp.prop}{
Optional: expected proportions of the observations which fall into the block clusters.
}
  \item{min.prop}{
Minimum proportion of observations in each cluster.
}
  \item{noise.type}{
Either \code{"s"} or \code{"I"}. The noise type \code{"s"} means that \code{n*noise} observations are drawn at random from the block-diagonal matrix. Then for these observations the category labels are permuted at random. \code{"I"}  adds noise in form of a random sample from the independence matrix with the same marginal totals as the block matrix.
}
\item{dimnames}{
A list of 2: The first entry defines the variable labels (default: A,B,C,...) and the second entry defines the category labels (default 1:k).
}
}
\details{
Not a very sophisticated way of generating random arrays but it serves for tests and illustrations of the other functions.
}
\value{
A simulated data array.
}
\examples{
A <- arsim(1000, c(12,12), 3, shuffle = FALSE)
fluctile(A)

A <- arsim(1000, c(12,12), 3, shuffle = FALSE, dimnames = list(NULL,letters))
dimnames(A)

A <- arsim(4000, c(11,7,5), 3, shuffle = TRUE, dimnames = list(0:2,letters))
dimnames(A)

\dontrun{
A2<- arsim(1000, c(12,12,12), 3, shuffle = FALSE)
fluctile3d(A2, shape ="oct")
}

}


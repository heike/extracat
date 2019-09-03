\name{setcover}
\alias{setcover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
greedy setcover optimisation
}
\description{
This function takes an indicator matrix with rows representing objects and columns representing sets and computes a minimal redundancy free set using the greedy setcover optimization algorithm. The aim is to find a minimal set of clusters which covers all objects (or a minimum proportion \code{rat}).

Alternatively the number of clusters \code{k} can be specified. Then the problem becomes a maximum covergae problem. Both versions also permit weights such as frequencies (weighted setcover/maximum coverage).
}
\usage{
setcover(x, k = NULL, rat = 1, s = NULL, w = NULL, check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The indicator matrix.
}
  \item{k}{
An optional number of clusters.
}
 \item{rat}{
The minimum proportion of objects that is to be covered by the cluster set. If weights are specified in \code{w} then those are respected.
}
\item{s}{
If weights are specified but not all objects are covered by one of the sets it can be necessary to specify the total weight in order to compute a sensible ratio.
}
  \item{w}{
Optional weights per object.
}
  \item{check}{
Whether or not to check for redundancies.
}
}

\value{
The indices of the clusters in the minimal redundancy-free set.
The result is not always the globally optiomal solution since the algorithm is greedy.
}

\author{
Alexander Pilhoefer
}
\note{
This is written supporting the GSAC algorithm.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{gsac}
}
\examples{
# compute 100 clusterings with 24 clusters each:
sc <- scale(olives[,3:10])
km100 <- as.data.frame(replicate(100, kmeans(sc,centers = 24)$cluster))

# convert to indicator matrix
I100 <- idat(km100)

# select from all clusters a minimum set:
scover <- setcover(as.matrix(I100))


cdata <- subtable(
	as.data.frame(cbind(olives[,1:2],
	I100[,scover])),1:(length(scover)+2))
scpcp(cdata,sel="Area")



}

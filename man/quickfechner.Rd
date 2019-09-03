\name{quickfechner}
\alias{quickfechner}

\title{
fechnerian scaling
}
\description{
This function computes a fechnerian distance matrix from either a similarity matrix or a dissimilarity matrix. In addition to the basic procedure which looks for the shortest paths between the objects in the dissimilarity matrix a second approach is offered which connects similarities in a multiplicative manner. 
}
\usage{
quickfechner(x, x.type = "diss", scale ="-", path.op = "+", sym.op = "+", 
rescale = FALSE, exclude.zero = FALSE, check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A similarity or dissimilarity matrix.
}
  \item{x.type}{
The type of the matrix (\code{'sim'} or \code{'diss'}). 
}
  \item{scale}{
Either divide the similarities by the diagonal entries 
(\code{'div'}, \code{'/'},\code{'*'},\code{'exp'},\code{'expected'},\code{'mult'},\code{'multiplicative'}) 
or subtract the diagonal entries in the dissimilarity matrix
  (\code{'-'}, \code{'+'}, \code{'add'},\code{'additive'})
}
  \item{path.op}{
Whether to use the similarities to find multiplicative paths (\code{'*'},\code{'exp'},\code{'expected'},\code{'mult'},\code{'multiplicative'}) or to use the dissimilarities and find additive paths (\code{'+'},\code{'add'},\code{'additive'},\code{'max'},\code{'maximum'}). 
Note that similarity matrices are simply converted to dissimilarity matrices by subtraction from 1. Other transitions such as 2M/(1+M) are not yet implemented and have to be done by hand.
}
  \item{sym.op}{
This sets the function which is used to ensure symmetry. \code{"min"} uses the minumum value, \code{"+"}, \code{"sum"} or \code{"mean"} use the sum. \code{"none"}, \code{NA} or \code{FALSE} stand for no operation and hence the resulting matrix will not necessarily be symmetric.
}
  \item{rescale}{
Whether or not the original diagonal will be used for a correction of the results.
}
  \item{exclude.zero}{
If \code{TRUE} zero-entries are not considered in the updating algorithm.
}
  \item{check}{
Whether or not to check for regular minimality or maximality.
}
}
\details{
The algorithm first computes a dissimilarity matrix with a zero-diagonal. Then it iteratively tries to find shorter paths between the items.
}
\value{
The Fechnerian distance matrix.
}

\author{
Alexander Pilhoefer
}

\examples{
data(olives)
#not a distance matrix, but a similarity matrix in some sense
cx <- 1-abs(cor(olives[-c(1,2,11)]))

cx2 <- quickfechner(cx)

rownames(cx2) <- names(olives)[-c(1,2,11)]
plot(hclust(as.dist(cx2)))

dm <- matrix(runif(100),10,10)
dm <- dm+t(dm)
diag(dm) <- 0
dm2 <- quickfechner(dm)

dmS <- 1-dm/max(dm)
dmS2 <- quickfechner(dmS, x.type="sim", path.op = "*")

\dontrun{
# check triangular inequality:
	extracat:::trinq(dm)
	extracat:::trinq(dm2)
	extracat:::trinq(dmS2)
}
}


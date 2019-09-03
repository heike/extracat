\name{combcl}
\alias{combcl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine categorical variables from cfluctile and cflcl
}
\description{
Combines variables obtained via \link[extracat]{cfcl} and \link[extracat]{cfluctile} to a single \code{factor} variable with one level per block-cluster and one level for the rest.
}
\usage{
combcl(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix, table or data.frame. All variables should have the same number of categories.}
}


\value{
A factor variable with 1 level per diagonal element and 1 level for the rest.}

\examples{

a <- arsim(2000, c(12,17),5, noise=0.2,shuffle = FALSE)
cfa <- cfluctile(a)

da <- as.data.frame(a)
clusters <- cfcl( da, ll = cfa)

dev.new()
fluctile(xtabs(da$Freq~clusters[,1] + clusters[,2]))

table(combcl(clusters))

}

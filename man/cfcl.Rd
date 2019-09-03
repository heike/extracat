\name{cfcl}
\alias{cfcl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract clusters from cfluctile
}
\description{
Extract clusters from cfluctile
}
\usage{
cfcl(x, y = NULL, ll)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector or dataframe.
}
  \item{y}{
if \code{x} is a vector, \code{y} needs to be specified.
}
  \item{ll}{
The list with the names of the levels which are combined.
}
}

\value{
A 2-column dataframe with the cluster factors.
}
\seealso{
\link[extracat]{cfluctile}
}
\examples{
a <- arsim(2000, c(12,17),5, noise=0.2,shuffle = FALSE)
cfa <- cfluctile(a)

da <- as.data.frame(a)
clusters <- cfcl( da, ll = cfa)

dev.new()
fluctile(xtabs(da$Freq~clusters[,1] + clusters[,2]))

table(combcl(clusters))
}


\name{dcorMVtable}
\alias{dcorMVtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Distance Correlation for two sets of variables
}
\description{
Computes the distances within two sets of variables and the corresponding distance correlation.
}
\usage{
dcorMVtable(x, ind = 1, method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A contingency table of class \code{table}.
}
  \item{ind}{
The indices for the first set of variables. The second set consists of all remaining variables.
}
  \item{method}{
The method for \link[stats]{dist}
}
}

\value{
The distance correlation between 0 and 1 for the distances from the two sets of variables.
}

\author{
Alexander Pilhoefer
}
\note{
This code has not been tested thoroughly and may still contain errors.
}


\seealso{
\link[extracat]{dcorMVdata}, \link[extracat]{wdcor}, \link[extracat]{approx.dcor}
}
\examples{
\dontrun{
	A2 <- arsim(2000,c(8,9),5,0.1)
	A2 <- optile(A2, iter=100)
	BCI(A2)
	wdcor(A2)
	
	p1 <- runif(11)+0.1
	p1 <- p1/sum(p1)
	A2b <- apply(A2,1:2,function(z) rmultinom(1,z,p1))
	
	# now the first variable is roughly independent from the other two:
	
	dcorMVtable(as.table(A2b),ind = 1)
	
	# here the third variable is NOT independent from the others:
	A3 <- arsim(2000,c(8,9,11),5,0.1)
	A3 <- optile(A3, iter=100)
	BCI(A3)
	dcorMVtable(A3,ind = 3)
	
	
}
}


\name{dcorMVdata}
\alias{dcorMVdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Distance Correlation for two sets of variables
}
\description{
Computes the distances within two sets of variables and the corresponding distance correlation.
}
\usage{
dcorMVdata(x, ind = 1, method = "euclidean", approx = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The \code{data.frame} which should only contain non-factor variables. For factor variables use \link[stats]{xtabs} in combination with \link[extracat]{dcorMVtable}.
}
  \item{ind}{
The indices for the first set of variables. The second set consists of all remaining variables.
}
  \item{method}{
The method for \link[stats]{dist}.
}
  \item{approx}{
\code{FALSE} for no approximation via binning or an integer value for the number of bins.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{dcorMVtable}, \link[extracat]{wdcor}, \link[extracat]{approx.dcor}
}
\examples{
\dontrun{
	so <- scale(olives[,3:8])
	dcorMVdata(so,ind=1)
	
	p1 <- princomp(so)
	so1 <- cbind(so,p1$scores[,1])
	so2 <- cbind(so,p1$scores[,2])
	so12 <- cbind(so,p1$scores[,1:2])
	
	dcorMVdata(so1,ind=7)
	dcorMVdata(so2,ind=7)
	dcorMVdata(so12,ind=7:8)
	# how about principal curves?
}
}


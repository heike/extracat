\name{cmat}
\alias{cmat}

\title{
pairwise association matrix
}
\description{
Computes pairwise BCI values via \link[extracat]{qBCI}. 
}
\usage{
cmat(x, sort = TRUE, crit = BCI, k = 5, iter = 20, 
p = NULL,  jitter = TRUE, freqvar = NULL, diag = NULL,
fun = "BCC", foreign = NULL)
}

\arguments{
  \item{x}{
A \code{data.frame} with \code{factor} variables or numeric variables which will be transformed to ordinal interval variables via \link[base]{cut}.
The breakpoints are quantiles of the variables such that for each pair of numeric variables the expected number of observations in each combination of intervals is at least \code{k}. }
\item{sort}{
Whether or not to sort the pairwise tables via \link[extracat]{optile}.
}
  \item{crit}{
The criterion function, e.g. \code{kendalls},  \code{BCI},  \code{WBCI} or  \code{wdcor}.
}
  \item{k}{
The minimum expected number for each cell after quantile binning. See also \link[extracat]{qBCI}.
}
  \item{iter}{
An \link[extracat]{optile} parameter.
}
  \item{p}{
The quantile distance. See \link[extracat]{qBCI}.
}
  \item{jitter}{
Whether or not to use jittering in order to avoid ties. This is equivalent to a random assignment of ranks to observations with the same value.
}
  \item{freqvar}{
Optional weights, e.g. a frequency variable.
}
  \item{diag}{
An optional value for the diagonal. Avoids unnecessary function calls for the diagonal elements. E.g. \code{diag = 0} for \code{crit = BCI} or \code{diag = 1} for \code{crit = kendalls} makes sense.
}
  \item{fun}{
See \link[extracat]{optile}.
}
  \item{foreign}{
See \link[extracat]{optile}.
}
}
\details{
Uses pairwise complete cases only!
}
\value{
A symmetric matrix.
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{qBCI}, See \link[extracat]{wdcor}.
}
\examples{
\dontrun{
m1 <- cmat(olives)
fluctile(1 - m1,shape="o")
}
}

\name{ahist}
\alias{ahist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histogram using active bins
}
\description{
A standard histogram using \link[extracat]{getbw} to compute the binwidth and breakpoints. 
}
\usage{
ahist(x, k = NULL, m = NULL, fun = "qplot", col = "grey", ival = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector.
}
  \item{k}{
The desired number of active bins. A bin is active if it contains at least \code{min_n} observations. The default is
\code{k <- 1 + 2*ceiling(log(N)/log(2))}.
}
  \item{m}{
The minimum number of observations necessary for a bin to count as an active bin. Defaults to
\code{m = max(log(N/10)/log(10),1)}.
}
  \item{fun}{
Either \code{"qplot"} or \code{"hist"}.
}
  \item{col}{
The color for the bars.
}
  \item{ival}{
If this is set to a numeric value in \code{(0,1)} then \code{x} is trimmed according to
\code{innerval(x, p = ival)}.
}
}

\value{
The \code{ggplot} object.
}

\author{
Alexander Pilhoefer
}
\note{
This is purely experimental at this time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{getbw}, \link[extracat]{cutbw}
}
\examples{
ahist(rnorm(100))
ahist(rnorm(1000))
ahist(rnorm(10000))


ahist(rexp(100))
ahist(rexp(1000))
ahist(rexp(10000))


\dontrun{
ahist(rcauchy(1000))
ahist(rcauchy(1000), ival = 0.95)

x <- c(rnorm(400),rnorm(200, mean=6))
ahist(x)

x <- c(rnorm(400),rnorm(200, mean=16))
ahist(x)


x <- c(rnorm(400),rnorm(200, mean=32))
ahist(x)
}
}


\name{cutbw}
\alias{cutbw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Active binning
}
\description{
Uses \link[base]{cut} with breakpoints derived by \link[extracat]{getbw}.
}
\usage{
cutbw(x, k = NULL, min_n = NULL, warn = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric variable.
}
  \item{k}{
The desired number of active bins. A bin is active if it contains at least \code{min_n} observations. The default is
\code{k <- 1 + 2*ceiling(log(N)/log(2))}.
}
  \item{min_n}{
The minimum number of observations necessary for a bin to count as an active bin. Defaults to
\code{min_n = max(log(N/10)/log(10),1)}.
}
  \item{warn}{
Whether or not to print a warning if the desired number of bins is not possible.
}
}

\value{
An ordinal factor variable.
}

\author{
Alexander Pilhoefer
}
\note{
Experimental.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{getbw}, \link[extracat]{ahist}
}
\examples{
y<-cutbw(c(rnorm(200),rnorm(100,mean=8)),k = 30, min_n = 1)
barplot(table(y))
}


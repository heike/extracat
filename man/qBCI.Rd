\name{qBCI}
\alias{qBCI.data.frame}
\alias{qBCI.default}
\alias{qBCI}
\title{
Quantile BCI
}
\description{
Bins numeric variables according to their quantiles and computes the Bertin Classification Index \link[extracat]{BCI}.
The \code{data.frame} method computes the multivariate qBCI and not the pairwise values (c.f. \link[extracat]{cmat}).
}
\usage{
qBCI(x, ...)
\method{qBCI}{default}(x, y, p = NULL, k = 5, iter=20, \dots)
\method{qBCI}{data.frame}(x,p = NULL, k = 5, sort = TRUE, iter=20, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector (in this case \code{y} needs to be specified) or a data.frame with numeric or factor variables.
}
\item{y}{
A numeric vector.
}
\item{p}{
A percentage to use for the quantiles sequence. See details.
}
\item{k}{
A minimum expected number of observations in each cell after the binning.
}
\item{sort}{
Whether or not to compute the BCI for the optimized tables or not. If not, \link[extracat]{kendalls} is usually a better alternative.
}
\item{iter}{
An \link[extracat]{optile} parameter.
}
  \item{\dots}{
dots
}
}
\details{
The breakpoints for the binning are the data quantiles according to equidistant probabilities \code{seq(0,1,p)} where \code{p} is minimal under the condition that each cell has an expected number of observations of at least \code{k}.
}
\value{
A value between 0 and 1.
}

\author{
Alexander Pilhoefer
}

\seealso{
\link[extracat]{BCI},\link[extracat]{kendalls}, \link[extracat]{wdcor}, \link[extracat]{cmat}
}
\examples{
\dontrun{
qBCI(rnorm(100),runif(100))


# non-functional relationship:
x1 <- runif(500,0,10)
x2 <- runif(500,0,10)
y1 <- x1+rnorm(500,sd=1) 
y2 <- 10-x2+rnorm(500,sd=1) 

x <- c(x1,x2)
y <- c(y1,y2)

plot(x,y, pch = 19)

wdcor(x,y)
1 - qBCI(x,y)


y1 <- x1+rnorm(500,sd=0.1) 
y2 <- 10-x2+rnorm(500,sd=0.1) 

x <- c(x1,x2)
y <- c(y1,y2)

plot(x,y, pch = 19)

wdcor(x,y)
1 - qBCI(x,y)

# or a quadratic curve:
test <- sapply(seq(0,4,0.2),function(s){
x <- runif(200,-1,1)
y <- 5+12*x^2+rnorm(200,sd=s)
return(c(cor(x,y),
wdcor(x,y),
1 - qBCI(x,y)))
})


plot(test[3,],type="l", ylim=c(-0.2,1))
lines(test[1,], col = 2, lwd = 2)

lines(test[2,], col = 3, lwd = 2)

}

}

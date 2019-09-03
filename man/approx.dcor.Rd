\name{approx.dcor}
\alias{approx.dcor}
\title{
Distance Correlation Approximation
}
\description{
Computes the distance correlation for two variables using an approximation based on binning and \link[extracat]{wdcor.table}.
The approximation underestimates the true value by a small error depending on the number of bins. (In simulations with the default of 50 bins the average error was about 0.001.) }
\usage{
approx.dcor(x, y, n = 50, ep = 1, bin = "eq")
}
\arguments{
  \item{x}{
A numeric vector.
}
  \item{y}{
A numeric vector.
}
  \item{n}{
The number of bins per variable.
}
 
\item{ep}{
The euclidean distances are taken to the power of \code{ep}.
}
\item{bin}{
Either \code{"eq"} or \code{"q"} for equidistant breakpoints or quantile breakpoints.
}
}

\value{
The correlation value which is between 0 and 1.
}
\references{
Szekely, G. J. Rizzo, M. L. and Bakirov, N. K. (2007). "Measuring and testing independence by correlation of distances", Annals of Statistics, 35/6, 2769-2794
}
\author{
Alexander Pilhoefer
}




\examples{
\dontrun{

# The straightforward way of approximating the distance correlation fails:
# for instance the computation of dcor for a random sample with 4000 observations
# takes pretty long but drawing samples of 500, 1000 or even 2000 observations
# leads to relatively big errors.
# The approximation via approx.dcor is very fast and for 
# n = 50 or n=100 the results are very close to the truth

require(energy)
x<- rnorm(4000,mean=10,sd=3)
y <- rnorm(1,sd=0.01)*(x-10)^3 + rnorm(1,sd=0.1)*(x-10)^2
 + rnorm(1)*(x-10)+rnorm(4000,sd=4)

system.time(dd <- dcor(x,y))
system.time(dd0 <- wdcor(x,y))[[3]]
dd - dd0


system.time(da100 <- approx.dcor(x,y,100))[[3]]
da100-dd0
	
# For a smaller sample size we can try out how good the approximation really is:
test<-replicate(100,{
N <- 1000
x<- rnorm(N,mean=10,sd=3)
y <- rnorm(1,sd=0.01)*(x-10)^3 + rnorm(1,sd=0.1)*(x-10)^2 
y <- y + rnorm(1)*(x-10)+rnorm(N,sd=4)

dd <- wdcor(x,y)
dd25 <- approx.dcor(x,y,25)
dd50 <- approx.dcor(x,y,50)
dd100 <- approx.dcor(x,y,100)
dd75 <- approx.dcor(x,y,75)

dd25c <- approx.dcor(x,y,25,correct = TRUE)
dd50c <- approx.dcor(x,y,50,correct = TRUE)
dd100c <- approx.dcor(x,y,100,correct = TRUE)
dd75c <- approx.dcor(x,y,75,correct = TRUE)
c(2*dd, dd25, dd50, dd75, dd100, dd25c, dd50c, dd75c, dd100c)-dd
})

rm<-apply(test,1,mean)

plot( seq(25,100,25), rm[2:5],type="l",
 ylim= c(min(rm),abs(min(rm))), xlab = "No. of bins per axis",ylab = "error")
points( seq(25,100,25), rm[2:5],pch=19)
points( seq(25,100,25), rm[6:9],type="l", col=2)
points( seq(25,100,25), rm[6:9],pch=19,col=2)	
abline(h=0,lwd=3)
legend( 25,abs(min(rm)),legend=c("raw value after binning","corrected value"),
col=1:2,lwd=3)
}
}


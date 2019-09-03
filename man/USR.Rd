\name{USR}
\alias{USR}
\docType{data}
\title{
MovieLens USER data
}
\description{
The user data from the MovieLens 1M data.
}
\usage{data(USR)}
\format{
  A data frame with 6040 observations on the following 25 variables.
  \describe{
    \item{\code{UserID}}{a numeric vector}
    \item{\code{Gender}}{a factor with levels \code{F} \code{M}}
    \item{\code{Age}}{a numeric vector}
    \item{\code{Occupation}}{a numeric vector}
    \item{\code{Zip.code}}{a factor}
    \item{\code{occupation}}{a factor}
    \item{\code{zip1}}{a numeric vector}
    \item{\code{zip2}}{a numeric vector}
    \item{\code{zip12}}{a numeric vector}
    \item{\code{UserVotes}}{a numeric vector}
    \item{\code{meanUserRat}}{a numeric vector}
    \item{\code{medianUserRat}}{a numeric vector}
    \item{\code{sdUserRat}}{a numeric vector}
    \item{\code{minTs}}{a numeric vector}
    \item{\code{AvgRat}}{a numeric vector}
    \item{\code{meanDiffRat}}{a numeric vector}
    \item{\code{meanAbsDiffRat}}{a numeric vector}
    \item{\code{Pct.Action}}{a numeric vector}
    \item{\code{Pct.Adventure}}{a numeric vector}
    \item{\code{Pct.Animation}}{a numeric vector}
    \item{\code{Pct.Children.s}}{a numeric vector}
    \item{\code{Pct.Fantasy}}{a numeric vector}
    \item{\code{Pct.Horror}}{a numeric vector}
    \item{\code{Pct.Sci.Fi}}{a numeric vector}
    \item{\code{Pct.Comedy}}{a numeric vector}
  }
}
\source{
MovieLens.org
}
\examples{
data(USR)
## maybe str(USR) ; plot(USR) ...
}
\keyword{datasets}

\name{CPScluster}
\alias{CPScluster}
\docType{data}
\title{
Clusterings for the US Current Population Survey.
}
\description{
Different hierarchical clusterings and k-means clusterings as well as a model-based clustering have been applied to several financial variables for a random sample of ten thousand observations.
}
\usage{data(CPScluster)}
\format{
  A data frame with 10000 observations on the following 39 variables.
  \describe{
    \item{\code{Age}}{a numeric vector}
    \item{\code{Sex}}{a factor with levels \code{female} \code{male}}
    \item{\code{Race}}{a factor with levels \code{Black} \code{White}}
    \item{\code{Ethnic}}{a factor}
    \item{\code{Marital.Status}}{a factor}
    \item{\code{Kind.of.Family}}{a factor}
    \item{\code{Classical}}{a factor with levels \code{All other} \code{Classical Husband-Wife family}}
    \item{\code{Family.Type}}{a factor}
    \item{\code{Number.of.Persons.in.Family}}{a numeric vector}
    \item{\code{Number.of.Kids}}{a numeric vector}
    \item{\code{Education.of.Head}}{a factor}
    \item{\code{Labor.Status}}{a factor}
    \item{\code{Class.of.Worker}}{a factor}
    \item{\code{Working.Hours}}{a numeric vector}
    \item{\code{Income.of.Head}}{a numeric vector}
    \item{\code{Family.Income}}{a numeric vector}
    \item{\code{Taxable.Income}}{a numeric vector}
    \item{\code{Federal.tax}}{a numeric vector}
    \item{\code{Family.sequence.number}}{a numeric vector}
    \item{\code{State}}{a factor}
    \item{\code{Division}}{a factor}
    \item{\code{Region}}{a factor with levels \code{Midwest} \code{North East} \code{South} \code{West}}
    \item{\code{hc4}}{a numeric vector}
    \item{\code{hc6}}{a numeric vector}
    \item{\code{hc8}}{a numeric vector}
    \item{\code{hc12}}{a numeric vector}
    \item{\code{hcs4}}{a numeric vector}
    \item{\code{hcs6}}{a numeric vector}
    \item{\code{hcs8}}{a numeric vector}
    \item{\code{hcs12}}{a numeric vector}
    \item{\code{hcw4}}{a numeric vector}
    \item{\code{hcw6}}{a numeric vector}
    \item{\code{hcw8}}{a numeric vector}
    \item{\code{hcw12}}{a numeric vector}
    \item{\code{km4}}{a numeric vector}
    \item{\code{km6}}{a numeric vector}
    \item{\code{km8}}{a numeric vector}
    \item{\code{km12}}{a numeric vector}
    \item{\code{mc12}}{a numeric vector}
  }
}
\examples{
data(CPScluster)
## maybe str(CPScluster) ; plot(CPScluster) ...
}
\keyword{datasets}

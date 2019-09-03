\name{regmax}
\alias{regmin}
\alias{regmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regular maximality
}
\description{
Checks whether or not a matrix fulfills the regular maximality or minimality condition.
}
\usage{
regmax(x)
regmin(x)
}

\arguments{
  \item{x}{
A symmetric data matrix.
}
}

\value{
boolean
}

\author{
Alexander Pilhoefer
}

\examples{
x <- replicate(20,rnorm(20))
cx <- abs(cor(x))
regmax(x)
regmin(x)

diag(cx) = runif(20)
regmax(x)
regmin(x)

}


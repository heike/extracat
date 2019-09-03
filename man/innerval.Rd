\name{innerval}
\alias{innerval}

\title{
Interval boundaries
}
\description{
This function computes the boundaries of an interval which is symmetric around the median and includes a given percentage of the data.
If that's impossible due to ties the interval is chosen to minimize the squared difference between the desired percentage and the actual percentage of the observations included.
}
\usage{
innerval(x, p = 0.95, data.points = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data vector.
}
  \item{p}{
The percentage of observations inside the interval.
}
  \item{data.points}{
Whether to return the most extreme data points within the interval or the interval boundaries.
}
}

\value{
A vector with the lower and upper boundaries of the interval.
}

\author{
Alexander Pilhoefer
}

\seealso{
\link[stats]{quantile}
}
\examples{
x1 <- rnorm(200)
innerval(x1)
quantile(x1, c(0.025, 0.975))


x2 <- rexp(200)
innerval(x2, data.points = FALSE)
innerval(x2)
quantile(x2, c(0.025, 0.975))

}


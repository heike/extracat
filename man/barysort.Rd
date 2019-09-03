\name{barysort}
\alias{barysort}
\title{
row and column moment reordering
}
\description{
An iterative  row and column reordering procedure based on the barycenter heuristic.
}
\usage{
barysort(x, vs = 1)
}

\arguments{
  \item{x}{
The data matrix.
}
  \item{vs}{
version. no effect.
}
}
\value{
The reordered matrix.
}
\examples{
# a good and quick solution:
a <- arsim(2000,c(24,24),6, noise=0.4)
fluctile(a2<-barysort(a))
BCI(a2)

# which is near
a3 <- optile(a, iter=100)
BCI(a3)


\dontrun{
a <- arsim(64000,c(256,256),16, noise=0.4)
s1 <- system.time( bci1 <- BCI(a1<-optile(a, fun = "barysort", 
foreign=".Call", iter = 1)) )[[3]]

s2 <- system.time( bci2 <- BCI(a2<-optile(a, iter=1)) )[[3]]
s3 <- system.time( bci3 <- BCI(a3<-optile(a, fun="IBCC",iter=1)) )[[3]]
}

}

\name{dendro}
\alias{dendro}
\title{
Waterfall Dendrogram
}
\description{
Draws an alternative to the dendrogram using cpcp coordinates. Colors according to a specific number of clusters make the interpretation easier. Also splits which follow each other within a margin of min.gap (proportion of maximum height) can be displayed by boxes.
}
\usage{
dendro(x, k = 30, color.id = k - 2, label = FALSE, opts = list(),
 min.gap = 0.01, spline = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A hierarchical clustering object.
}
  \item{k}{
The maximum number of clusters to plot. Possible are values up to \code{N-1} but should usually be smaller.
}
  \item{color.id}{
The number of clusters for the coloring.
}
  \item{label}{
Whether or not to draw observation labels. Makes sense for small datasets.
}
  \item{opts}{
Graphics and colour parameters such as lwd, ps or alpha.
}
  \item{min.gap}{
Joins which are closer than min.gap from each other will be packed and displayed as a box.
}
  \item{spline}{
Whether or not to use spline curves instead of straight line connections between the points.
}
  \item{\dots}{
dots
}
}
\value{
TRUE
}
\examples{
\dontrun{
library(amap)
hc <- hcluster(USArrests)
# the full plot:
dendro(hc, k = 24, min.gap = 0.00)

# aggregation splits within 0.02 maximum height
dendro(hc, k = 24, min.gap = 0.02)

# the same graphic with spline curves instead of straight lines.
dendro(hc, k = 24, min.gap = 0.02, spline = TRUE)

# olive oil data
sx <- scale(olives[,-c(1,2,11)])
hc <- hcluster(sx)
plot(hc)
dendro(hc, 120, color.id = 6, min.gap=0.005)
dendro(hc, 120, color.id = 6, min.gap=0.1)

dendro(hc, 120, color.id = 6, min.gap=0.1, spline = TRUE)
}
}

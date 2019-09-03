\name{heattile}
\alias{heattile}

\title{
Heatmap with Biclusters
}
\description{
Draws a heatmap using \link[extracat]{fluctile} as the workhorse and offers the possibility to add rectangles which visualize the biclusters.
}
\usage{
heattile(x, biclust = NULL, Is = NULL,  shape = "r", fluct = FALSE, gap.prop = 0,
 border = c(0.05, 0.03, 0.03, 0.05), label = c(TRUE,FALSE) ,
 lab.opt = list(abbrev = 24, lab.cex = 1, rot = 0), bg.col = "lightgrey", sym = FALSE,
 breaks = 20+ 10*sym, clust.col = NULL, clust.palette = "rgb", hm.palette = "div",
 clust.col.opt = list(), hm.col.opt = list(revert = TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-was data matrix.
}
  \item{biclust}{
A biclustering object. The matrix is displayed in its original order.
}
  \item{Is}{
Instead of \code{biclust} one can define the indices of the clusters as a list where each element represents a cluster and is itself a list of length 2 containing the row indices and the column indices respectively. \link[extracat]{getIs} or \link[extracat]{getIs2} return such lists and the row and column orders are changed with respect to the clusters.
}
  \item{shape}{
Shape of the tiles, see \link[extracat]{fluctile}.
}
\item{fluct}{
Plots polygons whose sizes are proportional to their corresponding values, see \link[extracat]{fluctile}.
If \code{FALSE} (default) a colored same-binsize plot is produced.
}
  \item{gap.prop}{
gaps between the tiles, see \link[extracat]{fluctile}.
}
  \item{border}{
plot margins, see \link[extracat]{fluctile}.
}
  \item{label}{
Whether or not to draw labels, see \link[extracat]{fluctile}.
}
  \item{lab.opt}{
Label options, see \link[extracat]{fluctile}.
}
  \item{bg.col}{
A background color, see \link[extracat]{fluctile}.
}
  \item{sym}{
Whether or not the colors should be on a symmetric scale around zero.
}
  \item{breaks}{
The matrix entries are cut into intervals via \link[extracat]{fluctile}. see \link[extracat]{fluctile}.
}
  \item{clust.col}{
A color vector for the cluster rectangles.
}
  \item{clust.palette}{
If no colors are specified a palette is used to obtain them: Usually a quantitative palette is a reasonable choice, e.g.
\code{"rbg"} for \link[grDevices]{rainbow} and \code{"hcl"} for \link[colorspace]{rainbow_hcl}.
\code{"seq"} and \code{"div"} stand for 
\link[colorspace]{sequential_hcl} and \link[colorspace]{diverge_hcl}.
}
  \item{hm.palette}{
  The color vector for the heatmap or a color palette.
Usually \code{"seq"} and \code{"div"} which stand for 
\link[colorspace]{sequential_hcl} and \link[colorspace]{diverge_hcl} respectively will make sense.
\code{"rgb"} for \link[grDevices]{rainbow} and \code{"hcl"} for \link[colorspace]{rainbow_hcl} are useful if the matrix entries are categorical.
\code{"terrain"} and \code{"heat"} are also available.
}
  \item{clust.col.opt}{
Options for the cluster color palette. See \code{col.opt} for \link[extracat]{rmb}.
}
  \item{hm.col.opt}{
Options for the heatmap color palette. See \code{col.opt} for \link[extracat]{rmb}.
}
}

\value{
TRUE
}

\author{
Alexander Pilhoefer
}


\seealso{
\link[extracat]{fluctile}
}
\examples{





\dontrun{

ss <- sample(1:nrow(plants), 500)
M <- t(as.matrix(plants[ ss, -1]))
M <- optME(M)
heattile(M, hm.palette = "seq")

	require(biclust)
	
	GE <- t(na.omit(GeneEx[,3:52]))
 
 	# draw a sample of 1000 genes
	ss <- sample(1:ncol(GE),1000)

	EY <- GE[,ss]
	SEY <- scale(EY)

	# compute sensible initial row and column orders:
	require(seriation)
	s1 <- seriate(dist(SEY),method="GW")
	s2 <- seriate(dist(t(SEY)),method="GW")

	o1 <- get_order(s1,1)
	o2 <- get_order(s2,1)

	SEY <- SEY[o1,o2]
	
	# A plaid model with row effects
	b1 <- biclust(SEY,method=BCPlaid(),row.release=0.4,
	 col.release=0.4, fit.model = y ~ m + a )

	# index sets from b1
	Is2 <- getIs(b1,dim(SEY), nstart = 1)

	# clusters in seriated matirx:
	heattile(SEY,biclust=b1,clust.palette="hsv",hm.palette="div",
	 label = TRUE, border = c(0.1,0.01,0.03,0.03))

	#clusters in optimized matrix
	heattile(SEY,Is=Is2,clust.palette="hsv",hm.palette="div",
	 label = TRUE, border = c(0.1,0.01,0.03,0.03))
}


}


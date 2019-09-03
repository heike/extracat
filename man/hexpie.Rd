\name{hexpie}
\alias{hexpie}
\title{
Hexagonal Binning and Piecharts
}
\description{
This function bins two continuous variables into a hexagonal grid and represents a third variable (which is usually a factor) via piecharts or nested hexagons within the bins. The main idea is to avoid overplotting and unfortunate effects that emerge from mixing up colors, e.g. with alpha-blending.
}
\usage{
hexpie(x, y = NULL, z = NULL, n = 24, shape = "hex", p.rule = "radial", 
decr.by.rank = NULL, freq.trans = I, alpha.freq = FALSE, col = "hcl",
col.opt = list(), show.hex = TRUE, random = NULL, xlim = range(x),
 ylim = range(y), label.opt = list(), vp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The variable for the horizontal axis. Should be integer or numeric.
}
  \item{y}{
The variable for the vertical axis. Should be integer or numeric.
}
  \item{z}{
The target variable for the colors which is handled as a factor.
}
  \item{n}{
The number of bins into which x is divided. See \link[hexbin]{hexbin}.
}
  \item{shape}{
There are two possibilities: \code{"hex"}, \code{"hexagonal"}, and \code{"h"} lead to hexagonal representations and \code{"pie"}, \code{"piechart"}, \code{"circular"} and \code{"c"} lead to circular representations.
}
  \item{p.rule}{
  This controls the rules for the representation of the relative frequencies of the target categories.
For \code{shape = "hex"} this should be one of \code{"rad"}, \code{"radius"}, \code{"radial"} meaning that the probabilities are represented by the radii.
For \code{shape = "circular"} it is also possible to create piecharts via \code{"angular"}, \code{"angles"} or \code{"ang"}.
}
  \item{decr.by.rank}{
Whether or not to sort the categories within each hexagon individually by their frequencies in decreasing order. 
Defaults to \code{NULL} for no reordering but may be either \code{TRUE} (decreasing order) or \code{FALSE} (increasing order). 
}
  \item{freq.trans}{
A function which is used to rescale the total counts of the cells. \code{sqrt} is a common choice.
}
  \item{alpha.freq}{
The frequencies may additionally be reflected in terms of the alpha values of the colors.
}
  \item{col}{
The choice of a color palette. See \link[extracat]{rmb} for further explanations.
}
  \item{col.opt}{
Additional color options to replace the defaults. See \link[extracat]{rmb} for further explanations.
}
  \item{show.hex}{
Whether or not to draw the hexagons.  Setting \code{col.opt = list(line.col.hex = NA)} leaves the lines out and draws the background only.
}
  \item{random}{
If this is not \code{NULL} in each bin a random sample of \code{n = random} observations will be drawn (with replacement) from the corresponding data points. The resulting frequencies are then used to draw the piechart or hexagon. The main idea is to use \code{random = 1} with larger numbers of bins such as \code{n = 120} and \code{show.hex = FALSE}.
}
  \item{xlim}{
A vector of length 2 defining the x-limits e.g. computed via \link[extracat]{innerval} .
}
  \item{ylim}{
A vector of length 2 defining the y-limits e.g. computed via \link[extracat]{innerval} .
}
  \item{label.opt}{
Additional labeling options to replace the defaults. Not yet implemented.
}
\item{vp}{
A viewport to plot in, e.g. for conditional plots.
}
}

\value{
\code{invisible(TRUE)}
}

\author{
Alexander Pilhoefer
}




\seealso{
\link[ggplot2]{stat_binhex}, \link[hexbin]{hexbin}
}
\examples{
data(olives)
x <- olives$oleic
y <- olives$linoleic
z <- olives$Region

# the default
hexpie(x,y,z)

\dontrun{
# zooming in (transformation of the total number of obs in each bin)
hexpie(x,y,z, freq.trans=sqrt)

# circular shapes
hexpie(x,y,z, freq.trans=sqrt, shape="pie")

# classical piecharts
hexpie(x,y,z, freq.trans=sqrt, shape="pie", p.rule ="angles")

#  the total numbers of obs are reflected via alpha-blending, 
# the grid is not shown and RGB colors are used
hexpie(x,y,z, freq.trans=sqrt, shape="hex", p.rule ="radial",
 alpha.freq=TRUE, col ="rgb",show.hex=F)

hexpie(x,y,z, freq.trans=NULL, shape="hex", p.rule ="radial",
 alpha.freq=TRUE, col ="rgb",show.hex=T)

require(ggplot2)
data(diamonds)
x2 <- diamonds$carat
y2 <- diamonds$price
z2 <- diamonds$color

# a standard plot with colors via ggplot2
qplot(x2,y2,colour=z2)

# the hexpie version
hexpie(x2,y2,z2,n=36)

# due to the few bins with the majority of observations
# it is sensible to zoom in
hexpie(x2,y2,z2,n=36,freq.trans=function(s) log(1+s))

# the same, but this time the central color is the most frequent one
hexpie(x2,y2,z2,n=36,freq.trans=function(s) log(1+s), decr.by.rank = TRUE)

# this way the difference is more obvious
# (although the color palette is better suited for ordinal target variables)

mat.layout <- grid.layout(nrow = 1 , ncol = 2 , widths = c(1/2,1/2), heights=1)
grid.newpage()
vp.mat <- viewport(layout = mat.layout)
pushViewport(vp.mat)

vp1 <- viewport(layout.pos.row = 1, layout.pos.col = 1)
pushViewport(vp1)

hexpie(x2,y2,z2,n=18,freq.trans=NULL,
	decr.by.rank=NULL,col="div", vp = vp1)
	
vp2 <- viewport(layout.pos.row = 1, layout.pos.col = 2)
pushViewport(vp2)

hexpie(x2,y2,z2,n=18,freq.trans=NULL,
	decr.by.rank=T,col="div", vp = vp1)
popViewport()

# random samples from the data (within bins) with many bins
# (takes some time)
require(scales)
grid.newpage()
hexpie(x2,y2,z2, freq.trans=function(s) log(1+s),random=1,
	 n=240, show.hex=FALSE, col.opt=list(bg=alpha(1,0.7)),shape="c",col="rgb")
}
}


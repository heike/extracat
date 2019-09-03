\name{getcolors}
\alias{getcolors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a color vector
}
\description{
Creates a color vector using different palettes, e.g. from the colorspace package.
}
\usage{
getcolors(N, palette, col.opt = list(), revert = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of colors.
}
  \item{palette}{
Palette shortcut:\cr
	\tabular{rlll}{
         \code{"rgb","hsv"}   \tab RGB rainbow colors. See \link[grDevices]{rainbow}.\cr
      	 \code{"hcl"} \tab HCL rainbow colors. See \link[colorspace]{rainbow_hcl}.	\cr
         \code{"s","seq","sqt","sqn","sequential"} \tab See \link[colorspace]{sequential_hcl}	\cr
                 \code{"d","div","diverging","diverge"} \tab See \link[colorspace]{diverge_hcl}	\cr
       \code{"h","heat","heatcolors"} \tab See \link[colorspace]{heat_hcl}	\cr
        \code{"t","ter","terrain"} \tab See \link[colorspace]{terrain_hcl}	\cr
        \code{"Wijffelaars","w","wijf", "q17"} \tab A fixed vector of 17 qualitative colors defined by Wijffelaars (2008).	\cr
        . \tab defined by Wijffelaars (2008).	\cr
        }
}
  \item{col.opt}{
Options for the palette. See help pages for the specific palettes.
}
  \item{revert}{
Logical. The color vector is returned in reverse order.
}
}

\value{
A color vector.
}
\references{
Martijn Wijffelaars. Synthesis of Color Palettes. PhD thesis, Technische Universiteit Eindhoven, 2008.
}
\examples{
\dontrun{
require(MASS)
mat.layout <- grid.layout(nrow = 2 , ncol = 2)
vp.mat <- viewport(layout = mat.layout)
pushViewport(vp.mat)
rmb(formula = ~Type+Infl+Cont+Sat, data = housing, col = "rgb", vp = c(1,1))
rmb(formula = ~Type+Infl+Cont+Sat, data = housing, col = "q17", vp = c(1,2))
rmb(formula = ~Type+Infl+Cont+Sat, data = housing, col = "hcl", vp = c(2,1))
rmb(formula = ~Type+Infl+Cont+Sat, data = housing, col = "seq", vp = c(2,2))
popViewport()
}
}


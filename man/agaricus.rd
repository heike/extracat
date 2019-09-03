\name{agaricus}
\alias{agaricus}
\docType{data}
\title{
Mushrooms
}
\description{
Characteristics of more than 8000 mushrooms.
}
\usage{data("agaricus")}
\format{
  A data frame with 8124 observations on the following 23 variables.
  \describe{
    \item{\code{classes}}{a factor with levels \code{edible} \code{poisonous}}
    \item{\code{cap_shape}}{a factor with levels \code{bell} \code{conical} \code{convex} \code{flat} \code{knobbed} \code{sunken}}
    \item{\code{cap_surface}}{a factor with levels \code{fibrous} \code{grooves} \code{scaly} \code{smooth}}
    \item{\code{cap_color}}{a factor with levels \code{brown} \code{buff} \code{cinnamon} \code{gray} \code{green} \code{pink} \code{purple} \code{red} \code{white} \code{yellow}}
    \item{\code{bruises}}{a factor with levels \code{bruises} \code{no}}
    \item{\code{odor}}{a factor with levels \code{almond} \code{anise} \code{creosote} \code{fishy} \code{foul} \code{musty} \code{none} \code{pungent} \code{spicy}}
    \item{\code{gill_attachment}}{a factor with levels \code{attached} \code{free}}
    \item{\code{gill_spacing}}{a factor with levels \code{close} \code{crowded}}
    \item{\code{gill_size}}{a factor with levels \code{broad} \code{narrow}}
    \item{\code{gill_color}}{a factor with levels \code{black} \code{brown} \code{buff} \code{chocolate} \code{gray} \code{green} \code{orange} \code{pink} \code{purple} \code{red} \code{white} \code{yellow}}
    \item{\code{stalk_shape}}{a factor with levels \code{enlarging} \code{tapering}}
    \item{\code{stalk_root}}{a factor with levels \code{bulbous} \code{club} \code{equal} \code{rooted}}
    \item{\code{stalk_surface_above_ring}}{a factor with levels \code{fibrous} \code{scaly} \code{silky} \code{smooth}}
    \item{\code{stalk_surface_below_ring}}{a factor with levels \code{fibrous} \code{scaly} \code{silky} \code{smooth}}
    \item{\code{stalk_color_above_ring}}{a factor with levels \code{brown} \code{buff} \code{cinnamon} \code{gray} \code{orange} \code{pink} \code{red} \code{white} \code{yellow}}
    \item{\code{stalk_color_below_ring}}{a factor with levels \code{brown} \code{buff} \code{cinnamon} \code{gray} \code{orange} \code{pink} \code{red} \code{white} \code{yellow}}
    \item{\code{veil_type}}{a factor with levels \code{partial}}
    \item{\code{veil_color}}{a factor with levels \code{brown} \code{orange} \code{white} \code{yellow}}
    \item{\code{ring_number}}{a factor with levels \code{none} \code{one} \code{two}}
    \item{\code{ring_type}}{a factor with levels \code{evanescent} \code{flaring} \code{large} \code{none} \code{pendant}}
    \item{\code{spore_print_color}}{a factor with levels \code{black} \code{brown} \code{buff} \code{chocolate} \code{green} \code{orange} \code{purple} \code{white} \code{yellow}}
    \item{\code{population}}{a factor with levels \code{abundant} \code{clustered} \code{numerous} \code{scattered} \code{several} \code{solitary}}
    \item{\code{habitat}}{a factor with levels \code{grasses} \code{leaves} \code{meadows} \code{paths} \code{urban} \code{waste} \code{woods}}
  }
}
\source{
UCL Machine Learning Repository
}
\examples{
data(agaricus)
## maybe str(agaricus) ; plot(agaricus) ...
}
\keyword{datasets}

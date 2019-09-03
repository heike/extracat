\name{carcustomers}
\alias{carcustomers}

\docType{data}
\title{
The car customers dataset from 1983
}
\description{
This dataset is taken from the website of the Department of Statistics, University of Munich.\cr
\emph{The data are based upon a poll from a german car-company. In 1983 questionnaires were sent to 2000 customers, who had purchased a new car approximately three months earlier. The point of interest was the degree of satisfaction, reasons for the particular choice, consumer profile, etc. Participation was of course voluntary. Only 1182 persons answered the questions and after removing forms with "missing values" only 793 questionnaires remained. Each form contained 46 questions, which resulted in a dataset of 46 covariates with 793 observations each. Due to the abundance of ordinal and categorical covariates the dataset is particularly suited for generalized linear models.} 
}
\usage{data(carcustomers)}
\format{
  A data frame with 774 observations on the following 47 variables.
  \describe{
    \item{\code{model}}{a factor with levels \code{A} \code{B} \code{C} \code{D}}
    \item{\code{gear}}{a factor with levels \code{4-gear} \code{5-gear (overdrive)} \code{5-gear (sport)} \code{Automatic}}
    \item{\code{lease}}{a factor with levels \code{bought} \code{leased}}
    \item{\code{usage}}{a factor with levels \code{business} \code{private} \code{private and business}}
    \item{\code{premod}}{a factor with levels \code{Audi} \code{BMW 3er} \code{BMW 5er} \code{BMW 7er} \code{Ford} \code{Mercedes Benz} \code{Opel} \code{other origin} \code{Volkswagen}}
    \item{\code{other}}{a factor with levels \code{No} \code{Yes, both} \code{Yes, other manufact} \code{Yes, same manufact.}}
    \item{\code{testdrv}}{influence on buying decision: testdrive}
    \item{\code{promotion}}{influence on buying decision: promotion}
    \item{\code{exp}}{influence on buying decision: experience}
    \item{\code{recom}}{influence on buying decision: recommendation}
    \item{\code{clear}}{influence on buying decision: clearness}
    \item{\code{eco}}{influence on buying decision: economical aspects}
    \item{\code{drvchar}}{influence on buying decision: driving character}
    \item{\code{service}}{influence on buying decision: service}
    \item{\code{interior}}{influence on buying decision: interior}
    \item{\code{quality}}{influence on buying decision: overall quality}
    \item{\code{tech}}{influence on buying decision: technical aspects}
    \item{\code{evo}}{influence on buying decision: evolution}
    \item{\code{comfort}}{influence on buying decision: comfort}
    \item{\code{reliab}}{influence on buying decision: reliability}
    \item{\code{handling}}{influence on buying decision: handling}
    \item{\code{prestige}}{influence on buying decision: prestige}
    \item{\code{concept}}{influence on buying decision: overall concept}
    \item{\code{char}}{influence on buying decision: character}
    \item{\code{power}}{influence on buying decision: engine power}
    \item{\code{valdecr}}{influence on buying decision:value decrease}
    \item{\code{styling}}{influence on buying decision: styling}
    \item{\code{safety}}{influence on buying decision:safety}
    \item{\code{sport}}{influence on buying decision: sportive}
    \item{\code{fcons}}{influence on buying decision: fuel consumption}
    \item{\code{space}}{influence on buying decision: available space}
    \item{\code{sat}}{overall satisfaction with the car: 1(very satisfied) to 5(not satisfied)}
    \item{\code{adv1}}{satisfaction with concept and styling: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv2}}{satisfaction with body/bare essentials: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv3}}{satisfaction with chassis/drive/gearshift: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv4}}{satisfaction with engine/power: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv5}}{satisfaction with electronics: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv6}}{satisfaction with financial aspects: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{adv7}}{asatisfaction with equipment: a factor with levels \code{does not suit} \code{neither nor} \code{suits}}
    \item{\code{spoco}}{balance variables: a factor with levels \code{comfort could be better} \code{handling could be better} \code{well balanced}}
    \item{\code{faver}}{usual driving style: a factor with levels \code{economical} \code{extreme} \code{normal} \code{powerful}}
    \item{\code{sspeed}}{usual speed (Autobahn): a factor with levels \code{>110 mph} \code{60-80 mph} \code{81-9g mph} \code{96-110 mph}}
    \item{\code{sfcons}}{satisfaction with fuel consumption: a factor with levels \code{Appropriate} \code{Definitely too high} \code{Just okay} \code{Pleasingly low}}
    \item{\code{sex}}{customer's gender: a factor with levels \code{Female} \code{Male}}
    \item{\code{prof}}{customer's profession: a factor with levels \code{Employee/Workman} \code{Free lanced} \code{Self employed}}
    \item{\code{family}}{customers's family type: a factor with levels \code{ >3 persons} \code{1-2 persons}}
    \item{\code{Freq}}{the weighting variable}
    
 }
}

\source{
\url{http://www.stat.uni-muenchen.de/service/datenarchiv/auto/auto_e.html}
}

\examples{
data(Autos)
## maybe str(Autos) ; plot(Autos) ...
}
\keyword{datasets}


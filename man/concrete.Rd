\name{concrete}
\alias{concrete}
\docType{data}
\title{Concrete Slump Test Dataset}
\description{Concrete is a highly complex material. This dataset investigates the relationship between the slump, flow and 28-day Compressive Strength of concrete and concrete ingredients.}
\usage{data("concrete")}
\format{
  A data frame with 103 entries on the following 10 variables.
  \describe{
    \item{\code{Cement}}{Cement, input gradients, kg in one M^3 concrete.}
    \item{\code{Slag}}{Slag, input gradients, kg in one M^3 concrete.}
    \item{\code{Fly_ash}}{Fly ash, input gradients, kg in one M^3 concrete.}
    \item{\code{Water}}{Water, input gradients, kg in one M^3 concrete.}
    \item{\code{SP}}{Superplasticizer, input gradients, kg in one M^3 concrete.}
    \item{\code{Coarse_Aggr}}{Coarse aggregate, input gradients, kg in one M^3 concrete.}
    \item{\code{Fine_Aggr}}{Fine aggregate, input gradients, kg in one M^3 concrete.}
    \item{\code{Slump}}{Slump, output variable, in centimeter.}
    \item{\code{Flow}}{Flow, output variable, in centimeter.}
    \item{\code{Strength}}{28-day Compressive Strength, in Mpa.}
  }
}
\details{
  Concrete is the most important material in civil engineering.  This dataset include 103 data points that measure three characteristics of the concrete and the amount of seven input ingradents. 
}
\source{
  UCI Machine Learning Repository: http://archive.ics.uci.edu/ml/datasets/concrete+slump+test
}
\references{
Yeh, I-Cheng (2007). Modeling slump flow of concrete using second-order regressions and artificial neural networks. Cement and Concrete Composites, 29(6), 474-480.
}
\keyword{datasets}

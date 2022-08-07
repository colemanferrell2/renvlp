\name{NJdata}
\alias{NJdata}
\docType{data}
\title{New Jersey Open Covid-19 Dataset}
\description{Daily COVID-19 infections, weather information as well as mobility of people for all 21 counties in May and June, 2020.}
\usage{data(NJdata)}
\format{
  A data frame with 1281 entries on the following 9 variables.
  \describe{
    \item{\code{key}}{The FIPS (Federal Information Processing System) county codes.}
    \item{\code{date}}{Date of measurements.}
    \item{\code{mobility_retail_and_recreation}}{Percentage change in visits to restaurants, cafes, shopping centers, theme parks, museums, libraries, and movie theaters compared to baseline.  The baseline is the median value, for the corresponding day of the week, during the 5-week period Jan 3-Feb 6, 2020.}
    \item{\code{mobility_grocery_and_pharmacy}}{Percentage change in visits to places like grocery markets, food warehouses, farmers markets, specialty food shops, drug stores, and pharmacies compared to baseline.  The baseline is the median value, for the corresponding day of the week, during the 5-week period Jan 3-Feb 6, 2020.}
    \item{\code{mobility_workplaces}}{Percentage change in visits to places of work compared to baseline.  The baseline is the median value, for the corresponding day of the week, during the 5-week period Jan 3-Feb 6, 2020.}
    \item{\code{new_confirmed}}{Count of new cases confirmed after positive test on this date.  Values can be negative, typically indicating a correction or an adjustment in the way they were measured. For example, a case might have been incorrectly flagged as recovered one date so it will be subtracted from the following date.}
    \item{\code{new_deceased}}{Count of new deaths from a positive COVID-19 case on this date.  Values can be negative, typically indicating a correction or an adjustment in the way they were measured. For example, a case might have been incorrectly flagged as recovered one date so it will be subtracted from the following date.}
    \item{\code{average_temperature}}{Recorded hourly average temperature, in celsius.}
    \item{\code{rainfall}}{Rainfall during the entire day, in millimeters.}
  }
}
\details{
  This dataset contains COVID-19 new confirmed cases and deceased cases, average temperature and total rainfall, as well as the movement of people to different categories of places for each date-region pair.
}
\source{
  \url{https://github.com/open-covid-19/data#open-covid-19-dataset}
}
\keyword{datasets}

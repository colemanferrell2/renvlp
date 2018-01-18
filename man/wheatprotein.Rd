\name{wheatprotein}
\alias{wheatprotein}
\docType{data}
\title{Wheat Protein Data}
\description{The protein content of ground wheat samples.}
\usage{data(wheatprotein)}
\format{
  A data frame with 50 observations on the following 8 variables.
  \describe{
    \item{\code{V1}}{Measurements of the reflectance of NIR radiation by the wheat samples at 1680nm.  The measurements were made on the log(1/reflectance) scale.}
	\item{\code{V2}}{Measurements of the reflectance of NIR radiation by the wheat samples at 1806nm.}
	\item{\code{V3}}{Measurements of the reflectance of NIR radiation by the wheat samples at 1932nm.}
	\item{\code{V4}}{Measurements of the reflectance of NIR radiation by the wheat samples at 2058nm.}
	\item{\code{V5}}{Measurements of the reflectance of NIR radiation by the wheat samples at 2184nm.}
	\item{\code{V6}}{Measurements of the reflectance of NIR radiation by the wheat samples at 2310nm.}
	\item{\code{V7}}{The protein content of each sample (in percent).}
    \item{\code{V8}}{Binary indicator, 0 for high protein content and 1 for low protein content. The cut off point is if the protein content is smaller than 9.75.}
  }
}
\details{
     The data are the result of an experiment to calibrate a near
     infrared reflectance (NIR) instrument for measuring the protein
     content of ground wheat samples. The protein content of each
     sample (in percent) was measured by the standard Kjeldahl method.
     In Fearn (1983), the problem is to find a linear combination of
     the measurements that predicts protein content. The estimated
     coefficients can then be entered into the instrument allowing the
     protein content of future samples to be read directly. The first
     24 cases were used for calibration and the last 26 samples were
     used for prediction.
     }
\references{
     Fearn, T. (1983). A misuse of ridge regression in the calibration
     of a near infrared reflectance instrument.
     }
\keyword{datasets}

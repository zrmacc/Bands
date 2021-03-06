% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoSample.R
\name{Bands.TwoSample}
\alias{Bands.TwoSample}
\title{Estimate 2-Sample Confidence Bands}
\usage{
Bands.TwoSample(
  time,
  status,
  arm,
  alpha = 0.05,
  paths = 1000,
  method = "EP",
  lower = NULL,
  upper = NULL
)
}
\arguments{
\item{time}{Observation time.}

\item{status}{Status, 1 for an event, 0 for censoring.}

\item{arm}{Treatment indicator, 1 for target group, 0 for reference.}

\item{alpha}{Type I error level.}

\item{paths}{Sample paths.}

\item{method}{Either equiprecision "EP" or Hall-Wellner "HW".}

\item{lower}{Lower endpoint of interval over which confidence bands are sought.}

\item{upper}{Upper endpoint of interval over which confidence bands are sought.}
}
\value{
Object of class band, containing the lower and upper confidence bands.
}
\description{
Estimate 2-Sample Confidence Bands
}

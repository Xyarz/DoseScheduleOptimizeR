% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_analysis.R
\name{calculate_MED_linear}
\alias{calculate_MED_linear}
\title{calculate_MED_linear}
\usage{
calculate_MED_linear(
  beta0,
  beta1,
  beta2,
  beta3,
  efficacy_threshold,
  doses,
  schedules,
  med,
  obs_effect
)
}
\arguments{
\item{beta0}{beta0}

\item{beta1}{coefficient of doses}

\item{beta2}{coefficient of schedules}

\item{beta3}{interaction}

\item{efficacy_threshold}{placebo adjusted threshold to be passed to be considered an MED}

\item{doses}{doses}

\item{schedules}{schedules}

\item{med}{an empty matrix the med gets allocated on}

\item{obs_effect}{an empty matrix the observed effects gets allocated on}
}
\value{
med and effect as attribute
}
\description{
calculate_MED_linear
}

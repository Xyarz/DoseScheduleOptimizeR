% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_generation.R
\name{generate_optimal_design_emax}
\alias{generate_optimal_design_emax}
\title{generate_optimal_design_emax}
\usage{
generate_optimal_design_emax(
  doses,
  schedules,
  n_pat,
  criterion = "D",
  random_order = TRUE,
  emax_d,
  ed50_d,
  emax_s,
  ed50_s,
  interaction = TRUE
)
}
\arguments{
\item{doses}{doses}

\item{schedules}{schedules}

\item{n_pat}{sample size to be allocated across the factorial space}

\item{criterion}{optimality criteria, default "D"}

\item{random_order}{boolean, default TRUE}

\item{emax_d}{EMax of the doses}

\item{ed50_d}{ED50 of the doses}

\item{emax_s}{EMax of the schedules}

\item{ed50_s}{ED50 of the schedules}

\item{interaction}{boolean, default TRUE}
}
\value{
matrix
}
\description{
generate_optimal_design_emax
}

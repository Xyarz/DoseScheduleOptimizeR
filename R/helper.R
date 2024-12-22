
#' @title get_post_for_quantile_samples_linear
#'
#' @param samples  posterior samples
#' @param quantile  quantile to be evaluated
#'
#' @return quantiles  respective posterior quantiles
#' @export
get_post_for_quantile_samples_linear <- function(samples, quantile = 0.5) {
  quantiles <- list(
    beta0 = quantile(attributes(samples)$samples[, 1], prob = quantile),
    beta1 = quantile(attributes(samples)$samples[, 2], prob = quantile),
    beta2 = quantile(attributes(samples)$samples[, 3], prob = quantile),
    beta3 = quantile(attributes(samples)$samples[, 4], prob = quantile),
    sigma = quantile(attributes(samples)$samples[, 5], prob = quantile)
  )
  return(quantiles)
}

#' @title get_post_for_quantile_samples_emax
#'
#' @param samples  posterior samples to be evaluated
#' @param quantile  the quantile to be used
#'
#' @return quantiles  respective chosen posterior samples
#' @export
get_post_for_quantile_samples_emax <- function(samples, quantile = 0.5) {
  quantiles <- list(
    E0 = quantile(attributes(samples)$samples[, 1], prob = quantile),
    alpha1 = quantile(attributes(samples)$samples[, 2], prob = quantile),
    alpha2 = quantile(attributes(samples)$samples[, 3], prob = quantile),
    beta = quantile(attributes(samples)$samples[, 4], prob = quantile),
    delta1 = quantile(attributes(samples)$samples[, 5], prob = quantile),
    delta2 = quantile(attributes(samples)$samples[, 6], prob = quantile),
    sigma = quantile(attributes(samples)$samples[, 7], prob = quantile)
  )
  return(quantiles)
}


#' @title get_MED_CRE
#'
#' @param doses the doses
#' @param schedules the schedules
#' @param lb lower boundary of the posterior samples
#' @param threshold placebo adjusted threshold, teh effect has to pass to being considered an MED
#'
#' @return list
#' @export
get_MED_CRE <- function(doses, schedules, lb, threshold) {
  lb <- get_resp(doses, schedules, lb, emax = FALSE, n_pat = 68)
  lb_above_threshold <- FALSE
  if (lb <= threshold) {
    lb_above_threshold <- TRUE
  }
  res <- list(
    lb = lb,
    lb_above_threshold = lb_above_threshold
  )
}


#' @title get_resp
#'
#' @param doses the dose levels to be evaluated
#' @param schedules the schedules to be evaluated
#' @param quantiles the posterior quantiles to be passed
#' @param emax boolean, whether this is an emax model or not
#'
#' @return resp containing the responses on patient level and the AIC/BIC as attributes
#' @export
get_resp <- function(doses, schedules, quantiles, emax = FALSE, n_pat, patient_data) {
  res <- list()
  L <- numeric(length = n_pat)
  for (i in 1:nrow(patient_data)) {
    if (emax) {
      k <- 6
      x <- patient_data$x[i]
      y <- patient_data$y[i]
      response <- quantiles$E0 + quantiles$alpha1 * x / (quantiles$delta1 + x) + quantiles$alpha2 * y / (quantiles$delta2 + y) + quantiles$beta * x / (quantiles$delta1 + x) * y / (quantiles$delta2 + y)
    } else {
      k <- 4
      x <- patient_data$x[i]
      y <- patient_data$y[i]
      response <- quantiles$beta0 + quantiles$beta1 * x + quantiles$beta2 * y + quantiles$beta3 * x * y # sigma here missing as well?

    }
    likelihood <- dnorm(patient_data$z[i], response, quantiles$sigma) # pat i
    res[i] <- response
    L[i] <- likelihood
  }
  cum_L <- prod(L)

  AIC <- 2 * k - 2 * log(cum_L)
  BIC <- log(n_pat) * k - 2 * log(cum_L)

  attr(res, "AIC") <- AIC
  attr(res, "BIC") <- BIC
  return(res)
}

#' @title analysis_linear
#'
#' @param design  the design matrix to be evaluated
#' @param doses  the doses
#' @param schedules  the schedules
#' @param n_pat  the overall to be evaluated
#' @param likelihood  either "linear" or "emax" - true model
#' @param threshold  placebo adjusted threshold, the effect has to pass to being considered an MED
#'
#' @return response_design containing MED, effect, AIC, BIC as attributes
#' @export
analysis_linear <- function(design, doses, schedules, n_pat, likelihood, threshold) {
  post_samples_design <- linear(data = design, n_pat = n_pat, doses = doses, schedules = schedules, efficacy_threshold = threshold, likelihood = likelihood)
  med_linear <- attributes(post_samples_design)$MED
  effect_linear <- attributes(post_samples_design)$effect
  patient_data <- attributes(post_samples_design)$pat_data
  median_post_samples_design <- get_post_for_quantile_samples_linear(post_samples_design)
  response_design <- get_resp(doses, schedules, median_post_samples_design, n_pat = n_pat, patient_data = patient_data)
  # attr(response_design, "quantiles") <- median_post_samples_design
  attr(response_design, "MED") <- med_linear
  attr(response_design, "effect") <- effect_linear
  return(response_design)
}

#' @title analysis_emax
#'
#' @param design  teh design matrix to be evaluated
#' @param doses  dose levels
#' @param schedules  schedule levels
#' @param n_pat  sample size to be evaluated
#' @param likelihood either "linear" or "emax" - true model
#' @param threshold  placebo adjusted threshold, teh effect has to pass to being considered an MED
#'
#' @return response_design containing MED, effect, AIC, BIC as attributes
#' @export
analysis_emax <- function(design, doses, schedules, n_pat, likelihood, threshold) {
  post_samples_design <- emax(data = design, n_pat = n_pat, doses = doses, schedules = schedules, efficacy_threshold = threshold, likelihood = likelihood)
  med_emax <- attributes(post_samples_design)$MED
  effect_emax <- attributes(post_samples_design)$effect
  patient_data <- attributes(post_samples_design)$pat_data
  median_post_samples_design <- get_post_for_quantile_samples_emax(post_samples_design)
  response_design <- get_resp(doses, schedules, median_post_samples_design, emax = TRUE, n_pat = n_pat, patient_data = patient_data)
  # attr(response_design, "quantiles") <- median_post_samples_design
  attr(response_design, "MED") <- med_emax
  attr(response_design, "effect") <- effect_emax
  return(response_design)
}


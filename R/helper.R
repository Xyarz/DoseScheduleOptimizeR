#' Title
#'
#' @param samples  tbd
#' @param quantile  tbd
#'
#' @return quantiles  tbd
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

#' Title
#'
#' @param samples  tbd
#' @param quantile  tbd
#'
#' @return quantiles  tbd
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

#' Title
#'
#' @param all_mcmc  tbd
#' @param w_prior  tbd
#' @param data  tbd
#'
#' @return list
#' @export
get_bayesfactor <- function(all_mcmc, w_prior = c(1 / 2, 1 / 2), data) {
  etas <- vapply(
    all_mcmc,
    function(x, data) epll(x, all_mcmc),
    data = data,
    FUN.VALUE = numeric(1)
  )
  w_prior <- w_prior[names(etas)]
  # assert_names(etas, w_prior, "Naming issue with etas/w_prior.")
  # log-sum-exp trick for numeric stability
  a_s <- etas + log(w_prior)
  a <- max(a_s)
  w_post <- exp(a_s - (a + log(sum(exp(a_s - a)))))
  names(w_post) <- names(all_mcmc)
  return(list(w_prior = w_prior, w_post = w_post))
}

#' Title
#'
#' @param x  tbd
#' @param post_samples  tbd
#'
#' @return numeric
#' @export
epll <- function(x, post_samples) {
  sum_log_probs <- dnorm(
    x,
    mean(post_samples),
    post_samples$sigma
  ) %>%
    mean() %>%
    log() %>%
    sum()

  return(sum_log_probs)
}

#' Title
#'
#' @param epll  tbd
#' @param q  tbd
#'
#' @return eta
#' @export
get_eta <- function(epll, q) {
  eta <- epll - q / 2
  return(eta)
}


#' Title
#'
#' @param doses tbd
#' @param schedules tbd
#' @param lb tbd
#' @param threshold tbd
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
#' @param doses tbd
#' @param schedules tbd
#' @param quantiles tbd
#' @param emax tbd
#'
#' @return resp
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

#' Title
#'
#' @param design  tbd
#' @param doses  tbd
#' @param schedules  tbd
#' @param n_pat  tbd
#' @param likelihood  tbd
#' @param threshold  tbd
#'
#' @return design
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

#' Title
#'
#' @param design  tbd
#' @param doses  tbd
#' @param schedules  tbd
#' @param n_pat  tbd
#' @param likelihood  tbd
#' @param threshold  tbd
#'
#' @return design
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


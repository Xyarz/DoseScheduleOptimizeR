
#' @title gen_likelihood_data_linear
#'
#' @param data_mat matrix
#' @param initial_beta0 prior of beta0
#' @param initial_beta1 prior of coefficient for doses
#' @param initial_beta2 prior of coefficient for schedukles
#' @param initial_beta3 prior for interaction
#' @param initial_sigma_mu prior for mean of epsilon
#' @param initial_sigma_sd prior of sigma
#'
#' @return res a list containing patient level data
#' @export
gen_likelihood_data_linear <- function(
    data_mat,
    beta0,
    beta1,
    beta2,
    beta3,
    sig_mu,
    sig_sd,
    n_pat
    ) {

  data <- as.data.frame(as.table(data_mat))
  colnames(data) <- c("x", "y", "n")
  # Create patient level data
  patient_data <- data_frame(
    x = rep(NA, n_pat),
    y = rep(NA, n_pat),
    z = rep(NA, n_pat)
  )
  iter <- 1
  for(i in 1:nrow(data)) {
    if(iter <= n_pat) {
      if(!data$n[i] == 0) {
        for(j in 1:data$n[i]) {
          x <- as.numeric(as.character(data$x[i]))
          y <- as.numeric(as.character(data$y[i]))
          patient_data$x[iter] <- x
          patient_data$y[iter] <- y
          patient_data$z[iter] <- beta0 + beta1 * x + beta2 * y + beta3 * x * y + rnorm(n = 1, mean = sig_mu, sd = sig_sd)
          iter <- iter + 1
        }
      }
    }
  }

  res <- list(
    patient_data = patient_data,
    data = data
  )
  return(res)
}

#' @title gen_likelihood_data_emax
#'
#' @param data_mat  matrix
#' @param initial_E0 prior of the E0
#' @param initial_alpha1 prior of the doses EMax
#' @param initial_alpha2 prior of the schedules EMax
#' @param initial_delta1 prior of the doses ED50
#' @param initial_delta2 prior of the schedules ED50
#' @param initial_beta interaction
#' @param initial_sigma_mu mean of epsilon
#' @param initial_sigma_sd sigma
#'
#' @return res a list containing patient level data
#' @export
gen_likelihood_data_emax <- function(
    data_mat,
    E0,
    alpha1,
    alpha2,
    delta1,
    delta2,
    beta,
    sig_mu,
    sig_sd,
    n_pat
  ) {

  data <- as.data.frame(as.table(data_mat))
  colnames(data) <- c("x", "y", "n")
  # Create patient level data
  patient_data <- data_frame(
    x = rep(NA, n_pat),
    y = rep(NA, n_pat),
    z = rep(NA, n_pat)
  )

  iter <- 1
  for(i in 1:nrow(data)) {
    if(iter <= n_pat) {
      if(!data$n[i] == 0) {
        for(j in 1:data$n[i]) {
          x <- as.numeric(as.character(data$x[i]))
          y <- as.numeric(as.character(data$y[i]))
          patient_data$x[iter] <- x
          patient_data$y[iter] <- y
          patient_data$z[iter] <- E0 + alpha1 * x / (delta1 + x) + alpha2 * y / (delta2 + y) + beta * x / (delta1 + x) * y / (delta2 + y) + rnorm(n = 1, mean = sig_mu, sd = sig_sd)
          iter <- iter + 1
        }
      }
    }
  }
  res <- list(
    patient_data = patient_data,
    data = data
  )
  return(res)
}


#' @title linear
#'
#' @param data a matrix
#' @param n_chain number of chains to be used for the JAGS model, default 3
#' @param n_burn number of burn samples to be used for the JAGS model, default 1000
#' @param n_sample number of samples to be used for the JAGS model, default 1000
#' @param n_thin thinning rate to be used for the JAGS model, default 5
#' @param doses doses
#' @param schedules schedules
#' @param n_pat sample size
#' @param efficacy_threshold placebo adjusted threshold to be passed to be considered an MED
#' @param likelihood_data likelihood patient level data
#'
#' @return res, containing the patient level responses, as attributes also the MED, effect, all samples and patient level data
#' @export
linear <- function(
    data = n_matrix,
    doses,
    schedules,
    n_chain = 3,
    n_burn = 1000,
    n_sample = 1000,
    n_thin = 5,
    n_pat,
    efficacy_threshold,
    likelihood_data
    ) {
  # likelihood data fixed as param - done
  # estimate _ effect should be captured - done
  # patient lvl data mean/sd output - likelihood - done
  # data mat output - done
  # amount of patients at med - doable
  # pass true med based on likelihood params - done
  # success criterion being above certain threshold -> any med being above certain probability e.g. 50% -> TRUE for success - doable
  # interim size 50% starting point - done
  data <- likelihood_data
  patient_data <- data$patient_data
  data <- data$data
  model_string <- "
model {
  for (i in 1:N) {
      z[i] ~ dnorm(mu[i], prec)
      mu[i] <- beta0 + beta1*x[i] + beta2*y[i] + beta3*x[i]*y[i]
  }
  beta0 ~ dnorm(0, 1.0E-6)
  beta1 ~ dnorm(0, 1.0E-6)
  beta2 ~ dnorm(0, 1.0E-6)
  beta3 ~ dnorm(0, 1.0E-6)
  prec ~ dgamma(0.001, 0.001)
  sigma <- sqrt(1 / prec)
}
"
  # data prep for JAGS model
  jags_data <- list(
    z = patient_data$z,
    x = patient_data$x,
    y = patient_data$y,
    N = nrow(patient_data)
  )


  # model init
  model <- jags.model(textConnection(model_string), data = jags_data, n.chains = 3, quiet = TRUE)

  # define extract params
  params <- c("beta0", "beta1", "beta2", "beta3", "sigma")

  # Burn-in
  update(model, n.iter = n_burn, progress.bar = "none")

  # retrieve post samples
  out <- coda.samples(model, variable.names = params, n.iter = n_sample * n_thin)[[1]] %>%
    as.data.frame()

  med <- matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(med) <- doses
  colnames(med) <- schedules

  obs_effect <- matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(obs_effect) <- doses
  colnames(obs_effect) <- schedules

  res <- list()

  for (i in 1:nrow(out)) {
    beta0 <- out$beta0[i]
    beta1 <- out$beta1[i]
    beta2 <- out$beta2[i]
    beta3 <- out$beta3[i]
    sigma <- out$sigma[i]

    response <- beta0 + beta1 * patient_data$x + beta2 * patient_data$y + beta3 * patient_data$x * patient_data$y
    med <- calculate_MED_linear(beta0, beta1, beta2, beta3, efficacy_threshold, med = med, doses = doses, schedules = schedules, obs_effect = obs_effect)
    obs_effect <- attributes(med)$effect

    res[[i]] <- list(
      response = response
    )
  }

  attr(res, 'samples') <- out
  attr(res, 'pat_data') <- patient_data
  attr(res, 'MED') <- med
  attr(res, 'effect') <- obs_effect
  return(res)
}

#' @title calculate_MED_emax
#'
#' @param E0  E0
#' @param alpha1  EMax of doses
#' @param alpha2  EMax of schedules
#' @param beta  interaction
#' @param delta1  ED50 of doses
#' @param delta2  ED50 of schedules
#' @param efficacy_threshold  placebo adjusted threshold to be passed to be considered an MED
#' @param doses  doses
#' @param schedules  schedules
#' @param med  an empty matrix the med gets allocated on
#' @param obs_effect  an empty matrix the observed effects gets allocated on
#'
#' @return med and effect as attribute
#' @export
calculate_MED_emax <- function(E0, alpha1, alpha2, beta, delta1, delta2, efficacy_threshold, doses,
                          schedules, med, obs_effect) {

  placebo <- E0 + alpha1 * 0 / (delta1 + 0) + alpha2 * 1 / (delta2 + 1) + beta * 0 / (delta1 + 0) * 1 / (delta2 + 1)
  dose_range <- doses
  schedule_range <- schedules



  for (dose in dose_range) {
    for (schedule in schedule_range){
      efficacy <- E0 + alpha1 * dose / (delta1 + dose) + alpha2 * schedule / (delta2 + schedule) + beta * dose / (delta1 + dose) * schedule / (delta2 + schedule)
      effect <- efficacy - placebo
      if (effect > efficacy_threshold) {
        obs_effect[which(doses == dose), which(schedules == schedule)] <- obs_effect[which(doses == dose), which(schedules == schedule)] + effect
      }
    }
  }
  vec <- as.vector(unlist(obs_effect))
  if(any(vec > 0)) {
    min_val_gt_zero <- min(vec[vec > 0])
    pos <- which(obs_effect == min_val_gt_zero)
    med[pos] <- med[pos] + 1
  }
  attr(med, "effect") <- obs_effect
  return(med)
}

#' @title calculate_MED_linear
#'
#' @param beta0  beta0
#' @param beta1  coefficient of doses
#' @param beta2  coefficient of schedules
#' @param beta3  interaction
#' @param efficacy_threshold  placebo adjusted threshold to be passed to be considered an MED
#' @param doses  doses
#' @param schedules  schedules
#' @param med  an empty matrix the med gets allocated on
#' @param obs_effect  an empty matrix the observed effects gets allocated on
#'
#' @return med and effect as attribute
#' @export
calculate_MED_linear <- function(beta0, beta1, beta2, beta3, efficacy_threshold, doses,
                               schedules, med, obs_effect) {

  placebo <- beta0 + beta1 * 0 + beta2 * 1 + beta3 * 0 * 1
  # plot_res   <- 1e2
  # dose_range <- seq(from       = min(doses),
  #                 to         = max(doses),
  #                 length.out = plot_res)
  # schedule_range <- seq(from = min(schedules),
  #                 to         = max(schedules),
  #                 length.out = plot_res)
  dose_range <- doses
  schedule_range <- schedules

  for (dose in dose_range) {
    for (schedule in schedule_range){
      efficacy <- beta0 + beta1 * dose + beta2 * schedule + beta3 * dose * schedule
      effect <- efficacy - placebo
      if (effect > efficacy_threshold) {
        # med[which(doses == dose), which(schedules == schedule)] <- med[which(doses == dose), which(schedules == schedule)] + 1
        obs_effect[which(doses == dose), which(schedules == schedule)] <- obs_effect[which(doses == dose), which(schedules == schedule)] + effect
      }
    }
  }
  vec <- as.vector(unlist(obs_effect))
  if(any(vec > 0)) {
    min_val_gt_zero <- min(vec[vec > 0])
    pos <- which(obs_effect == min_val_gt_zero)
    med[pos] <- med[pos] + 1
  }

  attr(med, "effect") <- obs_effect
  return(med)
}

#' @title emax
#'
#' @param data  matrix
#' @param doses  doses
#' @param schedules  schedules
#' @param n_chain number of chains to be used for the JAGS model, default 3
#' @param n_burn number of burn samples to be used for the JAGS model, default 1000
#' @param n_sample number of samples to be used for the JAGS model, default 1000
#' @param n_thin thinning rate to be used for the JAGS model, default 5
#' @param n_pat  sample size
#' @param efficacy_threshold  placebo adjusted threshold to be passed to be considered an MED
#' @param likelihood_data  likelihood patient level data
#'
#' @return res, containing the patient level responses, as attributes also the MED, effect, all samples and patient level data
#' @export
emax <- function(
    data = n_matrix,
    doses,
    schedules,
    n_chain = 3,
    n_burn = 1000,
    n_sample = 1000,
    n_thin = 5,
    n_pat,
    efficacy_threshold,
    likelihood_data
) {

  data <- likelihood_data
  patient_data <- data$patient_data
  data <- data$data

  model_string <- paste0("
model {
    for (i in 1:N) {
      z[i] ~ dnorm(mu[i], prec)
      mu[i] <- E0 + alpha1*x[i]/(delta1 + x[i]) + alpha2*y[i]/(delta2 + y[i]) + beta*x[i]/(delta1 + x[i])*y[i]/(delta2 + y[i])
  }

  E0 ~ dnorm(0, 1.0E-6)
  alpha1 ~ dnorm(0, 1.0E-6)
  alpha2 ~ dnorm(0, 1.0E-6)
  beta ~ dnorm(0, 1.0E-6)
  delta1 ~ dunif(0, ",max(doses)," )
  delta2 ~ dunif(0, ",max(schedules),")
  prec ~ dgamma(0.001, 0.001)
  sigma <- sqrt(1 / prec)
}
")
# data prep for JAGS model
jags_data <- list(
  z = patient_data$z,
  x = patient_data$x,
  y = patient_data$y,
  N = nrow(patient_data)
)


# model init
model <- jags.model(textConnection(model_string), data = jags_data, n.chains = 3, quiet = TRUE)

# define extract params
params <- c("E0", "alpha1", "alpha2", "beta", "delta1", "delta2", "sigma")

# Burn-in
update(model, n.iter = n_burn, progress.bar = "none")

# retrieve post samples
out <- coda.samples(model, variable.names = params, n.iter = n_sample * n_thin)[[1]] %>%
  as.data.frame()
res <- list()



med <- matrix(0, nrow = length(doses), ncol = length(schedules))
rownames(med) <- doses
colnames(med) <- schedules

obs_effect <- matrix(0, nrow = length(doses), ncol = length(schedules))
rownames(obs_effect) <- doses
colnames(obs_effect) <- schedules

for (i in 1:nrow(out)) {
  E0 <- out$E0[i]
  alpha1 <- out$alpha1[i]
  alpha2 <- out$alpha2[i]
  beta <- out$beta[i]
  delta1 <- out$delta1[i]
  delta2 <- out$delta2[i]
  sigma <- out$sigma[i]

  response <- E0 + alpha1 * patient_data$x / (delta1 + patient_data$x) + alpha2 * patient_data$y / (delta2 + patient_data$y) + beta * patient_data$x / (delta1 + patient_data$x) * patient_data$y / (delta2 + patient_data$y)
  med <- calculate_MED_emax(E0, alpha1, alpha2, beta, delta1, delta2, efficacy_threshold, med = med, doses = doses, schedules = schedules, obs_effect = obs_effect)
  obs_effect <- attributes(med)$effect

  res[[i]] <- list(
    response = response
  )
}



# MED <- MED/length(out)

attr(res, 'samples') <- out
attr(res, 'pat_data') <- patient_data
attr(res, 'MED') <- med
attr(res, 'effect') <- obs_effect
return(res)
}

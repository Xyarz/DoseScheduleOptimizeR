simulate_final_trial <- function(
    scenario_label = "null",
    doses = c(0, 4, 5, 6, 9),
    schedules = c(1, 3, 4),
    threshold = 8,
    alloc_mat = 1,
    n_pat = 200,
    interim = TRUE,
    likelihood = "linear",
    E0 = -12.8,
    alpha1 = 10, # positive 15 null 0
    alpha2 = 10,
    delta1 = 3,
    delta2 = 1,
    beta = 0,
    sig_mu = 1,
    sig_sd = 1,

    beta0 = -12.8,
    beta1 = 0.1319648, # positive 0.3 null 0
    beta2 = 0.1319648,
    beta3 = 0
) {

  bool <- interim

  mats <- list()
  if(interim == TRUE) {

    # full factorial
    full_factorial <- generate_full_factorial_design(
      doses = doses,
      schedules = schedules,
      n_pat = n_pat
    )

    likelihood_data <- NULL
    med_mat <- matrix(0, nrow = length(doses), ncol = length(schedules))
    rownames(med_mat) <- doses
    colnames(med_mat) <- schedules

    obs_effect <- matrix(0, nrow = length(doses), ncol = length(schedules))
    rownames(obs_effect) <- doses
    colnames(obs_effect) <- schedules

    if(likelihood == "linear") {
      likelihood_data <- gen_likelihood_data_linear(
        data_mat = full_factorial,
        n_pat = n_pat,
        beta0 = beta0,
        beta1 = beta1,
        beta2 = beta2,
        beta3 = beta3,
        sig_mu = sig_mu,
        sig_sd = sig_sd
        )
      med_mat <- calculate_MED_linear(
        beta0 = beta0,
        beta1 = beta1,
        beta2 = beta2,
        beta3 = beta3,
        efficacy_threshold = threshold,
        doses = doses,
        schedules = schedules,
        med = med_mat,
        obs_effect = obs_effect
      )
    } else if(likelihood == "emax") {
      likelihood_data <- gen_likelihood_data_emax(
        data_mat = full_factorial,
        n_pat = n_pat,
        E0 = E0,
        alpha1 = alpha1,
        alpha2 = alpha2,
        delta1 = delta1,
        delta2 = delta2,
        beta = beta,
        sig_mu = sig_mu,
        sig_sd = sig_sd
        )
      med_mat <- calculate_MED_emax(
        E0 = E0,
        alpha1 = alpha1,
        alpha2 = alpha2,
        beta = beta,
        delta1 = delta1,
        delta2 = delta2,
        efficacy_threshold = threshold,
        doses = doses,
        schedules = schedules,
        med = med_mat,
        obs_effect = obs_effect
      )
    }

    n_pat_interim <- n_pat * 0.5
      matrix <- readRDS(paste0("data/alloc_mat",alloc_mat,".RDS"))

      d_opt_design_linear <- as.matrix(matrix[[1]][[1]])
      d_opt_design_emax <- as.matrix(matrix[[1]][[2]])
      i_opt_design_linear <- as.matrix(matrix[[1]][[3]])
      i_opt_design_emax <- as.matrix(matrix[[1]][[4]])


    comp_design <- generate_custom_corner_mid(n_pat = n_pat_interim, doses = doses, schedules = schedules)
    # comp_alt_design <- generate_custom_corner_mid(n_pat = n_pat)


    full_factorial_interim <- generate_full_factorial_design(
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim
    )


    # linear ------------------------------------------------------------------


    res_ff_linear <- analysis_linear(
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_linear <- analysis_linear(
      design = comp_design, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_linear <- analysis_linear(
      design = d_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_emax <- analysis_linear(
      design = d_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_linear <- analysis_linear(
      design = i_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_emax <- analysis_linear(
      design = i_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_ff_interim <- analysis_linear(
      design = full_factorial_interim, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    AIC_linear_ff <- attributes(res_ff_linear)$AIC
    BIC_linear_ff <- attributes(res_ff_linear)$BIC
    AIC_linear_comp <- attributes(res_comp_linear)$AIC
    BIC_linear_comp <- attributes(res_comp_linear)$BIC

    AIC_linear_d_opt_design_linear <- attributes(res_linear_d_opt_design_linear)$AIC
    BIC_linear_d_opt_design_linear <- attributes(res_linear_d_opt_design_linear)$BIC
    AIC_linear_d_opt_design_emax <- attributes(res_linear_d_opt_design_emax)$AIC
    BIC_linear_d_opt_design_emax <- attributes(res_linear_d_opt_design_emax)$BIC
    AIC_linear_i_opt_design_linear <- attributes(res_linear_i_opt_design_linear)$AIC
    BIC_linear_i_opt_design_linear <- attributes(res_linear_i_opt_design_linear)$BIC
    AIC_linear_i_opt_design_emax <- attributes(res_linear_i_opt_design_emax)$AIC
    BIC_linear_i_opt_design_emax <- attributes(res_linear_i_opt_design_emax)$BIC
    AIC_linear_ff_interim <- attributes(res_linear_ff_interim)$AIC
    BIC_linear_ff_interim <- attributes(res_linear_ff_interim)$BIC

    # emax --------------------------------------------------------------------


    res_ff_emax <- analysis_emax(
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_emax <- analysis_emax(
      design = comp_design, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_linear <- analysis_emax(
      design = d_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_emax <- analysis_emax(
      design = d_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_linear <- analysis_emax(
      design = i_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_emax <- analysis_emax(
      design = i_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_ff_interim <- analysis_emax(
      design = full_factorial_interim, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_interim,
      likelihood = likelihood_data,
      threshold = threshold
    )


    AIC_emax_ff <- attributes(res_ff_emax)$AIC
    BIC_emax_ff <- attributes(res_ff_emax)$BIC
    AIC_emax_comp <- attributes(res_comp_emax)$AIC
    BIC_emax_comp <- attributes(res_comp_emax)$BIC

    AIC_emax_d_opt_design_linear <- attributes(res_emax_d_opt_design_linear)$AIC
    BIC_emax_d_opt_design_linear <- attributes(res_emax_d_opt_design_linear)$BIC
    AIC_emax_d_opt_design_emax <- attributes(res_emax_d_opt_design_emax)$AIC
    BIC_emax_d_opt_design_emax <- attributes(res_emax_d_opt_design_emax)$BIC
    AIC_emax_i_opt_design_linear <- attributes(res_emax_i_opt_design_linear)$AIC
    BIC_emax_i_opt_design_linear <- attributes(res_emax_i_opt_design_linear)$BIC
    AIC_emax_i_opt_design_emax <- attributes(res_emax_i_opt_design_emax)$AIC
    BIC_emax_i_opt_design_emax <- attributes(res_emax_i_opt_design_emax)$BIC

    AIC_emax_ff_interim <- attributes(res_emax_ff_interim)$BIC
    BIC_emax_ff_interim <- attributes(res_emax_ff_interim)$BIC



    out_linear <- list(
      ff_MED = attributes(res_ff_linear)$MED,
      comp_MED = attributes(res_comp_linear)$MED,
      ff_effect = attributes(res_ff_emax)$effect,
      comp_effect = attributes(res_comp_emax)$effect,
      AIC_ff = AIC_linear_ff,
      AIC_comp = AIC_linear_comp,
      BIC_ff = BIC_linear_ff,
      BIC_comp = BIC_linear_comp,

      d_opt_linear_MED = attributes(res_linear_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_linear_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_linear_d_opt_design_linear,
      d_opt_linear_BIC = BIC_linear_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_linear_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_linear_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_linear_d_opt_design_emax,
      d_opt_emax_BIC = BIC_linear_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_linear_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_linear_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_linear_i_opt_design_linear,
      i_opt_linear_BIC = BIC_linear_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_linear_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_linear_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_linear_i_opt_design_emax,
      i_opt_emax_BIC = BIC_linear_i_opt_design_emax,

      ff_interim_MED = attributes(res_linear_ff_interim)$MED,
      ff_interim_effect = attributes(res_linear_ff_interim)$effect,
      ff_interim_AIC = AIC_linear_ff_interim,
      ff_interim_BIC = BIC_linear_ff_interim

    )

    out_emax <- list(
      ff_MED = attributes(res_ff_emax)$MED,
      comp_MED = attributes(res_comp_emax)$MED,
      ff_effect = attributes(res_ff_emax)$effect,
      comp_effect = attributes(res_comp_emax)$effect,
      AIC_ff = AIC_emax_ff,
      AIC_comp = AIC_emax_comp,
      BIC_ff = BIC_emax_ff,
      BIC_comp = BIC_emax_comp,

      d_opt_linear_MED = attributes(res_emax_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_emax_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_emax_d_opt_design_linear,
      d_opt_linear_BIC = BIC_emax_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_emax_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_emax_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_emax_d_opt_design_emax,
      d_opt_emax_BIC = BIC_emax_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_emax_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_emax_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_emax_i_opt_design_linear,
      i_opt_linear_BIC = BIC_emax_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_emax_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_emax_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_emax_i_opt_design_emax,
      i_opt_emax_BIC = BIC_emax_i_opt_design_emax,

      ff_interim_MED = attributes(res_emax_ff_interim)$MED,
      ff_interim_effect = attributes(res_emax_ff_interim)$effect,
      ff_interim_AIC = AIC_emax_ff_interim,
      ff_interim_BIC = BIC_emax_ff_interim
    )
    interim <- tibble(
      scenario_label = scenario_label,
      linear = list(out_linear),
      emax = list(out_emax)
    )


    n_pat_final <- n_pat * 0.5

    comp_design_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$comp_MED)
    rownames(comp_design_linear) <- doses
    colnames(comp_design_linear) <- schedules
    comp_design_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$comp_MED)
    rownames(comp_design_emax) <- doses
    colnames(comp_design_emax) <- schedules

    comp_linear_d_opt_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$d_opt_linear_MED)
    rownames(comp_linear_d_opt_linear) <- doses
    colnames(comp_linear_d_opt_linear) <- schedules
    comp_linear_d_opt_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$d_opt_linear_MED)
    rownames(comp_linear_d_opt_emax) <- doses
    colnames(comp_linear_d_opt_emax) <- schedules

    comp_emax_d_opt_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$d_opt_emax_MED)
    rownames(comp_emax_d_opt_linear) <- doses
    colnames(comp_emax_d_opt_linear) <- schedules
    comp_emax_d_opt_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$d_opt_emax_MED)
    rownames(comp_emax_d_opt_emax) <- doses
    colnames(comp_emax_d_opt_emax) <- schedules

    comp_linear_i_opt_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$i_opt_linear_MED)
    rownames(comp_linear_i_opt_linear) <- doses
    colnames(comp_linear_i_opt_linear) <- schedules
    comp_linear_i_opt_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$i_opt_linear_MED)
    rownames(comp_linear_i_opt_emax) <- doses
    colnames(comp_linear_i_opt_emax) <- schedules

    comp_emax_i_opt_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$i_opt_emax_MED)
    rownames(comp_emax_i_opt_linear) <- doses
    colnames(comp_emax_i_opt_linear) <- schedules
    comp_emax_i_opt_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$i_opt_emax_MED)
    rownames(comp_emax_i_opt_emax) <- doses
    colnames(comp_emax_i_opt_emax) <- schedules

    ff_interim_linear <- generate_custom_design(n_pat_final, post_probs = interim$linear[[1]]$ff_interim_MED)
    rownames(ff_interim_linear) <- doses
    colnames(ff_interim_linear) <- schedules
    ff_interim_emax <- generate_custom_design(n_pat_final, post_probs = interim$emax[[1]]$ff_interim_MED)
    rownames(ff_interim_emax) <- doses
    colnames(ff_interim_emax) <- schedules

    # full_factorial_linear <- generate_custom_design(n_pat, post_probs = interim$linear[[1]]$ff_MED)
    # rownames(full_factorial_linear) <- doses
    # colnames(full_factorial_linear) <- schedules
    # full_factorial_emax <- generate_custom_design(n_pat, post_probs = interim$emax[[1]]$ff_MED)
    # rownames(full_factorial_emax) <- doses
    # colnames(full_factorial_emax) <- schedules

    res_ff_linear <- analysis_linear(
      # design = full_factorial_linear,
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_linear <- analysis_linear(
      design = comp_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_linear <- analysis_linear(
      design = comp_linear_d_opt_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_emax <- analysis_linear(
      design = comp_emax_d_opt_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_linear <- analysis_linear(
      design = comp_linear_i_opt_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_emax <- analysis_linear(
      design = comp_emax_i_opt_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_ff_interim <- analysis_linear(
      design = ff_interim_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_ff_emax <- analysis_emax(
      # design = full_factorial_emax,
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_emax <- analysis_emax(
      design = comp_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_linear <- analysis_emax(
      design = comp_linear_d_opt_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_emax <- analysis_emax(
      design = comp_emax_d_opt_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_linear <- analysis_emax(
      design = comp_linear_i_opt_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_emax <- analysis_emax(
      design = comp_emax_i_opt_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )
    res_emax_ff_interim <- analysis_emax(
      design = ff_interim_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat_final,
      likelihood = likelihood_data,
      threshold = threshold
    )

    mats[["ff_mat"]] <- full_factorial
    mats[["comp_mat"]] <- comp_design
    # mats[["ff_mat_linear"]] <- full_factorial_linear
    # mats[["ff_mat_emax"]] <- full_factorial_emax
    mats[["comp_mat_linear"]] <- comp_design_linear
    mats[["comp_mat_emax"]] <- comp_design_emax
    mats[["comp_linear_d_opt_linear"]] <- comp_linear_d_opt_linear
    mats[["comp_emax_d_opt_linear"]] <- comp_emax_d_opt_linear
    mats[["comp_linear_i_opt_linear"]] <- comp_linear_i_opt_linear
    mats[["comp_emax_i_opt_linear"]] <- comp_emax_i_opt_linear
    mats[["comp_linear_d_opt_emax"]] <- comp_linear_d_opt_emax
    mats[["comp_emax_d_opt_emax"]] <- comp_emax_d_opt_emax
    mats[["comp_linear_i_opt_emax"]] <- comp_linear_i_opt_emax
    mats[["comp_emax_i_opt_emax"]] <- comp_emax_i_opt_emax
    mats[["ff_interim_linear"]] <- ff_interim_linear
    mats[["ff_interim_emax"]] <- ff_interim_emax
  } else {
    comp_design <- generate_custom_corner_mid(n_pat = n_pat, doses = doses, schedules = schedules)
    if(interim) {
      matrix <- readRDS(paste0("data/alloc_mat",alloc_mat,".RDS"))
      d_opt_design_linear <- as.matrix(matrix[[2]][[1]])
      d_opt_design_emax <- as.matrix(matrix[[2]][[2]])
      i_opt_design_linear <- as.matrix(matrix[[2]][[3]])
      i_opt_design_emax <- as.matrix(matrix[[2]][[4]])
    } else {
      matrix <- readRDS(paste0("data/alloc_mat_no_interim",alloc_mat,".RDS"))
      d_opt_design_linear <- as.matrix(matrix[[1]])
      d_opt_design_emax <- as.matrix(matrix[[2]])
      i_opt_design_linear <- as.matrix(matrix[[3]])
      i_opt_design_emax <- as.matrix(matrix[[4]])
    }





    # full factorial
    full_factorial <- generate_full_factorial_design(
      doses = doses,
      schedules = schedules,
      n_pat = n_pat
    )

    likelihood_data <- NULL
    med_mat <- matrix(0, nrow = length(doses), ncol = length(schedules))
    rownames(med_mat) <- doses
    colnames(med_mat) <- schedules

    obs_effect <- matrix(0, nrow = length(doses), ncol = length(schedules))
    rownames(obs_effect) <- doses
    colnames(obs_effect) <- schedules

    if(likelihood == "linear") {
      likelihood_data <- gen_likelihood_data_linear(
        data_mat = full_factorial,
        n_pat = n_pat,
        beta0 = beta0,
        beta1 = beta1,
        beta2 = beta2,
        beta3 = beta3,
        sig_mu = sig_mu,
        sig_sd = sig_sd
        )
      med_mat <- calculate_MED_linear(
        beta0 = beta0,
        beta1 = beta1,
        beta2 = beta2,
        beta3 = beta3,
        efficacy_threshold = threshold,
        doses = doses,
        schedules = schedules,
        med = med_mat,
        obs_effect = obs_effect
      )
    } else if(likelihood == "emax") {
      likelihood_data <- gen_likelihood_data_emax(
        data_mat = full_factorial,
        n_pat = n_pat,
        E0 = E0,
        alpha1 = alpha1,
        alpha2 = alpha2,
        delta1 = delta1,
        delta2 = delta2,
        beta = beta,
        sig_mu = sig_mu,
        sig_sd = sig_sd
        )
      med_mat <- calculate_MED_emax(
        E0 = E0,
        alpha1 = alpha1,
        alpha2 = alpha2,
        beta = beta,
        delta1 = delta1,
        delta2 = delta2,
        efficacy_threshold = threshold,
        doses = doses,
        schedules = schedules,
        med = med_mat,
        obs_effect = obs_effect
      )
    }

    res_ff_linear <- analysis_linear(
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_linear <- analysis_linear(
      design = comp_design, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_linear <- analysis_linear(
      design = d_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_d_opt_design_emax <- analysis_linear(
      design = d_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_linear <- analysis_linear(
      design = i_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_linear_i_opt_design_emax <- analysis_linear(
      design = i_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_ff_emax <- analysis_emax(
      design = full_factorial,
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_comp_emax <- analysis_emax(
      design = comp_design, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_linear <- analysis_emax(
      design = d_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_d_opt_design_emax <- analysis_emax(
      design = d_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_linear <- analysis_emax(
      design = i_opt_design_linear, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    res_emax_i_opt_design_emax <- analysis_emax(
      design = i_opt_design_emax, #$mat if optimal
      doses = doses,
      schedules = schedules,
      n_pat = n_pat,
      likelihood = likelihood_data,
      threshold = threshold
    )

    mats[["ff_mat"]] <- full_factorial
    mats[["comp_mat"]] <- comp_design
    interim <- NULL
  }

  # linear ------------------------------------------------------------------

  AIC_linear_ff <- attributes(res_ff_linear)$AIC
  AIC_linear_comp <- attributes(res_ff_linear)$BIC
  BIC_linear_ff <- attributes(res_comp_linear)$AIC
  BIC_linear_comp <- attributes(res_comp_linear)$BIC

  AIC_linear_d_opt_design_linear <- attributes(res_linear_d_opt_design_linear)$AIC
  BIC_linear_d_opt_design_linear <- attributes(res_linear_d_opt_design_linear)$BIC
  AIC_linear_d_opt_design_emax <- attributes(res_linear_d_opt_design_emax)$AIC
  BIC_linear_d_opt_design_emax <- attributes(res_linear_d_opt_design_emax)$BIC
  AIC_linear_i_opt_design_linear <- attributes(res_linear_i_opt_design_linear)$AIC
  BIC_linear_i_opt_design_linear <- attributes(res_linear_i_opt_design_linear)$BIC
  AIC_linear_i_opt_design_emax <- attributes(res_linear_i_opt_design_emax)$AIC
  BIC_linear_i_opt_design_emax <- attributes(res_linear_i_opt_design_emax)$BIC

  # emax --------------------------------------------------------------------

  AIC_emax_ff <- attributes(res_ff_emax)$AIC
  AIC_emax_comp <- attributes(res_ff_emax)$BIC
  BIC_emax_ff <- attributes(res_comp_emax)$AIC
  BIC_emax_comp <- attributes(res_comp_emax)$BIC

  AIC_emax_d_opt_design_linear <- attributes(res_emax_d_opt_design_linear)$AIC
  BIC_emax_d_opt_design_linear <- attributes(res_emax_d_opt_design_linear)$BIC
  AIC_emax_d_opt_design_emax <- attributes(res_emax_d_opt_design_emax)$AIC
  BIC_emax_d_opt_design_emax <- attributes(res_emax_d_opt_design_emax)$BIC
  AIC_emax_i_opt_design_linear <- attributes(res_emax_i_opt_design_linear)$AIC
  BIC_emax_i_opt_design_linear <- attributes(res_emax_i_opt_design_linear)$BIC
  AIC_emax_i_opt_design_emax <- attributes(res_emax_i_opt_design_emax)$AIC
  BIC_emax_i_opt_design_emax <- attributes(res_emax_i_opt_design_emax)$BIC

  if(bool == TRUE) {
    out_linear <- list(
      ff_MED = attributes(res_ff_linear)$MED,
      comp_MED = attributes(res_comp_linear)$MED,
      ff_effect = attributes(res_ff_linear)$effect,
      comp_effect = attributes(res_comp_linear)$effect,
      AIC_ff = AIC_linear_ff,
      AIC_comp = AIC_linear_comp,
      BIC_ff = BIC_linear_ff,
      BIC_comp = BIC_linear_comp,

      d_opt_linear_MED = attributes(res_linear_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_linear_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_linear_d_opt_design_linear,
      d_opt_linear_BIC = BIC_linear_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_linear_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_linear_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_linear_d_opt_design_emax,
      d_opt_emax_BIC = BIC_linear_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_linear_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_linear_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_linear_i_opt_design_linear,
      i_opt_linear_BIC = BIC_linear_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_linear_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_linear_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_linear_i_opt_design_emax,
      i_opt_emax_BIC = BIC_linear_i_opt_design_emax,
      ff_interim_MED = attributes(res_linear_ff_interim)$MED,
      ff_interim_effect = attributes(res_linear_ff_interim)$effect,
      ff_interim_AIC = attributes(res_linear_ff_interim)$AIC,
      ff_interim_BIC = attributes(res_linear_ff_interim)$BIC
    )

    out_emax <- list(
      ff_MED = attributes(res_ff_emax)$MED,
      comp_MED = attributes(res_comp_emax)$MED,
      ff_effect = attributes(res_ff_emax)$effect,
      comp_effect = attributes(res_comp_emax)$effect,
      AIC_ff = AIC_emax_ff,
      AIC_comp = AIC_emax_comp,
      BIC_ff = BIC_emax_ff,
      BIC_comp = BIC_emax_comp,

      d_opt_linear_MED = attributes(res_emax_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_emax_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_emax_d_opt_design_linear,
      d_opt_linear_BIC = BIC_emax_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_emax_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_emax_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_emax_d_opt_design_emax,
      d_opt_emax_BIC = BIC_emax_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_emax_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_emax_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_emax_i_opt_design_linear,
      i_opt_linear_BIC = BIC_emax_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_emax_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_emax_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_emax_i_opt_design_emax,
      i_opt_emax_BIC = BIC_emax_i_opt_design_emax,
      ff_interim_MED = attributes(res_emax_ff_interim)$MED,
      ff_interim_effect = attributes(res_emax_ff_interim)$effect,
      ff_interim_AIC = attributes(res_emax_ff_interim)$AIC,
      ff_interim_BIC = attributes(res_emax_ff_interim)$BIC
    )

  } else {
    out_linear <- list(
      ff_MED = attributes(res_ff_linear)$MED,
      comp_MED = attributes(res_comp_linear)$MED,
      ff_effect = attributes(res_ff_linear)$effect,
      comp_effect = attributes(res_comp_linear)$effect,
      AIC_ff = AIC_linear_ff,
      AIC_comp = AIC_linear_comp,
      BIC_ff = BIC_linear_ff,
      BIC_comp = BIC_linear_comp,

      d_opt_linear_MED = attributes(res_linear_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_linear_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_linear_d_opt_design_linear,
      d_opt_linear_BIC = BIC_linear_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_linear_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_linear_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_linear_d_opt_design_emax,
      d_opt_emax_BIC = BIC_linear_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_linear_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_linear_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_linear_i_opt_design_linear,
      i_opt_linear_BIC = BIC_linear_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_linear_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_linear_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_linear_i_opt_design_emax,
      i_opt_emax_BIC = BIC_linear_i_opt_design_emax
    )

    out_emax <- list(
      ff_MED = attributes(res_ff_emax)$MED,
      comp_MED = attributes(res_comp_emax)$MED,
      ff_effect = attributes(res_ff_emax)$effect,
      comp_effect = attributes(res_comp_emax)$effect,
      AIC_ff = AIC_emax_ff,
      AIC_comp = AIC_emax_comp,
      BIC_ff = BIC_emax_ff,
      BIC_comp = BIC_emax_comp,

      d_opt_linear_MED = attributes(res_emax_d_opt_design_linear)$MED,
      d_opt_linear_effect = attributes(res_emax_d_opt_design_linear)$effect,
      d_opt_linear_AIC = AIC_emax_d_opt_design_linear,
      d_opt_linear_BIC = BIC_emax_d_opt_design_linear,

      d_opt_emax_MED = attributes(res_emax_d_opt_design_emax)$MED,
      d_opt_emax_effect = attributes(res_emax_d_opt_design_emax)$effect,
      d_opt_emax_AIC = AIC_emax_d_opt_design_emax,
      d_opt_emax_BIC = BIC_emax_d_opt_design_emax,

      i_opt_linear_MED = attributes(res_emax_i_opt_design_linear)$MED,
      i_opt_linear_effect = attributes(res_emax_i_opt_design_linear)$effect,
      i_opt_linear_AIC = AIC_emax_i_opt_design_linear,
      i_opt_linear_BIC = BIC_emax_i_opt_design_linear,

      i_opt_emax_MED = attributes(res_emax_i_opt_design_emax)$MED,
      i_opt_emax_effect = attributes(res_emax_i_opt_design_emax)$effect,
      i_opt_emax_AIC = AIC_emax_i_opt_design_emax,
      i_opt_emax_BIC = BIC_emax_i_opt_design_emax
    )

  }

  df <- likelihood_data$patient_data

  result_mean <- df %>%
    group_by(x, y) %>%
    summarise(mean_z = mean(z, na.rm = TRUE))

  # Spread the data into a wide format to create the matrix
  mean_z_matrix <- result_mean %>%
    spread(key = y, value = mean_z)

  # Convert to a matrix
  mean_z_matrix <- as.matrix(mean_z_matrix)

  # Set row names to x values
  rownames(mean_z_matrix) <- mean_z_matrix[, 1]
  obs_mean <- mean_z_matrix[, -1]

  result_sd <- df %>%
    group_by(x, y) %>%
    summarise(mean_z = sd(z, na.rm = TRUE))

  # Spread the data into a wide format to create the matrix
  mean_z_matrix <- result_sd %>%
    spread(key = y, value = mean_z)

  # Convert to a matrix
  mean_z_matrix <- as.matrix(mean_z_matrix)

  # Set row names to x values
  rownames(mean_z_matrix) <- mean_z_matrix[, 1]
  obs_sd <- mean_z_matrix[, -1]


  out <- tibble(
    scenario_label = list(scenario_label),
    linear = list(out_linear),
    emax = list(out_emax),
    obs_mean = list(obs_mean),
    obs_sd = list(obs_sd),
    # interim = list(interim),
    mats = list(mats),
    true_med = list(med_mat),
    true_effect = list(attributes(med_mat)$effect)
  )

  return(out)
}

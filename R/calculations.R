# RME ---------------------------------------------------------------------


#' @title get_effect_true
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param var  variable to be looked at, numeric
#'
#' @return sum numeric value of the effect
#' @export
get_effect_true <- function(results, scenario, var) {
  sum <- numeric(length(results[[scenario]][[var]][[1]]))

  for(i in 1:length(results[[scenario]]$linear)) {
    if(i == 1) {
      sum <- results[[scenario]][[var]][[i]]
    } else {
      sum <- sum + results[[scenario]][[var]][[i]]
    }
  }

  return(sum)
}

#' @title get_effect
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param likelihood  true model, either "linear" or "emax"
#' @param var  numeric value of the value to be looked at
#'
#' @return sum
#' @export
get_effect <- function(results, scenario, likelihood, var) {
  sum <- numeric(length(results[[scenario]][[likelihood]]))
  for(i in 1:length(results[[scenario]][[likelihood]])) {
    pos <- which(results[[scenario]][[likelihood]][[i]][[var]] / 5000 > 8)[1]
    true_med_effect <- results[[scenario]][["true_effect"]][[i]][pos]
    obs_effect <- results[[scenario]][[likelihood]][[i]][[var]][pos] / 5000
    sum[i] <- (obs_effect - true_med_effect)^2

  }

  # sum_1 <- sum / length(results[[scenario]][[likelihood]])
  return(sum)
}

#' @title calc_MSE_RMSE_across_scenarios
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param interim  boolean, whether interim was conducted or not
#'
#' @return list containg MSE and RMSE
#' @export
calc_MSE_RMSE_across_scenarios <- function(results, scenario, interim) {
  observed_effect_1_lin <- get_effect(results, scenario, "linear", 3)
  observed_effect_2_lin <- get_effect(results, scenario, "linear", 4)
  observed_effect_3_lin <- get_effect(results, scenario, "linear", 10)
  observed_effect_4_lin <- get_effect(results, scenario, "linear", 14)
  observed_effect_5_lin <- get_effect(results, scenario, "linear", 18)
  observed_effect_6_lin <- get_effect(results, scenario, "linear", 22)

  observed_effect_1_ema <- get_effect(results, scenario, "emax", 3)
  observed_effect_2_ema <- get_effect(results, scenario, "emax", 4)
  observed_effect_3_ema <- get_effect(results, scenario, "emax", 10)
  observed_effect_4_ema <- get_effect(results, scenario, "emax", 14)
  observed_effect_5_ema <- get_effect(results, scenario, "emax", 18)
  observed_effect_6_ema <- get_effect(results, scenario, "emax", 22)

  if(interim) {
    observed_effect_7_lin <- get_effect(results, scenario, "linear", 26)
    observed_effect_7_ema <- get_effect(results, scenario, "emax", 26)

    df <- list(
      MSE = c(
        mean((observed_effect_1_lin), na.rm = T),
        mean((observed_effect_2_lin), na.rm = T),
        mean((observed_effect_3_lin), na.rm = T),
        mean((observed_effect_4_lin), na.rm = T),
        mean((observed_effect_5_lin), na.rm = T),
        mean((observed_effect_6_lin), na.rm = T),
        mean((observed_effect_7_lin), na.rm = T),
        mean((observed_effect_1_ema), na.rm = T),
        mean((observed_effect_2_ema), na.rm = T),
        mean((observed_effect_3_ema), na.rm = T),
        mean((observed_effect_4_ema), na.rm = T),
        mean((observed_effect_5_ema), na.rm = T),
        mean((observed_effect_6_ema), na.rm = T),
        mean((observed_effect_7_ema), na.rm = T)
      ),
      RMSE = c(
        sqrt(mean((observed_effect_1_lin), na.rm = T)),
        sqrt(mean((observed_effect_2_lin), na.rm = T)),
        sqrt(mean((observed_effect_3_lin), na.rm = T)),
        sqrt(mean((observed_effect_4_lin), na.rm = T)),
        sqrt(mean((observed_effect_5_lin), na.rm = T)),
        sqrt(mean((observed_effect_6_lin), na.rm = T)),
        sqrt(mean((observed_effect_7_lin), na.rm = T)),
        sqrt(mean((observed_effect_1_ema), na.rm = T)),
        sqrt(mean((observed_effect_2_ema), na.rm = T)),
        sqrt(mean((observed_effect_3_ema), na.rm = T)),
        sqrt(mean((observed_effect_4_ema), na.rm = T)),
        sqrt(mean((observed_effect_5_ema), na.rm = T)),
        sqrt(mean((observed_effect_6_ema), na.rm = T)),
        sqrt(mean((observed_effect_7_ema), na.rm = T))
      )
    )
  } else {
    df <- list(
      MSE = c(
        mean((observed_effect_1_lin), na.rm = T),
        mean((observed_effect_2_lin), na.rm = T),
        mean((observed_effect_3_lin), na.rm = T),
        mean((observed_effect_4_lin), na.rm = T),
        mean((observed_effect_5_lin), na.rm = T),
        mean((observed_effect_6_lin), na.rm = T),
        mean((observed_effect_7_lin), na.rm = T),
        mean((observed_effect_1_ema), na.rm = T),
        mean((observed_effect_2_ema), na.rm = T),
        mean((observed_effect_3_ema), na.rm = T),
        mean((observed_effect_4_ema), na.rm = T),
        mean((observed_effect_5_ema), na.rm = T),
        mean((observed_effect_6_ema), na.rm = T)
      ),
      RMSE = c(
        sqrt(mean((observed_effect_1_lin), na.rm = T)),
        sqrt(mean((observed_effect_2_lin), na.rm = T)),
        sqrt(mean((observed_effect_3_lin), na.rm = T)),
        sqrt(mean((observed_effect_4_lin), na.rm = T)),
        sqrt(mean((observed_effect_5_lin), na.rm = T)),
        sqrt(mean((observed_effect_6_lin), na.rm = T)),
        sqrt(mean((observed_effect_7_lin), na.rm = T)),
        sqrt(mean((observed_effect_1_ema), na.rm = T)),
        sqrt(mean((observed_effect_2_ema), na.rm = T)),
        sqrt(mean((observed_effect_3_ema), na.rm = T)),
        sqrt(mean((observed_effect_4_ema), na.rm = T)),
        sqrt(mean((observed_effect_5_ema), na.rm = T)),
        sqrt(mean((observed_effect_6_ema), na.rm = T))
      )
    )
  }


  return(df)
}

# pats at med -------------------------------------------------------------

#' @title get_alloc
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param mat  matrix to be looked at
#'
#' @return sum
#' @export
get_alloc <- function(results, scenario, mat) {
  sum <- numeric(length(results[[scenario]][["linear"]]))

  for(i in 1:length(results[[scenario]]$linear)) {
    true_med <- results[[scenario]][["true_med"]][[i]]
    true_med_pos <- which(true_med > 0)[1]
    if(is.na(true_med_pos)) {
      sum[i] <- NA
    }else{
      sum[i] <- results[[scenario]]$mats[[i]][[mat]][[true_med_pos]]
    }
  }

  sum_1 <- mean(sum)
  return(sum_1)
}

#' @title get_MED_alloc
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param mat  matrix to be looked at
#'
#' @return sum
#' @export
get_MED_alloc <- function(results, scenario, mat) {

  ff_alloc <- generate_full_factorial_design(
    doses = results[[scenario]]$doses[[1]],
    schedules = results[[scenario]]$schedules[[1]],
    n_pat = results[[scenario]]$n_pat[[1]] * 0.5
  )

  alpha1 <- results[[scenario]]$alpha1[[1]]
  alpha2 <- results[[scenario]]$alpha2[[1]]
  delta1 <- results[[scenario]]$delta1[[1]]
  delta2 <- results[[scenario]]$delta2[[1]]
  n_pat <- results[[scenario]]$n_pat[[1]] * 0.5
  d_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules, n_pat = n_pat)$mat
  d_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,
                                                    emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat
  i_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I")$mat
  i_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I",
                                                    emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat




  sum_1 = get_alloc(results, scenario, 1)
  sum_2 = get_alloc(results, scenario, 2)
  sum_3 = get_alloc(results, scenario, 3)
  sum_4 = get_alloc(results, scenario, 4)
  sum_5 = get_alloc(results, scenario, 5)
  sum_6 = get_alloc(results, scenario, 6)
  sum_7 = get_alloc(results, scenario, 7)
  sum_8 = get_alloc(results, scenario, 8)
  sum_9 = get_alloc(results, scenario, 9)
  sum_10 = get_alloc(results, scenario, 10)
  sum_11 = get_alloc(results, scenario, 11)
  sum_12 = get_alloc(results, scenario, 12)
  sum_13 = get_alloc(results, scenario, 13)
  sum_14 = get_alloc(results, scenario, 14)


  true_med <- results[[scenario]][["true_med"]][[1]]
  pos_true_med <- which(true_med > 0)[1]
  sum_mat1_d_interim = d_opt_design_linear[pos_true_med]
  sum_mat2_d_interim = d_opt_design_emax[pos_true_med]
  sum_mat1_i_interim = i_opt_design_linear[pos_true_med]
  sum_mat2_i_interim = i_opt_design_emax[pos_true_med]

  sum_values_ff_mat = sum_1
  sum_values_comp_mat = sum_2
  sum_values_comp_mat_linear = sum_3 + sum_values_comp_mat
  sum_values_comp_mat_emax = sum_4 + sum_values_comp_mat
  sum_values_comp_linear_d_opt_linear = sum_5 + sum_mat1_d_interim
  sum_values_comp_emax_d_opt_linear = sum_6 + sum_mat2_d_interim
  sum_values_comp_linear_i_opt_linear = sum_7 + sum_mat1_i_interim
  sum_values_comp_emax_i_opt_linear = sum_8 + sum_mat2_i_interim
  sum_values_comp_linear_d_opt_emax = sum_9 + sum_mat1_d_interim
  sum_values_comp_emax_d_opt_emax = sum_10 + sum_mat2_d_interim
  sum_values_comp_linear_i_opt_emax = sum_11 + sum_mat1_i_interim
  sum_values_comp_emax_i_opt_emax = sum_12 + sum_mat2_i_interim
  sum_values_ff_interim_linear = sum_13 + ff_alloc[pos_true_med]
  sum_values_ff_interim_emax = sum_14 + ff_alloc[pos_true_med]
  res <- list(
    sum_values_ff_mat,
    sum_values_comp_mat_linear,
    sum_values_comp_linear_d_opt_linear,
    sum_values_comp_emax_d_opt_linear,
    sum_values_comp_linear_i_opt_linear,
    sum_values_comp_emax_i_opt_linear,
    sum_values_ff_interim_linear,
    sum_values_comp_mat_emax,
    sum_values_comp_linear_d_opt_emax,
    sum_values_comp_emax_d_opt_emax,
    sum_values_comp_linear_i_opt_emax,
    sum_values_comp_emax_i_opt_emax,
    sum_values_ff_interim_emax
  )

  return(res)
}

#' @title get_MED_alloc_no_interim
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param mat  matrix to be looked at
#'
#' @return sum
#' @export
get_MED_alloc_no_interim <- function(results, scenario, mat) {

  ff_alloc <- generate_full_factorial_design(
    doses = results[[scenario]]$doses[[1]],
    schedules = results[[scenario]]$schedules[[1]],
    n_pat = results[[scenario]]$n_pat[[1]]
  )
  true_med <- results[[scenario]][["true_med"]][[1]]
  pos_true_med <- which(true_med > 0)[1]

  alpha1 <- results[[scenario]]$alpha1[[1]]
  alpha2 <- results[[scenario]]$alpha2[[1]]
  delta1 <- results[[scenario]]$delta1[[1]]
  delta2 <- results[[scenario]]$delta2[[1]]
  n_pat <- results[[scenario]]$n_pat[[1]]
  d_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules, n_pat = n_pat)$mat
  d_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,
                                                    emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat
  i_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I")$mat
  i_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I",
                                                    emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat


  sum_1 = get_alloc(results, scenario, 1)
  sum_2 = get_alloc(results, scenario, 2)
  sum_mat1_d = d_opt_design_linear[pos_true_med]
  sum_mat2_d = d_opt_design_emax[pos_true_med]
  sum_mat1_i = i_opt_design_linear[pos_true_med]
  sum_mat2_i = i_opt_design_emax[pos_true_med]

  sum_values_ff_mat = sum_1
  sum_values_comp_mat = sum_2
  sum_values_comp_d_opt_linear = sum_mat1_d
  sum_values_comp_d_opt_emax = sum_mat2_d
  sum_values_comp_i_opt_linear = sum_mat1_i
  sum_values_comp_i_opt_emax = sum_mat2_i

  res <- list(
    sum_values_ff_mat,
    sum_values_comp_mat,
    sum_values_comp_d_opt_linear,
    sum_values_comp_d_opt_emax,
    sum_values_comp_i_opt_linear,
    sum_values_comp_i_opt_emax
  )

  return(res)
}



# correct med -------------------------------------------------------------

#' @title get_correct_med_est
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param likelihood  tbd
#' @param var  tbd
#'
#' @return tbd
#' @export
get_correct_med_est <- function(results, scenario, likelihood, var) {
  sum <- numeric(length(results[[scenario]][[likelihood]]))

  for(i in 1:length(results[[scenario]][[likelihood]])) {
    if(i == 1) {
      sum <- as.numeric(any(results[[scenario]][[likelihood]][[i]][[var]][which(results[[scenario]]$true_med[[1]] > 0)[1]] > 0))
    } else {
      sum <- sum +  as.numeric(any(results[[scenario]][[likelihood]][[i]][[var]][which(results[[scenario]]$true_med[[1]] > 0)[1]] > 0))
    }
  }

  sum_1 <- (sum / length(results[[scenario]][[likelihood]]))
  return(sum_1)
}

#' @title get_all_correct_med_est
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param interim  boolean, indicating whether an interim analysis has been conducted
#'
#' @return list
#' @export
get_all_correct_med_est <- function(results, scenario, interim){
  med_1 = get_correct_med_est(results, scenario, "linear", 1)
  med_2 = get_correct_med_est(results, scenario, "linear", 2)
  med_3 = get_correct_med_est(results, scenario, "linear", 9)
  med_4 = get_correct_med_est(results, scenario, "linear", 13)
  med_5 = get_correct_med_est(results, scenario, "linear", 17)
  med_6 = get_correct_med_est(results, scenario, "linear", 21)
  med_7 = get_correct_med_est(results, scenario, "emax", 1)
  med_8 = get_correct_med_est(results, scenario, "emax", 2)
  med_9 = get_correct_med_est(results, scenario, "emax", 9)
  med_10 = get_correct_med_est(results, scenario, "emax", 13)
  med_11 = get_correct_med_est(results, scenario, "emax", 17)
  med_12 = get_correct_med_est(results, scenario, "emax", 21)

  if(interim) {
    med_13 = get_correct_med_est(results, scenario, "linear", 25)
    med_14 = get_correct_med_est(results, scenario, "emax", 25)

    list(
      med_1, med_2, med_3, med_4, med_5, med_6, med_13,
      med_7, med_8, med_9, med_10, med_11, med_12, med_14
    )
  } else {
    list(
      med_1, med_2, med_3, med_4, med_5, med_6,
      med_7, med_8, med_9, med_10, med_11, med_12
    )
  }

}


# gonogo ------------------------------------------------------------------

#' @title get_GoNoGo
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param GoNoGo_Threshold  gonogo threshold to be passed, e.g., 0.5
#' @param interim  boolean, indicating whether an interim analysis has been conducted
#'
#' @return list
#' @export
get_GoNoGo <- function(results, scenario, GoNoGo_Threshold, interim) {

  if(interim) {
    success_linear_ff <- numeric(length(results[[1]]$linear))
    success_linear_comp <- numeric(length(results[[1]]$linear))
    success_linear_d_opt_linear <- numeric(length(results[[1]]$linear))
    success_linear_d_opt_emax <- numeric(length(results[[1]]$linear))
    success_linear_i_opt_linear <- numeric(length(results[[1]]$linear))
    success_linear_i_opt_emax <- numeric(length(results[[1]]$linear))
    success_linear_ff_interim <- numeric(length(results[[1]]$linear))

    success_emax_ff <- numeric(length(results[[1]]$linear))
    success_emax_comp <- numeric(length(results[[1]]$linear))
    success_emax_d_opt_linear <- numeric(length(results[[1]]$linear))
    success_emax_d_opt_emax <- numeric(length(results[[1]]$linear))
    success_emax_i_opt_linear <- numeric(length(results[[1]]$linear))
    success_emax_i_opt_emax <- numeric(length(results[[1]]$linear))
    success_emax_ff_interim <- numeric(length(results[[1]]$linear))

    for(i in 1:length(results[[1]]$linear)) {
      res_lin <- results[[scenario]]$linear[[i]]
      res_ema <- results[[scenario]]$emax[[i]]


      success_linear_ff[i] <-  ifelse(any((res_lin[[1]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_comp[i] <-  ifelse(any((res_lin[[2]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_d_opt_linear[i] <-  ifelse(any((res_lin[[9]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_d_opt_emax[i] <-  ifelse(any((res_lin[[13]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_i_opt_linear[i] <-  ifelse(any((res_lin[[17]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_i_opt_emax[i] <-  ifelse(any((res_lin[[21]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_ff_interim[i] <-  ifelse(any((res_lin[[25]] / 5000) > GoNoGo_Threshold), 1, 0)

      success_emax_ff[i] <-  ifelse(any((res_ema[[1]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_comp[i] <-  ifelse(any((res_ema[[2]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_d_opt_linear[i] <-  ifelse(any((res_ema[[9]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_d_opt_emax[i] <-  ifelse(any((res_ema[[13]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_i_opt_linear[i] <-  ifelse(any((res_ema[[17]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_i_opt_emax[i] <-  ifelse(any((res_ema[[21]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_ff_interim[i] <-  ifelse(any((res_ema[[25]] / 5000) > GoNoGo_Threshold), 1, 0)
    }

    success_linear_ff <- sum(success_linear_ff) / length(results[[1]]$linear)
    success_linear_comp <-  sum(success_linear_comp) / length(results[[1]]$linear)
    success_linear_d_opt_linear <-  sum(success_linear_d_opt_linear) / length(results[[1]]$linear)
    success_linear_d_opt_emax <-  sum(success_linear_d_opt_emax) / length(results[[1]]$linear)
    success_linear_i_opt_linear <-  sum(success_linear_i_opt_linear) / length(results[[1]]$linear)
    success_linear_i_opt_emax <-  sum(success_linear_i_opt_emax) / length(results[[1]]$linear)
    success_linear_ff_interim <-  sum(success_linear_ff_interim) / length(results[[1]]$linear)

    success_emax_ff <-  sum(success_emax_ff) / length(results[[1]]$linear)
    success_emax_comp <-  sum(success_emax_comp) / length(results[[1]]$linear)
    success_emax_d_opt_linear <-  sum(success_emax_d_opt_linear) / length(results[[1]]$linear)
    success_emax_d_opt_emax <-  sum(success_emax_d_opt_emax) / length(results[[1]]$linear)
    success_emax_i_opt_linear <-  sum(success_emax_i_opt_linear) / length(results[[1]]$linear)
    success_emax_i_opt_emax <-  sum(success_emax_i_opt_emax) / length(results[[1]]$linear)
    success_emax_ff_interim <-  sum(success_emax_ff_interim) / length(results[[1]]$linear)

    res <- list(
      success_linear_ff, success_linear_comp, success_linear_d_opt_linear, success_linear_d_opt_emax, success_linear_i_opt_linear, success_linear_i_opt_emax,
      success_linear_ff_interim,
      success_emax_ff, success_emax_comp, success_emax_d_opt_linear, success_emax_d_opt_emax, success_emax_i_opt_linear, success_emax_i_opt_emax,
      success_emax_ff_interim
    )
  } else {
    success_linear_ff <- numeric(length(results[[1]]$linear))
    success_linear_comp <- numeric(length(results[[1]]$linear))
    success_linear_d_opt_linear <- numeric(length(results[[1]]$linear))
    success_linear_d_opt_emax <- numeric(length(results[[1]]$linear))
    success_linear_i_opt_linear <- numeric(length(results[[1]]$linear))
    success_linear_i_opt_emax <- numeric(length(results[[1]]$linear))


    success_emax_ff <- numeric(length(results[[1]]$linear))
    success_emax_comp <- numeric(length(results[[1]]$linear))
    success_emax_d_opt_linear <- numeric(length(results[[1]]$linear))
    success_emax_d_opt_emax <- numeric(length(results[[1]]$linear))
    success_emax_i_opt_linear <- numeric(length(results[[1]]$linear))
    success_emax_i_opt_emax <- numeric(length(results[[1]]$linear))


    for(i in 1:length(results[[1]]$linear)) {
      res_lin <- results[[scenario]]$linear[[i]]
      res_ema <- results[[scenario]]$emax[[i]]

      success_linear_ff[i] <-  ifelse(any((res_lin[[1]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_comp[i] <-  ifelse(any((res_lin[[2]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_d_opt_linear[i] <-  ifelse(any((res_lin[[9]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_d_opt_emax[i] <-  ifelse(any((res_lin[[13]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_i_opt_linear[i] <-  ifelse(any((res_lin[[17]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_linear_i_opt_emax[i] <-  ifelse(any((res_lin[[21]] / 5000) > GoNoGo_Threshold), 1, 0)


      success_emax_ff[i] <-  ifelse(any((res_ema[[1]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_comp[i] <-  ifelse(any((res_ema[[2]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_d_opt_linear[i] <-  ifelse(any((res_ema[[9]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_d_opt_emax[i] <-  ifelse(any((res_ema[[13]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_i_opt_linear[i] <-  ifelse(any((res_ema[[17]] / 5000) > GoNoGo_Threshold), 1, 0)
      success_emax_i_opt_emax[i] <-  ifelse(any((res_ema[[21]] / 5000) > GoNoGo_Threshold), 1, 0)

    }

    success_linear_ff <- sum(success_linear_ff) / length(results[[1]]$linear)
    success_linear_comp <-  sum(success_linear_comp) / length(results[[1]]$linear)
    success_linear_d_opt_linear <-  sum(success_linear_d_opt_linear) / length(results[[1]]$linear)
    success_linear_d_opt_emax <-  sum(success_linear_d_opt_emax) / length(results[[1]]$linear)
    success_linear_i_opt_linear <-  sum(success_linear_i_opt_linear) / length(results[[1]]$linear)
    success_linear_i_opt_emax <-  sum(success_linear_i_opt_emax) / length(results[[1]]$linear)


    success_emax_ff <-  sum(success_emax_ff) / length(results[[1]]$linear)
    success_emax_comp <-  sum(success_emax_comp) / length(results[[1]]$linear)
    success_emax_d_opt_linear <-  sum(success_emax_d_opt_linear) / length(results[[1]]$linear)
    success_emax_d_opt_emax <-  sum(success_emax_d_opt_emax) / length(results[[1]]$linear)
    success_emax_i_opt_linear <-  sum(success_emax_i_opt_linear) / length(results[[1]]$linear)
    success_emax_i_opt_emax <-  sum(success_emax_i_opt_emax) / length(results[[1]]$linear)


    res <- list(
      success_linear_ff, success_linear_comp, success_linear_d_opt_linear, success_linear_d_opt_emax, success_linear_i_opt_linear, success_linear_i_opt_emax,
      success_emax_ff, success_emax_comp, success_emax_d_opt_linear, success_emax_d_opt_emax, success_emax_i_opt_linear, success_emax_i_opt_emax
    )
  }


  return(res)
}

# AIC ---------------------------------------------------------------------


#' @title get_ABIC_score
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param abic  location of AIC/BIC values
#' @param med  location of MED
#' @param method  for which option sum should get calulated. Options are:
#' "gonogo", "corr_med", "mse", "rmse", "pat_alloc"
#' @param GoNoGo_Threshold  gonogo threshold to be passed, e.g., 0.5
#'
#' @return sum
#' @export
get_ABIC_score <- function(results, scenario, abic, med, method, GoNoGo_Threshold) {
  mat <- results[[scenario]]$alloc_mat[[1]]
  interim <- results[[scenario]]$interim[[1]]
  likelihood <- results[[scenario]]$likelihood[[1]]
  sum <- numeric(length(results[[1]]$linear))
  if(likelihood == "linear") {
    other <- "emax"
  } else {
    other <- "linear"
  }
  sel_model_vec <- character(length(results[[scenario]][[likelihood]]))
  for(i in 1:length(results[[scenario]][[likelihood]])) {
    model_sel <- results[[scenario]][[likelihood]][[i]][[abic]] < results[[scenario]][[other]][[i]][[abic]]
    if(model_sel) {
      sel_model_vec[i] <- likelihood
    } else {
      sel_model_vec[i] <- other
    }
  }


  if(method == "gonogo") {
    for(i in 1:length(results[[scenario]][[likelihood]])) {
      sel_model <- sel_model_vec[i]
      sum[i] <- ifelse(any((results[[scenario]][[sel_model]][[i]][[med]] / 5000) > GoNoGo_Threshold), 1, 0)
    }
    sum_1 <- sum(sum) / length(results[[1]]$linear)
  } else if(method == "corr_med") {
    for(i in 1:length(results[[scenario]][[likelihood]])) {
      sel_model <- sel_model_vec[i]
      true_med <- results[[scenario]][["true_med"]][[i]]
      true_med_pos <- which(true_med > 0)[1]
      if(i == 1) {
        sum <- as.numeric(results[[scenario]][[sel_model]][[i]][[med]][true_med_pos] > 0)
      } else {
        sum <- sum + as.numeric(results[[scenario]][[sel_model]][[i]][[med]][true_med_pos] > 0)
      }
    }
    sum_1 <- sum / length(results[[scenario]][[likelihood]])
  } else if(method == "mse") {
    for(i in 1:length(results[[scenario]][[likelihood]])) {
      sel_model <- sel_model_vec[i]
      true_med <- results[[scenario]][["true_med"]][[i]]
      true_med_pos <- which(true_med > 0)[1]
      if(med == 1 | med == 2) {
        med_effect <- med + 2
      } else {
        med_effect <- med + 1
      }

      pos <- which(results[[scenario]][[likelihood]][[i]][[med_effect]] / 5000 > 8)[1]
      true_med_effect <- results[[scenario]][["true_effect"]][[i]][pos]
      obs_effect <- results[[scenario]][[likelihood]][[i]][[med_effect]][pos] / 5000
      sum[i] <- (obs_effect - true_med_effect)^2

    }
    sum_1 <- mean(sum)
  } else if(method == "rmse") {
    for(i in 1:length(results[[scenario]][[likelihood]])) {
      sel_model <- sel_model_vec[i]
      true_med <- results[[scenario]][["true_med"]][[i]]
      true_med_pos <- which(true_med > 0)[1]
      if(med == 1 | med == 2) {
        med_effect <- med + 2
      } else {
        med_effect <- med + 1
      }
      pos <- which(results[[scenario]][[likelihood]][[i]][[med_effect]] / 5000 > 8)[1]
      true_med_effect <- results[[scenario]][["true_effect"]][[i]][pos]
      obs_effect <- results[[scenario]][[likelihood]][[i]][[med_effect]][pos] / 5000
      sum[i] <- (obs_effect - true_med_effect)^2

    }
    sum_1 <- sqrt(mean(sum))
  } else if(method == "pat_alloc") {
    if(interim) {
      ff_alloc <- generate_full_factorial_design(
        doses = results[[scenario]]$doses[[1]],
        schedules = results[[scenario]]$schedules[[1]],
        n_pat = results[[scenario]]$n_pat[[1]] * 0.5
      )
      alpha1 <- results[[scenario]]$alpha1[[1]]
      alpha2 <- results[[scenario]]$alpha2[[1]]
      delta1 <- results[[scenario]]$delta1[[1]]
      delta2 <- results[[scenario]]$delta2[[1]]
      n_pat <- results[[scenario]]$n_pat[[1]] * 0.5
      d_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules, n_pat = n_pat)$mat
      d_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,
                                                        emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat
      i_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I")$mat
      i_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I",
                                                        emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat



      true_med <- results[[scenario]][["true_med"]][[1]]
      pos_true_med <- which(true_med > 0)[1]
      sum_mat1_d_interim = d_opt_design_linear[pos_true_med]
      sum_mat2_d_interim = d_opt_design_emax[pos_true_med]
      sum_mat1_i_interim = i_opt_design_linear[pos_true_med]
      sum_mat2_i_interim = i_opt_design_emax[pos_true_med]
      sum <- numeric(length(results[[scenario]][["linear"]]))

      for(i in 1:length(results[[scenario]][[likelihood]])) {
        sel_model <- sel_model_vec[i]
        if(sel_model == "linear") {

          if(med == 1) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[1]][[pos_true_med]]
            }

          } else if (med == 2) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[2]][[pos_true_med]] +  results[[scenario]]$mats[[i]][[3]][[pos_true_med]]
            }

          }else if (med == 9) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[5]][[pos_true_med]] + sum_mat1_d_interim
            }

          }else if (med == 13) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[6]][[pos_true_med]] + sum_mat1_d_interim
            }

          }else if (med == 17) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[7]][[pos_true_med]] + sum_mat1_i_interim
            }

          }else if (med == 21) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[8]][[pos_true_med]] + sum_mat1_i_interim
            }

          }else if (med == 25) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[13]][[pos_true_med]] + ff_alloc[pos_true_med]
            }

          }
        } else {
          if(med == 1) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[1]][[pos_true_med]]
            }

          } else if (med == 2) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[2]][[pos_true_med]] +  results[[scenario]]$mats[[i]][[4]][[pos_true_med]]
            }

          }else if (med == 9) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[9]][[pos_true_med]] + sum_mat2_d_interim
            }

          }else if (med == 13) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[10]][[pos_true_med]] + sum_mat2_d_interim
            }

          }else if (med == 17) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[11]][[pos_true_med]] + sum_mat2_i_interim
            }

          }else if (med == 21) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[12]][[pos_true_med]] + sum_mat2_i_interim
            }

          }else if (med == 25) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[14]][[pos_true_med]] + ff_alloc[pos_true_med]
            }

          }
        }

      }
    } else {
      true_med <- results[[scenario]][["true_med"]][[1]]
      pos_true_med <- which(true_med > 0)[1]

      alpha1 <- results[[scenario]]$alpha1[[1]]
      alpha2 <- results[[scenario]]$alpha2[[1]]
      delta1 <- results[[scenario]]$delta1[[1]]
      delta2 <- results[[scenario]]$delta2[[1]]
      n_pat <- results[[scenario]]$n_pat[[1]]
      d_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules, n_pat = n_pat)$mat
      d_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,
                                                        emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat
      i_opt_design_linear <- generate_optimal_design_linear(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I")$mat
      i_opt_design_emax <- generate_optimal_design_emax(doses = doses, schedules = schedules,n_pat = n_pat,criterion = "I",
                                                        emax_d = alpha1, ed50_d = delta1, emax_s = alpha2, ed50_s = delta2)$mat


      sum_mat1_d = d_opt_design_linear[pos_true_med]
      sum_mat2_d = d_opt_design_emax[pos_true_med]
      sum_mat1_i = i_opt_design_linear[pos_true_med]
      sum_mat2_i = i_opt_design_emax[pos_true_med]
      sum <- numeric(length(results[[scenario]][["linear"]]))

      for(i in 1:length(results[[scenario]][[likelihood]])) {
        sel_model <- sel_model_vec[i]
        if(sel_model == "linear") {
          if(med == 1) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[1]][[pos_true_med]]
            }

          } else if (med == 2) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[2]][[pos_true_med]]
            }

          }else if (med == 9) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat1_d
            }

          }else if (med == 13) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat1_d
            }

          }else if (med == 17) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat1_i
            }

          }else if (med == 21) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat1_i
            }

          }
        } else {
          if(med == 1) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[1]][[pos_true_med]]
            }

          } else if (med == 2) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- results[[scenario]]$mats[[i]][[2]][[pos_true_med]]
            }

          }else if (med == 9) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat2_d
            }

          }else if (med == 13) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat2_d
            }

          }else if (med == 17) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat2_i
            }

          }else if (med == 21) {
            if(is.na(pos_true_med)) {
              sum[i] <- NA
            }else{
              sum[i] <- sum_mat2_i
            }

          }
        }


      }
    }

    sum_1 <- mean(sum)
  }




  return(sum_1)
}

#' @title get_AIC_BIC_scores_across
#'
#' @param results  results dataset
#' @param scenario numeric value of scenario
#' @param interim  boolean, indicating whether an interim analysis has been conducted
#' @param gonogo_threshold  gonogo threshold to be passed, e.g., 0.5
#'
#' @return dataframe
#' @export
get_AIC_BIC_scores_across <- function(results, scenario, interim, gonogo_threshold) {
  if(interim) {
    med_1_A_gonogo <- get_ABIC_score(results, scenario, 5, 1,  "gonogo", gonogo_threshold)
    med_1_A_corr_med <- get_ABIC_score(results, scenario, 5, 1,  "corr_med", gonogo_threshold)
    med_1_A_mse <- get_ABIC_score(results, scenario, 5, 1,  "mse", gonogo_threshold)
    med_1_A_rmse <- get_ABIC_score(results, scenario, 5, 1,  "rmse", gonogo_threshold)
    med_1_A_pat_alloc<- get_ABIC_score(results, scenario, 5, 1,  "pat_alloc", gonogo_threshold)
    med_1_B_gonogo <- get_ABIC_score(results, scenario, 7, 1,  "gonogo", gonogo_threshold)
    med_1_B_corr_med <- get_ABIC_score(results, scenario, 7, 1,  "corr_med", gonogo_threshold)
    med_1_B_mse <- get_ABIC_score(results, scenario, 7, 1,  "mse", gonogo_threshold)
    med_1_B_rmse <- get_ABIC_score(results, scenario, 7, 1,  "rmse", gonogo_threshold)
    med_1_B_pat_alloc <- get_ABIC_score(results, scenario, 7, 1,  "pat_alloc", gonogo_threshold)

    med_2_A_gonogo <- get_ABIC_score(results, scenario, 6, 2,  "gonogo", gonogo_threshold)
    med_2_A_corr_med <- get_ABIC_score(results, scenario, 6, 2,  "corr_med", gonogo_threshold)
    med_2_A_mse <- get_ABIC_score(results, scenario, 6, 2,  "mse", gonogo_threshold)
    med_2_A_rmse <- get_ABIC_score(results, scenario, 6, 2,  "rmse", gonogo_threshold)
    med_2_A_pat_alloc<- get_ABIC_score(results, scenario, 6, 2,  "pat_alloc", gonogo_threshold)
    med_2_B_gonogo <- get_ABIC_score(results, scenario, 8, 2,  "gonogo", gonogo_threshold)
    med_2_B_corr_med <- get_ABIC_score(results, scenario, 8, 2,  "corr_med", gonogo_threshold)
    med_2_B_mse <- get_ABIC_score(results, scenario, 8, 2,  "mse", gonogo_threshold)
    med_2_B_rmse <- get_ABIC_score(results, scenario, 8, 2,  "rmse", gonogo_threshold)
    med_2_B_pat_alloc <- get_ABIC_score(results, scenario, 8, 2,  "pat_alloc", gonogo_threshold)

    med_3_A_gonogo <- get_ABIC_score(results, scenario, 11, 9,  "gonogo", gonogo_threshold)
    med_3_A_corr_med <- get_ABIC_score(results, scenario, 11, 9,  "corr_med", gonogo_threshold)
    med_3_A_mse <- get_ABIC_score(results, scenario, 11, 9,  "mse", gonogo_threshold)
    med_3_A_rmse <- get_ABIC_score(results, scenario, 11, 9,  "rmse", gonogo_threshold)
    med_3_A_pat_alloc<- get_ABIC_score(results, scenario, 11, 9,  "pat_alloc", gonogo_threshold)
    med_3_B_gonogo <- get_ABIC_score(results, scenario, 12, 9,  "gonogo", gonogo_threshold)
    med_3_B_corr_med <- get_ABIC_score(results, scenario, 12, 9,  "corr_med", gonogo_threshold)
    med_3_B_mse <- get_ABIC_score(results, scenario, 12, 9,  "mse", gonogo_threshold)
    med_3_B_rmse <- get_ABIC_score(results, scenario, 12, 9,  "rmse", gonogo_threshold)
    med_3_B_pat_alloc <- get_ABIC_score(results, scenario, 12, 9,  "pat_alloc", gonogo_threshold)

    med_4_A_gonogo <- get_ABIC_score(results, scenario, 15, 13,  "gonogo", gonogo_threshold)
    med_4_A_corr_med <- get_ABIC_score(results, scenario, 15, 13,  "corr_med", gonogo_threshold)
    med_4_A_mse <- get_ABIC_score(results, scenario, 15, 13,  "mse", gonogo_threshold)
    med_4_A_rmse <- get_ABIC_score(results, scenario, 15, 13,  "rmse", gonogo_threshold)
    med_4_A_pat_alloc<- get_ABIC_score(results, scenario, 15, 13,  "pat_alloc", gonogo_threshold)
    med_4_B_gonogo <- get_ABIC_score(results, scenario, 16, 13,  "gonogo", gonogo_threshold)
    med_4_B_corr_med <- get_ABIC_score(results, scenario, 16, 13,  "corr_med", gonogo_threshold)
    med_4_B_mse <- get_ABIC_score(results, scenario, 16, 13,  "mse", gonogo_threshold)
    med_4_B_rmse <- get_ABIC_score(results, scenario, 16, 13,  "rmse", gonogo_threshold)
    med_4_B_pat_alloc <- get_ABIC_score(results, scenario, 16, 13,  "pat_alloc", gonogo_threshold)

    med_5_A_gonogo <- get_ABIC_score(results, scenario, 19, 17,  "gonogo", gonogo_threshold)
    med_5_A_corr_med <- get_ABIC_score(results, scenario, 19, 17,  "corr_med", gonogo_threshold)
    med_5_A_mse <- get_ABIC_score(results, scenario, 19, 17,  "mse", gonogo_threshold)
    med_5_A_rmse <- get_ABIC_score(results, scenario, 19, 17,  "rmse", gonogo_threshold)
    med_5_A_pat_alloc<- get_ABIC_score(results, scenario, 19, 17,  "pat_alloc", gonogo_threshold)
    med_5_B_gonogo <- get_ABIC_score(results, scenario, 20, 17,  "gonogo", gonogo_threshold)
    med_5_B_corr_med <- get_ABIC_score(results, scenario, 20, 17,  "corr_med", gonogo_threshold)
    med_5_B_mse <- get_ABIC_score(results, scenario, 20, 17,  "mse", gonogo_threshold)
    med_5_B_rmse <- get_ABIC_score(results, scenario, 20, 17,  "rmse", gonogo_threshold)
    med_5_B_pat_alloc <- get_ABIC_score(results, scenario, 20, 17,  "pat_alloc", gonogo_threshold)

    med_6_A_gonogo <- get_ABIC_score(results, scenario, 23, 21,  "gonogo", gonogo_threshold)
    med_6_A_corr_med <- get_ABIC_score(results, scenario, 23, 21,  "corr_med", gonogo_threshold)
    med_6_A_mse <- get_ABIC_score(results, scenario, 23, 21,  "mse", gonogo_threshold)
    med_6_A_rmse <- get_ABIC_score(results, scenario, 23, 21,  "rmse", gonogo_threshold)
    med_6_A_pat_alloc<- get_ABIC_score(results, scenario, 23, 21,  "pat_alloc", gonogo_threshold)
    med_6_B_gonogo <- get_ABIC_score(results, scenario, 24, 21,  "gonogo", gonogo_threshold)
    med_6_B_corr_med <- get_ABIC_score(results, scenario, 24, 21,  "corr_med", gonogo_threshold)
    med_6_B_mse <- get_ABIC_score(results, scenario, 24, 21,  "mse", gonogo_threshold)
    med_6_B_rmse <- get_ABIC_score(results, scenario, 24, 21,  "rmse", gonogo_threshold)
    med_6_B_pat_alloc <- get_ABIC_score(results, scenario, 24, 21,  "pat_alloc", gonogo_threshold)

    med_7_A_gonogo <- get_ABIC_score(results, scenario, 27, 25,  "gonogo", gonogo_threshold)
    med_7_A_corr_med <- get_ABIC_score(results, scenario, 27, 25,  "corr_med", gonogo_threshold)
    med_7_A_mse <- get_ABIC_score(results, scenario, 27, 25,  "mse", gonogo_threshold)
    med_7_A_rmse <- get_ABIC_score(results, scenario, 27, 25,  "rmse", gonogo_threshold)
    med_7_A_pat_alloc<- get_ABIC_score(results, scenario, 27, 25,  "pat_alloc", gonogo_threshold)
    med_7_B_gonogo <- get_ABIC_score(results, scenario, 28, 25,  "gonogo", gonogo_threshold)
    med_7_B_corr_med <- get_ABIC_score(results, scenario, 28, 25,  "corr_med", gonogo_threshold)
    med_7_B_mse <- get_ABIC_score(results, scenario, 28, 25,  "mse", gonogo_threshold)
    med_7_B_rmse <- get_ABIC_score(results, scenario, 28, 25,  "rmse", gonogo_threshold)
    med_7_B_pat_alloc <- get_ABIC_score(results, scenario, 28, 25,  "pat_alloc", gonogo_threshold)


    df <- data.frame(
      AIC = c(
        med_1_A_gonogo,
        med_1_A_corr_med,
        med_1_A_mse,
        med_1_A_rmse,
        med_1_A_pat_alloc,
        med_2_A_gonogo,
        med_2_A_corr_med,
        med_2_A_mse,
        med_2_A_rmse,
        med_2_A_pat_alloc,
        med_3_A_gonogo,
        med_3_A_corr_med,
        med_3_A_mse,
        med_3_A_rmse,
        med_3_A_pat_alloc,
        med_4_A_gonogo,
        med_4_A_corr_med,
        med_4_A_mse,
        med_4_A_rmse,
        med_4_A_pat_alloc,
        med_5_A_gonogo,
        med_5_A_corr_med,
        med_5_A_mse,
        med_5_A_rmse,
        med_5_A_pat_alloc,
        med_6_A_gonogo,
        med_6_A_corr_med,
        med_6_A_mse,
        med_6_A_rmse,
        med_6_A_pat_alloc,
        med_7_A_gonogo,
        med_7_A_corr_med,
        med_7_A_mse,
        med_7_A_rmse,
        med_7_A_pat_alloc
      ),
      BIC = c(
        med_1_B_gonogo,
        med_1_B_corr_med,
        med_1_B_mse,
        med_1_B_rmse,
        med_1_B_pat_alloc,
        med_2_B_gonogo,
        med_2_B_corr_med,
        med_2_B_mse,
        med_2_B_rmse,
        med_2_B_pat_alloc,
        med_3_B_gonogo,
        med_3_B_corr_med,
        med_3_B_mse,
        med_3_B_rmse,
        med_3_B_pat_alloc,
        med_4_B_gonogo,
        med_4_B_corr_med,
        med_4_B_mse,
        med_4_B_rmse,
        med_4_B_pat_alloc,
        med_5_B_gonogo,
        med_5_B_corr_med,
        med_5_B_mse,
        med_5_B_rmse,
        med_5_B_pat_alloc,
        med_6_B_gonogo,
        med_6_B_corr_med,
        med_6_B_mse,
        med_6_B_rmse,
        med_6_B_pat_alloc,
        med_7_B_gonogo,
        med_7_B_corr_med,
        med_7_B_mse,
        med_7_B_rmse,
        med_7_B_pat_alloc
      )
    )
  } else {
    med_1_A_gonogo <- get_ABIC_score(results, scenario, 5, 1,  "gonogo", gonogo_threshold)
    med_1_A_corr_med <- get_ABIC_score(results, scenario, 5, 1,  "corr_med", gonogo_threshold)
    med_1_A_mse <- get_ABIC_score(results, scenario, 5, 1,  "mse", gonogo_threshold)
    med_1_A_rmse <- get_ABIC_score(results, scenario, 5, 1,  "rmse", gonogo_threshold)
    med_1_A_pat_alloc<- get_ABIC_score(results, scenario, 5, 1,  "pat_alloc", gonogo_threshold)
    med_1_B_gonogo <- get_ABIC_score(results, scenario, 7, 1,  "gonogo", gonogo_threshold)
    med_1_B_corr_med <- get_ABIC_score(results, scenario, 7, 1,  "corr_med", gonogo_threshold)
    med_1_B_mse <- get_ABIC_score(results, scenario, 7, 1,  "mse", gonogo_threshold)
    med_1_B_rmse <- get_ABIC_score(results, scenario, 7, 1,  "rmse", gonogo_threshold)
    med_1_B_pat_alloc <- get_ABIC_score(results, scenario, 7, 1,  "pat_alloc", gonogo_threshold)

    med_2_A_gonogo <- get_ABIC_score(results, scenario, 6, 2,  "gonogo", gonogo_threshold)
    med_2_A_corr_med <- get_ABIC_score(results, scenario, 6, 2,  "corr_med", gonogo_threshold)
    med_2_A_mse <- get_ABIC_score(results, scenario, 6, 2,  "mse", gonogo_threshold)
    med_2_A_rmse <- get_ABIC_score(results, scenario, 6, 2,  "rmse", gonogo_threshold)
    med_2_A_pat_alloc<- get_ABIC_score(results, scenario, 6, 2,  "pat_alloc", gonogo_threshold)
    med_2_B_gonogo <- get_ABIC_score(results, scenario, 8, 2,  "gonogo", gonogo_threshold)
    med_2_B_corr_med <- get_ABIC_score(results, scenario, 8, 2,  "corr_med", gonogo_threshold)
    med_2_B_mse <- get_ABIC_score(results, scenario, 8, 2,  "mse", gonogo_threshold)
    med_2_B_rmse <- get_ABIC_score(results, scenario, 8, 2,  "rmse", gonogo_threshold)
    med_2_B_pat_alloc <- get_ABIC_score(results, scenario, 8, 2,  "pat_alloc", gonogo_threshold)

    med_3_A_gonogo <- get_ABIC_score(results, scenario, 11, 9,  "gonogo", gonogo_threshold)
    med_3_A_corr_med <- get_ABIC_score(results, scenario, 11, 9,  "corr_med", gonogo_threshold)
    med_3_A_mse <- get_ABIC_score(results, scenario, 11, 9,  "mse", gonogo_threshold)
    med_3_A_rmse <- get_ABIC_score(results, scenario, 11, 9,  "rmse", gonogo_threshold)
    med_3_A_pat_alloc<- get_ABIC_score(results, scenario, 11, 9,  "pat_alloc", gonogo_threshold)
    med_3_B_gonogo <- get_ABIC_score(results, scenario, 12, 9,  "gonogo", gonogo_threshold)
    med_3_B_corr_med <- get_ABIC_score(results, scenario, 12, 9,  "corr_med", gonogo_threshold)
    med_3_B_mse <- get_ABIC_score(results, scenario, 12, 9,  "mse", gonogo_threshold)
    med_3_B_rmse <- get_ABIC_score(results, scenario, 12, 9,  "rmse", gonogo_threshold)
    med_3_B_pat_alloc <- get_ABIC_score(results, scenario, 12, 9,  "pat_alloc", gonogo_threshold)

    med_4_A_gonogo <- get_ABIC_score(results, scenario, 15, 13,  "gonogo", gonogo_threshold)
    med_4_A_corr_med <- get_ABIC_score(results, scenario, 15, 13,  "corr_med", gonogo_threshold)
    med_4_A_mse <- get_ABIC_score(results, scenario, 15, 13,  "mse", gonogo_threshold)
    med_4_A_rmse <- get_ABIC_score(results, scenario, 15, 13,  "rmse", gonogo_threshold)
    med_4_A_pat_alloc<- get_ABIC_score(results, scenario, 15, 13,  "pat_alloc", gonogo_threshold)
    med_4_B_gonogo <- get_ABIC_score(results, scenario, 16, 13,  "gonogo", gonogo_threshold)
    med_4_B_corr_med <- get_ABIC_score(results, scenario, 16, 13,  "corr_med", gonogo_threshold)
    med_4_B_mse <- get_ABIC_score(results, scenario, 16, 13,  "mse", gonogo_threshold)
    med_4_B_rmse <- get_ABIC_score(results, scenario, 16, 13,  "rmse", gonogo_threshold)
    med_4_B_pat_alloc <- get_ABIC_score(results, scenario, 16, 13,  "pat_alloc", gonogo_threshold)

    med_5_A_gonogo <- get_ABIC_score(results, scenario, 19, 17,  "gonogo", gonogo_threshold)
    med_5_A_corr_med <- get_ABIC_score(results, scenario, 19, 17,  "corr_med", gonogo_threshold)
    med_5_A_mse <- get_ABIC_score(results, scenario, 19, 17,  "mse", gonogo_threshold)
    med_5_A_rmse <- get_ABIC_score(results, scenario, 19, 17,  "rmse", gonogo_threshold)
    med_5_A_pat_alloc<- get_ABIC_score(results, scenario, 19, 17,  "pat_alloc", gonogo_threshold)
    med_5_B_gonogo <- get_ABIC_score(results, scenario, 20, 17,  "gonogo", gonogo_threshold)
    med_5_B_corr_med <- get_ABIC_score(results, scenario, 20, 17,  "corr_med", gonogo_threshold)
    med_5_B_mse <- get_ABIC_score(results, scenario, 20, 17,  "mse", gonogo_threshold)
    med_5_B_rmse <- get_ABIC_score(results, scenario, 20, 17,  "rmse", gonogo_threshold)
    med_5_B_pat_alloc <- get_ABIC_score(results, scenario, 20, 17,  "pat_alloc", gonogo_threshold)

    med_6_A_gonogo <- get_ABIC_score(results, scenario, 23, 21,  "gonogo", gonogo_threshold)
    med_6_A_corr_med <- get_ABIC_score(results, scenario, 23, 21,  "corr_med", gonogo_threshold)
    med_6_A_mse <- get_ABIC_score(results, scenario, 23, 21,  "mse", gonogo_threshold)
    med_6_A_rmse <- get_ABIC_score(results, scenario, 23, 21,  "rmse", gonogo_threshold)
    med_6_A_pat_alloc<- get_ABIC_score(results, scenario, 23, 21,  "pat_alloc", gonogo_threshold)
    med_6_B_gonogo <- get_ABIC_score(results, scenario, 24, 21,  "gonogo", gonogo_threshold)
    med_6_B_corr_med <- get_ABIC_score(results, scenario, 24, 21,  "corr_med", gonogo_threshold)
    med_6_B_mse <- get_ABIC_score(results, scenario, 24, 21,  "mse", gonogo_threshold)
    med_6_B_rmse <- get_ABIC_score(results, scenario, 24, 21,  "rmse", gonogo_threshold)
    med_6_B_pat_alloc <- get_ABIC_score(results, scenario, 24, 21,  "pat_alloc", gonogo_threshold)

    df <- data.frame(
      AIC = c(
        med_1_A_gonogo,
        med_1_A_corr_med,
        med_1_A_mse,
        med_1_A_rmse,
        med_1_A_pat_alloc,
        med_2_A_gonogo,
        med_2_A_corr_med,
        med_2_A_mse,
        med_2_A_rmse,
        med_2_A_pat_alloc,
        med_3_A_gonogo,
        med_3_A_corr_med,
        med_3_A_mse,
        med_3_A_rmse,
        med_3_A_pat_alloc,
        med_4_A_gonogo,
        med_4_A_corr_med,
        med_4_A_mse,
        med_4_A_rmse,
        med_4_A_pat_alloc,
        med_5_A_gonogo,
        med_5_A_corr_med,
        med_5_A_mse,
        med_5_A_rmse,
        med_5_A_pat_alloc,
        med_6_A_gonogo,
        med_6_A_corr_med,
        med_6_A_mse,
        med_6_A_rmse,
        med_6_A_pat_alloc
      ),
      BIC = c(
        med_1_B_gonogo,
        med_1_B_corr_med,
        med_1_B_mse,
        med_1_B_rmse,
        med_1_B_pat_alloc,
        med_2_B_gonogo,
        med_2_B_corr_med,
        med_2_B_mse,
        med_2_B_rmse,
        med_2_B_pat_alloc,
        med_3_B_gonogo,
        med_3_B_corr_med,
        med_3_B_mse,
        med_3_B_rmse,
        med_3_B_pat_alloc,
        med_4_B_gonogo,
        med_4_B_corr_med,
        med_4_B_mse,
        med_4_B_rmse,
        med_4_B_pat_alloc,
        med_5_B_gonogo,
        med_5_B_corr_med,
        med_5_B_mse,
        med_5_B_rmse,
        med_5_B_pat_alloc,
        med_6_B_gonogo,
        med_6_B_corr_med,
        med_6_B_mse,
        med_6_B_rmse,
        med_6_B_pat_alloc
      )
    )
  }


  return(df)
}

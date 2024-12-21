
#  Table of parameters for targets
scenarios_table <- tibble::tribble(
  ~scenario_label, ~doses, ~schedules, ~threshold, ~alloc_mat, ~n_pat, ~interim, ~likelihood, ~E0, ~alpha1, ~alpha2, ~delta1, ~delta2, ~beta, ~sig_mu, ~sig_sd, ~beta0, ~beta1, ~beta2, ~beta3,
  "viable_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "viable_no_interim_interaction_0_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_no_interim_interaction_0_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_interim_interaction_0_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_interim_interaction_0_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,

  "viable_no_interim_full_interaction", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_no_interim_full_interaction", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_interim_full_interaction", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_interim_full_interaction", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,

  "viable_no_interim_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_no_interim_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "viable_no_interim_interaction_0_5_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_no_interim_interaction_0_5_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_interim_interaction_0_5_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,
  "viable_interim_interaction_0_5_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0.5, 1, 1, -12.8, 0.9, 0.9, 0.5,

  "viable_no_interim_full_interaction_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_no_interim_full_interaction_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_interim_full_interaction_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,
  "viable_interim_full_interaction_equal", c(0, 5, 9), c(1, 3, 4), 8, 2, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 1, 1, 1, -12.8, 0.9, 0.9, 1,

  "viable_no_interim_low", c(0, 5, 9), c(1, 3, 4), 8, 3, 50, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_no_interim_low", c(0, 5, 9), c(1, 3, 4), 8, 3, 50, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_low", c(0, 5, 9), c(1, 3, 4), 8, 3, 50, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_low", c(0, 5, 9), c(1, 3, 4), 8, 3, 50, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "viable_no_interim_high", c(0, 5, 9), c(1, 3, 4), 8, 4, 400, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_no_interim_high", c(0, 5, 9), c(1, 3, 4), 8, 4, 400, FALSE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_high", c(0, 5, 9), c(1, 3, 4), 8, 4, 400, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "viable_interim_high", c(0, 5, 9), c(1, 3, 4), 8, 4, 400, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "null_equal", c(0, 5, 9), c(1, 3, 4), 8, 5, 200, TRUE, "linear",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,
  "null_equal", c(0, 5, 9), c(1, 3, 4), 8, 5, 200, TRUE, "emax",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,

  "null", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 6, 200, TRUE, "linear",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,
  "null", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 6, 200, TRUE, "emax",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,

  "null_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 6, 200, FALSE, "linear",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,
  "null_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 6, 200, FALSE, "emax",  -12.8, 0, 0, 3, 1, 0, 1, 1, -12.8, 0, 0, 0,

  "positive_equal", c(0, 5, 9), c(1, 3, 4), 8, 7, 200, TRUE, "linear",  -12.8, 15, 15, 3, 1, 0, 1, 1, -12.8, 1.35, 1.35, 0,
  "positive_equal", c(0, 5, 9), c(1, 3, 4), 8, 7, 200, TRUE, "emax",  -12.8, 15, 15, 3, 1, 0, 1, 1, -12.8, 1.35, 1.35, 0,

  "positive", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 8, 200, TRUE, "linear",  -12.8, 15, 15, 3, 1, 0, 1, 1, -12.8, 1.35, 1.35, 0,
  "positive", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 8, 200, TRUE, "emax",  -12.8, 15, 15, 3, 1, 0, 1, 1, -12.8, 1.35, 1.35, 0,

  "slight_positive", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 9, 200, TRUE, "linear",  -12.8, 12, 12, 3, 1, 0, 1, 1, -12.8, 1.25, 1.25, 0,
  "slight_positive", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 9, 200, TRUE, "emax",  -12.8, 12, 12, 3, 1, 0, 1, 1, -12.8, 1.25, 1.25, 0,

  "slight_positive_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 9, 200, FALSE, "linear",  -12.8, 12, 12, 3, 1, 0, 1, 1, -12.8, 1.25, 1.25, 0,
  "slight_positive_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 9, 200, FALSE, "emax",  -12.8, 12, 12, 3, 1, 0, 1, 1, -12.8, 1.25, 1.25, 0,

  "negative", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 10, 200, TRUE, "linear",  -12.8, 6, 6, 3, 1, 0, 1, 1, -12.8, 0.5, 0.5, 0,
  "negative", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 10, 200, TRUE, "emax",  -12.8, 6, 6, 3, 1, 0, 1, 1, -12.8, 0.5, 0.5, 0,

  "negative_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 10, 200, FALSE, "linear",  -12.8, 6, 6, 3, 1, 0, 1, 1, -12.8, 0.5, 0.5, 0,
  "negative_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 10, 200, FALSE, "emax",  -12.8, 6, 6, 3, 1, 0, 1, 1, -12.8, 0.5, 0.5, 0,

  "less_negative", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 13, 200, TRUE, "linear",  -12.8, 8, 8, 3, 1, 0, 1, 1, -12.8, 0.7, 0.7, 0,
  "less_negative", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 13, 200, TRUE, "emax",  -12.8, 8, 8, 3, 1, 0, 1, 1, -12.8, 0.7, 0.7, 0,

  "less_negative_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 13, 200, FALSE, "linear",  -12.8, 8, 8, 3, 1, 0, 1, 1, -12.8, 0.7, 0.7, 0,
  "less_negative_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 13, 200, FALSE, "emax",  -12.8, 8, 8, 3, 1, 0, 1, 1, -12.8, 0.7, 0.7, 0,

  "low_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.01, -12.8, 0.9, 0.9, 0,
  "low_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.01, -12.8, 0.9, 0.9, 0,

  "low_var_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.01, -12.8, 0.9, 0.9, 0,
  "low_var_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.01, -12.8, 0.9, 0.9, 0,

  "mid_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.1, -12.8, 0.9, 0.9, 0,
  "mid_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.1, -12.8, 0.9, 0.9, 0,

  "mid_var_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.1, -12.8, 0.9, 0.9, 0,
  "mid_var_no_interim", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.1, -12.8, 0.9, 0.9, 0,

  "plateau", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4), 8, 11, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "plateau", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4), 8, 11, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "plateau_high", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4, 6, 9), 8, 12, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
  "plateau_high", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4, 6, 9), 8, 12, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

  "negative_EMAX7", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "linear",  -12.8, 7, 7, 3, 1, 0, 1, 1, -12.8, 0.6, 0.6, 0,
  "negative_EMAX7", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "emax",  -12.8, 7, 7, 3, 1, 0, 1, 1, -12.8, 0.6, 0.6, 0,

  "negative_EMAX6_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "linear",  -12.8, 6.5, 6.5, 3, 1, 0, 1, 1, -12.8, 0.65, 0.65, 0,
  "negative_EMAX6_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "emax",  -12.8, 6.5, 6.5, 3, 1, 0, 1, 1, -12.8, 0.65, 0.65, 0
)

# scenario generating function - add each scenario here
#' Title
#'
#' @param scenario_label  tbd
#'
#' @return df
#' @export
get_scenario <- function(scenario_label) {
  if (scenario_label == "viable_no_interim") {
    return(viable_no_interim$scenario)
  }
  if (scenario_label == "viable_interim") {
    return(viable_interim$scenario)
  }
  if (scenario_label == "viable_no_interim_interaction_0_5") {
    return(viable_no_interim_interaction_0_5$scenario)
  }
  if (scenario_label == "viable_interim_interaction_0_5") {
    return(viable_interim_interaction_0_5$scenario)
  }
  if (scenario_label == "viable_no_interim_full_interaction") {
    return(viable_no_interim_full_interaction$scenario)
  }
  if (scenario_label == "viable_interim_full_interaction") {
    return(viable_interim_full_interaction$scenario)
  }
  if (scenario_label == "viable_no_interim_equal") {
    return(viable_no_interim$scenario)
  }
  if (scenario_label == "viable_interim_equal") {
    return(viable_interim$scenario)
  }
  if (scenario_label == "viable_no_interim_interaction_0_5_equal") {
    return(viable_no_interim_interaction_0_5$scenario)
  }
  if (scenario_label == "viable_interim_interaction_0_5_equal") {
    return(viable_interim_interaction_0_5$scenario)
  }
  if (scenario_label == "viable_no_interim_full_interaction_equal") {
    return(viable_no_interim_full_interaction$scenario)
  }
  if (scenario_label == "viable_interim_full_interaction_equal") {
    return(viable_interim_full_interaction$scenario)
  }
  if (scenario_label == "viable_no_interim_low") {
    return(viable_no_interim$scenario)
  }
  if (scenario_label == "viable_interim_low") {
    return(viable_interim$scenario)
  }
  if (scenario_label == "viable_no_interim_high") {
    return(viable_no_interim$scenario)
  }
  if (scenario_label == "viable_interim_high") {
    return(viable_interim$scenario)
  }
  if (scenario_label == "null_equal") {
    return(null_equal$scenario)
  }
  if (scenario_label == "null") {
    return(null$scenario)
  }
  if (scenario_label == "positive_equal") {
    return(positive_equal$scenario)
  }
  if (scenario_label == "positive") {
    return(positive$scenario)
  }
  if (scenario_label == "low_var") {
    return(low_var$scenario)
  }
  if (scenario_label == "mid_var") {
    return(mid_var$scenario)
  }
  if (scenario_label == "plateau") {
    return(plateau$scenario)
  }
  if (scenario_label == "plateau_high") {
    return(plateau_high$scenario)
  }
  if (scenario_label == "low_var_no_interim") {
    return(low_var_no_interim$scenario)
  }
  if (scenario_label == "mid_var_no_interim") {
    return(mid_var_no_interim$scenario)
  }
  if (scenario_label == "slight_positive") {
    return(slight_positive$scenario)
  }
  if (scenario_label == "negative") {
    return(negative$scenario)
  }
  if (scenario_label == "negative_no_interim") {
    return(negative_no_interim$scenario)
  }
  if (scenario_label == "slight_positive_no_interim") {
    return(slight_positive_no_interim$scenario)
  }
  if (scenario_label == "null_no_interim") {
    return(null_no_interim$scenario)
  }
  if (scenario_label == "negative_EMAX7") {
    return(negative_EMAX7$scenario)
  }
  if (scenario_label == "negative_EMAX6_5") {
    return(negative_EMAX6_5$scenario)
  }

}

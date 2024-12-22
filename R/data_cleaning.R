# tar_load(trialsims)

#' @title calculate_columns
#'
#' @param scenario_row the row of the scenario table
#' @param trialsims  the result dataset - potentially a targets object
#'
#' @return result dataframe containing all relevent information for that scenario
#' @export
calculate_columns <- function(scenario_row, trialsims) {
  scenario_label <- scenario_row$scenario_label
  likelihood <- scenario_row$likelihood

  result <- trialsims %>%
    dplyr::filter(unlist(.$scenario_label) == !!scenario_label, unlist(.$likelihood) == !!likelihood) %>%
    dplyr::mutate(
      true_med == .$true_med
    )

  # Add scenario_label and rho as the first two columns
  result <- result %>%
    mutate(scenario_label = scenario_label, likelihood = likelihood)

  return(result)
}

# results <- lapply(seq_len(nrow(scenarios_table)), function(i) {
#   calculate_columns(scenarios_table[i, ], trialsims)
# })

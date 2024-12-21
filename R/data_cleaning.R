# tar_load(trialsims)

#' Title
#'
#' @param scenario_row tbd
#' @param trialsims  tbd
#'
#' @return tbd
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

library(shiny)
library(gt)
library(ggplot2)
library(reshape2)
library(gridExtra)

ui <- fluidPage(
  # Application title
  titlePanel("Dose Schedule Finding Simulation"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "scen",
        label = "Scenario selection:",
        choices = 1:93,
        selected = 4
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(id = "tabset",
                  tabPanel(
                    "Scenarios",
                    value = "scenarios",
                    gt::gt_output("scenarios"),
                  ), tabPanel(
                    "Tabular Overview",
                    gt::gt_output("table")
                  ),
                  tabPanel(
                    "MED Probabilities",
                    h3("Linear"),
                    plotOutput("props_1"),
                    h3("EMAX"),
                    plotOutput("props_2")
                  ),
                  tabPanel(
                    "MED Effects",
                    h3("Linear"),
                    plotOutput("effects_1"),
                    h3("EMAX"),
                    plotOutput("effects_2")
                  )
      )
    )
  )
)

server <- function(input, output, session) {

  get_MED <- function(results, scenario, likelihood, var) {

    sum <- numeric(length(results[[likelihood]][[scenario]][[var]]))

    for(i in 1:length(results[[likelihood]][[scenario]])) {
      if(i == 1) {
        sum <- results[[likelihood]][[scenario]][[i]][[var]]
      } else {
        sum <- sum + (results[[likelihood]][[scenario]][[i]][[var]])
      }
    }

    sum_1 <- sum / length(results[[likelihood]][[scenario]])
    return(sum_1 / 5000)
  }

  get_effect <- function(results, scenario, likelihood, var) {
    sum <- numeric(length(results[[likelihood]][[scenario]][[var]]))

    for(i in 1:length(results[[likelihood]][[scenario]])) {
      if(i == 1) {
        sum <- results[[likelihood]][[scenario]][[i]][[var]]
      } else {
        sum <- sum + (results[[likelihood]][[scenario]][[i]][[var]])
      }
    }

    sum_1 <- sum / length(results[[likelihood]][[scenario]])
    return(sum_1 / 5000)
  }
  output$scenarios <- gt::render_gt({
    scenarios_table <- tribble(
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

      "plateau", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4), 8, 18, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "plateau", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4), 8, 18, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

      "plateau_high", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4, 6, 9), 8, 19, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "plateau_high", c(0, 4, 5, 6, 9, 10, 12), c(1, 3, 4, 6, 9), 8, 19, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

      "negative_EMAX7", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "linear",  -12.8, 7, 7, 3, 1, 0, 1, 1, -12.8, 0.6, 0.6, 0,
      "negative_EMAX7", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 14, 200, TRUE, "emax",  -12.8, 7, 7, 3, 1, 0, 1, 1, -12.8, 0.6, 0.6, 0,

      "negative_EMAX6_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 15, 200, TRUE, "linear",  -12.8, 6.5, 6.5, 3, 1, 0, 1, 1, -12.8, 0.65, 0.65, 0,
      "negative_EMAX6_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 15, 200, TRUE, "emax",  -12.8, 6.5, 6.5, 3, 1, 0, 1, 1, -12.8, 0.65, 0.65, 0,

      "mid_var5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.02, -12.8, 0.9, 0.9, 0,
      "mid_var5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.02, -12.8, 0.9, 0.9, 0,

      "mid_var2_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear",  -12.8, 10, 10, 3, 1, 0, 1, 0.15, -12.8, 0.9, 0.9, 0,
      "mid_var2_5", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.15, -12.8, 0.9, 0.9, 0,

      "mid_var1_8", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.18, -12.8, 0.9, 0.9, 0,
      "mid_var1_6", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.16, -12.8, 0.9, 0.9, 0,
      "mid_var1_4", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax",  -12.8, 10, 10, 3, 1, 0, 1, 0.14, -12.8, 0.9, 0.9, 0,

      "viable_no_interim_100", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 16, 100, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_no_interim_100", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 16, 100, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_interim_100", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 16, 100, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_interim_100", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 16, 100, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

      "viable_no_interim_150", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 17, 150, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_no_interim_150", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 17, 150, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_interim_150", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 17, 150, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,
      "viable_interim_150", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 17, 150, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1, -12.8, 0.9, 0.9, 0,

      "viable_no_interim_high_var2", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1.5, -12.8, 0.9, 0.9, 0,
      "viable_no_interim_high_var2", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1.5, -12.8, 0.9, 0.9, 0,
      "viable_interim_high_var2", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 1.5, -12.8, 0.9, 0.9, 0,
      "viable_interim_high_var2", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 1.5, -12.8, 0.9, 0.9, 0,

      "viable_no_interim_high_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 2, -12.8, 0.9, 0.9, 0,
      "viable_no_interim_high_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, FALSE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 2, -12.8, 0.9, 0.9, 0,
      "viable_interim_high_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "linear", -12.8, 10, 10, 3, 1, 0, 1, 2, -12.8, 0.9, 0.9, 0,
      "viable_interim_high_var", c(0, 4, 5, 6, 9), c(1, 3, 4), 8, 1, 200, TRUE, "emax", -12.8, 10, 10, 3, 1, 0, 1, 2, -12.8, 0.9, 0.9, 0
    )


    gt(scenarios_table)
  })

  data <- reactive({
    results <- readRDS("inst/data/results_short.RDS")
    df <- readRDS("inst/data/table.RDS")
    df_no_interim <- readRDS("inst/data/table_no_interim.RDS")

    list(
      results = results,
      df = df,
      df_no_interim = df_no_interim
    )
  })

  output$table <- gt::render_gt({
    data <- data()
    results <- data$results
    df  <- data$df
    df_no_interim <- data$df_no_interim

    if(input$scen %in% which(scenarios_table$interim == TRUE)) {

      df <- df %>%
        dplyr::filter(scenario == input$scen)
      df <- df %>%
        dplyr::select(
          -c(scenario, scenario_name, effect, interim, epsilon_sd, true_model)
        )
      number <- as.numeric(input$scen)
      gt(df) %>%
        tab_header(
          title = paste0("Dose Schedule Finding Scenario: ", input$scen),
          subtitle = paste0(
            "Sample Size: ", results$n_pat[[number]],
            ", Doses: ", paste(results$doses[[number]], collapse = ","),
            ", Schedules: ", paste(results$schedules[[number]], collapse = ","),
            ", E_max: ", results$alpha1[[number]],
            ", Beta2 ", results$beta2[[number]],
            ", Interaction ", results$beta3[[number]],
            ", Interim ", results$interim[[number]],
            ", Allocation Matrix ", results$alloc_mat[[number]],
            ", Epsilon Sigma ", results$sig_sd[[number]]
          )
        ) %>%
        tab_spanner(
          label = "Go/No-Go",
          columns = c(gng_lin, gng_ema, gng_aic, gng_bic)
        ) %>%
        tab_spanner(
          label = "Correct MED Estimation",
          columns = c(corr_med_lin, corr_med_ema, corr_med_aic, corr_med_bic)
        ) %>%
        tab_spanner(
          label = "MSE",
          columns = c(mse_lin, mse_ema, mse_aic, mse_bic)
        ) %>%
        tab_spanner(
          label = "RMSE",
          columns = c(rmse_lin, rmse_ema, rmse_aic, rmse_bic)
        ) %>%
        tab_spanner(
          label = "Patient allocation to MED",
          columns = c(pat_alloc_lin, pat_alloc_ema, pat_alloc_aic, pat_alloc_bic)
        ) %>%
        fmt_number(decimals = 5, columns = c(mse_lin, mse_ema, mse_aic, mse_bic,
                                             rmse_lin, rmse_ema, rmse_aic, rmse_bic)) %>%
        cols_label(
          design = "Design",
          gng_lin = "Linear",
          gng_ema = "EMAX",
          gng_aic = "AIC",
          gng_bic = "BIC",
          corr_med_lin = "Linear",
          corr_med_ema = "EMAX",
          corr_med_aic = "AIC",
          corr_med_bic = "BIC",
          mse_lin = "Linear",
          mse_ema = "EMAX",
          mse_aic = "AIC",
          mse_bic = "BIC",
          rmse_lin = "Linear",
          rmse_ema  = "EMAX",
          rmse_aic = "AIC",
          rmse_bic = "BIC",
          pat_alloc_lin = "Linear",
          pat_alloc_ema = "EMAX",
          pat_alloc_aic = "AIC",
          pat_alloc_bic = "BIC"
        )

    } else {

      df <- df_no_interim %>%
        dplyr::filter(scenario == input$scen)
      df <- df %>%
        dplyr::select(
          -c(scenario, scenario_name, effect, interim, epsilon_sd, true_model)
        )
      number <- as.numeric(input$scen)
      gt(df) %>%
        tab_header(
          title = paste0("Dose Schedule Finding Scenario: ", input$scen),
          subtitle = paste0(
            "Sample Size: ", results$n_pat[[number]],
            ", Doses: ", paste(results$doses[[number]], collapse = ","),
            ", Schedules: ", paste(results$schedules[[number]], collapse = ","),
            ", E_max: ", results$alpha1[[number]],
            ", Beta2 ", results$beta2[[number]],
            ", Interaction ", results$beta3[[number]],
            ", Interim ", results$interim[[number]],
            ", Allocation Matrix ", results$alloc_mat[[number]],
            ", Epsilon Sigma ", results$sig_sd[[number]]
          )
        ) %>%
        tab_spanner(
          label = "Go/No-Go",
          columns = c(gng_lin, gng_ema, gng_aic, gng_bic)
        ) %>%
        tab_spanner(
          label = "Correct MED Estimation",
          columns = c(corr_med_lin, corr_med_ema, corr_med_aic, corr_med_bic)
        ) %>%
        tab_spanner(
          label = "MSE",
          columns = c(mse_lin, mse_ema, mse_aic, mse_bic)
        ) %>%
        tab_spanner(
          label = "RMSE",
          columns = c(rmse_lin, rmse_ema, rmse_aic, rmse_bic)
        ) %>%
        tab_spanner(
          label = "Patient allocation to MED",
          columns = c(pat_alloc, pat_alloc_aic, pat_alloc_bic)
        ) %>%
        fmt_number(decimals = 5, columns = c(mse_lin, mse_ema, mse_aic, mse_bic,
                                             rmse_lin, rmse_ema, rmse_aic, rmse_bic)) %>%
        cols_label(
          design = "Design",
          gng_lin = "Linear",
          gng_ema = "EMAX",
          gng_aic = "AIC",
          gng_bic = "BIC",
          corr_med_lin = "Linear",
          corr_med_ema = "EMAX",
          corr_med_aic = "AIC",
          corr_med_bic = "BIC",
          mse_lin = "Linear",
          mse_ema = "EMAX",
          mse_aic = "AIC",
          mse_bic = "BIC",
          rmse_lin = "Linear",
          rmse_ema  = "EMAX",
          rmse_aic = "AIC",
          rmse_bic = "BIC",
          pat_alloc = "Fixed",
          pat_alloc_aic = "AIC",
          pat_alloc_bic = "BIC"
        )

    }


  })

  output$props_1 <- renderPlot({
    data <- data()
    results <- data$results
    scenario = as.numeric(input$scen)
    doses <- results$doses[[scenario]]
    schedules <- results$schedules[[scenario]]




    if(input$scen %in% which(scenarios_table$interim == TRUE)) {
      observed_med_1_lin <- get_MED(results, scenario, "linear", 1)
      observed_med_2_lin <- get_MED(results, scenario, "linear", 2)
      observed_med_3_lin <- get_MED(results, scenario, "linear", 9)
      observed_med_4_lin <- get_MED(results, scenario, "linear", 13)
      observed_med_5_lin <- get_MED(results, scenario, "linear", 17)
      observed_med_6_lin <- get_MED(results, scenario, "linear", 21)
      observed_med_7_lin <- get_MED(results, scenario, "linear", 25)

      observed_med_1_ema <- get_MED(results, scenario, "emax", 1)
      observed_med_2_ema <- get_MED(results, scenario, "emax", 2)
      observed_med_3_ema <- get_MED(results, scenario, "emax", 9)
      observed_med_4_ema <- get_MED(results, scenario, "emax", 13)
      observed_med_5_ema <- get_MED(results, scenario, "emax", 17)
      observed_med_6_ema <- get_MED(results, scenario, "emax", 21)
      observed_med_7_ema <- get_MED(results, scenario, "emax", 25)


      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)
      mat13 <- matrix(observed_med_7_lin, nrow = 3, byrow = TRUE)
      mat14 <- matrix(observed_med_7_ema, nrow = 3, byrow = TRUE)
      mat15 <- matrix(results$true_med[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat13,
        mat14,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat13_melt <- melt(mat13)
      mat14_melt <- melt(mat14)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(
        mat_melt,
        mat2_melt,
        mat3_melt,
        mat4_melt,
        mat5_melt,
        mat6_melt
      )


      # Create a list to store the plots
      plots <- list()

      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Probability") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }
    } else {
      observed_med_1_lin <- get_MED(results, scenario, "linear", 1)
      observed_med_2_lin <- get_MED(results, scenario, "linear", 2)
      observed_med_3_lin <- get_MED(results, scenario, "linear", 9)
      observed_med_4_lin <- get_MED(results, scenario, "linear", 13)
      observed_med_5_lin <- get_MED(results, scenario, "linear", 17)
      observed_med_6_lin <- get_MED(results, scenario, "linear", 21)


      observed_med_1_ema <- get_MED(results, scenario, "emax", 1)
      observed_med_2_ema <- get_MED(results, scenario, "emax", 2)
      observed_med_3_ema <- get_MED(results, scenario, "emax", 9)
      observed_med_4_ema <- get_MED(results, scenario, "emax", 13)
      observed_med_5_ema <- get_MED(results, scenario, "emax", 17)
      observed_med_6_ema <- get_MED(results, scenario, "emax", 21)



      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)

      mat15 <- matrix(results$true_med[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(
        mat_melt,
        mat2_melt,
        mat3_melt,
        mat4_melt,
        mat5_melt,
        mat6_melt
      )


      # Create a list to store the plots
      plots <- list()

      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Probability") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }
    }



    # Arrange the plots side by side
    do.call(grid.arrange, c(plots, ncol = 3))

  })

  output$props_2 <- renderPlot({
    data <- data()
    results <- data$results
    scenario = as.numeric(input$scen)
    doses <- results$doses[[scenario]]
    schedules <- results$schedules[[scenario]]



    if(input$scen %in% which(scenarios_table$interim == TRUE)) {
      observed_med_1_lin <- get_MED(results, scenario, "linear", 1)
      observed_med_2_lin <- get_MED(results, scenario, "linear", 2)
      observed_med_3_lin <- get_MED(results, scenario, "linear", 9)
      observed_med_4_lin <- get_MED(results, scenario, "linear", 13)
      observed_med_5_lin <- get_MED(results, scenario, "linear", 17)
      observed_med_6_lin <- get_MED(results, scenario, "linear", 21)
      observed_med_7_lin <- get_MED(results, scenario, "linear", 25)

      observed_med_1_ema <- get_MED(results, scenario, "emax", 1)
      observed_med_2_ema <- get_MED(results, scenario, "emax", 2)
      observed_med_3_ema <- get_MED(results, scenario, "emax", 9)
      observed_med_4_ema <- get_MED(results, scenario, "emax", 13)
      observed_med_5_ema <- get_MED(results, scenario, "emax", 17)
      observed_med_6_ema <- get_MED(results, scenario, "emax", 21)
      observed_med_7_ema <- get_MED(results, scenario, "emax", 25)


      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)
      mat13 <- matrix(observed_med_7_lin, nrow = 3, byrow = TRUE)
      mat14 <- matrix(observed_med_7_ema, nrow = 3, byrow = TRUE)
      mat15 <- matrix(results$true_med[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat13,
        mat14,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat13_melt <- melt(mat13)
      mat14_melt <- melt(mat14)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(
        mat7_melt,
        mat8_melt,
        mat9_melt,
        mat10_melt,
        mat11_melt,
        mat12_melt,
        mat14_melt,
        mat15_melt
      )


      # Create a list to store the plots
      plots <- list()
      vec <- c("Full Factorial, Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
               "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim")
      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim", "True MED")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Probability") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

      # Arrange the plots side by side
      do.call(grid.arrange, c(plots, ncol = 3))

    } else {
      observed_med_1_lin <- get_MED(results, scenario, "linear", 1)
      observed_med_2_lin <- get_MED(results, scenario, "linear", 2)
      observed_med_3_lin <- get_MED(results, scenario, "linear", 9)
      observed_med_4_lin <- get_MED(results, scenario, "linear", 13)
      observed_med_5_lin <- get_MED(results, scenario, "linear", 17)
      observed_med_6_lin <- get_MED(results, scenario, "linear", 21)

      observed_med_1_ema <- get_MED(results, scenario, "emax", 1)
      observed_med_2_ema <- get_MED(results, scenario, "emax", 2)
      observed_med_3_ema <- get_MED(results, scenario, "emax", 9)
      observed_med_4_ema <- get_MED(results, scenario, "emax", 13)
      observed_med_5_ema <- get_MED(results, scenario, "emax", 17)
      observed_med_6_ema <- get_MED(results, scenario, "emax", 21)

      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)
      mat15 <- matrix(results$true_med[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(
        mat7_melt,
        mat8_melt,
        mat9_melt,
        mat10_melt,
        mat11_melt,
        mat12_melt,
        mat15_melt
      )


      # Create a list to store the plots
      plots <- list()
      vec <- c("Full Factorial, Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
               "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim")
      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "True MED")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Probability") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

      # Arrange the plots side by side
      do.call(grid.arrange, c(plots, ncol = 3))

    }


  })


  output$effects_1 <- renderPlot({
    data <- data()
    results <- data$results
    scenario = as.numeric(input$scen)
    doses <- results$doses[[scenario]]
    schedules <- results$schedules[[scenario]]

    if(input$scen %in% which(scenarios_table$interim == TRUE)) {

      observed_med_1_lin <- get_effect(results, scenario, "linear", 3)
      observed_med_2_lin <- get_effect(results, scenario, "linear", 4)
      observed_med_3_lin <- get_effect(results, scenario, "linear", 10)
      observed_med_4_lin <- get_effect(results, scenario, "linear", 14)
      observed_med_5_lin <- get_effect(results, scenario, "linear", 18)
      observed_med_6_lin <- get_effect(results, scenario, "linear", 22)
      observed_med_7_lin <- get_effect(results, scenario, "linear", 26)

      observed_med_1_ema <- get_effect(results, scenario, "emax", 3)
      observed_med_2_ema <- get_effect(results, scenario, "emax", 4)
      observed_med_3_ema <- get_effect(results, scenario, "emax", 10)
      observed_med_4_ema <- get_effect(results, scenario, "emax", 14)
      observed_med_5_ema <- get_effect(results, scenario, "emax", 18)
      observed_med_6_ema <- get_effect(results, scenario, "emax", 22)
      observed_med_7_ema <- get_effect(results, scenario, "emax", 26)


      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)
      mat13 <- matrix(observed_med_7_lin, nrow = 3, byrow = TRUE)
      mat14 <- matrix(observed_med_7_ema, nrow = 3, byrow = TRUE)
      mat15 <- matrix(results$true_effect[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat13,
        mat14,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat13_melt <- melt(mat13)
      mat14_melt <- melt(mat14)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(
        mat_melt,
        mat2_melt,
        mat3_melt,
        mat4_melt,
        mat5_melt,
        mat6_melt,
        mat7_melt
      )


      # Create a list to store the plots
      plots <- list()

      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Effect") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

    } else {

      observed_med_1_lin <- get_effect(results, scenario, "linear", 3)
      observed_med_2_lin <- get_effect(results, scenario, "linear", 4)
      observed_med_3_lin <- get_effect(results, scenario, "linear", 10)
      observed_med_4_lin <- get_effect(results, scenario, "linear", 14)
      observed_med_5_lin <- get_effect(results, scenario, "linear", 18)
      observed_med_6_lin <- get_effect(results, scenario, "linear", 22)


      observed_med_1_ema <- get_effect(results, scenario, "emax", 3)
      observed_med_2_ema <- get_effect(results, scenario, "emax", 4)
      observed_med_3_ema <- get_effect(results, scenario, "emax", 10)
      observed_med_4_ema <- get_effect(results, scenario, "emax", 14)
      observed_med_5_ema <- get_effect(results, scenario, "emax", 18)
      observed_med_6_ema <- get_effect(results, scenario, "emax", 22)



      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)


      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)


      # Create a list of the melted matrices
      melted_mats <- list(
        mat_melt,
        mat2_melt,
        mat3_melt,
        mat4_melt,
        mat5_melt,
        mat6_melt
      )


      # Create a list to store the plots
      plots <- list()

      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base")
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Effect") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

    }

    do.call(grid.arrange, c(plots, ncol = 3))


  })

  output$effects_2 <- renderPlot({
    data <- data()
    results <- data$results
    scenario = as.numeric(input$scen)
    doses <- results$doses[[scenario]]
    schedules <- results$schedules[[scenario]]


    if(input$scen %in% which(scenarios_table$interim == TRUE)) {


      observed_med_1_lin <- get_effect(results, scenario, "linear", 3)
      observed_med_2_lin <- get_effect(results, scenario, "linear", 4)
      observed_med_3_lin <- get_effect(results, scenario, "linear", 10)
      observed_med_4_lin <- get_effect(results, scenario, "linear", 14)
      observed_med_5_lin <- get_effect(results, scenario, "linear", 18)
      observed_med_6_lin <- get_effect(results, scenario, "linear", 22)
      observed_med_7_lin <- get_effect(results, scenario, "linear", 26)

      observed_med_1_ema <- get_effect(results, scenario, "emax", 3)
      observed_med_2_ema <- get_effect(results, scenario, "emax", 4)
      observed_med_3_ema <- get_effect(results, scenario, "emax", 10)
      observed_med_4_ema <- get_effect(results, scenario, "emax", 14)
      observed_med_5_ema <- get_effect(results, scenario, "emax", 18)
      observed_med_6_ema <- get_effect(results, scenario, "emax", 22)
      observed_med_7_ema <- get_effect(results, scenario, "emax", 26)


      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)
      mat13 <- matrix(observed_med_7_lin, nrow = 3, byrow = TRUE)
      mat14 <- matrix(observed_med_7_ema, nrow = 3, byrow = TRUE)
      mat15 <- matrix(results$true_effect[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,
        mat12,
        mat13,
        mat14,
        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)
      mat13_melt <- melt(mat13)
      mat14_melt <- melt(mat14)
      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(

        mat8_melt,
        mat9_melt,
        mat10_melt,
        mat11_melt,
        mat12_melt,
        mat13_melt,
        mat14_melt,
        mat15_melt
      )


      # Create a list to store the plots
      plots <- list()

      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base", "Full Factorial with Interim",
                 paste0("True MED Effect ", results$true_effect[[scenario]][which(results$true_med[[scenario]]>0)]))
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Effect") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

    } else {


      observed_med_1_lin <- get_effect(results, scenario, "linear", 3)
      observed_med_2_lin <- get_effect(results, scenario, "linear", 4)
      observed_med_3_lin <- get_effect(results, scenario, "linear", 10)
      observed_med_4_lin <- get_effect(results, scenario, "linear", 14)
      observed_med_5_lin <- get_effect(results, scenario, "linear", 18)
      observed_med_6_lin <- get_effect(results, scenario, "linear", 22)


      observed_med_1_ema <- get_effect(results, scenario, "emax", 3)
      observed_med_2_ema <- get_effect(results, scenario, "emax", 4)
      observed_med_3_ema <- get_effect(results, scenario, "emax", 10)
      observed_med_4_ema <- get_effect(results, scenario, "emax", 14)
      observed_med_5_ema <- get_effect(results, scenario, "emax", 18)
      observed_med_6_ema <- get_effect(results, scenario, "emax", 22)



      mat <- matrix(observed_med_1_lin, nrow = 3, byrow = TRUE)
      mat2 <- matrix(observed_med_2_lin, nrow = 3, byrow = TRUE)
      mat3 <- matrix(observed_med_3_lin, nrow = 3, byrow = TRUE)
      mat4 <- matrix(observed_med_4_lin, nrow = 3, byrow = TRUE)
      mat5 <- matrix(observed_med_5_lin, nrow = 3, byrow = TRUE)
      mat6 <- matrix(observed_med_6_lin, nrow = 3, byrow = TRUE)
      mat7 <- matrix(observed_med_1_ema, nrow = 3, byrow = TRUE)
      mat8 <- matrix(observed_med_2_ema, nrow = 3, byrow = TRUE)
      mat9 <- matrix(observed_med_3_ema, nrow = 3, byrow = TRUE)
      mat10 <- matrix(observed_med_4_ema, nrow = 3, byrow = TRUE)
      mat11 <- matrix(observed_med_5_ema, nrow = 3, byrow = TRUE)
      mat12 <- matrix(observed_med_6_ema, nrow = 3, byrow = TRUE)

      mat15 <- matrix(results$true_effect[[scenario]], nrow = 3, byrow = TRUE)

      mats <- list(
        mat,
        mat2,
        mat3,
        mat4,
        mat5,
        mat6,
        mat7,
        mat8,
        mat9,
        mat10,
        mat11,

        mat15
      )

      for(i in 1:length(mats)) {
        colnames(mats[[i]]) <- doses
        rownames(mats[[i]]) <- schedules
      }


      # Combine the matrices into a list
      # mat_list <- list(mat, mat2, mat3)

      # Create a new matrix that combines all three matrices
      combined_mat <- do.call(cbind, mats)


      # Melt the matrices into long format for ggplot
      mat_melt <- melt(mat)
      mat2_melt <- melt(mat2)
      mat3_melt <- melt(mat3)
      mat4_melt <- melt(mat4)
      mat5_melt <- melt(mat5)
      mat6_melt <- melt(mat6)
      mat7_melt <- melt(mat7)
      mat8_melt <- melt(mat8)
      mat9_melt <- melt(mat9)
      mat10_melt <- melt(mat10)
      mat11_melt <- melt(mat11)
      mat12_melt <- melt(mat12)

      mat15_melt <- melt(mat15)

      # Create a list of the melted matrices
      melted_mats <- list(

        mat7_melt,
        mat8_melt,
        mat9_melt,
        mat10_melt,
        mat11_melt,
        mat12_melt,
        mat15_melt
      )


      # Create a list to store the plots
      plots <- list()

      # Create a heatmap for each matrix and store it in the list
      for (i in 1:length(melted_mats)) {
        melted_mats[[i]]$Var1 <- factor(melted_mats[[i]]$Var1)
        melted_mats[[i]]$Var2 <- factor(melted_mats[[i]]$Var2)
        vec <- c("Full Factorial", "Custom Corner Mid", "D-Optimitality with Linear Base", "D-Optimitality with EMAX Base",
                 "I-Optimitality with Linear Base", "I-Optimitality with EMAX Base",
                 paste0("True MED Effect ", results$true_effect[[scenario]][which(results$true_med[[scenario]]>0)]))
        # Get the minimum and maximum values in the data
        min_val <- min(melted_mats[[i]]$value[melted_mats[[i]]$value!=0], na.rm = TRUE)
        max_val <- max(melted_mats[[i]]$value, na.rm = TRUE)
        p <- ggplot(data = melted_mats[[i]], aes(x = Var1, y = Var2, fill = value)) +
          geom_tile() +
          geom_text(aes(label = round(value, 1)), size = 3) +
          scale_fill_gradientn(colors = colorRampPalette(c("white", "orange"))(20)) +
          theme_minimal() +
          ggtitle(paste("Heatmap of Design", vec[i])) +
          labs(x = "Schedules", y = "Doses", fill = "Effect") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_x_discrete(breaks = unique(melted_mats[[i]]$Var1)) +
          scale_y_discrete(breaks = unique(melted_mats[[i]]$Var2))

        plots[[i]] <- p
      }

    }

    # Arrange the plots side by side
    do.call(grid.arrange, c(plots, ncol = 3))


  })


}

shinyApp(ui, server)

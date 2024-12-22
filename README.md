<!-- badges: start -->

[![R-CMD-check](https://github.com/Xyarz/DoseScheduleOptimizeR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Xyarz/DoseScheduleOptimizeR/actions/workflows/R-CMD-check.yaml) <!-- badges: end -->

# DoseScheduleOptimizeR

Dose Schedule Optimization - Simulation based

## Installation

``` r
install.packages("remotes")
remotes::install_github("https://github.com/xyarz/DoseScheduleOptimizeR")
```

## Usage of ShinyApp

``` r
library(DoseScheduleOptimizeR)

DoseScheduleOptimizeR::runApp()
```

The aim of the application is to provide the option to look into all 93 scenarios created in this research across tables and visualizations provided in this thesis. The data and tables for this have been attached to the R package inside the inst/data folder in RDS format and the heatmap visualizations of the MED probabilities and the respective effects are calculated based on the scenario selection made by the user, that can be found inside a drop-down menu on the left side. In the following there are a couple screenshots, displaying the capabilities and tools available inside the application.

Firstly, there is the Landing page, if you start the application, one will find as aforementioned a drop-down selection for the scenario to be looked at on the left side and in the middle a tabular overview of the parameters chosen for each of the 93 available scenarios.

![Landing Page](inst/img/landing_page.png)

Secondly, one can see the tabular overview of the performance measurements for the chosen scenario.

![Tabular Overview](inst/img/table.png)

Here, the user can investigate the MED probabilities on the factorial space across designs and evaluation model.

![MED Probabilities](inst/img/props.png)

Lastly, Here, the user can investigate the placebo adjusted MED effects that surpassed the pre-defined threshold on the factorial space across designs and evaluation model.

![MED Effects](inst/img/effects.png)

## Further Information

## Initialize
library(ferx)
library(vpc)
library(dplyr)
library(tidyr)
library(ggplot2)

ex <- ferx_example("warfarin")

## Fit
result <- ferx_fit(
  model = ex$model, 
  data = ex$data, 
  method = "gn", ## Gauss-Newton
  covariance = FALSE
)

result_foce <- ferx_fit(
  model = ex$model, 
  data = ex$data, 
  method = "gn" ## Gauss-Newton
)

result_saem <- ferx_fit(
  model = ex$model, 
  data = ex$data, 
  method = "saem"
)

## Basic GOF
result$sdtab |>
  tidyr::pivot_longer(cols = c("PRED", "IPRED")) |>
  ggplot(aes(x = value, y = DV)) +
    geom_point() +
    facet_wrap(~name)
result$sdtab |>
  tidyr::pivot_longer(cols = c("CWRES", "DV")) |>
  ggplot(aes(x = PRED, y = value)) +
    geom_point() +
    facet_wrap(~name, scales = "free") +
    ylab("")

# Simulate and create VPC
sim <- ferx_simulate(
  ex$model, 
  ex$data, 
  n_sim = 200, 
  seed = 42,
  fit = result
)
obs <- read.csv(file = ex$data)
vpc(
  obs = obs |> 
    mutate(DV = as.numeric(DV)), 
  sim = sim, 
  obs_cols = list(), 
  sim_cols = list(dv = "DV_SIM")
)

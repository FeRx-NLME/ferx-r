## Initialize
library(vpc)
library(dplyr)

model <- "warfarin.ferx"
data <- "../data/warfarin.csv"

# model <- "two_cpt_oral_cov.ferx"
# data <- "../data/two_cpt_oral_cov.csv"

## Fit
result <- ferx_fit(
  model = model, 
  data = data, 
  method = "focei"
)
# result$theta    # named vector of estimates
# result$omega    # BSV covariance matrix
# result$sdtab    # data.frame with PRED, IPRED, CWRES, IWRES

# Simulate and create VPC
sim <- ferx_simulate(model, data, n_sim = 100, seed = 42)
obs <- read.csv(file = data)
vpc(
  obs = obs |> 
    mutate(DV = as.numeric(DV)), 
  sim = sim, 
  obs_cols = list(), 
  sim_cols = list(dv = "DV_SIM")
)

# Population predictions
preds <- ferx_predict(
  model = model, 
  data = data
)

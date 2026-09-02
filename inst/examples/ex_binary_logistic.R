library(ferx)

# Binary / logistic endpoint ([binary_model], ferx-core #760 Track C).
#
# A fixed-effects (n_eta = 0) logistic regression: the 0/1 outcome DV on CMT 3 is
# Bernoulli with log-odds TH0 + THX * X + THT * TIME. With no random effects this
# is exactly base-R glm(DV ~ X + TIME, family = binomial); the bundled data was
# generated from TH0 = -0.4, THX = 0.9, THT = 0.5.
ex <- ferx_example("binary_logistic")

fit <- ferx_fit(ex$model, ex$data, method = "focei")
print(fit)

# Simulate: ferx draws a fresh 0/1 outcome per record (Bernoulli at the fitted
# probability) and returns it in DV_SIM on the binary CMT (3). Before ferx-r #271
# these rows came back as DV_SIM = NA with no CMT column, so a simulated binary
# outcome was indistinguishable from a PK row that failed to predict.
sim <- ferx_simulate(ex$model, ex$data, n_sim = 20L, seed = 1L, fit = fit)

bin <- sim[sim$CMT == 3, ]
stopifnot(
  all(bin$DV_SIM %in% c(0, 1)),   # simulated outcomes are 0/1, never NA
  !anyNA(bin$DV_SIM)
)

# Simulated event rate per record time should track the fitted probability.
cat("simulated P(DV = 1) by TIME:\n")
print(tapply(bin$DV_SIM, bin$TIME, mean))

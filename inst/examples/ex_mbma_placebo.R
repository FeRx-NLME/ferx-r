library(ferx)

# Model-based meta-analysis (MBMA) of arm-level summary data.
#
# Six synthetic studies, each with a placebo arm and 2-3 dose arms, report the
# arm mean of a symptom score (lower = better) at week 0 and three later
# visits. The model separates a study-by-visit placebo time course
# (PLACEBO[STUDY, TIME], one level per cell) from an Emax dose-response that
# builds up over time, with between-study variability on the baseline (omega)
# and between-arm variability (kappa).
#
# The data are synthetic and seeded: mbma_placebo/simulate_dataset.R (in the
# same examples folder) generates them in closed form with base R and shapes
# them with ferx_mbma_data(). Published MBMA data such as the naproxen dataset
# ferx-core fits is licensed CC BY-NC, which is why it is not bundled here.
ex <- ferx_example("mbma_placebo")

d <- utils::read.csv(ex$data)
head(d)
# 88 rows: one per arm and visit. ID and STUDY are the study, OCC the arm,
# DV the arm mean, SE its standard error and NARM the arm size.

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# -- The placebo level block --------------------------------------------------
# One row per (STUDY, TIME) cell in the data: 24 cells. With
# contrast = sum_to_zero_within the levels sum to zero within each study, so
# one level per study (6 in all) is derived as minus the sum of the others.
# Those have no theta of their own, so theta_name is NA; they are not
# separate estimates.
fit$theta_levels
sum(is.na(fit$theta_levels$theta_name))

# -- Truth versus estimate ----------------------------------------------------
# TVE0 is not the week-0 baseline. With sum_to_zero_within, TVE0 + ETA_E0 is
# each study's MEAN placebo level over its own visits, and the PLACEBO levels
# are deviations from it. The generator's nominal 50 is therefore not the
# truth to compare against. A study's mean placebo level is
# 50 + eta_s + (its mean placebo deviation), and the truth TVE0 estimates is the
# realised mean of that over the 6 studies, which simulate_dataset.R prints:
# 48.1949.
truth <- c(TVE0 = 48.1949, EMAX = 12, ED50 = 40, ET50 = 3)
est <- fit$theta[names(truth)]
se <- fit$se_theta[names(truth)]
data.frame(
  truth = truth,
  estimate = round(est, 3),
  lower95 = round(est - 1.96 * se, 3),
  upper95 = round(est + 1.96 * se, 3)
)

# -- Residual and between-arm variability -------------------------------------
# The error model is DV ~ additive(ADD_ERR) weight = SE with sigma FIXed to 1:
# inverse-variance weighting, where each arm mean counts 1 / SE^2. Sigma is
# not estimated. DV, PRED and the simulations below stay on the natural scale.
#
# KAPPA_ARM uses weight = NARM, so its estimate is the unweighted gamma^2 (the
# generator used 100): an arm of NARM patients has between-arm SD
# sqrt(gamma^2 / NARM). Read it as a variance. KAPPA_ARM is additive, so
# print(fit) shows its SD, sqrt(gamma^2), "at weight 1" - the SD of a
# one-patient arm - and not a CV%.
fit$kappa_param_types
gamma2 <- fit$omega_iov[1, 1]
gamma2
# The arm SD at the typical (median) arm size, as print(fit) shows it.
fit$kappa_weight_typical
sqrt(gamma2 / fit$kappa_weight_typical)
fit$shrinkage_kappa                    # a fraction, not a percent

# -- Visual predictive check --------------------------------------------------
# Simulate at the fit's own studies and visits. The placebo levels exist only
# for the (study, visit) cells the fit saw, so this is the grid to simulate on.
n_sim <- 200
sim <- ferx_simulate(ex$model, ex$data, n_sim = n_sim, fit = fit)
# Each replicate holds the data rows in file order, so the arm's dose can be
# carried across by position.
sim$DOSE <- rep(d$DOSE, times = n_sim)
sim$ARM <- ifelse(sim$DOSE == 0, "placebo", "active")
obs <- transform(d, ARM = ifelse(DOSE == 0, "placebo", "active"))
# Simulated 5th / 50th / 95th percentiles of the arm means, by visit ...
stats::aggregate(DV_SIM ~ ARM + TIME, sim, stats::quantile,
                 probs = c(0.05, 0.5, 0.95))
# ... and the observed medians to set beside them.
stats::aggregate(DV ~ ARM + TIME, obs, stats::median)

# A denser grid is refused, by design: a visit the fit never saw has no
# placebo level, and the engine will not invent one (it does not default the
# level to 0). The same holds for a new study.
dense <- d[d$STUDY == 1 & d$OCC == 1, ][rep(1, 3), ]
dense$TIME <- c(1, 3, 6)
dense_path <- tempfile(fileext = ".csv")
utils::write.csv(dense, dense_path, row.names = FALSE)
refused <- tryCatch(ferx_simulate(ex$model, dense_path, n_sim = 1, fit = fit),
                    error = function(e) e)
stopifnot(inherits(refused, "error"))
cat(conditionMessage(refused), "\n")

# To pool visits into windows (week 1-3 as one level, say), number the cells
# yourself and use the counted form theta PLACEBO[N] with PLACEBO[LEVEL_IDX]:
# see ?ferx_mbma_data, arguments `levels` and `index`.

# Population predictions at the fitted cells.
pred <- ferx_predict(ex$model, ex$data, fit = fit)
head(pred)

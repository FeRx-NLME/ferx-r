library(ferx)

# Nested variability levels below the subject (ferx-core #1082).
#
# 12 sponsors each ran 2 studies with 2 arms, and every arm reports its mean
# response at 5 visits. A sponsor is the subject (ID). Random effects sit at
# three levels:
#
#   ETA_SP  once per sponsor              (omega)
#   K_ST    once per study, in a sponsor  (kappa ... at STUDY)
#   K_ARM   once per arm, in a study      (kappa ... at ARM weight = NARM)
#
# The data are synthetic: ferx-core's NONMEM anchor nested_levels_3, MIT
# licensed. NONMEM 7.6.0 fits the same model, written as occasion-style ETA
# slots with $OMEGA SAME, to OFV 107.04357106 (FOCEI). The model's initial
# estimates are NONMEM's initial estimates, not its converged values, so the
# fit below has to find the optimum itself.
ex <- ferx_example("nested_levels")

d <- utils::read.csv(ex$data)
head(d)
# 240 rows: one per arm and visit. ARMID, STUDYID and SLOT are the NONMEM
# coding of the same design and are not read by the ferx model.

# -- The levels ----------------------------------------------------------------
# `level STUDY by STUDY in ID` makes a study a group of rows that share a STUDY
# value inside one sponsor; `level ARM by ARM in STUDY` does the same for an
# arm inside a study. The labels are local: every sponsor has studies 1 and 2,
# and every study has arms 1 and 2. They are nested by key, so the 12 x 2 = 24
# studies and 24 x 2 = 48 arms are told apart by (ID, STUDY) and
# (ID, STUDY, ARM), never by STUDY or ARM alone.
sort(unique(d$STUDY))
sort(unique(d$ARM))
c(sponsors = nrow(unique(d["ID"])),
  studies  = nrow(unique(d[c("ID", "STUDY")])),
  arms     = nrow(unique(d[c("ID", "STUDY", "ARM")])))

fit <- ferx_fit(ex$model, ex$data)
print(fit)
# The OMEGA_IOV rows name each kappa with its level: "K_ST at STUDY" and
# "K_ARM at ARM", in the estimates table and in SHRINKAGE alike.

# -- Against NONMEM ------------------------------------------------------------
# NONMEM's estimates (nested_levels_3.ext) and its OFV at that optimum
# (nested_levels_3_eval). Sigma is the residual variance; ferx reports an
# additive sigma as its SD, so it is squared here.
nonmem <- c(OFV = 107.04357106, TVE0 = 10.1818, TVEMAX = 22.0140,
            ED50 = 33.4048, ET50 = 3.20150, ETA_SP = 0.0699388,
            K_ST = 0.0904362, K_ARM = 3.75478, ADD = 0.255164)
ferx <- c(OFV = fit$ofv, fit$theta, ETA_SP = fit$omega[1, 1],
          diag(fit$omega_iov), ADD = fit$sigma[[1]]^2)
data.frame(nonmem = nonmem, ferx = signif(ferx[names(nonmem)], 9))

# -- Each kappa's level --------------------------------------------------------
fit$kappa_levels

# One row per group of the kappa's own level: 24 study estimates and 48 arm
# estimates. GROUP numbers the groups inside a sponsor in data order.
by_level <- fit$ebe_kappas_by_level
names(by_level)
head(by_level$STUDY)
head(by_level$ARM)
sapply(by_level, nrow)

# Both arms of a study share its K_ST: the sdtab carries each row's own study
# and arm estimate.
head(fit$sdtab[c("ID", "TIME", "DV", "IPRED", "K_ST", "K_ARM")], 10)

# -- The arm-size weight -------------------------------------------------------
# `weight = NARM` scales the between-arm variance by arm size. The K_ARM
# estimate is the unweighted gamma^2: an arm of NARM patients has between-arm
# SD sqrt(gamma^2 / NARM) on the log scale of EMAX.
gamma2 <- fit$omega_iov["K_ARM", "K_ARM"]
narm <- sort(unique(d$NARM))
data.frame(NARM = narm, between_arm_sd = round(sqrt(gamma2 / narm), 4))

# Shrinkage is pooled over each kappa's own groups: 24 studies, 48 arms.
fit$shrinkage_kappa

# -- Simulation ----------------------------------------------------------------
# A simulation draws ETA_SP once per sponsor, K_ST once per study and K_ARM
# once per arm, on the fitted design. The observed arm means at the last visit
# against the simulated 5th-95th percentiles:
sim <- ferx_simulate(ex$model, ex$data, n_sim = 200, seed = 1, fit = fit)
last <- d$TIME == max(d$TIME)
sim_last <- sim$DV_SIM[sim$TIME == max(d$TIME)]
rbind(observed  = stats::quantile(d$DV[last], c(0.05, 0.5, 0.95)),
      simulated = stats::quantile(sim_last, c(0.05, 0.5, 0.95)))

# -- Save and reload -----------------------------------------------------------
# The levels and the per-level estimates survive a .fitrx round trip.
f <- tempfile(fileext = ".fitrx")
ferx_save_fit(fit, f)
back <- ferx_load_fit(f)
identical(back$kappa_levels, fit$kappa_levels)
identical(back$ebe_kappas_by_level, fit$ebe_kappas_by_level)
unlink(f)

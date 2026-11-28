# Writes mixture_tte.csv: mixture_iv.csv plus, per subject, one competing-risks
# pair of TTE rows (CMT 2 = cause A, CMT 3 = cause B) at the first event time,
# administratively censored at 24. Run from tests/testthat/fixtures.
set.seed(675)
pk <- utils::read.csv("mixture_iv.csv")
ids <- unique(pk$ID)
t_a <- stats::rexp(length(ids), 0.05 * exp(stats::rnorm(length(ids), sd = 0.3)))
t_b <- stats::rexp(length(ids), 0.03)
t <- pmin(t_a, t_b, 24)
tte <- do.call(rbind, lapply(seq_along(ids), function(i) {
  data.frame(ID = ids[i], TIME = round(t[i], 4),
             DV = c(t_a[i] == t[i], t_b[i] == t[i]) * 1,
             AMT = 0, EVID = 0, CMT = 2:3, WT = pk$WT[pk$ID == ids[i]][1])
}))
out <- rbind(pk, tte)
out <- out[order(out$ID, out$TIME, out$CMT), ]
utils::write.csv(out, "mixture_tte.csv", row.names = FALSE, quote = FALSE)

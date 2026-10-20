- **New bundled example `nested_levels`**: random effects at three nested
  levels below the subject ([ferx-core #1082](https://github.com/FeRx-NLME/ferx-core/issues/1082),
  #PR). 12 sponsors (`ID`) each ran 2 studies of 2 arms, observed at 5 visits;
  `level STUDY by STUDY in ID` and `level ARM by ARM in STUDY` declare the
  levels, `kappa K_ST ... at STUDY` and `kappa K_ARM ... at ARM weight = NARM`
  put a kappa at each. `ex_nested_levels.R` fits it with FOCEI from NONMEM's
  initial estimates to NONMEM 7.6.0's optimum (OFV 107.04357), and shows
  `fit$kappa_levels`, `fit$ebe_kappas_by_level`, the arm-size weight, a
  simulation on the fitted design and a `.fitrx` round trip. The data are
  ferx-core's synthetic, MIT-licensed NONMEM anchor `nested_levels_3`.

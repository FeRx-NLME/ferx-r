# Tests for ferx_vpc() (#787): the continuous VPC / pcVPC adapter over
# vpc::vpc(). Each block names the mutation of R/ferx_vpc.R it must catch.

skip_if_not_installed("vpc")

# A fit of a bundled example, cached per (example, ID style, data edit).
# `ids = "text"` rewrites every subject ID to "S<id>", so each row is also
# checked on IDs that are not numbers (#786: an ID is the subject's text).
# `edit` is a function of the data frame (all columns character) returning
# the edited data; `tag` names it in the cache.
vpc_fit <- local({
  cache <- list()
  function(name, ids = c("numeric", "text"), edit = NULL, tag = "",
           bloq_method = NULL) {
    ids <- match.arg(ids)
    key <- paste(name, ids, tag, bloq_method %||% "", sep = "/")
    if (is.null(cache[[key]])) {
      ex <- ferx_example(name)
      data <- ex$data
      if (ids == "text" || !is.null(edit)) {
        d <- utils::read.csv(ex$data, colClasses = "character",
                             check.names = FALSE)
        if (ids == "text") d$ID <- paste0("S", d$ID)
        if (!is.null(edit)) d <- edit(d)
        data <- tempfile(fileext = ".csv")
        utils::write.csv(d, data, row.names = FALSE, quote = FALSE)
      }
      args <- list(ex$model, data, verbose = FALSE, covariance = FALSE,
                   settings = list(maxiter = 30L))
      if (!is.null(bloq_method)) args$bloq_method <- bloq_method
      cache[[key]] <<- suppressWarnings(do.call(ferx_fit, args))
    }
    cache[[key]]
  }
})

# The Gaussian rows of a fit's sdtab, as ferx_vpc() reads them.
vpc_obs <- function(fit) {
  fit$sdtab[.fitrx_sdtab_blocks(fit)$gaussian, , drop = FALSE]
}

vpc_err <- function(expr) {
  tryCatch({
    expr
    NULL
  }, error = function(e) e)
}

for (ids in c("numeric", "text")) {

  test_that(paste0("plain VPC: every replicate simulates exactly the fit's observations (", ids, " IDs)"), {
    # Mutation: the match replaced by a positional rep() of the sim rows.
    fit <- vpc_fit("warfarin", ids)
    obs <- vpc_obs(fit)
    db <- ferx_vpc(fit, n_sim = 20L, vpcdb = TRUE)
    expect_equal(nrow(db$obs), nrow(obs))
    expect_equal(nrow(db$sim), 20L * nrow(obs))
    expect_type(db$sim$id, "character")
    expect_false(anyNA(db$sim$id))
    for (r in c(1L, 7L, 20L)) {
      s <- db$sim[db$sim$sim == r, ]
      expect_identical(as.character(s$id), as.character(obs$ID))
      expect_equal(s$idv, obs$TIME)
    }
    if (ids == "text") expect_true(all(startsWith(db$obs$id, "S")))
  })

  test_that(paste0("data selection: the rows the fit ignored are not simulated (", ids, " IDs)"), {
    # Mutations: a positional join; the occurrence index dropped from the key
    # (see the duplicate-time block below for the latter).
    fit <- vpc_fit("warfarin_data_selection", ids)
    obs <- vpc_obs(fit)
    expect_equal(nrow(obs), 108L)
    db <- ferx_vpc(fit, n_sim = 10L, vpcdb = TRUE)
    expect_equal(nrow(db$sim), 10L * 108L)
    pre <- if (ids == "text") "S" else ""
    ignored <- paste0(pre, c("7", "8"), "|", 120)
    expect_false(any(paste(db$sim$id, db$sim$idv, sep = "|") %in% ignored))
    s <- db$sim[db$sim$sim == 3L, ]
    expect_identical(as.character(s$id), as.character(obs$ID))
    expect_equal(s$idv, obs$TIME)
  })

  test_that(paste0("repeated times pair by occurrence, each sim row used once (", ids, " IDs)"), {
    # Mutation: the occurrence index dropped from the key - both observations
    # at the repeated time then pair with the first simulated row.
    dup <- function(d) {
      i <- which(d$ID %in% c("1", "S1") & d$EVID == "0")[[3L]]
      extra <- d[i, ]
      extra$DV <- format(as.numeric(d$DV[i]) * 1.1)
      rbind(d[seq_len(i), ], extra, d[-seq_len(i), ])
    }
    fit <- vpc_fit("warfarin", ids, edit = dup, tag = "dup")
    obs <- vpc_obs(fit)
    k <- paste(obs$ID, obs$TIME)
    rows <- which(k %in% k[duplicated(k)])
    expect_length(rows, 2L)
    db <- ferx_vpc(fit, n_sim = 5L, vpcdb = TRUE)
    for (r in 1:5) {
      v <- db$sim$dv[db$sim$sim == r][rows]
      expect_false(isTRUE(all.equal(v[[1L]], v[[2L]])))
    }
  })

  test_that(paste0("pcVPC: each sim row carries its observation's PRED, and vpc corrects by the closed form (", ids, " IDs)"), {
    # Mutations: PRED attached in obs order to sim rows that are not; PRED
    # taken from IPRED.
    fit <- vpc_fit("warfarin", ids)
    obs <- vpc_obs(fit)
    key_obs <- paste(obs$ID, obs$TIME)
    for (lb in c(0, 0.5)) {
      db <- ferx_vpc(fit, n_sim = 20L, pred_corr = TRUE,
                     pred_corr_lower_bnd = lb, vpcdb = TRUE)
      m <- match(paste(db$sim$id, db$sim$idv), key_obs)
      expect_false(anyNA(m))
      expect_equal(db$sim$pred, obs$PRED[m])
      o <- db$obs
      i <- o$bin == 1
      expect_gt(sum(i), 5L)
      mo <- match(paste(o$id[i], o$idv[i]), key_obs)
      pred <- obs$PRED[mo]
      pred_bin <- stats::median(pred)
      hand <- lb + (obs$DV[mo] - lb) * (pred_bin - lb) / (pred - lb)
      expect_equal(o$dv[i], hand, tolerance = 1e-12)
    }
  })

  test_that(paste0("censoring: CENS = 1 observations are blanked and lloq comes from them (", ids, " IDs)"), {
    # Mutation: the censoring blank removed (the 10 rows come back at DV = 2).
    for (bm in list(NULL, "drop")) {
      fit <- vpc_fit("warfarin_bloq", ids, bloq_method = bm)
      obs <- vpc_obs(fit)
      cens <- obs$CENS %in% 1
      expect_equal(sum(cens), 10L)
      db <- ferx_vpc(fit, n_sim = 5L, vpcdb = TRUE)
      expect_equal(nrow(db$obs), nrow(obs))
      expect_true(all(is.na(db$obs$dv[cens])))
      expect_false(anyNA(db$obs$dv[!cens]))
      expect_equal(db$lloq, 2)
    }
  })

  test_that(paste0("pcVPC with a limit is refused before simulating, naming it (", ids, " IDs)"), {
    # Mutation: the refusal removed (vpc 1.2.4 then fails in quantile()).
    fit <- vpc_fit("warfarin_bloq", ids)
    local_mocked_bindings(.ferx_vpc_simulate = function(...) stop("simulated"))
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L, pred_corr = TRUE))
    expect_s3_class(e, "ferx_vpc_error")
    expect_match(conditionMessage(e), "lloq = 2 (from the fit's censored", fixed = TRUE)
    expect_match(conditionMessage(e), "pred_corr = FALSE", fixed = TRUE)
    e <- vpc_err(ferx_vpc(vpc_fit("warfarin", ids), n_sim = 5L,
                          pred_corr = TRUE, uloq = 20))
    expect_s3_class(e, "ferx_vpc_error")
    expect_match(conditionMessage(e), "(uloq = 20):", fixed = TRUE)
  })

  test_that(paste0("several LLOQ values are refused unless lloq is given (", ids, " IDs)"), {
    # Mutation: the default silently takes one of the values.
    two_lloq <- function(d) {
      i <- d$CENS == "1" & d$ID %in% unique(d$ID[d$CENS == "1"])[[1L]]
      d$DV[i] <- "1.5"
      d
    }
    fit <- vpc_fit("warfarin_bloq", ids, edit = two_lloq, tag = "two_lloq")
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L))
    expect_s3_class(e, "ferx_vpc_error")
    expect_match(conditionMessage(e), "values (1.5, 2),", fixed = TRUE)
    expect_match(conditionMessage(e), "single `lloq`", fixed = TRUE)
    expect_match(conditionMessage(e), "`lloq =`", fixed = TRUE)
    db <- ferx_vpc(fit, n_sim = 5L, lloq = 2, vpcdb = TRUE)
    expect_equal(db$lloq, 2)
    expect_true(all(is.na(db$obs$dv[vpc_obs(fit)$CENS %in% 1])))
  })

  test_that(paste0("ULOQ: CENS = -1 observations are blanked and uloq comes from them (", ids, " IDs)"), {
    # Mutation: CENS = -1 never marked as censored (review round 1, #1).
    above_only <- function(d) {
      d$CENS[d$CENS == "1"] <- "0"
      i <- which(d$EVID == "0" & d$CENS == "0")[c(2L, 3L, 13L)]
      d$CENS[i] <- "-1"
      d$DV[i] <- "15"
      d
    }
    fit <- vpc_fit("warfarin_bloq", ids, edit = above_only, tag = "uloq")
    obs <- vpc_obs(fit)
    above <- obs$CENS %in% -1
    expect_equal(sum(above), 3L)
    expect_false(any(obs$CENS %in% 1))
    db <- ferx_vpc(fit, n_sim = 5L, vpcdb = TRUE)
    expect_true(all(is.na(db$obs$dv[above])))
    expect_false(anyNA(db$obs$dv[!above]))
    expect_equal(db$uloq, 15)
    expect_null(db$lloq)
  })

  test_that(paste0("censored on both sides: refused before simulating; a given limit replaces both defaults (", ids, " IDs)"), {
    # Mutations: the both-sides refusal removed (vpc's own stop after the
    # simulation); a given lloq still defaulting uloq.
    both <- function(d) {
      i <- which(d$EVID == "0" & d$CENS == "0")[c(2L, 3L, 13L)]
      d$CENS[i] <- "-1"
      d$DV[i] <- "15"
      d
    }
    fit <- vpc_fit("warfarin_bloq", ids, edit = both, tag = "both")
    obs <- vpc_obs(fit)
    cens <- obs$CENS %in% c(1, -1)
    expect_equal(sum(cens), 13L)
    local({
      local_mocked_bindings(.ferx_vpc_simulate = function(...) stop("simulated"))
      e <- vpc_err(ferx_vpc(fit, n_sim = 5L))
      expect_s3_class(e, "ferx_vpc_error")
      expect_match(conditionMessage(e),
                   "(lloq = 2 from CENS = 1, uloq = 15 from CENS = -1)", fixed = TRUE)
      expect_match(conditionMessage(e), "Pass `lloq =` or `uloq =`", fixed = TRUE)
      e <- vpc_err(ferx_vpc(fit, n_sim = 5L, lloq = 2, uloq = 15))
      expect_s3_class(e, "ferx_vpc_error")
      expect_match(conditionMessage(e), "not both", fixed = TRUE)
    })
    db <- ferx_vpc(fit, n_sim = 5L, lloq = 2, vpcdb = TRUE)
    expect_equal(db$lloq, 2)
    expect_null(db$uloq)
    expect_true(all(is.na(db$obs$dv[cens])))
    expect_false(anyNA(db$obs$dv[!cens]))
  })

  test_that(paste0("several endpoints: cmt = NULL refused, cmt picks one in obs and sim (", ids, " IDs)"), {
    # Mutation: CMT dropped from the key (PK and PD share times, so CMT 3
    # observations then pair with CMT 2 simulations).
    fit <- vpc_fit("emax_pkpd", ids)
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L))
    expect_s3_class(e, "ferx_vpc_error")
    expect_match(conditionMessage(e), "CMT 2, 3", fixed = TRUE)
    expect_match(conditionMessage(e), "`cmt =`", fixed = TRUE)
    expect_no_match(conditionMessage(e), "invalid")
    for (c in c(2, 3)) {
      db <- ferx_vpc(fit, n_sim = 5L, cmt = c, vpcdb = TRUE)
      n <- sum(vpc_obs(fit)$CMT == c)
      expect_equal(nrow(db$obs), n)
      expect_equal(nrow(db$sim), 5L * n)
      expect_true(all(db$obs$CMT == c))
      expect_true(all(db$sim$CMT == c))
    }
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L, cmt = 7))
    expect_match(conditionMessage(e), "its continuous endpoints are CMT 2, 3.",
                 fixed = TRUE)
  })

  test_that(paste0("stratify = \"CMT\" with no cmt shows every continuous endpoint (", ids, " IDs)"), {
    # Mutation: the CMT filter applied regardless of stratify (review round 1, #3).
    fit <- vpc_fit("emax_pkpd", ids)
    obs <- vpc_obs(fit)
    db <- ferx_vpc(fit, n_sim = 5L, stratify = "CMT", vpcdb = TRUE)
    expect_equal(nrow(db$obs), nrow(obs))
    expect_equal(nrow(db$sim), 5L * nrow(obs))
    expect_equal(as.numeric(table(db$obs$CMT)[c("2", "3")]), c(15, 12))
    expect_setequal(as.character(unique(db$obs$strat)), c("2", "3"))
    expect_identical(as.character(db$sim$strat), rep(as.character(db$obs$strat), 5L))
  })

  test_that(paste0("joint PK-TTE: horizon is required and named; with it the continuous rows match (", ids, " IDs)"), {
    # Mutation: the E_SIMULATE_HORIZON rewrite removed (the engine's message
    # names `[simulation] horizon`, not ferx_vpc()'s argument).
    fit <- vpc_fit("pktte_joint", ids)
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L))
    expect_s3_class(e, "ferx_vpc_error")
    expect_identical(e$code, "E_SIMULATE_HORIZON")
    expect_match(conditionMessage(e), "`horizon =`", fixed = TRUE)
    expect_match(conditionMessage(e), "ferx_vpc()", fixed = TRUE)
    expect_no_match(conditionMessage(e), "ferx_vpc_tte", fixed = TRUE)
    db <- ferx_vpc(fit, n_sim = 5L, horizon = 200, vpcdb = TRUE)
    obs <- vpc_obs(fit)
    expect_equal(nrow(db$obs), 12L)
    expect_equal(nrow(db$sim), 5L * 12L)
    expect_true(all(db$sim$CMT == 2))
    # The TTE endpoint is refused as a `cmt`, pointing at ferx_vpc_tte().
    e <- vpc_err(ferx_vpc(fit, n_sim = 5L, horizon = 200, cmt = 3))
    expect_match(conditionMessage(e), "CMT 2", fixed = TRUE)
    expect_match(conditionMessage(e), "ferx_vpc_tte()", fixed = TRUE)
  })

  test_that(paste0("stratify: a data column reaches obs and sim, an MDV = 1 record included (", ids, " IDs)"), {
    # Mutation: the column attached to obs only.
    strat <- function(d) {
      d$WT_GROUP <- ifelse(as.integer(sub("^S", "", d$ID)) %% 2L == 0L,
                           "even", "odd")
      i <- which(d$ID %in% c("2", "S2") & d$EVID == "0")[[2L]]
      extra <- d[i, ]
      extra$TIME <- format(as.numeric(d$TIME[i]) + 0.25)
      extra$DV <- "."
      extra$MDV <- "1"
      extra$WT_GROUP <- "mdv"
      rbind(d[seq_len(i), ], extra, d[-seq_len(i), ])
    }
    fit <- vpc_fit("warfarin", ids, edit = strat, tag = "strat")
    db <- ferx_vpc(fit, n_sim = 5L, stratify = "WT_GROUP", vpcdb = TRUE)
    want <- ifelse(as.integer(sub("^S", "", db$obs$id)) %% 2L == 0L,
                   "even", "odd")
    expect_identical(db$obs$WT_GROUP, want)
    want_sim <- ifelse(as.integer(sub("^S", "", db$sim$id)) %% 2L == 0L,
                       "even", "odd")
    expect_identical(db$sim$WT_GROUP, want_sim)
    expect_false("mdv" %in% c(db$obs$WT_GROUP, db$sim$WT_GROUP))
    expect_setequal(unique(db$obs$strat), c("even", "odd"))
  })
}

test_that("a fit loaded from .fitrx gives the live fit's VPC data", {
  fit <- vpc_fit("warfarin")
  p <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, p)
  loaded <- suppressWarnings(ferx_load_fit(p))
  a <- ferx_vpc(fit, n_sim = 5L, pred_corr = TRUE, vpcdb = TRUE)
  b <- ferx_vpc(loaded, n_sim = 5L, pred_corr = TRUE, vpcdb = TRUE)
  expect_equal(as.data.frame(b$obs), as.data.frame(a$obs))
  expect_equal(as.data.frame(b$sim), as.data.frame(a$sim))
})

test_that("the default plot is vpc's ggplot", {
  p <- ferx_vpc(vpc_fit("warfarin"), n_sim = 5L)
  expect_true(inherits(p, "ggplot"))
})

# -- User-facing text, one row per input class --

test_that("vpc not installed: the message names the package and how to install it", {
  # Mutation: the vpc check removed.
  local_mocked_bindings(.ferx_has_vpc = function() FALSE)
  e <- vpc_err(ferx_vpc("not a fit"))
  expect_s3_class(e, "ferx_vpc_error")
  expect_match(conditionMessage(e), "needs the vpc package", fixed = TRUE)
  expect_match(conditionMessage(e), "install.packages(\"vpc\")", fixed = TRUE)
  expect_no_match(conditionMessage(e), "fit", fixed = TRUE)
})

test_that("not a ferx_fit: the message names ferx_fit() / ferx_load_fit()", {
  # Mutation: the inherits() check removed.
  e <- vpc_err(ferx_vpc(ferx_example("warfarin")$model))
  expect_s3_class(e, "ferx_vpc_error")
  expect_match(conditionMessage(e), "ferx_fit() or ferx_load_fit()", fixed = TRUE)
  expect_no_match(conditionMessage(e), "path")
})

test_that("no continuous rows: TTE and categorical VPCs are named", {
  # Mutation: the empty-Gaussian refusal removed.
  fit <- vpc_fit("tte_weibull")
  e <- vpc_err(ferx_vpc(fit, n_sim = 5L))
  expect_s3_class(e, "ferx_vpc_error")
  expect_match(conditionMessage(e), "ferx_vpc_tte() or ferx_vpc_cif()", fixed = TRUE)
  expect_match(conditionMessage(e), "vpc::vpc_cat()", fixed = TRUE)
  expect_no_match(conditionMessage(e), "empty")
})

test_that("an adapter-owned argument in ... is refused with its replacement", {
  # Mutation: the ... guard removed.
  fit <- vpc_fit("warfarin")
  e <- vpc_err(ferx_vpc(fit, obs = data.frame()))
  expect_s3_class(e, "ferx_vpc_error")
  expect_match(conditionMessage(e), "`obs`", fixed = TRUE)
  expect_match(conditionMessage(e), "`fit =`", fixed = TRUE)
  e <- vpc_err(ferx_vpc(fit, sim = data.frame()))
  expect_match(conditionMessage(e), "`n_sim =`", fixed = TRUE)
  e <- vpc_err(ferx_vpc(fit, sim_cols = list()))
  expect_match(conditionMessage(e), "`sim_cols`", fixed = TRUE)
  e <- vpc_err(ferx_vpc(fit, 5L, 1L, FALSE, NULL, NULL, NULL, NULL, NULL, 3))
  expect_match(conditionMessage(e), "must be named", fixed = TRUE)
})

test_that("a replicate that does not match is reported as a ferx bug with the counts and key", {
  # Mutation: the anyNA(idx) refusal removed (match returns NA rows).
  obs <- data.frame(ID = c("1", "1", "2"), TIME = c(1, 2, 1), CMT = 1)
  sim <- data.frame(SIM = rep(1:2, each = 3), ID = c("1", "1", "2", "1", "2", "2"),
                    TIME = c(1, 2, 1, 1, 1, 2), CMT = 1)
  e <- vpc_err(.ferx_vpc_match(obs, sim, 2L))
  expect_s3_class(e, "ferx_vpc_error")
  expect_match(conditionMessage(e), "replicate 2 matches 2 of the fit's 3 observations",
               fixed = TRUE)
  expect_match(conditionMessage(e), "(CMT, ID, TIME, occurrence)", fixed = TRUE)
  expect_match(conditionMessage(e), "ferx bug", fixed = TRUE)
  expect_no_match(conditionMessage(e), "change your data")
  # The same frames with replicate 2 complete match row for row.
  sim$ID[5:6] <- c("1", "2"); sim$TIME[5:6] <- c(2, 1)
  expect_identical(.ferx_vpc_match(obs, sim, 2L), 1:6)
})

test_that("a replicate holding the records in another order pairs by key, not by replicate 1's positions", {
  # Mutations: replicate 1's positions reused without the structure check; a
  # differing replicate refused instead of keyed on its own.
  obs <- data.frame(ID = c("1", "1", "2"), TIME = c(1, 2, 1), CMT = 1)
  sim <- data.frame(SIM = rep(1:2, each = 3), ID = c("1", "1", "2", "2", "1", "1"),
                    TIME = c(1, 2, 1, 1, 1, 2), CMT = 1)
  expect_identical(.ferx_vpc_match(obs, sim, 2L), c(1L, 2L, 3L, 5L, 6L, 4L))
})

test_that("the pairing equals a per-row key match on a randomized simulation", {
  # Reference: every row keyed by (SIM, CMT, ID, TIME, occurrence), as before
  # #796 finding 6. Replicates differ in length (dropped rows), order and
  # repeated times, so both branches of .ferx_vpc_match() run.
  set.seed(796)
  n <- 60L
  base <- data.frame(CMT = sample(2:3, n, TRUE), ID = as.character(sample(1:8, n, TRUE)),
                     TIME = sample(c(0.5, 1, 2, 4), n, TRUE))
  keep <- !duplicated(base) | seq_len(n) %% 7L == 0L
  base <- base[keep, ]
  obs <- base[sample(nrow(base), nrow(base) - 5L), ]
  obs <- obs[order(as.integer(rownames(obs))), ]
  sim <- do.call(rbind, lapply(1:6, function(r) {
    b <- base
    if (r %in% c(2L, 5L)) b <- b[sample(nrow(b)), ]
    if (r == 4L) b <- rbind(b, base[1:3, ])
    cbind(SIM = r, b)
  }))
  ref_key <- function(sim_, cmt, id, time) {
    k <- paste(sim_, cmt, id, sprintf("%.17g", time), sep = "|")
    paste(k, stats::ave(seq_along(k), k, FUN = seq_along), sep = "#")
  }
  k_sim <- ref_key(sim$SIM, sim$CMT, sim$ID, sim$TIME)
  k_obs <- unlist(lapply(1:6, function(r) ref_key(rep(r, nrow(obs)), obs$CMT, obs$ID, obs$TIME)))
  expect_identical(.ferx_vpc_match(obs, sim, 6L), match(k_obs, k_sim))
})

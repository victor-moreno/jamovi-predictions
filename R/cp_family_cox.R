# Time-to-event outcome (Cox proportional hazards) family. The outcome y is a
# two-column matrix (time, status 0/1). Predictions are risks of the event by
# a horizon h, and the "linear predictor" of this family is the log cumulative
# hazard at h: lp = log H0(h) + Xb, so risk = 1 − exp(−exp(lp)). log H0(h)
# (Efron's baseline estimator) is carried as the intercept of the coefficient
# vector; it does not enter the partial likelihood, and everything that only
# ranks the patients (concordance, D, slopes, SHAP values) ignores it.

# Distinct event times with the sums Efron's method needs, on a scale shifted
# by max(eta) for numerical safety
.cox_sets <- function(time, status, eta) {
    m <- max(eta)
    o <- order(time)
    t <- time[o]; s <- status[o]; e <- exp(eta[o] - m)
    rc <- rev(cumsum(rev(e)))
    ev <- s == 1
    ut <- unique(t[ev])
    sR <- rc[match(ut, t)]
    sD <- rowsum(e[ev], t[ev], reorder = TRUE)[, 1]
    d <- as.vector(table(factor(t[ev], levels = ut)))
    etaD <- rowsum(eta[o][ev], t[ev], reorder = TRUE)[, 1]
    list(ut = ut, sR = sR, sD = sD, d = d, etaD = etaD, shift = m)
}

# Expansion over tied events: d terms (SR − k/d SD), k = 0..d−1
.efron_terms <- function(rs) {
    g <- rep(seq_along(rs$d), rs$d)
    k <- sequence(rs$d) - 1
    list(g = g, den = rs$sR[g] - k / rs$d[g] * rs$sD[g])
}

# Partial log-likelihood (Efron ties), as coxph
.cox_ll <- function(y, lp) {
    rs <- .cox_sets(y[, 1], y[, 2], lp)
    if (!length(rs$d)) return(0)
    ef <- .efron_terms(rs)
    sum(rs$etaD) - sum(log(ef$den)) - sum(rs$d) * rs$shift
}

# Baseline cumulative hazard at h for linear predictor eta (Efron's estimator,
# as survfit.coxph after an Efron fit)
.cox_H0 <- function(time, status, eta, h) {
    rs <- .cox_sets(time, status, eta)
    if (!length(rs$d)) return(0)
    ef <- .efron_terms(rs)
    sum((1 / ef$den)[rs$ut[ef$g] <= h]) / exp(rs$shift)
}

cp_cox_outcome <- function(prep, timeName, statusName, eventLevel = "", horizon = 0) {
    tm <- prep$data[[timeName]]
    if (is.factor(tm)) tm <- suppressWarnings(as.numeric(as.character(tm)))
    tm <- as.numeric(tm)
    if (any(tm < 0)) .fail("The time variable '", timeName, "' has negative values.")
    st <- factor(prep$data[[statusName]], ordered = FALSE)
    st <- droplevels(st)
    lev <- levels(st)
    if (length(lev) < 2)
        .fail("The event variable '", statusName, "' has a single level in the complete cases.")
    if (is.null(eventLevel) || !nzchar(eventLevel)) eventLevel <- lev[2]
    if (!(eventLevel %in% lev))
        .fail("The event level '", eventLevel, "' does not occur in the complete cases.")
    status <- as.integer(st == eventLevel)
    if (sum(status) < 2) .fail("There are fewer than 2 events in the complete cases.")
    if (is.null(horizon) || !is.finite(horizon) || horizon <= 0) {
        horizon <- unname(stats::median(tm))
        auto <- TRUE
    } else auto <- FALSE
    if (!any(status == 1 & tm <= horizon))
        .fail("No event happens before the prediction horizon (", .num(horizon), "); choose a later horizon.")
    if (!any(tm > horizon))
        .fail("Nobody is followed beyond the prediction horizon (", .num(horizon), "); choose an earlier horizon.")
    prep$data[[timeName]] <- tm
    prep$spec$outcome <- list(name = statusName, time = timeName, levels = lev, event = eventLevel,
                              horizon = horizon, autoHorizon = auto)
    prep$y <- cbind(time = tm, status = status)
    prep
}

# Cox fit (Efron ties); coefficients prefixed by log H0(h), the vcov padded
# with a zero row and column for it
cp_cox_fit <- function(X, y, horizon) {
    w <- NULL
    k <- ncol(X)
    if (k == 0) {
        b <- numeric(0); V <- matrix(0, 0, 0); conv <- TRUE
    } else {
        ctrl <- survival::coxph.control()
        fit <- withCallingHandlers(
            survival::coxph.fit(X, survival::Surv(y[, 1], y[, 2]), strata = NULL, offset = NULL,
                                init = NULL, control = ctrl, weights = NULL, method = "efron",
                                rownames = NULL, resid = FALSE),
            warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") })
        b <- fit$coefficients
        V <- fit$var
        conv <- fit$iter < ctrl$iter.max && all(is.finite(b[!is.na(b)]))
    }
    bb <- b; bb[is.na(bb)] <- 0
    eta <- if (k) drop(X %*% bb) else rep(0, nrow(y))
    H0 <- .cox_H0(y[, 1], y[, 2], eta, horizon)
    beta <- c(`(Intercept)` = log(H0), stats::setNames(b, colnames(X)))
    attr(beta, "converged") <- conv && H0 > 0
    attr(beta, "separation") <- any(grepl("infinite", w))
    Vf <- matrix(0, k + 1, k + 1)
    if (k) Vf[-1, -1] <- V
    attr(beta, "vcov") <- Vf
    beta
}

# Kaplan–Meier survival at times `at` (left = TRUE: just before them), with
# Greenwood's variance
.km <- function(time, status, at, left = FALSE) {
    ut <- sort(unique(time[status == 1]))
    if (!length(ut)) return(list(surv = rep(1, length(at)), var = rep(0, length(at))))
    st <- sort(time)
    nrisk <- length(time) - findInterval(ut, st, left.open = TRUE)
    d <- tabulate(match(time[status == 1], ut), length(ut))
    S <- cumprod(1 - d / nrisk)
    gw <- cumsum(d / (nrisk * pmax(nrisk - d, 1)))
    idx <- findInterval(at, ut, left.open = left)
    list(surv = c(1, S)[idx + 1], var = (c(1, S)^2 * c(0, gw))[idx + 1])
}

# Inverse probability of censoring weights at horizon h: events by h weigh
# 1/G(T−), those still event-free at h weigh 1/G(h), the rest (censored
# before h) weigh 0
.ipcw <- function(time, status, h) {
    Gt <- .km(time, 1 - status, time, left = TRUE)$surv
    Gh <- .km(time, 1 - status, h)$surv
    list(case = ifelse(time <= h & status == 1, 1 / Gt, 0), ctrl = ifelse(time > h, 1 / Gh, 0))
}

# Cumulative/dynamic time-dependent AUC at h with IPCW (Uno 2007; Blanche 2013)
.auc_td <- function(r, time, status, h, wcase) {
    cs <- which(time <= h & status == 1); ct <- which(time > h)
    if (!length(cs) || !length(ct)) return(NA_real_)
    rc <- sort(r[ct])
    below <- findInterval(r[cs], rc, left.open = TRUE)
    eq <- findInterval(r[cs], rc) - below
    sum(wcase[cs] * (below + 0.5 * eq)) / (sum(wcase[cs]) * length(ct))
}

# One-predictor Cox slope (and its variance)
.cox_slope <- function(y, x) {
    f <- tryCatch(suppressWarnings(survival::coxph.fit(matrix(x), survival::Surv(y[, 1], y[, 2]), strata = NULL,
        offset = NULL, init = NULL, control = survival::coxph.control(), weights = NULL, method = "efron",
        rownames = NULL, resid = FALSE)), error = function(e) NULL)
    if (is.null(f) || !is.finite(f$coefficients[1])) return(c(NA_real_, NA_real_))
    c(unname(f$coefficients[1]), f$var[1, 1])
}

# Calibration at h (Austin, Harrell & van Klaveren 2020): a Cox model of the
# outcome on a restricted cubic spline (3 knots) of the cloglog predicted risk;
# its predicted risk at h is the "observed" risk. Returns a function giving it
# (with a pointwise CI) for any cloglog values.
.cox_calfit <- function(y, lp, h) {
    kn <- tryCatch(.rcs_knots(lp, 3, "lp"), error = function(e) NULL)
    basis <- function(x) if (is.null(kn)) matrix(x) else cbind(x, .rcs_basis(x, kn))
    B <- basis(lp)
    f <- tryCatch(suppressWarnings(survival::coxph.fit(B, survival::Surv(y[, 1], y[, 2]), strata = NULL,
        offset = NULL, init = NULL, control = survival::coxph.control(), weights = NULL, method = "efron",
        rownames = NULL, resid = FALSE)), error = function(e) NULL)
    if (is.null(f) || any(!is.finite(f$coefficients))) return(NULL)
    b <- f$coefficients
    H0 <- .cox_H0(y[, 1], y[, 2], drop(B %*% b), h)
    function(x, z = 0) {
        Bx <- basis(x)
        eta <- drop(Bx %*% b)
        se <- sqrt(pmax(0, rowSums((Bx %*% f$var) * Bx)))
        list(fit = 1 - exp(-H0 * exp(eta)), lo = 1 - exp(-H0 * exp(eta - z * se)),
             hi = 1 - exp(-H0 * exp(eta + z * se)))
    }
}

.cox_risk <- function(lp) 1 - exp(-exp(lp))

# R² measures from the partial likelihood (n = number of patients)
cp_cox_r2 <- function(y, lp) {
    n <- nrow(y)
    ll <- .cox_ll(y, lp); ll0 <- .cox_ll(y, rep(0, n))
    r2cs <- 1 - exp(-2 * (ll - ll0) / n)
    c(r2mf = 1 - ll / ll0, r2n = r2cs / (1 - exp(2 * ll0 / n)), r2cs = r2cs)
}

# Royston & Sauerbrei's D: Cox slope on the normal scores of the prognostic
# index, times kappa = sqrt(8/pi); R²_D = (D²/κ²) / (π²/6 + D²/κ²)
.royston_d <- function(y, lp) {
    n <- length(lp); kap <- sqrt(8 / pi)
    z <- stats::qnorm((rank(lp) - 3 / 8) / (n + 0.25)) / kap
    s <- .cox_slope(y, z)
    c(d = s[1], var = s[2])
}

# Every time-to-event measure at horizon h
cp_cox_metrics <- function(y, lp, h) {
    lp <- pmin(pmax(lp, -15), 5)
    time <- y[, 1]; status <- y[, 2]; n <- length(time)
    r <- .cox_risk(lp)
    S <- survival::Surv(time, status)
    cH <- survival::concordancefit(S, lp, reverse = TRUE, std.err = FALSE)$concordance
    cU <- survival::concordancefit(S, lp, reverse = TRUE, timewt = "n/G2", ymax = h, std.err = FALSE)$concordance
    w <- .ipcw(time, status, h)
    obs <- 1 - .km(time, status, h)$surv
    cal <- .cox_calfit(y, lp, h)
    err <- if (!is.null(cal)) abs(cal(lp)$fit - r) else rep(NA_real_, n)
    brier <- mean(w$case * (1 - r)^2 + w$ctrl * r^2)
    brier0 <- mean(w$case * (1 - obs)^2 + w$ctrl * obs^2)
    D <- .royston_d(y, lp)[["d"]]
    dk <- D^2 / (8 / pi)
    r2 <- cp_cox_r2(y, lp)
    c(c = unname(cH), uno = unname(cU),
      auc = .auc_td(r, time, status, h, w$case),
      d = D,
      oe = obs / mean(r),
      slope = .cox_slope(y, lp)[1],
      ici = mean(err), e50 = unname(stats::quantile(err, 0.5, na.rm = TRUE)),
      e90 = unname(stats::quantile(err, 0.9, na.rm = TRUE)),
      brier = brier, sbrier = 1 - brier / brier0,
      r2d = dk / (pi^2 / 6 + dk), r2n = r2[["r2n"]], r2cs = r2[["r2cs"]])
}

cp_cox_ci <- function(y, lp, h, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    out <- list()
    cf <- survival::concordancefit(survival::Surv(y[, 1], y[, 2]), lp, reverse = TRUE)
    out$c <- cf$concordance + c(-1, 1) * z * sqrt(cf$var[1])
    s <- .cox_slope(y, lp)
    out$slope <- s[1] + c(-1, 1) * z * sqrt(s[2])
    d <- .royston_d(y, lp)
    out$d <- d[["d"]] + c(-1, 1) * z * sqrt(d[["var"]])
    km <- .km(y[, 1], y[, 2], h)
    obs <- 1 - km$surv; r <- .cox_risk(lp)
    out$oe <- (obs / mean(r)) * exp(c(-1, 1) * z * sqrt(km$var) / obs)
    out
}

# Observed risk at h by Kaplan–Meier for the calibration plot: smoothed curve
# (with CI), risk groups (Greenwood CI) and a histogram of the predicted risks
cp_cox_calcurve <- function(y, p, h, groups = 10, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    lp <- log(-log(1 - p))
    grid <- seq(stats::quantile(p, 0.01), stats::quantile(p, 0.99), length.out = 100)
    cal <- .cox_calfit(y, lp, h)
    curve <- if (!is.null(cal)) {
        cc <- cal(log(-log(1 - grid)), z)
        data.frame(x = grid, y = cc$fit, lo = cc$lo, hi = cc$hi)
    } else data.frame(x = grid, y = NA_real_, lo = NA_real_, hi = NA_real_)
    curve <- curve[stats::complete.cases(curve$y), ]
    br <- unique(stats::quantile(p, seq(0, 1, length.out = groups + 1)))
    g <- cut(p, br, include.lowest = TRUE)
    grp <- do.call(rbind, lapply(split(seq_along(p), g), function(i) {
        km <- .km(y[i, 1], y[i, 2], h)
        o <- 1 - km$surv; se <- sqrt(km$var)
        data.frame(pred = mean(p[i]), obs = o, obs_lo = max(0, o - z * se), obs_hi = min(1, o + z * se),
                   n = length(i))
    }))
    edges <- seq(0, max(p), length.out = 51)
    cnt <- graphics::hist(p, breaks = edges, plot = FALSE)$counts
    list(curve = curve, groups = grp,
         hist = data.frame(xmin = utils::head(edges, -1), xmax = edges[-1], count = cnt))
}

# Net benefit at h (Vickers et al. 2008): the event risk among those above the
# threshold estimated by Kaplan–Meier, so censoring is accounted for
cp_net_benefit_surv <- function(y, p, thresholds, h) {
    n <- nrow(y)
    obsAll <- 1 - .km(y[, 1], y[, 2], h)$surv
    nb <- vapply(thresholds, function(t) {
        pos <- p >= t
        if (!any(pos)) return(0)
        ev <- 1 - .km(y[pos, 1], y[pos, 2], h)$surv
        mean(pos) * (ev - (1 - ev) * t / (1 - t))
    }, numeric(1))
    data.frame(threshold = thresholds, model = nb,
               all = obsAll - (1 - obsAll) * thresholds / (1 - thresholds), none = 0)
}

# Kaplan–Meier curves by groups of predicted risk (cut points from the
# development risks), with the mean predicted risk of each group at h
cp_km_groups <- function(y, p, cuts, h, maxPts = 200) {
    g <- cut(p, c(-Inf, cuts, Inf), labels = FALSE)
    ng <- length(cuts) + 1
    labs <- paste0(seq_len(ng), ifelse(seq_len(ng) == 1, " (lowest)", ifelse(seq_len(ng) == ng, " (highest)", "")))
    curves <- list(); pts <- list()
    for (k in sort(unique(g))) {
        i <- which(g == k)
        ut <- sort(unique(y[i, 1]))
        if (length(ut) > maxPts) ut <- ut[unique(round(seq(1, length(ut), length.out = maxPts)))]
        s <- .km(y[i, 1], y[i, 2], ut)$surv
        curves[[k]] <- data.frame(time = c(0, ut), surv = c(1, s), group = labs[k])
        pts[[k]] <- data.frame(time = h, surv = 1 - mean(p[i]), group = labs[k])
    }
    list(curves = do.call(rbind, curves), pred = do.call(rbind, pts))
}

# Minimum sample size for a time-to-event model (Riley et al. 2020, criteria
# as pmsampsize): shrinkage >= 0.9; optimism in Nagelkerke R² <= 0.05; and
# the precision of the overall risk at the horizon that this size gives
cp_cox_samplesize <- function(nparams, r2cs, rate, meanfup, horizon, shrinkage = 0.9, delta = 0.05) {
    nn <- 10000
    events <- ceiling(rate * meanfup * nn)
    lnull <- events * log(events / nn) - events
    r2max <- 1 - exp(2 * lnull / nn)
    n_for <- function(S) ceiling(nparams / ((S - 1) * log(1 - r2cs / S)))
    S2 <- r2cs / (r2cs + delta * r2max)
    n <- c(shrinkage = n_for(shrinkage), optimism = n_for(S2))
    py <- meanfup * max(n)
    cum <- 1 - exp(-rate * horizon)
    ci <- 1 - exp(-(rate + c(-1, 1) * 1.96 * sqrt(rate / py)) * horizon)
    list(n = n, r2max = r2max, S2 = S2, risk = cum, riskCI = ci)
}

.family_cox <- function(horizon = NULL) list(
    name = "cox",
    intercept = TRUE,
    showIntercept = FALSE,
    bounded = TRUE,
    horizon = horizon,
    fit = function(X, y) cp_cox_fit(X, y, horizon),
    vcov = function(beta) attr(beta, "vcov"),
    dfResid = function(beta) Inf,
    nullLoss = function(y) -.cox_ll(y, rep(0, nrow(y))),
    lp = cp_lp,
    metrics = function(y, lp) cp_cox_metrics(y, lp, horizon),
    ci = function(y, lp, conf) cp_cox_ci(y, lp, horizon, conf),
    risk = .cox_risk,
    smooth = function(y, p, newp) {
        cal <- .cox_calfit(y, log(-log(1 - p)), horizon)
        if (is.null(cal)) rep(NA_real_, length(newp)) else cal(log(-log(1 - newp)))$fit
    },
    r2 = function(y, lp, measure) cp_cox_r2(y, lp)[[measure]],
    loss = function(y, lp) -.cox_ll(y, lp),
    quick = function(y, lp) {
        r <- .cox_risk(pmin(pmax(lp, -15), 5)); w <- .ipcw(y[, 1], y[, 2], horizon)
        c(c = survival::concordancefit(survival::Surv(y[, 1], y[, 2]), lp, reverse = TRUE,
                                        std.err = FALSE)$concordance[[1]],
          brier = mean(w$case * (1 - r)^2 + w$ctrl * r^2))
    },
    permKeys = c("c", "brier"), permHigher = c(TRUE, FALSE),
    strata = function(y) y[, 2],
    exp = TRUE)

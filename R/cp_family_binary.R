# Binary outcome (logistic regression) family: outcome coding, fit, metrics.
# The validation engine and importance code call these through cp_family("binary").

# Outcome coding: y = 1 for the event level, 0 for every other level
cp_binary_outcome <- function(prep, eventLevel = "") {
    x <- factor(prep$data[[prep$spec$outcome$name]], ordered = FALSE)
    x <- droplevels(x)
    lev <- levels(x)
    if (length(lev) < 2)
        .fail(sprintf(.("The outcome '%s' has a single level in the complete cases."), prep$spec$outcome$name))
    if (is.null(eventLevel) || !nzchar(eventLevel)) eventLevel <- lev[2]
    if (!(eventLevel %in% lev))
        .fail(sprintf(.("The event level '%s' does not occur in the complete cases."), eventLevel))
    y <- as.integer(x == eventLevel)
    if (sum(y) < 2 || sum(1 - y) < 2)
        .fail(sprintf(.("The outcome needs at least 2 events and 2 non-events; there are %s and %s."), sum(y), sum(1 - y)))
    prep$spec$outcome$levels <- lev
    prep$spec$outcome$event <- eventLevel
    prep$y <- y
    prep
}

# Logistic fit by IRLS; returns coefficients (intercept first) and flags.
# Aliased columns get NA (treated as 0 when predicting).
cp_binary_fit <- function(X, y) {
    w <- NULL
    fit <- withCallingHandlers(
        stats::glm.fit(cbind(1, X), y, family = stats::binomial()),
        warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") })
    beta <- fit$coefficients
    names(beta) <- c("(Intercept)", colnames(X))
    attr(beta, "converged") <- isTRUE(fit$converged)
    attr(beta, "separation") <- any(grepl("fitted probabilities numerically 0 or 1", w))
    attr(beta, "fit") <- fit
    beta
}

cp_lp <- function(beta, X) {
    b <- beta; b[is.na(b)] <- 0
    drop(cbind(1, X) %*% b)
}

# Log-likelihood of y given risks p
.binary_ll <- function(y, p) {
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    sum(y * log(p) + (1 - y) * log(1 - p))
}

# AUC with DeLong's variance (placement values)
.auc_delong <- function(y, p) {
    p1 <- p[y == 1]; p0 <- p[y == 0]
    n1 <- length(p1); n0 <- length(p0)
    r <- rank(c(p1, p0))
    auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    # placement values: share of the other group ranked below (ties count half)
    r1 <- rank(p1); r0 <- rank(p0)
    v10 <- (r[seq_len(n1)] - r1) / n0
    v01 <- 1 - (r[n1 + seq_len(n0)] - r0) / n1
    se <- sqrt(stats::var(v10) / n1 + stats::var(v01) / n0)
    list(auc = auc, se = se)
}

# Smoothed calibration curve: loess of y on p (Austin & Steyerberg 2019)
.cal_smooth <- function(y, p, newp = p) {
    # approximate statistics: same curve, ~15x faster (matters in the bootstrap)
    ctrl <- stats::loess.control(statistics = "approximate", trace.hat = "approximate")
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, degree = 2, span = 0.75, control = ctrl)),
                    error = function(e) NULL)
    sm <- if (!is.null(fit)) suppressWarnings(stats::predict(fit, newdata = data.frame(p = newp))) else NULL
    if (is.null(sm) || anyNA(sm)) {
        lw <- stats::lowess(p, y, iter = 0)
        sm <- stats::approx(lw$x, lw$y, xout = newp, ties = mean, rule = 2)$y
    }
    pmin(pmax(sm, 0), 1)
}

# Calibration intercept (in the large, slope fixed at 1) and slope
.cal_glm <- function(y, lp) {
    f0 <- suppressWarnings(stats::glm.fit(matrix(1, length(y), 1), y, offset = lp,
                                          family = stats::binomial()))
    f1 <- suppressWarnings(stats::glm.fit(cbind(1, lp), y, family = stats::binomial()))
    # a fit that diverges (separated bootstrap models) gives no estimate
    list(citl = if (f0$converged) f0$coefficients[1] else NA_real_, citl_fit = f0,
         slope = if (f1$converged) f1$coefficients[2] else NA_real_, slope_fit = f1)
}

# All binary performance measures for outcome y and linear predictor lp.
# Fast enough to run on every bootstrap sample.
cp_binary_metrics <- function(y, lp) {
    # risks closer to 0 or 1 than 3e-7 add only numerical noise (overfitted
    # bootstrap models can produce them)
    lp <- pmin(pmax(lp, -15), 15)
    p <- stats::plogis(lp)
    n <- length(y); ybar <- mean(y)
    auc <- .auc_delong(y, p)$auc
    cg <- .cal_glm(y, lp)
    sm <- .cal_smooth(y, p)
    err <- abs(sm - p)
    ll <- .binary_ll(y, p)
    ll0 <- n * (ybar * log(ybar) + (1 - ybar) * log(1 - ybar))
    lr <- 2 * (ll - ll0)
    r2cs <- 1 - exp(-lr / n)
    r2max <- 1 - exp(2 * ll0 / n)
    brier <- mean((y - p)^2)
    c(auc = auc,
      dslope = mean(p[y == 1]) - mean(p[y == 0]),
      citl = unname(cg$citl),
      slope = unname(cg$slope),
      oe = sum(y) / sum(p),
      ici = mean(err),
      e50 = unname(stats::quantile(err, 0.5)),
      e90 = unname(stats::quantile(err, 0.9)),
      emax = max(err),
      brier = brier,
      sbrier = 1 - brier / (ybar * (1 - ybar)),
      r2n = r2cs / r2max,
      r2cs = r2cs,
      r2mf = 1 - ll / ll0)
}

# Confidence intervals where an analytic one exists; others stay NA
cp_binary_ci <- function(y, lp, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    p <- stats::plogis(lp)
    out <- list()
    a <- .auc_delong(y, p)
    # logit-transformed interval stays inside (0, 1)
    if (a$auc > 0 && a$auc < 1) {
        la <- stats::qlogis(a$auc); sel <- a$se / (a$auc * (1 - a$auc))
        out$auc <- stats::plogis(la + c(-1, 1) * z * sel)
    }
    cg <- .cal_glm(y, lp)
    wald <- function(fit, k) fit$coefficients[k] + c(-1, 1) * z * sqrt(diag(.glm_vcov(fit)))[k]
    out$citl <- unname(wald(cg$citl_fit, 1))
    out$slope <- unname(wald(cg$slope_fit, 2))
    O <- sum(y); E <- sum(p)
    out$oe <- (O / E) * exp(c(-1, 1) * z * sqrt((1 - mean(y)) / O))
    b <- (y - p)^2
    out$brier <- mean(b) + c(-1, 1) * z * stats::sd(b) / sqrt(length(b))
    out
}

# Metric rows: key, label, and whether "apparent" equals a definition (no CI)
cp_metric_defs <- function(family) {
    switch(family,
        binary = data.frame(
            key = c("auc", "dslope", "citl", "slope", "oe", "ici", "e50", "e90", "emax",
                    "brier", "sbrier", "r2n", "r2cs"),
            label = c(.("C-statistic (AUC)"), .("Discrimination slope"),
                      .("Calibration-in-the-large"), .("Calibration slope"), .("O/E ratio"),
                      .("ICI"), .("E50"), .("E90"), .("Emax"),
                      .("Brier score"), .("Scaled Brier score"), .("Nagelkerke R²"), .("Cox–Snell R²")),
            group = c(.("Discrimination"), .("Discrimination"),
                      .("Calibration"), .("Calibration"), .("Calibration"), .("Calibration"),
                      .("Calibration"), .("Calibration"), .("Calibration"),
                      .("Overall"), .("Overall"), .("Overall"), .("Overall")),
            stringsAsFactors = FALSE),
        cox = data.frame(
            key = c("c", "uno", "auc", "d", "oe", "slope", "ici", "e50", "e90",
                    "brier", "sbrier", "r2d", "r2n", "r2cs"),
            label = c(.("Harrell's C"), .("Uno's C (to the horizon)"), .("Time-dependent AUC (at the horizon)"),
                      .("Royston's D"), .("O/E ratio (at the horizon)"), .("Calibration slope"),
                      .("ICI"), .("E50"), .("E90"), .("Brier score (at the horizon)"),
                      .("Scaled Brier score (IPA)"), .("Royston's R²_D"), .("Nagelkerke R²"), .("Cox–Snell R²")),
            group = c(rep(.("Discrimination"), 4), rep(.("Calibration"), 5), rep(.("Overall"), 5)),
            stringsAsFactors = FALSE),
        gaussian = data.frame(
            key = c("r2", "rmse", "mae", "citl", "slope"),
            label = c(.("R²"), .("RMSE"), .("MAE"), .("Calibration-in-the-large"), .("Calibration slope")),
            group = c(.("Overall"), .("Overall"), .("Overall"), .("Calibration"), .("Calibration")),
            stringsAsFactors = FALSE),
        .fail(.("Unknown model family.")))
}

# Classification at a risk threshold, with Wilson intervals
cp_binary_classify <- function(y, p, threshold, conf = 0.95) {
    pos <- p >= threshold
    tp <- sum(pos & y == 1); fp <- sum(pos & y == 0)
    fn <- sum(!pos & y == 1); tn <- sum(!pos & y == 0)
    wilson <- function(x, n) {
        if (n == 0) return(c(NA_real_, NA_real_, NA_real_))
        z <- stats::qnorm(1 - (1 - conf) / 2)
        ph <- x / n
        centre <- (ph + z^2 / (2 * n)) / (1 + z^2 / n)
        half <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
        c(ph, centre - half, centre + half)
    }
    list(tp = tp, fp = fp, fn = fn, tn = tn,
         sens = wilson(tp, tp + fn), spec = wilson(tn, tn + fp),
         ppv = wilson(tp, tp + fp), npv = wilson(tn, tn + fn),
         acc = wilson(tp + tn, length(y)))
}

# Data-driven thresholds
cp_binary_threshold <- function(y, p, method, value) {
    switch(method,
        value = value,
        prevalence = mean(y),
        youden = {
            # sensitivity and specificity at every distinct threshold from cumulative counts
            # (p >= t is positive), linear after sorting
            u <- sort(unique(p))
            g <- match(p, u)
            tp <- rev(cumsum(rev(tabulate(g[y == 1], length(u)))))
            fp <- rev(cumsum(rev(tabulate(g[y == 0], length(u)))))
            j <- tp / sum(y == 1) + (1 - fp / sum(y == 0)) - 1
            u[which.max(j)]
        })
}

# Net benefit of the model and of treating everyone (Vickers & Elkin 2006)
cp_net_benefit <- function(y, p, thresholds) {
    n <- length(y); prev <- mean(y)
    nb <- vapply(thresholds, function(t) {
        pos <- p >= t
        sum(pos & y == 1) / n - sum(pos & y == 0) / n * t / (1 - t)
    }, numeric(1))
    data.frame(threshold = thresholds, model = nb,
               all = prev - (1 - prev) * thresholds / (1 - thresholds), none = 0)
}

# Minimum sample size for developing a binary-outcome model (Riley et al. 2020):
# shrinkage >= 0.9, optimism in Nagelkerke R2 <= 0.05, intercept within +-0.05
cp_binary_samplesize <- function(nparams, prevalence, r2cs, shrinkage = 0.9, delta = 0.05,
                                 margin = 0.05) {
    lnull <- prevalence * log(prevalence) + (1 - prevalence) * log(1 - prevalence)
    r2max <- 1 - exp(2 * lnull)
    n_for <- function(S) nparams / ((S - 1) * log(1 - r2cs / S))
    S2 <- r2cs / (r2cs + delta * r2max)
    n <- c(shrinkage = n_for(shrinkage),
           optimism = n_for(S2),
           intercept = (stats::qnorm(0.975) / margin)^2 * prevalence * (1 - prevalence))
    list(n = ceiling(n), S2 = S2, r2max = r2max)
}

# The family object used by the generic engine
cp_family <- function(family, horizon = NULL) {
    switch(family,
        binary = list(
            name = "binary",
            intercept = TRUE,
            fit = cp_binary_fit,
            vcov = function(beta) .glm_vcov(attr(beta, "fit")),
            nullLoss = function(y) { p <- mean(y); -length(y) * (p * log(p) + (1 - p) * log(1 - p)) },
            lp = cp_lp,
            metrics = cp_binary_metrics,
            risk = stats::plogis,
            # value function for Shapley decomposition: McFadden or Nagelkerke R2
            r2 = function(y, lp, measure) {
                m <- cp_binary_metrics_r2(y, lp)
                m[[measure]]
            },
            loss = function(y, lp) -.binary_ll(y, stats::plogis(lp)),
            # the measures permutation importance needs, without the loess smoother
            quick = function(y, lp) {
                p <- stats::plogis(lp)
                c(auc = .auc_delong(y, p)$auc, brier = mean((y - p)^2))
            },
            strata = function(y) y,
            # hooks shared with the other families (see cp_family_gaussian.R)
            ci = cp_binary_ci,
            smooth = .cal_smooth,
            dfResid = function(beta) Inf,
            permKeys = c("auc", "brier"), permHigher = c(TRUE, FALSE),
            exp = TRUE,
            bounded = TRUE),
        gaussian = .family_gaussian(),
        cox = .family_cox(horizon),
        .fail(.("Unknown model family.")))
}

# Cheaper R2 only (for the Shapley enumeration)
cp_binary_metrics_r2 <- function(y, lp) {
    p <- stats::plogis(lp)
    n <- length(y); ybar <- mean(y)
    ll <- .binary_ll(y, p)
    ll0 <- n * (ybar * log(ybar) + (1 - ybar) * log(1 - ybar))
    r2cs <- 1 - exp(-2 * (ll - ll0) / n)
    c(r2mf = 1 - ll / ll0, r2n = r2cs / (1 - exp(2 * ll0 / n)), r2cs = r2cs,
      auc = .auc_delong(y, p)$auc)
}

# Data for the calibration plot: smoothed curve with a pointwise CI, grouped
# observed proportions (Wilson CI) and a histogram of the predicted risks
cp_binary_calcurve <- function(y, p, groups = 10, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    grid <- seq(stats::quantile(p, 0.01), stats::quantile(p, 0.99), length.out = 100)
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, degree = 2, span = 0.75)),
                    error = function(e) NULL)
    curve <- if (!is.null(fit)) {
        pr <- suppressWarnings(stats::predict(fit, newdata = data.frame(p = grid), se = TRUE))
        data.frame(x = grid, y = pr$fit, lo = pr$fit - z * pr$se.fit, hi = pr$fit + z * pr$se.fit)
    } else {
        data.frame(x = grid, y = .cal_smooth(y, p, grid), lo = NA_real_, hi = NA_real_)
    }
    curve[, c("y", "lo", "hi")] <- lapply(curve[, c("y", "lo", "hi")], function(v) pmin(pmax(v, 0), 1))
    curve <- curve[stats::complete.cases(curve$y), ]
    br <- unique(stats::quantile(p, seq(0, 1, length.out = groups + 1)))
    g <- cut(p, br, include.lowest = TRUE)
    grp <- do.call(rbind, lapply(split(seq_along(p), g), function(i) {
        n <- length(i); x <- sum(y[i]); ph <- x / n
        centre <- (ph + z^2 / (2 * n)) / (1 + z^2 / n)
        half <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
        data.frame(pred = mean(p[i]), obs = ph, obs_lo = centre - half, obs_hi = centre + half, n = n)
    }))
    edges <- seq(0, max(p), length.out = 51)
    cnt <- graphics::hist(p, breaks = edges, plot = FALSE)$counts
    list(curve = curve, groups = grp,
         hist = data.frame(xmin = utils::head(edges, -1), xmax = edges[-1], count = cnt))
}

# Covariance of glm.fit coefficients (dispersion 1), NA for aliased columns
.glm_vcov <- function(fit) {
    k <- length(fit$coefficients); r <- fit$rank
    V <- matrix(NA_real_, k, k)
    piv <- fit$qr$pivot[seq_len(r)]
    V[piv, piv] <- chol2inv(fit$qr$qr[seq_len(r), seq_len(r), drop = FALSE])
    V
}

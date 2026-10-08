# Continuous outcome (linear regression) family: outcome, fit, metrics,
# calibration curve, sample size. Same interface as cp_family("binary").

cp_gaussian_outcome <- function(prep) {
    nm <- prep$spec$outcome$name
    y <- prep$data[[nm]]
    if (is.factor(y)) y <- suppressWarnings(as.numeric(as.character(y)))
    y <- as.numeric(y)
    if (anyNA(y)) .fail("The outcome '", nm, "' has values that are not numbers.")
    if (length(y) < 3 || stats::sd(y) == 0)
        .fail("The outcome '", nm, "' does not vary in the complete cases.")
    prep$data[[nm]] <- y
    prep$spec$outcome$mean <- mean(y)
    prep$spec$outcome$sd <- stats::sd(y)
    prep$y <- y
    prep
}

# Least squares; coefficients with the intercept first, NA for aliased columns
cp_gaussian_fit <- function(X, y) {
    fit <- stats::lm.fit(cbind(1, X), y)
    beta <- fit$coefficients
    names(beta) <- c("(Intercept)", colnames(X))
    attr(beta, "converged") <- TRUE
    attr(beta, "fit") <- fit
    beta
}

# Covariance of the coefficients (residual variance with n − p − 1 df)
.lm_vcov <- function(fit) {
    k <- length(fit$coefficients); r <- fit$rank
    s2 <- sum(fit$residuals^2) / fit$df.residual
    V <- matrix(NA_real_, k, k)
    piv <- fit$qr$pivot[seq_len(r)]
    V[piv, piv] <- s2 * chol2inv(fit$qr$qr[seq_len(r), seq_len(r), drop = FALSE])
    V
}

# Profile negative log-likelihood (σ at its maximum likelihood value), so that
# 2 × the difference of two fits is the likelihood ratio χ², n·log(SSE0/SSE1)
.gaussian_loss <- function(y, lp) length(y) / 2 * log(mean((y - lp)^2))

# Smoothed calibration: loess of the observed on the predicted outcome
.cal_smooth_cont <- function(y, p, newp = p) {
    ctrl <- stats::loess.control(statistics = "approximate", trace.hat = "approximate")
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, degree = 2, span = 0.75, control = ctrl)),
                    error = function(e) NULL)
    sm <- if (!is.null(fit)) suppressWarnings(stats::predict(fit, newdata = data.frame(p = newp))) else NULL
    if (is.null(sm) || anyNA(sm)) {
        lw <- stats::lowess(p, y, iter = 0)
        sm <- stats::approx(lw$x, lw$y, xout = newp, ties = mean, rule = 2)$y
    }
    sm
}

# Performance of predictions lp of outcome y. R² is 1 − SSE / SST around the
# mean of these data, so it is negative when the predictions do worse than
# that mean; calibration-in-the-large is the mean residual (observed −
# predicted), the slope that of the regression of y on the predictions.
cp_gaussian_metrics <- function(y, lp) {
    res <- y - lp
    sse <- sum(res^2); sst <- sum((y - mean(y))^2)
    slope <- if (stats::sd(lp) > 0) stats::cov(y, lp) / stats::var(lp) else NA_real_
    c(r2 = 1 - sse / sst,
      rmse = sqrt(mean(res^2)),
      mae = mean(abs(res)),
      citl = mean(res),
      slope = slope)
}

cp_gaussian_ci <- function(y, lp, conf = 0.95) {
    n <- length(y); res <- y - lp
    tq <- stats::qt(1 - (1 - conf) / 2, n - 1)
    out <- list()
    out$citl <- mean(res) + c(-1, 1) * tq * stats::sd(res) / sqrt(n)
    out$mae <- mean(abs(res)) + c(-1, 1) * tq * stats::sd(abs(res)) / sqrt(n)
    f <- stats::lm.fit(cbind(1, lp), y)
    se <- sqrt(diag(.lm_vcov(f)))[2]
    out$slope <- unname(f$coefficients[2] + c(-1, 1) * stats::qt(1 - (1 - conf) / 2, n - 2) * se)
    out
}

# Calibration plot data: observed against predicted, loess with a pointwise CI
cp_gaussian_calcurve <- function(y, p, conf = 0.95, maxPts = 2000) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    grid <- seq(stats::quantile(p, 0.01), stats::quantile(p, 0.99), length.out = 100)
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, degree = 2, span = 0.75)), error = function(e) NULL)
    curve <- if (!is.null(fit)) {
        pr <- suppressWarnings(stats::predict(fit, newdata = data.frame(p = grid), se = TRUE))
        data.frame(x = grid, y = pr$fit, lo = pr$fit - z * pr$se.fit, hi = pr$fit + z * pr$se.fit)
    } else data.frame(x = grid, y = .cal_smooth_cont(y, p, grid), lo = NA_real_, hi = NA_real_)
    curve <- curve[stats::complete.cases(curve$y), ]
    # points thinned deterministically (results state is size-limited)
    i <- if (length(y) > maxPts) round(seq(1, length(y), length.out = maxPts)) else seq_along(y)
    list(curve = curve, points = data.frame(pred = signif(p[i], 5), obs = signif(y[i], 5)))
}

# Minimum sample size for developing a model of a continuous outcome (Riley et
# al. 2019, Stat Med 38:1262; criteria as in pmsampsize): (1) expected
# shrinkage >= S, (2) optimism of the apparent R² <= delta, (3) residual SD
# within a multiplicative margin mmoe, (4) mean outcome (intercept) within mmoe
cp_gaussian_samplesize <- function(nparams, r2adj, intercept, sd, shrinkage = 0.9, delta = 0.05,
                                   mmoe = 1.1) {
    p <- nparams
    S <- function(n) 1 + (p - 2) / (n * log(1 - (r2adj * (n - p - 1) + p) / (n - 1)))
    n1 <- p + 2
    while (!(is.finite(S(n1)) && S(n1) >= shrinkage)) n1 <- n1 + 1
    n2 <- ceiling(1 + p * (1 - r2adj) / delta)
    # residual variance: the smallest df whose 95% interval of σ is within mmoe
    df3 <- 1
    while (max(sqrt(df3 / stats::qchisq(0.025, df3)), sqrt(stats::qchisq(0.975, df3) / df3)) > mmoe)
        df3 <- df3 + 1
    n3 <- df3 + p + 1
    n4 <- max(n1, n2, n3)
    if (is.finite(intercept) && intercept != 0) {
        half <- function(n) stats::qt(0.975, n - p - 1) * sqrt(sd^2 * (1 - r2adj) / n)
        while ((abs(intercept) + half(n4)) / abs(intercept) > mmoe) n4 <- n4 + 1
    }
    list(n = c(shrinkage = n1, optimism = n2, residual = n3, intercept = n4))
}

# Ordinary least squares family object (see cp_family)
.family_gaussian <- function() list(
    name = "gaussian",
    intercept = TRUE,
    fit = cp_gaussian_fit,
    vcov = function(beta) .lm_vcov(attr(beta, "fit")),
    # a penalised model has no lm fit (and no standard errors): normal quantiles
    dfResid = function(beta) attr(beta, "fit")$df.residual %||% Inf,
    nullLoss = function(y) .gaussian_loss(y, mean(y)),
    lp = cp_lp,
    metrics = cp_gaussian_metrics,
    ci = cp_gaussian_ci,
    risk = identity,
    smooth = .cal_smooth_cont,
    r2 = function(y, lp, measure) 1 - sum((y - lp)^2) / sum((y - mean(y))^2),
    loss = .gaussian_loss,
    quick = function(y, lp) c(r2 = 1 - sum((y - lp)^2) / sum((y - mean(y))^2),
                              rmse = sqrt(mean((y - lp)^2))),
    permKeys = c("r2", "rmse"), permHigher = c(TRUE, FALSE),
    # cross-validation folds balanced over the quintiles of the outcome
    strata = function(y) findInterval(y, unique(stats::quantile(y, c(0.2, 0.4, 0.6, 0.8)))),
    exp = FALSE)

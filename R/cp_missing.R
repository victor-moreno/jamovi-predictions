# Missing values: exploration, the missing-category (indicator) method, and
# imputation by chained equations with predictive mean matching (pure R).

# ---- exploration ---------------------------------------------------------

# Per variable: number and share missing, and whether missingness goes with
# the outcome: the share missing among the two outcome groups (binary; Cox:
# censored and event) with a chi-square test, or the outcome mean among rows
# observed and missing (continuous) with Welch's t test. group: 0/1 outcome
# or a numeric outcome; vars: the variables to describe.
cp_miss_explore <- function(data, vars, group, numericOutcome = FALSE) {
    rows <- lapply(vars, function(v) {
        m <- is.na(data[[v]])
        out <- list(var = v, n = sum(m), pct = mean(m), a = NA_real_, b = NA_real_, p = NA_real_)
        ok <- !is.na(group)
        if (sum(m) > 0 && sum(!m) > 0) {
            if (numericOutcome) {
                out$a <- mean(group[ok & !m]); out$b <- mean(group[ok & m])
                if (sum(ok & m) > 1 && sum(ok & !m) > 1)
                    out$p <- tryCatch(stats::t.test(group[ok & m], group[ok & !m])$p.value, error = function(e) NA_real_)
            } else {
                g <- group[ok]; mm <- m[ok]
                out$a <- mean(mm[g == 0]); out$b <- mean(mm[g == 1])
                tab <- table(factor(g, 0:1), factor(mm, c(FALSE, TRUE)))
                if (all(dim(tab) == 2) && all(rowSums(tab) > 0))
                    out$p <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value)
            }
        }
        out
    })
    do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# The most frequent combinations of missing variables
cp_miss_patterns <- function(data, vars, top = 10) {
    M <- is.na(data[, vars, drop = FALSE])
    # "(complete)": no variable missing (table() cannot index an empty name)
    key <- apply(M, 1, function(r) if (any(r)) paste(vars[r], collapse = ", ") else "(complete)")
    tab <- sort(table(key), decreasing = TRUE)
    k <- utils::head(names(tab), top)
    data.frame(pattern = k, n = as.integer(tab[k]), pct = as.numeric(tab[k]) / nrow(data),
               nvars = vapply(k, function(x) if (x != "(complete)") length(strsplit(x, ", ")[[1]]) else 0L, 1L),
               stringsAsFactors = FALSE)
}

# ---- missing category / indicator ---------------------------------------

# Learnt on the development data: factors with missing values get a
# "(missing)" level; numeric predictors with missing values are filled with
# their development median and get an indicator predictor "<name> (missing)".
cp_miss_category <- function(data, numVars, facVars) {
    num <- Filter(function(v) anyNA(data[[v]]), numVars)
    fac <- Filter(function(v) anyNA(data[[v]]), facVars)
    list(num = num, fac = fac, median = vapply(num, function(v) stats::median(data[[v]], na.rm = TRUE), 0),
         indicator = stats::setNames(paste(num, "(missing)"), num))
}

# Apply it to any data with the same predictor names (development or validation)
cp_miss_category_apply <- function(h, data) {
    for (v in h$num) {
        if (is.null(data[[v]])) next
        x <- data[[v]]
        if (is.factor(x)) x <- suppressWarnings(as.numeric(as.character(x)))
        data[[h$indicator[[v]]]] <- as.numeric(is.na(x))
        x[is.na(x)] <- h$median[[v]]
        data[[v]] <- x
    }
    for (v in h$fac) {
        if (is.null(data[[v]])) next
        x <- as.character(data[[v]])
        x[is.na(x)] <- "(missing)"
        lev <- setdiff(if (is.factor(data[[v]])) levels(data[[v]]) else sort(unique(x)), "(missing)")
        data[[v]] <- factor(x, c(lev, "(missing)"))
    }
    data
}

# ---- imputation by chained equations ----------------------------------

# Draw from N(mu, V)
.rmvn <- function(mu, V) {
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) {
        ev <- eigen(V, symmetric = TRUE)
        R <- t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(mu)))
    }
    drop(mu + crossprod(R, stats::rnorm(length(mu))))
}

# Design of the imputation model for v from the other columns (factors as dummies)
.imp_design <- function(df, preds) {
    X <- stats::model.matrix(~ ., data = df[, preds, drop = FALSE])
    X[, colSums(abs(X)) > 0 | colnames(X) == "(Intercept)", drop = FALSE]
}

# Predictive mean matching (type 1, k donors) with a draw of the regression parameters
.imp_pmm <- function(y, X, obs, k = 5) {
    Xo <- X[obs, , drop = FALSE]; yo <- y[obs]
    XtX <- crossprod(Xo); ridge <- 1e-5 * diag(diag(XtX), ncol(Xo))
    Vu <- solve(XtX + ridge)
    b <- drop(Vu %*% crossprod(Xo, yo))
    df <- max(1, sum(obs) - ncol(Xo))
    s2 <- sum((yo - Xo %*% b)^2) / stats::rchisq(1, df)
    bs <- .rmvn(b, s2 * Vu)
    yo_hat <- drop(Xo %*% b); ym_hat <- drop(X[!obs, , drop = FALSE] %*% bs)
    vapply(ym_hat, function(h) {
        d <- abs(yo_hat - h)
        idx <- order(d)[seq_len(min(k, length(d)))]
        yo[idx[sample.int(length(idx), 1)]]
    }, numeric(1))
}

# Logistic regression with a parameter draw (two levels)
.imp_logreg <- function(f, X, obs) {
    lev <- levels(f)
    yb <- as.integer(f[obs] == lev[2])
    fit <- suppressWarnings(stats::glm.fit(X[obs, , drop = FALSE], yb, family = stats::binomial()))
    b <- fit$coefficients; b[is.na(b)] <- 0
    V <- .glm_vcov(fit); V[is.na(V)] <- 0
    bs <- .rmvn(b, V + diag(1e-8, length(b)))
    p <- stats::plogis(drop(X[!obs, , drop = FALSE] %*% bs))
    factor(lev[1 + (stats::runif(length(p)) < p)], lev)
}

# Multinomial regression (more than two levels), categories drawn from the fitted probabilities
.imp_polyreg <- function(f, X, obs) {
    lev <- levels(f)
    # through a data frame with plain names (predict.multinom cannot use a matrix term)
    Z <- as.data.frame(X[, -1, drop = FALSE]); names(Z) <- paste0("z", seq_len(ncol(Z)))
    fit <- nnet::multinom(.f ~ ., data = cbind(.f = f[obs], Z[obs, , drop = FALSE]), trace = FALSE, MaxNWts = 10000)
    P <- stats::predict(fit, newdata = Z[!obs, , drop = FALSE], type = "probs")
    P <- matrix(P, nrow = sum(!obs))
    if (ncol(P) == 1) P <- cbind(1 - P, P)
    u <- stats::runif(nrow(P))
    factor(lev[1 + rowSums(u > t(apply(P, 1, cumsum)))], lev)
}

# m completed copies of df: the variables in impVars are imputed (numeric by
# predictive mean matching, factors by logistic or multinomial regression)
# from all other columns of df, which must be complete (e.g. the outcome).
# Fully conditional specification (van Buuren 2018), maxit iterations, each
# imputation an independent chain started from random observed values.
cp_impute <- function(df, impVars, m = 5, maxit = 10, seed = 1, donors = 5) {
    impVars <- Filter(function(v) anyNA(df[[v]]), impVars)
    others <- setdiff(names(df), impVars)
    if (any(vapply(others, function(v) anyNA(df[[v]]), TRUE)))
        stop("cp_impute: the auxiliary columns must be complete")
    miss <- lapply(stats::setNames(impVars, impVars), function(v) is.na(df[[v]]))
    # one seed per chain: a chain is the same whichever core runs it
    .cp_lapply(seq_len(m), function(i) .with_seed(seed + 7919L * (i - 1L), {
        d <- df
        for (v in impVars) {
            o <- !miss[[v]]
            d[[v]][!o] <- d[[v]][o][sample.int(sum(o), sum(!o), replace = TRUE)]
        }
        if (!length(impVars)) return(d)
        for (it in seq_len(maxit)) for (v in impVars) {
            preds <- setdiff(names(d), v)
            X <- .imp_design(d, preds)
            o <- !miss[[v]]
            x <- d[[v]]
            d[[v]][!o] <- if (is.factor(x)) {
                if (nlevels(droplevels(x[o])) < 2) x[o][1]
                else if (nlevels(x) == 2) .imp_logreg(droplevels(x), X, o) else .imp_polyreg(droplevels(x), X, o)
            } else .imp_pmm(as.numeric(x), X, o, donors)
        }
        d
    }))
}

# ---- pooling ---------------------------------------------------------------

# Rubin's rules for coefficient vectors and their covariances
cp_rubin <- function(betas, Vs) {
    m <- length(betas)
    B <- do.call(rbind, betas)
    qbar <- colMeans(B)
    W <- Reduce(`+`, Vs) / m
    Bv <- if (m > 1) stats::cov(B) else matrix(0, length(qbar), length(qbar))
    Tm <- W + (1 + 1 / m) * Bv
    list(est = qbar, vcov = Tm, within = W, between = Bv, m = m)
}

# Pooled Wald test (D1; Li, Raghunathan & Rubin 1991) that the coefficients
# idx are all zero, with its F reference distribution
cp_d1 <- function(pool, idx, dfcom = Inf) {
    k <- length(idx); m <- pool$m
    q <- pool$est[idx]
    W <- pool$within[idx, idx, drop = FALSE]; B <- pool$between[idx, idx, drop = FALSE]
    r <- if (m > 1) (1 + 1 / m) * sum(diag(B %*% solve(W))) / k else 0
    Tt <- (1 + r) * W
    D <- drop(t(q) %*% solve(Tt) %*% q) / k
    t <- k * (m - 1)
    v <- if (m == 1 || r == 0) Inf else if (t > 4) 4 + (t - 4) * (1 + (1 - 2 / t) / r)^2 else t * (1 + 1 / k) * (1 + 1 / r)^2 / 2
    c(stat = D, df1 = k, df2 = v, p = stats::pf(D, k, v, lower.tail = FALSE))
}

# Elastic net for the three families (Friedman, Hastie & Tibshirani 2010;
# Simon et al. 2011 for Cox), in pure R. The objective is glmnet's:
#   loss/n + lambda * sum_j pf_j ((1 - alpha)/2 b_j^2 + alpha |b_j|)
# with loss the negative log-likelihood (half the residual sum of squares for
# the linear model; the Breslow partial likelihood for Cox), the columns
# standardised (1/n variance) and the penalty factors rescaled to sum to the
# number of columns, as glmnet does. Penalty factor 0 (baseline covariates)
# leaves a column unpenalised.
#
# Each outer step replaces the loss by its quadratic approximation (weights w,
# working response z; one exact step for the linear model) and minimises it by
# coordinate descent on the p x p weighted Gram matrix, so a sweep costs
# O(p^2) whatever n.

# Working weights and response of the current linear predictor eta
.enet_work <- function(family, y, eta) {
    switch(family,
        gaussian = list(w = rep(1, length(y)), z = y),
        binary = {
            p <- stats::plogis(eta); w <- pmax(p * (1 - p), 1e-5)
            list(w = w, z = eta + (y - p) / w)
        },
        cox = {
            # Breslow partial likelihood: gradient and diagonal of the Hessian
            tm <- y[, 1]; d <- y[, 2]
            o <- order(tm); e <- exp(eta - max(eta))
            S <- rev(cumsum(rev(e[o])))                       # risk-set sums (sorted)
            first <- match(tm[o], tm[o])                     # tied times share a risk set
            S <- S[first]
            ev <- d[o] == 1
            # for each subject: sum over event times t_k <= t_i of d_k/S_k and d_k/S_k^2
            a1 <- cumsum(ifelse(ev, 1 / S, 0)); a2 <- cumsum(ifelse(ev, 1 / S^2, 0))
            last <- findInterval(tm[o], tm[o])                # last sorted index with time <= t
            A1 <- a1[last]; A2 <- a2[last]
            g <- d[o] - e[o] * A1
            w <- pmax(e[o] * A1 - e[o]^2 * A2, 1e-8)
            gg <- ww <- numeric(length(tm)); gg[o] <- g; ww[o] <- w
            list(w = ww, z = eta + gg / ww)
        })
}

.enet_loss <- function(family, y, eta) switch(family,
    gaussian = sum((y - eta)^2) / 2,
    binary = -sum(y * eta - log1p(exp(eta))),
    cox = {
        tm <- y[, 1]; d <- y[, 2]
        o <- order(tm); m <- max(eta); e <- exp(eta - m)
        S <- rev(cumsum(rev(e[o])))[match(tm[o], tm[o])]
        -sum((d[o] == 1) * (eta[o] - log(S) - m))
    })

# Minimiser of 0.5 b'Hb - g'b + sum(lam1 |b| + lam2 b^2 / 2), warm start b.
# Coordinate-descent sweeps find the active set and signs; for that set the
# optimum solves (H_AA + diag(lam2_A)) b_A = g_A - lam1_A sign(b_A) exactly,
# accepted when the signs hold and the zero coefficients satisfy their KKT
# condition |g_j - (Hb)_j| <= lam1_j. Unpenalised columns are always active.
.enet_cd <- function(H, g, b, lam1, lam2, maxit = 1000) {
    p <- length(b); hd <- diag(H); den <- hd + lam2
    free <- lam1 == 0
    Hb <- drop(H %*% b)
    for (iter in seq_len(maxit)) {
        for (j in seq_len(p)) {
            zj <- g[j] - Hb[j] + hd[j] * b[j]
            l <- lam1[j]
            bj <- if (zj > l) (zj - l) / den[j] else if (zj < -l) (zj + l) / den[j] else 0
            if (bj != b[j]) { Hb <- Hb + H[, j] * (bj - b[j]); b[j] <- bj }
        }
        A <- which(b != 0 | free)
        sgn <- sign(b[A]); sgn[free[A]] <- 0
        bA <- if (length(A)) tryCatch(solve(H[A, A, drop = FALSE] + diag(lam2[A], length(A)), g[A] - lam1[A] * sgn),
                                      error = function(e) NULL) else numeric(0)
        if (is.null(bA)) next
        pen <- !free[A]
        if (any(sign(bA[pen]) != sgn[pen])) next
        b2 <- numeric(p); b2[A] <- bA
        Hb2 <- if (length(A)) drop(H[, A, drop = FALSE] %*% bA) else numeric(p)
        out <- setdiff(seq_len(p), A)
        if (all(abs(g[out] - Hb2[out]) <= lam1[out] * (1 + 1e-9) + 1e-12)) return(b2)
        b <- b2; Hb <- Hb2
    }
    b
}

# Penalised fit on standardised columns for a decreasing lambda sequence
# (warm starts). Returns the coefficients on the original scale, one column per lambda.
cp_enet_path <- function(X, y, family, alpha = 0.5, pf = rep(1, ncol(X)), lambda = NULL,
                         nlambda = 100, ratio = NULL, tol = 1e-7) {
    n <- NROW(X); p <- ncol(X)
    pf <- pf * p / sum(pf)                                     # glmnet's rescaling
    mu <- colMeans(X); sdv <- sqrt(colMeans(sweep(X, 2, mu)^2)); sdv[sdv == 0] <- 1
    Z <- sweep(sweep(X, 2, mu), 2, sdv, "/")
    intercept <- family != "cox"
    # glmnet standardises a continuous outcome, which divides its ridge part by SD(y)
    ridgeScale <- if (family == "gaussian") 1 / sqrt(mean((y - mean(y))^2)) else 1
    fitStep <- function(b, b0, lam) {
        for (outer in 1:100) {
            eta <- b0 + drop(Z %*% b)
            wk <- .enet_work(family, y, eta)
            w <- wk$w / n
            ZW <- Z * w
            H <- crossprod(ZW, Z); g <- drop(crossprod(ZW, wk$z))
            if (intercept) {
                # centred on the weighted means (the intercept is not penalised),
                # algebraically: Z'WZ - sw xw xw', Z'Wz - sw xw zw
                sw <- sum(w); xw <- colSums(ZW) / sw; zw <- sum(w * wk$z) / sw
                H <- H - sw * tcrossprod(xw); g <- g - sw * xw * zw
            }
            bNew <- .enet_cd(H, g, b, lam * alpha * pf, lam * (1 - alpha) * pf * ridgeScale)
            b0New <- if (intercept) zw - sum(xw * bNew) else 0
            done <- max(abs(bNew - b)) < tol && abs(b0New - b0) < tol
            b <- bNew; b0 <- b0New
            if (family == "gaussian" || done) break
        }
        list(b = b, b0 = b0)
    }
    # the largest useful lambda: every penalised column at 0 (unpenalised ones fitted)
    st <- fitStep(rep(0, p), if (family == "binary") stats::qlogis(mean(y)) else if (family == "gaussian") mean(y) else 0, 1e10)
    stopEarly <- is.null(lambda)
    if (is.null(lambda)) {
        eta <- st$b0 + drop(Z %*% st$b)
        wk <- .enet_work(family, y, eta)
        grad <- abs(drop(crossprod(Z, wk$w * (wk$z - eta)))) / n
        pen <- pf > 0
        lmax <- max(grad[pen] / pf[pen]) / max(alpha, 1e-3)
        if (is.null(ratio)) ratio <- if (n > p) 1e-4 else 1e-2
        lambda <- exp(seq(log(lmax), log(lmax * ratio), length.out = nlambda))
    }
    B <- matrix(0, p, length(lambda)); B0 <- numeric(length(lambda))
    b <- st$b; b0 <- st$b0
    devPrev <- NA_real_; last <- length(lambda)
    for (k in seq_along(lambda)) {
        f <- fitStep(b, b0, lambda[k])
        b <- f$b; b0 <- f$b0
        B[, k] <- b / sdv
        B0[k] <- b0 - sum(mu * b / sdv)
        # a path computed here (not given) stops once the deviance hardly moves (as glmnet)
        if (stopEarly) {
            dev <- .enet_loss(family, y, b0 + drop(Z %*% b))
            if (k >= 5 && is.finite(devPrev) && abs(devPrev - dev) < 1e-5 * abs(devPrev)) { last <- k; break }
            devPrev <- dev
        }
    }
    keep <- seq_len(last)
    B <- B[, keep, drop = FALSE]; B0 <- B0[keep]; lambda <- lambda[keep]
    dimnames(B) <- list(colnames(X), NULL)
    list(beta = B, a0 = B0, lambda = lambda)
}

# Cross-validated deviance along the path (k folds, stratified); the Cox
# deviance is the grouped partial likelihood (Verweij & van Houwelingen 1993)
cp_enet_cv <- function(X, y, family, alpha, pf, nfolds = 10, seed = 1, strata = NULL) {
    full <- cp_enet_path(X, y, family, alpha, pf)
    lam <- full$lambda
    n <- NROW(X)
    fold <- .with_seed(seed, .folds(strata %||% rep(1, n), nfolds))
    # the folds are fitted possibly on several cores (not inside a resampled fit)
    perFold <- .cp_lapply(seq_len(nfolds), function(f) {
        tr <- fold != f
        ytr <- if (is.matrix(y)) y[tr, , drop = FALSE] else y[tr]
        fit <- cp_enet_path(X[tr, , drop = FALSE], ytr, family, alpha, pf, lambda = lam)
        dv <- vapply(seq_along(lam), function(k) {
            eta <- fit$a0[k] + drop(X %*% fit$beta[, k])
            if (family == "cox") 2 * (.enet_loss("cox", y, eta) - .enet_loss("cox", ytr, eta[tr])) else
                2 * .enet_loss(family, if (is.matrix(y)) y[!tr, , drop = FALSE] else y[!tr], eta[!tr])
        }, 0)
        list(dev = dv, w = if (family == "cox") sum(y[!tr, 2]) else sum(!tr))
    })
    dev <- do.call(rbind, lapply(perFold, `[[`, "dev")); wts <- vapply(perFold, `[[`, 0, "w")
    cvm <- colSums(dev) / sum(wts)
    cvsd <- sqrt(colSums(sweep(dev / wts, 2, cvm)^2 * wts) / sum(wts) / (nfolds - 1))
    imin <- which.min(cvm)
    i1se <- min(which(cvm <= cvm[imin] + cvsd[imin]))
    list(path = full, cvm = cvm, cvsd = cvsd, lambda = lam, imin = imin, i1se = i1se, nfolds = nfolds)
}

# The penalised model, in the form of the other fits: coefficients with the
# intercept first (log H0(h) for Cox), lambda chosen by CV unless given
cp_enet_fit <- function(X, y, fam, alpha = 0.5, pf = rep(1, ncol(X)), rule = "min", nfolds = 10, seed = 1,
                        lambda = NULL) {
    family <- fam$name
    cv <- NULL
    if (is.null(lambda)) {
        cv <- cp_enet_cv(X, y, family, alpha, pf, nfolds, seed, fam$strata(y))
        k <- if (rule == "1se") cv$i1se else cv$imin
        b <- cv$path$beta[, k]; a0 <- cv$path$a0[k]; lambda <- cv$lambda[k]
    } else {
        path <- cp_enet_path(X, y, family, alpha, pf, lambda = lambda)
        b <- path$beta[, 1]; a0 <- path$a0[1]
    }
    if (family == "cox") {
        a0 <- log(.cox_H0(y[, 1], y[, 2], drop(X %*% b), fam$horizon))
    }
    beta <- c(`(Intercept)` = a0, b)
    attr(beta, "converged") <- TRUE
    attr(beta, "penalized") <- TRUE
    attr(beta, "lambda") <- lambda
    attr(beta, "cv") <- cv
    beta
}

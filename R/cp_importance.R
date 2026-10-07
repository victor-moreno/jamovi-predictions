# Variable importance, family-agnostic. A "variable" is a predictor as the user
# entered it: all design columns of a factor or a spline count as one player.

# Shapley decomposition of an R2 measure (LMG; Lindeman, Merenda & Gold 1980;
# general dominance, Azen & Traxel 2009): each variable's average gain in R2
# over all orders of entry. Exact up to maxExact variables, otherwise estimated
# from nPerm random orders (a Monte Carlo Shapley value).
cp_shapley_r2 <- function(X, y, spec, fam, measure = "r2mf", maxExact = 10, nPerm = 200,
                          seed = 1) {
    vars <- names(spec$vars)
    m <- length(vars)
    cache <- new.env(hash = TRUE)
    value <- function(set) {
        if (!length(set)) return(0)
        key <- paste(sort(set), collapse = ",")
        if (!is.null(cache[[key]])) return(cache[[key]])
        cols <- spec$colVar %in% vars[set]
        Xs <- X[, cols, drop = FALSE]
        beta <- fam$fit(Xs, y)
        v <- fam$r2(y, fam$lp(beta, Xs), measure)
        assign(key, v, envir = cache)
        v
    }
    if (m <= maxExact) {
        phi <- numeric(m)
        w <- factorial(0:(m - 1)) * factorial((m - 1):0) / factorial(m)
        for (j in seq_len(m)) {
            others <- setdiff(seq_len(m), j)
            for (s in 0:(m - 1)) {
                subsets <- if (s == 0) list(integer(0)) else
                    if (length(others) == 1) list(others) else utils::combn(others, s, simplify = FALSE)
                for (S in subsets) phi[j] <- phi[j] + w[s + 1] * (value(c(S, j)) - value(S))
            }
        }
        exact <- TRUE
    } else {
        phi <- numeric(m)
        .with_seed(seed, {
            for (r in seq_len(nPerm)) {
                ord <- sample.int(m)
                prev <- 0
                for (i in seq_len(m)) {
                    cur <- value(ord[seq_len(i)])
                    phi[ord[i]] <- phi[ord[i]] + cur - prev
                    prev <- cur
                }
            }
        })
        phi <- phi / nPerm
        exact <- FALSE
    }
    total <- value(seq_len(m))
    out <- data.frame(var = vars, shapley = phi, share = phi / total,
                      rank = rank(-phi, ties.method = "min"), stringsAsFactors = FALSE)
    structure(out, total = total, exact = exact, nFits = length(ls(cache)))
}

# Permutation importance: loss of performance of the fitted model when one
# variable's values are shuffled (Breiman 2001; Fisher, Rudin & Dominici 2019).
# Measured on the data given (development or validation).
cp_perm_importance <- function(beta, spec, df, y, fam, nPerm = 20, seed = 1,
                               keys = c("auc", "brier"), higherBetter = c(TRUE, FALSE)) {
    base <- fam$quick(y, fam$lp(beta, cp_design(spec, df)))[keys]
    # positive = loss of performance, whichever direction the measure runs
    sgn <- ifelse(higherBetter, 1, -1)
    vars <- names(spec$vars)
    res <- .with_seed(seed, lapply(vars, function(v) {
        d <- vapply(seq_len(nPerm), function(r) {
            dp <- df
            dp[[v]] <- dp[[v]][sample.int(nrow(dp))]
            sgn * (base - fam$quick(y, fam$lp(beta, cp_design(spec, dp)))[keys])
        }, numeric(length(keys)))
        d <- matrix(d, nrow = length(keys))
        c(mean = rowMeans(d), sd = apply(d, 1, stats::sd))
    }))
    out <- data.frame(var = vars, do.call(rbind, res), stringsAsFactors = FALSE, check.names = FALSE)
    names(out) <- c("var", paste0("d_", keys), paste0("sd_", keys))
    attr(out, "base") <- base
    out
}

# SHAP values of an additive model on the linear-predictor scale (log-odds for
# logistic, log hazard for Cox). With main effects only, the exact
# interventional Shapley value of variable j for individual i is
#   phi_ij = sum over j's columns c of beta_c * (X_ic - mean of X_c in the background)
# (Aas, Jullum & Loland 2021), and the phi of an individual add up to its linear
# predictor minus the mean linear predictor of the background (development) data.
cp_shap <- function(beta, spec, X, Xbg = X, intercept = TRUE) {
    b <- beta; b[is.na(b)] <- 0
    if (intercept) b <- b[-1]
    centred <- sweep(X, 2, colMeans(Xbg))
    contrib <- sweep(centred, 2, b, `*`)
    vars <- names(spec$vars)
    phi <- vapply(vars, function(v) rowSums(contrib[, spec$colVar == v, drop = FALSE]),
                  numeric(nrow(X)))
    phi <- matrix(phi, nrow = nrow(X), dimnames = list(NULL, vars))
    attr(phi, "base") <- sum(c(if (intercept) beta[1] else 0, b * colMeans(Xbg)), na.rm = TRUE)
    phi
}

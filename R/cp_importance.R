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
        # a column enters when all its variables are in the coalition (an
        # interaction joins once both of its variables have)
        cols <- vapply(spec$colVars, function(vs) all(vs %in% vars[set]), TRUE)
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

# SHAP values of the model on the linear-predictor scale (log-odds for
# logistic, log hazard for Cox): exact interventional Shapley values
# (Aas, Jullum & Loland 2021), the background being the development data.
# A main-effect column c of variable j gives phi_ij += beta_c (X_ic - mean X_c).
# An interaction column beta_c a(x_A) b(x_B) is a two-player game with
#   v0 = beta_c mean(a b), vA = beta_c a_i mean(b), vB = beta_c mean(a) b_i, vAB = beta_c a_i b_i
# so phi_A = ((vA - v0) + (vAB - vB)) / 2 and phi_B = ((vB - v0) + (vAB - vA)) / 2.
# The phi of an individual add up to its linear predictor minus the mean
# linear predictor of the background.
cp_shap <- function(beta, spec, X, Xbg = X, intercept = TRUE) {
    b <- beta; b[is.na(b)] <- 0
    if (intercept) b <- b[-1]
    vars <- names(spec$vars)
    phi <- matrix(0, nrow(X), length(vars), dimnames = list(NULL, vars))
    for (c in seq_along(spec$columns)) {
        if (b[c] == 0) next
        ab <- spec$colParts[[c]]
        if (is.null(ab)) {
            v <- spec$colVars[[c]]
            phi[, v] <- phi[, v] + b[c] * (X[, c] - mean(Xbg[, c]))
        } else {
            vs <- spec$colVars[[c]]
            a <- X[, ab[1]]; bb <- X[, ab[2]]
            v0 <- b[c] * mean(Xbg[, c]); vA <- b[c] * a * mean(Xbg[, ab[2]])
            vB <- b[c] * mean(Xbg[, ab[1]]) * bb; vAB <- b[c] * a * bb
            phi[, vs[1]] <- phi[, vs[1]] + ((vA - v0) + (vAB - vB)) / 2
            phi[, vs[2]] <- phi[, vs[2]] + ((vB - v0) + (vAB - vA)) / 2
        }
    }
    attr(phi, "base") <- sum(c(if (intercept) beta[1] else 0, b * colMeans(Xbg)), na.rm = TRUE)
    phi
}

# Incremental value of predictors added to a baseline model: differences in
# performance with paired confidence intervals, and reclassification.

# Difference of two correlated AUCs (DeLong, DeLong & Clarke-Pearson 1988):
# AUC(p2) - AUC(p1) with its standard error from the placement values
cp_delta_auc <- function(y, p1, p2, conf = 0.95) {
    pv <- function(p) {
        p1 <- p[y == 1]; p0 <- p[y == 0]; n1 <- length(p1); n0 <- length(p0)
        r <- rank(c(p1, p0))
        list(auc = (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0),
             v10 = (r[seq_len(n1)] - rank(p1)) / n0, v01 = 1 - (r[n1 + seq_len(n0)] - rank(p0)) / n1)
    }
    a <- pv(p1); b <- pv(p2)
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    s10 <- stats::cov(cbind(a$v10, b$v10)); s01 <- stats::cov(cbind(a$v01, b$v01))
    S <- s10 / n1 + s01 / n0
    d <- b$auc - a$auc
    se <- sqrt(max(0, S[1, 1] + S[2, 2] - 2 * S[1, 2]))
    z <- stats::qnorm(1 - (1 - conf) / 2)
    c(est = d, lower = d - z * se, upper = d + z * se, p = 2 * stats::pnorm(-abs(d / se)))
}

# Difference of Harrell's C of two Cox linear predictors, with survival's
# covariance of the two concordances
cp_delta_c <- function(y, lp1, lp2, conf = 0.95) {
    cc <- survival::concordance(survival::Surv(y[, 1], y[, 2]) ~ lp1 + lp2, reverse = TRUE)
    d <- unname(cc$concordance[2] - cc$concordance[1])
    se <- sqrt(max(0, cc$var[1, 1] + cc$var[2, 2] - 2 * cc$var[1, 2]))
    z <- stats::qnorm(1 - (1 - conf) / 2)
    c(est = d, lower = d - z * se, upper = d + z * se, p = 2 * stats::pnorm(-abs(d / se)))
}

# Continuous net reclassification improvement (Pencina et al. 2011) and the
# integrated discrimination improvement (Pencina et al. 2008) of risks p2 over p1
cp_nri_idi <- function(y, p1, p2, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    e <- y == 1
    up <- p2 > p1; down <- p2 < p1
    ne <- sum(e); nn <- sum(!e)
    pue <- mean(up[e]); pde <- mean(down[e]); pun <- mean(up[!e]); pdn <- mean(down[!e])
    nri <- (pue - pde) + (pdn - pun)
    se <- sqrt((pue + pde - (pue - pde)^2) / ne + (pun + pdn - (pun - pdn)^2) / nn)
    dd <- p2 - p1
    idi <- mean(dd[e]) - mean(dd[!e])
    sei <- sqrt(stats::var(dd[e]) / ne + stats::var(dd[!e]) / nn)
    list(nri = c(est = nri, lower = nri - z * se, upper = nri + z * se, p = 2 * stats::pnorm(-abs(nri / se))),
         nriEvents = pue - pde, nriNonevents = pdn - pun,
         idi = c(est = idi, lower = idi - z * sei, upper = idi + z * sei, p = 2 * stats::pnorm(-abs(idi / sei))))
}

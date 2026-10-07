# Model updating in validation data (Steyerberg 2019, ch. 20; Vergouwe et al.
# 2017), family-agnostic: the original model, recalibration in the large (new
# intercept), recalibration (new intercept and slope of the linear predictor)
# and model revision (every coefficient re-estimated). The closed testing
# procedure tests the revised model against the original, then against
# recalibration in the large, then against recalibration, and keeps the first
# method that is not significantly worse than revision.

cp_update <- function(fam, X, y, lp, alpha = 0.05) {
    n <- NROW(X)
    a1 <- unname(fam$metrics(y, lp)[["citl"]])
    f2 <- fam$fit(matrix(lp, ncol = 1, dimnames = list(NULL, "lp")), y)
    br <- tryCatch(fam$fit(X, y), error = function(e) NULL)
    revOk <- !is.null(br) && isTRUE(attr(br, "converged"))
    methods <- list(
        original = list(a = 0, b = 1, lp = lp, k = 0),
        citl = list(a = a1, b = 1, lp = lp + a1, k = 1),
        recal = list(a = unname(f2[1]), b = unname(f2[2]), lp = fam$lp(f2, matrix(lp, ncol = 1)), k = 2),
        revision = if (revOk) list(a = NA_real_, b = NA_real_, lp = fam$lp(br, X),
                                   k = sum(!is.na(br)), beta = br) else NULL)
    methods <- Filter(Negate(is.null), methods)
    for (m in names(methods)) methods[[m]]$loss <- fam$loss(y, methods[[m]]$lp)

    chosen <- NA_character_
    if (revOk) {
        rev <- methods$revision
        dfResid <- n - rev$k
        for (m in c("original", "citl", "recal")) {
            t <- .nested_test(fam, methods[[m]]$loss, rev$loss, rev$k - methods[[m]]$k, n, dfResid)
            methods[[m]]$test <- t
            if (is.na(chosen) && t[["p"]] >= alpha) chosen <- m
        }
        if (is.na(chosen)) chosen <- "revision"
    }
    list(methods = methods, chosen = chosen)
}

# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container: the apparent AUC and the
# coefficients of a logistic model on MASS::Pima.tr against glm.
d <- MASS::Pima.tr
r <- clinpred::cpLogistic(data = d, dep = "type", eventLevel = "Yes",
                          covs = c("glu", "bmi", "ped"), coefs = TRUE)
co <- r$coefTable$asDF
fit <- glm(type ~ glu + bmi + ped, family = binomial, data = d)
stopifnot(isTRUE(all.equal(co$est, unname(coef(fit)), tolerance = 1e-6)))
p <- fitted(fit); y <- d$type == "Yes"
auc <- mean(outer(p[y], p[!y], ">") + 0.5 * outer(p[y], p[!y], "=="))
stopifnot(isTRUE(all.equal(r$perfDev$asDF$apparent[1], auc, tolerance = 1e-8)))
cat(sprintf("   smoke test passed: AUC = %.4f, glu = %.4f\n", auc, co$est[2]))

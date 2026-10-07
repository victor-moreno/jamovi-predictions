# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container: each analysis against the base
# fit it wraps (glm, lm, coxph) on example data.
d <- MASS::Pima.tr
r <- clinpred::cpLogistic(data = d, dep = "type", eventLevel = "Yes",
                          covs = c("glu", "bmi", "ped"), coefs = TRUE)
co <- r$coefTable$asDF
fit <- glm(type ~ glu + bmi + ped, family = binomial, data = d)
stopifnot(isTRUE(all.equal(co$est, unname(coef(fit)), tolerance = 1e-6)))
p <- fitted(fit); y <- d$type == "Yes"
auc <- mean(outer(p[y], p[!y], ">") + 0.5 * outer(p[y], p[!y], "=="))
stopifnot(isTRUE(all.equal(r$perfDev$asDF$apparent[1], auc, tolerance = 1e-8)))
cat(sprintf("   smoke test passed: logistic AUC = %.4f, glu = %.4f\n", auc, co$est[2]))

r <- clinpred::cpLinear(data = d, dep = "glu", covs = c("bmi", "ped"), factors = "type", coefs = TRUE)
m <- lm(glu ~ bmi + ped + type, d)
stopifnot(isTRUE(all.equal(r$coefTable$asDF$est, unname(coef(m)), tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(r$perfDev$asDF$apparent[1], summary(m)$r.squared, tolerance = 1e-8)))
cat(sprintf("   smoke test passed: linear R² = %.4f\n", summary(m)$r.squared))

g <- survival::gbsg
g$status <- factor(g$status)
r <- clinpred::cpCox(data = g, time = "rfstime", status = "status", eventLevel = "1",
                     covs = c("age", "nodes"), horizon = 1826, coefs = TRUE)
f <- survival::coxph(survival::Surv(rfstime, status == "1") ~ age + nodes, g)
stopifnot(isTRUE(all.equal(r$coefTable$asDF$est, unname(coef(f)), tolerance = 1e-6)))
stopifnot(isTRUE(all.equal(r$perfDev$asDF$apparent[1], unname(survival::concordance(f)$concordance), tolerance = 1e-8)))
cat(sprintf("   smoke test passed: Cox C = %.4f\n", r$perfDev$asDF$apparent[1]))

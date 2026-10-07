# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container. The runner has already put
# jamovi's base library and the installed module first on .libPaths() and
# attached the module; the module name is in Sys.getenv("MODULE") in docker.
#
# Call one analysis the way jamovi would and stop() on a wrong answer, e.g.:
#
#   r <- mymodule::myanalysis(data = data.frame(x = c(1, 2, 3)), vars = "x")
#   stopifnot(isTRUE(all.equal(r$table$asDF$mean[1], 2)))
#   cat("   smoke test passed\n")
cat("   no smoke test defined yet (edit tools/smoke.R)\n")

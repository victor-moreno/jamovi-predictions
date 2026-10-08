# Opt-in multicore computation (macOS and Linux: forked workers from base R's
# parallel package; Windows always sequential). The analysis sets the number of
# cores and a checkpoint for the duration of a run; every random draw is made
# in the main process before the work is shared out, so results do not depend
# on the number of cores. Workers never fork again (nested calls run
# sequentially). The checkpoint runs between chunks, so jamovi can abandon a
# long computation when the options change.

.cp_par <- new.env()
.cp_par$cores <- 1L
.cp_par$tick <- NULL

.cp_cores <- function() {
    if (.Platform$OS.type != "unix") return(1L)
    max(1L, as.integer(.cp_par$cores %||% 1L))
}

# lapply over X, in chunks, on .cp_cores() forked workers; an error in a
# worker is raised again here
.cp_lapply <- function(X, FUN) {
    n <- length(X)
    out <- vector("list", n)
    if (!n) return(out)
    cores <- .cp_cores()
    chunk <- if (cores > 1) cores * 4L else max(1L, ceiling(n / 20))
    for (start in seq(1L, n, by = chunk)) {
        idx <- start:min(n, start + chunk - 1L)
        res <- if (cores > 1) {
            parallel::mclapply(X[idx], function(x) {
                .cp_par$cores <- 1L; .cp_par$tick <- NULL
                tryCatch(FUN(x), error = function(e) e)
            }, mc.cores = min(cores, length(idx)), mc.preschedule = TRUE)
        } else lapply(X[idx], FUN)
        if (cores > 1) for (r in res) if (inherits(r, "error")) stop(r)
        out[idx] <- res
        if (!is.null(.cp_par$tick)) .cp_par$tick()
    }
    out
}

# Run expr with the cores and checkpoint of an analysis
.cp_with_cores <- function(cores, tick, expr) {
    old <- list(cores = .cp_par$cores, tick = .cp_par$tick)
    on.exit({ .cp_par$cores <- old$cores; .cp_par$tick <- old$tick })
    .cp_par$cores <- cores
    .cp_par$tick <- tick
    expr
}

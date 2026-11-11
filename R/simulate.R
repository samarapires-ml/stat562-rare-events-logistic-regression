library(parallel)

# Generate simulated data
gen_data <- function(n1, n0, p0, OR) {
  K <- length(n1)
  p0 <- rep_len(p0, K)
  p1 <- plogis(qlogis(p0) + log(OR))
  
  list(
    x1 = rbinom(K, n1, p1),
    n1 = n1,
    x0 = rbinom(K, n0, p0),
    n0 = n0
  )
}

# Apply all methods
fit_all <- function(s) {
  do.call(rbind, lapply(names(methods), function(nm) {
    r <- tryCatch(
      suppressWarnings(methods[[nm]](s$x1, s$n1, s$x0, s$n0)),
      error = function(e) c(logOR = NA, lo = NA, hi = NA, p = NA)
    )
    
    data.frame(
      method = nm,
      logOR = r[["logOR"]],
      lo = r[["lo"]],
      hi = r[["hi"]],
      p = r[["p"]],
      error = is.na(r[["logOR"]])
    )
  }))
}

# Run simulation
run_sim <- function(R, OR, p0, n1, n0, seed, cores = 1) {
  if (.Platform$OS.type == "windows") cores <- 1
  
  one_rep <- function(r) {
    set.seed(seed * 100000 + r)
    
    s <- gen_data(n1, n0, p0, OR)
    res <- fit_all(s)
    
    res$rep <- r
    res$zero_any <- sum(s$x1 == 0 | s$x0 == 0)
    res$zero_double <- sum(s$x1 == 0 & s$x0 == 0)
    
    res
  }
  
  out <- if (cores > 1) mclapply(seq_len(R), one_rep, mc.cores = cores)
  else lapply(seq_len(R), one_rep)
  
  out <- do.call(rbind, out)
  out$true_logOR <- log(OR)
  
  out
}
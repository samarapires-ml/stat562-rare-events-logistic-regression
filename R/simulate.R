library(parallel)

# Generate simulated data
gen_data <- function(n1, n0, p0, OR, p0_sd = 0) {
  K <- length(n1)
  p0 <- rep_len(p0, K)
  
  if (p0_sd > 0) p0 <- plogis(qlogis(p0) + rnorm(K, 0, p0_sd))
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
run_sim <- function(R, OR, p0, n1, n0, seed, cores = 1, p0_sd = 0) {
  force(OR)
  force(p0)
  force(n1)
  force(n0)
  force(seed)
  force(p0_sd)
  
  one_rep <- function(r) {
    set.seed(seed * 100000 + r)
    
    s <- gen_data(n1, n0, p0, OR, p0_sd)
    res <- fit_all(s)
    
    res$rep <- r
    res$zero_any <- sum(s$x1 == 0 | s$x0 == 0)
    res$zero_double <- sum(s$x1 == 0 & s$x0 == 0)
    
    res
  }
  
  if (cores <= 1) {
    out <- lapply(seq_len(R), one_rep)
  } else if (.Platform$OS.type == "windows") {
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterCall(cl, setwd, getwd())
    clusterEvalQ(cl, {
      suppressMessages(source("R/methods.R"))
      source("R/simulate.R")
    })
    
    out <- parLapply(cl, seq_len(R), one_rep)
  } else {
    out <- mclapply(seq_len(R), one_rep, mc.cores = cores)
  }
  
  out <- do.call(rbind, out)
  out$true_logOR <- log(OR)
  
  out
}

# Summarise simulation results
summarise_sim <- function(sim, alpha = 0.05) {
  do.call(rbind, lapply(split(sim, sim$method), function(s) {
    truth <- s$true_logOR[1]
    ok <- !s$error
    est <- s$logOR[ok]
    n <- sum(ok)
    
    cover <- mean(s$lo[ok] <= truth & truth <= s$hi[ok])
    reject <- mean(s$p[ok] < alpha)
    
    data.frame(
      method = s$method[1],
      n_ok = n,
      n_failed = sum(!ok),
      bias = mean(est) - truth,
      bias_mcse = sd(est) / sqrt(n),
      emp_se = sd(est),
      coverage = cover,
      cover_mcse = sqrt(cover * (1 - cover) / n),
      reject = reject,
      reject_mcse = sqrt(reject * (1 - reject) / n)
    )
  }))
}
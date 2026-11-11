## =====================================================================
## Conditional Firth estimator of a common odds ratio
## for stratified 2x2 tables (e.g. rare-event meta-analysis of trials)
##
## STAT 562 group project
##
## Model: within trial k, logit P(event) = alpha_k + beta * treatment.
## Conditioning on each trial's total number of events m_k removes alpha_k:
## the number of treatment-arm events X_k then follows a noncentral
## hypergeometric distribution that depends only on beta (the log OR):
##
##   P(X_k = x | m_k) = C(n1k, x) C(n0k, m_k - x) exp(beta x) / normalising sum
##
## Conditional log-likelihood:  l_c(beta) = sum_k [ x_k beta - log S_k(beta) ]
## Conditional information:     I_c(beta) = sum_k Var_beta(X_k | m_k)
## Firth (Jeffreys) penalty:    l*(beta)  = l_c(beta) + 0.5 * log I_c(beta)
##
## Trials with m_k = 0 (double-zero) or with all patients having the event
## carry no conditional information and drop out automatically.
## =====================================================================


## ---- Per-trial building blocks -------------------------------------
## Returns the mean and variance of X_k | m_k, and log S_k, at a given beta.
trial_moments <- function(beta, n1, n0, m) {
  u   <- max(0, m - n0):min(n1, m)                      # possible values of X_k
  lw  <- lchoose(n1, u) + lchoose(n0, m - u) + beta * u # log weights
  mx  <- max(lw)
  w   <- exp(lw - mx)                                   # stabilised weights
  S   <- sum(w)
  p   <- w / S
  mu  <- sum(u * p)
  v   <- sum((u - mu)^2 * p)
  list(logS = log(S) + mx, mean = mu, var = v)
}

## ---- Penalised (or unpenalised) conditional log-likelihood -----------
cond_loglik <- function(beta, x, n1, n0, penalty = TRUE) {
  m  <- x + attr(x, "x0")
  ll <- 0; info <- 0
  for (k in seq_along(x)) {
    tm   <- trial_moments(beta, n1[k], n0[k], m[k])
    ll   <- ll + x[k] * beta - tm$logS
    info <- info + tm$var
  }
  if (penalty) ll + 0.5 * log(info) else ll
}


## ---- Main function -----------------------------------------------------
## x1, n1 : events and sample size in the treatment arm of each trial
## x0, n0 : events and sample size in the control arm of each trial
## penalty: TRUE = conditional Firth; FALSE = ordinary conditional MLE
## level  : confidence level for the profile (penalised) likelihood interval
cond_firth <- function(x1, n1, x0, n0, penalty = TRUE, level = 0.95,
                       bound = 15) {

  m   <- x1 + x0
  ## keep only informative trials: at least one event and at least one non-event
  keep <- m > 0 & m < (n1 + n0)
  x1 <- x1[keep]; n1 <- n1[keep]; x0 <- x0[keep]; n0 <- n0[keep]
  if (length(x1) == 0) stop("No informative trials (all double-zero).")
  x <- x1; attr(x, "x0") <- x0

  f   <- function(b) cond_loglik(b, x, n1, n0, penalty)
  opt <- optimize(f, interval = c(-bound, bound), maximum = TRUE,
                  tol = 1e-10)
  bhat <- opt$maximum
  lmax <- opt$objective

  ## profile (penalised) likelihood confidence interval
  crit <- qchisq(level, df = 1)
  g    <- function(b) 2 * (lmax - f(b)) - crit
  lo <- if (g(-bound) > 0) uniroot(g, c(-bound, bhat))$root else -Inf
  hi <- if (g( bound) > 0) uniroot(g, c(bhat,  bound))$root else  Inf

  ## (penalised) likelihood-ratio test of H0: beta = 0
  lrt <- 2 * (lmax - f(0))
  p   <- pchisq(lrt, df = 1, lower.tail = FALSE)

  ## Wald-type SE from the conditional information (for reference only)
  info <- sum(mapply(function(a, b, c)
    trial_moments(bhat, a, b, c)$var, n1, n0, x1 + x0))

  structure(list(
    method   = if (penalty) "Conditional Firth" else "Conditional MLE",
    logOR    = bhat,
    OR       = exp(bhat),
    se       = 1 / sqrt(info),
    ci.OR    = exp(c(lo, hi)),
    LRT      = lrt,
    p.value  = p,
    k.used   = sum(keep),
    k.dropped= sum(!keep),
    at.bound = abs(bhat) > bound - 1e-3),
    class = "cond_firth")
}

print.cond_firth <- function(x, ...) {
  cat(x$method, "\n")
  cat(sprintf("  Trials used: %d (dropped as uninformative: %d)\n",
              x$k.used, x$k.dropped))
  cat(sprintf("  OR = %.3f   95%% CI (profile) = [%.3f, %.3f]\n",
              x$OR, x$ci.OR[1], x$ci.OR[2]))
  cat(sprintf("  log OR = %.4f (SE %.4f)   LRT = %.3f, p = %.4f\n",
              x$logOR, x$se, x$LRT, x$p.value))
  if (x$at.bound) cat("  WARNING: estimate at search boundary (likely infinite).\n")
  invisible(x)
}

library(metafor)
library(logistf)
source("R/conditional_firth.R")

out <- function(est, lo, hi, p) {
  c(logOR = unname(est), lo = unname(lo),
    hi = unname(hi), p = unname(p))
}

# Mantel-Haenszel without continuity correction
m_mh <- function(x1, n1, x0, n0) {
  m <- rma.mh(
    ai = x1, n1i = n1, ci = x0, n2i = n0,
    measure = "OR",
    add = c(0.5, 0), to = c("only0", "none"),
    drop00 = c(TRUE, TRUE)
  )
  
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# Mantel-Haenszel with continuity correction
m_mh05 <- function(x1, n1, x0, n0) {
  m <- rma.mh(
    ai = x1, n1i = n1, ci = x0, n2i = n0,
    measure = "OR",
    add = c(0.5, 0.5), to = c("only0", "only0"),
    drop00 = c(TRUE, TRUE)
  )
  
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# Peto method
m_peto <- function(x1, n1, x0, n0) {
  m <- rma.peto(ai = x1, n1i = n1, ci = x0, n2i = n0)
  
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# Unconditional Firth logistic regression
m_ufirth <- function(x1, n1, x0, n0) {
  K <- length(x1)
  
  long <- data.frame(
    trial = factor(rep(seq_len(K), 4)),
    treat = rep(c(1, 1, 0, 0), each = K),
    y = rep(c(1, 0, 1, 0), each = K),
    w = c(x1, n1 - x1, x0, n0 - x0)
  )
  
  long <- long[long$w > 0, ]
  long$trial <- droplevels(long$trial)
  
  m <- logistf(y ~ treat + trial, data = long, weights = w)
  i <- which(names(coef(m)) == "treat")
  
  out(coef(m)[i], m$ci.lower[i], m$ci.upper[i], m$prob[i])
}

# Conditional Firth method
m_cfirth <- function(x1, n1, x0, n0) {
  m <- cond_firth(x1, n1, x0, n0)
  
  out(m$logOR, log(m$ci.OR[1]), log(m$ci.OR[2]), m$p.value)
}

methods <- list(
  MH = m_mh,
  MH05 = m_mh05,
  Peto = m_peto,
  UFirth = m_ufirth,
  CFirth = m_cfirth
)
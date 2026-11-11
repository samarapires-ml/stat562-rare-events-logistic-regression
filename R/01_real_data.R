## 01_real_data.R
## Output: results/real_data_results.csv & results/real_data_summary.csv





source("R/methods.R")
dir.create("results", showWarnings = FALSE)




## Load the data
if (requireNamespace("metadat", quietly = TRUE) &&
    "dat.nissen2007" %in% data(package = "metadat")$results[, "Item"]) {
  d <- metadat::dat.nissen2007
} else {
  load("data/dat.nissen2007.rda"); d <- dat.nissen2007
}

outcomes <- c(heart_attack = "infarction", cv_death = "death")




## Describe the data
ratio <- d$treat.total / d$cont.total
summary_tab <- do.call(rbind, lapply(names(outcomes), function(o) {
  x1 <- d[[paste0("treat.", outcomes[o])]]
  x0 <- d[[paste0("cont.",  outcomes[o])]]
  data.frame(
    outcome            = o,
    trials             = nrow(d),
    events_treat       = sum(x1),
    events_control     = sum(x0),
    risk_treat         = sum(x1) / sum(d$treat.total),
    risk_control       = sum(x0) / sum(d$cont.total),
    trials_zero_any    = sum(x1 == 0 | x0 == 0),   # zero in at least one arm
    trials_zero_double = sum(x1 == 0 & x0 == 0),   # zero in both arms
    trials_ratio_gt1.5 = sum(ratio > 1.5)          # allocation > 1.5 : 1
  )
}))
cat("\n===== DATA SUMMARY =====\n")
print(summary_tab, digits = 3, row.names = FALSE)
cat(sprintf("\nAllocation ratio (treatment : control): median %.2f, range %.2f to %.2f\n",
            median(ratio), min(ratio), max(ratio)))






## Fit every method to each outcome
results <- do.call(rbind, lapply(names(outcomes), function(o) {
  x1 <- d[[paste0("treat.", outcomes[o])]]
  x0 <- d[[paste0("cont.",  outcomes[o])]]
  do.call(rbind, lapply(names(methods), function(nm) {
    r <- suppressWarnings(methods[[nm]](x1, d$treat.total, x0, d$cont.total))
    data.frame(outcome = o, method = nm,
               OR = exp(r[["logOR"]]), CI_lower = exp(r[["lo"]]),
               CI_upper = exp(r[["hi"]]), p = r[["p"]])
  }))
}))

cat("\n===== ALL METHODS ON THE REAL DATA =====\n")
print(results, digits = 3, row.names = FALSE)

write.csv(summary_tab, "results/real_data_summary.csv", row.names = FALSE)
write.csv(results,     "results/real_data_results.csv", row.names = FALSE)
cat("\nSaved to results/real_data_summary.csv and results/real_data_results.csv\n")

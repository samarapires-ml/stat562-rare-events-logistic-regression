## 02_pilot.R
## Output: results/pilot_raw.rds, results/pilot_summary.csv




source("R/methods.R")
source("R/simulate.R")
dir.create("results", showWarnings = FALSE)

if (requireNamespace("metadat", quietly = TRUE) &&
    "dat.nissen2007" %in% data(package = "metadat")$results[, "Item"]) {
  d <- metadat::dat.nissen2007
} else {
  load("data/dat.nissen2007.rda"); d <- dat.nissen2007
}






## Settings
R     <- 500       # replicates per scenario (use 100 first for a quick test)
CORES <- max(1, parallel::detectCores() - 2)   # all cores but two, so the computer stays usable

## Control-group risk: the pooled control rate in the real data
p0_mi  <- sum(d$cont.infarction) / sum(d$cont.total)   # about 0.0059
p0_cvd <- sum(d$cont.death)      / sum(d$cont.total)   # about 0.0018

pilot <- list(
  heart_attack = run_sim(R, OR = 1, p0 = p0_mi,  n1 = d$treat.total,
                         n0 = d$cont.total, seed = 101, cores = CORES),
  cv_death     = run_sim(R, OR = 1, p0 = p0_cvd, n1 = d$treat.total,
                         n0 = d$cont.total, seed = 102, cores = CORES)
)
saveRDS(pilot, "results/pilot_raw.rds")







## Summary
summ <- do.call(rbind, lapply(names(pilot), function(o)
  cbind(outcome = o, summarise_sim(pilot[[o]]))))

cat("\n===== PILOT: true OR = 1,", R, "replicates per outcome =====\n")
cat("reject = false-positive rate (target 0.05); coverage target 0.95\n\n")
print(summ, digits = 3, row.names = FALSE)

## average number of zero-event trials per simulated data set
for (o in names(pilot)) {
  z <- unique(pilot[[o]][, c("rep", "zero_any", "zero_double")])
  cat(sprintf("\n%s: on average %.1f trials with a zero arm, %.1f double-zero",
              o, mean(z$zero_any), mean(z$zero_double)))
}
cat("\n")

write.csv(summ, "results/pilot_summary.csv", row.names = FALSE)
cat("\nSaved to results/pilot_raw.rds and results/pilot_summary.csv\n")

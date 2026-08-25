# R/12b_break_robustness.R
# -------------------------------------------------------------------------
# Robustness of the per-pixel confidence to the class-break scheme.
# Addresses the concern that "0.0% of pixels below 0.5" is partly an artefact
# of the wide equal-width breaks. It recomputes the confidence factor P(modal)
# under three break schemes and reports how the headline numbers move.
#
# Mirrors 12_hazard_posterior.R exactly (same weighted-sum overlay, same
# findInterval cut, same 8000 posterior draws), so the baseline row reproduces
# your reported 0.775 mean / 0.0% below 0.5. Only the breaks change.
#
# SAFETY: writes ONLY output/tables/break_robustness.csv. It does not touch
# conf_*.rds, hazard_*, smaa_*, or any other pipeline output.
#
# Run from the project root after step 11 (Bayesian) has been run:
#     Rscript R/12b_break_robustness.R
# or in RStudio:  source("R/12b_break_robustness.R")
# -------------------------------------------------------------------------

cat("\n--- 12b_break_robustness: confidence under alternative class breaks ---\n")
suppressPackageStartupMessages({ library(terra) })

if (!file.exists("output/bayes_weights.rds")) stop("Missing output/bayes_weights.rds -- run R/11_bayesian_ahp.R first.")
if (!file.exists("output/reclass_stack.rds")) stop("Missing output/reclass_stack.rds -- run R/03_reclassify.R first.")

W  <- readRDS("output/bayes_weights.rds")            # M x 5 posterior weight draws
rs <- unwrap(readRDS("output/reclass_stack.rds"))    # 5 reclassified factor layers
M  <- nrow(W)

# valid pixels as a matrix (same as 12_hazard_posterior.R)
vals <- values(rs)
ok   <- stats::complete.cases(vals)
Xok  <- vals[ok, , drop = FALSE]
n_ok <- nrow(Xok)
cat(sprintf("%d valid pixels, %d posterior draws.\n", n_ok, M))

# mean posterior overlay, used only to place the quantile breaks
mean_score <- as.numeric(Xok %*% colMeans(W))

# ---- three break schemes ----------------------------------------------------
brk_ew5  <- seq(1, 3, length.out = 6)                                  # baseline
brk_q5   <- as.numeric(quantile(mean_score, probs = seq(0, 1, 0.2)))   # quintiles
brk_q5[1] <- brk_q5[1] - 1e-9                                          # include lowest
brk_ew10 <- seq(1, 3, length.out = 11)                                 # finer (10 class)
schemes <- list(
  "equal-width 5 (baseline)" = brk_ew5,
  "quantile 5 (dense boundaries)" = brk_q5,
  "equal-width 10 (finer)" = brk_ew10
)

cut_to_class <- function(score, brk) {
  ncl <- length(brk) - 1L
  cl <- findInterval(score, brk, rightmost.closed = TRUE)
  cl[cl < 1] <- 1L; cl[cl > ncl] <- ncl
  as.integer(cl)
}

# ---- accumulate class counts for all schemes in one pass over the draws ----
counts <- lapply(schemes, function(b) matrix(0L, nrow = n_ok, ncol = length(b) - 1L))
cat("Streaming draws"); pb <- txtProgressBar(min = 0, max = M, style = 3)
for (m in 1:M) {
  score <- as.numeric(Xok %*% W[m, ])
  for (nm in names(schemes)) {
    cl <- cut_to_class(score, schemes[[nm]])
    counts[[nm]][cbind(seq_len(n_ok), cl)] <- counts[[nm]][cbind(seq_len(n_ok), cl)] + 1L
  }
  setTxtProgressBar(pb, m)
}
close(pb)

# ---- summarise each scheme --------------------------------------------------
summ <- data.frame()
for (nm in names(schemes)) {
  prob   <- counts[[nm]] / M
  ncl    <- ncol(prob)
  pmodal <- prob[cbind(seq_len(n_ok), max.col(prob, ties.method = "first"))]
  ent    <- -rowSums(ifelse(prob > 0, prob * log(prob), 0)) / log(ncl)
  summ <- rbind(summ, data.frame(
    break_scheme = nm,
    n_classes    = ncl,
    mean_Pmodal  = round(mean(pmodal), 3),
    median_Pmodal = round(median(pmodal), 3),
    pct_above_0.8 = round(100 * mean(pmodal > 0.8), 1),
    pct_below_0.5 = round(100 * mean(pmodal < 0.5), 1),
    mean_norm_entropy = round(mean(ent), 3)
  ))
}
cat("\n\nConfidence under alternative break schemes:\n"); print(summ, row.names = FALSE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(summ, "output/tables/break_robustness.csv", row.names = FALSE)
cat("\nSaved ONLY: output/tables/break_robustness.csv  (no other file touched).\n")
cat("The 'equal-width 5 (baseline)' row should reproduce your reported 0.775 / 0.0% below 0.5.\n\n")

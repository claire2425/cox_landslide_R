# R/08_oat_sensitivity.R
# Step 8: OAT (one-at-a-time) sensitivity check.
#
# What it does: Finds out which pairwise judgements actually matter.
# I'm tweaking them one by one (+/- 5, 10, 20%) to see how much the final weights move.
# This gives me the high-leverage judgements to watch out for in the MC stage.
#
# Inputs: output/ahp_5factor_baseline.rds
# Outputs: CSVs and the tornado plot.
# Safe to re-run anytime.

cat("\n--- 08_oat_sensitivity: one-at-a-time sensitivity of AHP weights ---\n")

if (!exists("ahp_weights")) source("R/00_functions.R")

## --- load baseline submatrix -----------------------------------------------
rds_path<- "output/ahp_5factor_baseline.rds"
if (!file.exists(rds_path)) {
  stop("Missing ", rds_path, " -- run R/01b_build_ahp_matrix.R first.")
}
base<- readRDS(rds_path)
A0<- base$matrix
w0<- ahp_weights(A0)
cr0<- ahp_cr(A0)$CR
n<- nrow(A0)
fac<- rownames(A0)

cat("\nBaseline weights:\n"); print(round(w0, 4))
cat(sprintf("Baseline CR = %.4f\n", cr0))

## --- perturbation grid -----------------------------------------------------
deltas<- c(-0.20, -0.10, -0.05, 0.05, 0.10, 0.20)

## helper: rebuild a reciprocal matrix after scaling one upper-triangle entry
perturb<- function(A, i, j, d) {
  Ap<- A
  Ap[i, j]<- A[i, j] * (1 + d)
  Ap[j, i]<- 1 / Ap[i, j]      # keep reciprocity
  Ap
}

## --- run OAT over every upper-triangle entry -------------------------------
rows<- list()
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    entry<- paste(fac[i], fac[j], sep = "~")
    for (d in deltas) {
      Ap<- perturb(A0, i, j, d)
      wp<- ahp_weights(Ap)
      crp<- ahp_cr(Ap)$CR
      rows[[length(rows) + 1]]<- data.frame(
        entry        = entry,
        i_factor     = fac[i],
        j_factor     = fac[j],
        delta_pct    = d * 100,
        CR           = round(crp, 4),
        # weight changes (perturbed - baseline) for each factor
        d_Slope      = round(wp["Slope"]     - w0["Slope"],     4),
        d_LandCover  = round(wp["LandCover"] - w0["LandCover"], 4),
        d_NDVI       = round(wp["NDVI"]      - w0["NDVI"],       4),
        d_SoilType   = round(wp[fac[grep("Soil", fac)]] - w0[fac[grep("Soil", fac)]], 4),
        d_DistRoad   = round(wp["DistRoad"]  - w0["DistRoad"],  4),
        # scalar leverage = max absolute weight change across all factors
        max_abs_shift = round(max(abs(wp - w0)), 4),
        stringsAsFactors = FALSE
      )
    }
  }
}
oat<- do.call(rbind, rows)
rownames(oat)<- NULL

## --- rank entries by leverage (largest weight move at +/-20%) --------------
ranked<- aggregate(max_abs_shift ~ entry, data = oat, FUN = max)
ranked<- ranked[order(-ranked$max_abs_shift), ]
rownames(ranked)<- NULL

cat("\nEntries ranked by leverage (max |weight shift| at +/-20%):\n")
print(ranked)

most<- ranked$entry[1]
cat(sprintf("\nMost sensitive judgement: %s (shifts a weight by up to %.3f at +/-20%%)\n",
            most, ranked$max_abs_shift[1]))

## --- write tables ----------------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(oat,    "output/tables/oat_sensitivity.csv",        row.names = FALSE)
write.csv(ranked, "output/tables/oat_sensitivity_ranked.csv", row.names = FALSE)
cat("\nSaved: output/tables/oat_sensitivity.csv, oat_sensitivity_ranked.csv\n")

## --- tornado plot (signed, left = negative / right = positive perturbation) -
# Classic OAT tornado that keeps DIRECTION: each judgement gets bars to the LEFT
# (negative perturbation) and to the RIGHT (positive perturbation). Three shades
# show the +/-5/10/20% magnitudes, so the plot reveals BOTH the leverage ranking
# AND that larger perturbations move the weights further, in opposite directions.
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

mags<- c(5, 10, 20)
# leverage magnitude per (entry, |delta|): the max |weight shift| across factors
lev<- function(entry, m) {
  v<- oat[oat$entry == entry & abs(oat$delta_pct) == m, "max_abs_shift"]
  if (length(v)) max(v) else 0
}
entries<- unique(oat$entry)
entries<- entries[order(sapply(entries, lev, m = 20))]   # smallest at bottom
ne<- length(entries)

# build signed matrices: negative side (left) and positive side (right)
neg<- t(sapply(entries, function(e) -sapply(mags, lev, entry = e)))  # ne x 3, negative
pos<- t(sapply(entries, function(e)  sapply(mags, lev, entry = e)))  # ne x 3, positive

png("output/figures/oat_tornado.png", width = 1150, height = 760, res = 130)
op<- par(mar = c(4.2, 9, 3, 2))
xmax<- max(pos) * 1.12
plot(NA, xlim = c(-xmax, xmax), ylim = c(0.5, ne + 0.5),
     yaxt = "n", xlab = "Weight shift (left = -perturbation, right = +perturbation)",
     ylab = "", main = "OAT sensitivity: leverage of each pairwise judgement at +/-5/10/20%")
abline(v = 0, col = "grey40")

# colours: light -> dark for 5/10/20%
neg_cols<- c("cadetblue2", "cadetblue3", "cadetblue4")   # greens, left
pos_cols<- c("orange", "darkorange1", "darkorange3")   # oranges, right
bh<- 0.26   # bar height per magnitude band

for (k in seq_along(entries)) {
  yc<- k
  # draw the three magnitude bands as nested bars (20% widest behind, 5% in front)
  ord3<- order(mags, decreasing = TRUE)   # draw 20 first (back), then 10, then 5
  for (b in ord3) {
    rect(neg[k, b], yc - bh, 0, yc + bh, col = neg_cols[b], border = NA)
    rect(0, yc - bh, pos[k, b], yc + bh, col = pos_cols[b], border = NA)
  }
}
axis(2, at = seq_len(ne), labels = entries, las = 1, cex.axis = 0.8)
legend("bottomright", bty = "n", cex = 0.8, title = "perturbation",
       fill = c(neg_cols[3], neg_cols[2], neg_cols[1], pos_cols[1], pos_cols[2], pos_cols[3]),
       legend = c("-20%", "-10%", "-5%", "+5%", "+10%", "+20%"), ncol = 2)
par(op); dev.off()
cat("Saved: output/figures/oat_tornado.png\n")

cat("\n--- end 08_oat_sensitivity ---\n\n")

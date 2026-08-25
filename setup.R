# setup.R
# One-off package installation + sanity check. Run once on a fresh machine.
# Set the working directory to the repository root before running
# (in RStudio use "Open Project", or setwd() to the folder that contains
#  this setup.R and the R/ folder).

cat("=== Cox's Bazar LSM project: setup ===\n\n")

required<- c(
  "terra",      # raster operations
  "sf",         # vector / shapefile operations
  "tmap",       # thematic GIS maps (v3 or v4 both supported by 06_make_map.R)
  "ggplot2",    # non-map plotting
  "ggspatial",  # scale bars / north arrows on ggplot maps
  "cowplot",    # combining plots into panels
  "pROC",       # ROC / AUC validation (used in 05_validate.R)
  "rstan",      # Bayesian AHP (Stan / NUTS)
  "posterior",  # posterior draw handling
  "geodata",    # GADM administrative boundaries
  "knitr"       # rendering the .Rmd reports
)

installed<- rownames(installed.packages())
to_install<- setdiff(required, installed)
if (length(to_install) > 0) {
  cat("Installing missing packages:\n  ", paste(to_install, collapse = ", "), "\n\n")
  install.packages(to_install)
} else {
  cat("All required packages already installed.\n\n")
}

cat("Loaded packages and versions:\n")
for (pkg in required) {
  ok<- suppressPackageStartupMessages(
    requireNamespace(pkg, quietly = TRUE))
  cat(sprintf("  %-12s  %s\n", pkg,
              if (ok) as.character(packageVersion(pkg)) else "NOT INSTALLED"))
}

cat("\nR version:", R.version.string, "\n")
if (getRversion() < "4.1.0")
  warning("R >= 4.1.0 recommended.")

# Project-root sanity check
expected<- c("setup.R", "README.md", "R")
missing<- expected[!file.exists(expected)]
if (length(missing) > 0) {
  warning(sprintf(
    "Working directory may not be the project root. Missing: %s\nUse setwd() or RStudio 'Open Project'.",
    paste(missing, collapse = ", ")))
} else {
  cat("\nWorking directory looks correct.\n")
}

cat("\n=== Setup complete ===\n")
cat("Next: source('scripts/run_all.R')\n") #if wanna run it all at a time

# scripts/run_all.R
# Master script: runs the full Cox's Bazar LSM pipeline in order, stopping
# cleanly at the first error. Run from the project ROOT:
#
#     source("setup.R")            # one-off package install
#     source("scripts/run_all.R")  # full pipeline
#
# PIPELINE = BASELINE (01-07) + UNCERTAINTY CHAIN (08-12)
#
# BASELINE
#   01_ahp_pairwise      -> output/ahp_weights.rds (+ matrix, diagnostics)
#   01b_build_ahp_matrix -> output/ahp_5factor_baseline.rds, ahp_construction_trace.csv
#                           (traceable Ahmed Table 3 -> 5-factor submatrix.
#                            NOTE: the soil factor is relabelled SoilPerm ->
#                            SoilType here as a declared PROXY; see the header
#                            of 01b and proxy_variable_methodology.md.)
#   02_load_data         -> output/rasters_raw.rds, inventory.rds, boundary.rds
#   03_reclassify        -> output/reclass_stack.rds
#   [GUARD] factor-name alignment check (matrix rownames == stack names)
#   07_verify_inputs (A) -> CHECK 1: code->class mapping vs .vat.dbf
#   04_hazard_overlay    -> output/hazard_class.rds, hazard_continuous.rds
#   05_validate          -> output/tables/validation_summary.csv
#   07_verify_inputs (B) -> CHECK 2: reclass direction (AUC)
#   06_make_map          -> output/figures/*.png
#
# UNCERTAINTY CHAIN
#   08_oat_sensitivity   -> RQ1:  one-at-a-time leverage   (reads 01b)
#   09_monte_carlo       -> RQ2:  weight-spread under noise (reads 01b)
#   10_bootstrap         -> RQ2:  AUC/lift sampling CIs     (reads 04 + inventory)
#   11_bayesian_ahp      -> RQ3a: weight POSTERIOR + tau    (reads 01b)
#   11b_prior_sensitivity-> RQ3a: is the tau finding robust to its prior?
#                                  (reads 01b; Stan-only; gated by RUN_PRIOR_SENS)
#   12_hazard_posterior  -> RQ3b: per-pixel confidence map  (reads 11 + 03/04)
#   13_smaa              -> RQ3c: per-pixel SMAA acceptability + central weights
#                                  (reads 11 + 03/04; gated by RUN_SMAA)
#
# DEPENDENCIES THAT FIX THE ORDER
#   08, 09, 11 only need the AHP matrix (01b) -> placed right after the guard.
#   11b refits the same model under different priors -> right after 11.
#   10 needs the baseline hazard (04) + inventory (02) -> after 05.
#   12 needs the posterior (11) AND the reclass stack (03) + class cuts (04).
#   13 needs the posterior (11) AND the reclass stack + the SAME class cuts as 12
#      -> right after 12.
#
# 07 is run TWICE on purpose (see CHECK 1 / CHECK 2 above).
#
# GIS-DEPENDENT STEPS: 02, 03, 04, 06, 10, 12, 13 touch rasters (terra/sf). If you
# run this on a machine without the GIS stack, set RUN_GIS<- FALSE to run only
# the matrix-side steps (01, 01b, 08, 09, 11, 11b); the raster steps are skipped
# with a clear message. 12 and 13 are skipped automatically when RUN_GIS is FALSE
# because they need the reclass stack.
#
# COMPLETE PROJECT FILE INVENTORY (so no file looks unused)
#
# Of the 20 source files (18 .R + 2 .stan), this orchestrator runs 15 directly.
# The remaining files are not "missing" from the pipeline — they have distinct
# roles and are deliberately not in the step list, as explained below.
#
# (A) MAIN PIPELINE — run automatically by this script (15 files):
#       01_ahp_pairwise, 01b_build_ahp_matrix, 02_load_data, 03_reclassify,
#       04_hazard_overlay, 05_validate, 06_make_map, 07_verify_inputs (x2),
#       08_oat_sensitivity, 09_monte_carlo, 10_bootstrap, 11_bayesian_ahp,
#       11b_prior_sensitivity, 12_hazard_posterior, 13_smaa.
#       (11b gated by RUN_PRIOR_SENS; 13 gated by RUN_SMAA — both default TRUE.)
#
# (B) SHARED FUNCTION LIBRARY — loaded, not executed as a step (1 file):
#       00_functions.R holds the reusable helpers (reclass_codes, ahp_weights,
#       align_to, ...). Every pipeline script sources it at its top, so the
#       whole pipeline depends on it; listing it as a step would merely redefine
#       the same functions. It is therefore loaded indirectly, never run alone.
#
# (C) ON-DEMAND DIAGNOSTIC / COMPARISON TOOLS — run individually when writing
#     up, not on every pipeline run (3 files):
#       03b_slope_flat_comparison : gentle-slope "include vs exclude" AUC table
#                                   (supports the appendix justification).
#       05b_direction_diagnostic  : tests whether any factor's reclass DIRECTION
#                                   is reversed (this is how the NDVI inversion
#                                   was diagnosed).
#       diagnostic_test           : quantifies the AUC cost of keeping the
#                                   physically-principled settings vs alternatives
#                                   (evidence for, e.g., rejecting the SoilType
#                                   flip). Each tool sources whatever earlier
#                                   steps it needs, so it runs standalone after a
#                                   normal pipeline pass. Re-running these on
#                                   every pipeline pass would waste time, so they
#                                   are intentionally excluded here.
#
# (D) STAN MODELS — read as data, not sourced (2 files):
#       R/stan/bayes_ahp.stan        : the Bayesian AHP model, read by 11 when
#                                      ENGINE = "stan" (rstan::stan(file = ...)).
#       R/stan/bayes_ahp_priors.stan : the SAME model with prior SDs passed in as
#                                      data, read by 11b for the prior sweep.
#       Neither is sourced as R code.
#
# Also present, outside this orchestrator's remit:
#   setup.R                    : one-off package installation (run once, first).

## ---- options --------------------------------------------------------------
RUN_GIS<- TRUE   # FALSE = skip raster steps (02,03,04,06,10,12,13)
RUN_UNCERTAINTY<- TRUE   # FALSE = baseline only (stop after 06)
RUN_PRIOR_SENS<- TRUE   # 11b: prior-sensitivity sweep on tau (Stan-only, ~3x an 11 run)
RUN_SMAA<- TRUE   # 13: SMAA per-pixel acceptability + central weights (needs GIS)

GIS_STEPS<- c("R/02_load_data.R", "R/03_reclassify.R", "R/04_hazard_overlay.R",
               "R/06_make_map.R", "R/10_bootstrap.R", "R/12_hazard_posterior.R",
               "R/13_smaa.R")

## ---- runner ---------------------------------------------------------------
run_step<- function(path) {
  if (!RUN_GIS && path %in% GIS_STEPS) {
    cat(sprintf("\n--- SKIP %s (RUN_GIS = FALSE) ---\n", path)); return(invisible(FALSE))
  }
  cat(sprintf("\n=================== %s ===================\n", path))
  ok<- tryCatch({ source(path, local = new.env()); TRUE },
                 error = function(e) {
                   cat(sprintf("\n!! ERROR in %s:\n   %s\n", path, conditionMessage(e)))
                   FALSE
                 })
  if (!ok) stop(sprintf("Pipeline stopped at %s. Fix the error above and re-run.", path),
                call. = FALSE)
  invisible(TRUE)
}

## ---- factor-alignment guard ----------------------------------------------
# Catches the SoilType/SoilPerm class of error (matrix label vs raster label)
# the MOMENT it appears, instead of letting it surface only at step 11/12.
check_factor_alignment<- function() {
  cat("\n=================== [GUARD] factor-name alignment ===================\n")
  if (!RUN_GIS) { cat("RUN_GIS = FALSE -> reclass stack not built; skipping guard.\n"); return(invisible()) }
  m_path<- "output/ahp_5factor_baseline.rds"
  s_path<- "output/reclass_stack.rds"
  if (!file.exists(m_path) || !file.exists(s_path)) {
    cat("Cannot run guard yet (missing matrix or stack rds).\n"); return(invisible())
  }
  A<- readRDS(m_path)$matrix
  rs<- readRDS(s_path)
  mn<- sort(rownames(A)); sn<- sort(names(rs))
  if (identical(mn, sn)) {
    cat("OK: matrix factors == reclass-stack factors:\n   ",
        paste(mn, collapse = ", "), "\n")
  } else {
    cat("MISMATCH between AHP matrix and reclass stack:\n")
    cat("  in matrix only :", paste(setdiff(mn, sn), collapse = ", "), "\n")
    cat("  in stack only  :", paste(setdiff(sn, mn), collapse = ", "), "\n")
    stop("Factor names do not align. Most likely the soil factor: the matrix ",
         "should carry the PROXY label 'SoilType' (re-run R/01b_build_ahp_matrix.R ",
         "after the proxy fix), and the reclass stack (03) should name it 'SoilType' ",
         "too. Reconcile the names before continuing — otherwise weights get applied ",
         "to the wrong rasters silently.", call. = FALSE)
  }
}

## ---- baseline steps -------------------------------------------------------
baseline_pre<- c(
  "R/01_ahp_pairwise.R",
  "R/01b_build_ahp_matrix.R",   # proxy: SoilPerm -> SoilType (declared)
  "R/02_load_data.R",
  "R/03_reclassify.R"
)
baseline_post<- c(
  "R/07_verify_inputs.R",   # CHECK 1 runs; CHECK 2 will SKIP (no hazard yet)
  "R/04_hazard_overlay.R",
  "R/05_validate.R",
  "R/07_verify_inputs.R",   # CHECK 2 now runs (stack + inventory exist)
  "R/06_make_map.R"
)

## ---- uncertainty chain ----------------------------------------------------
# Split around the baseline because 10, 12 and 13 need baseline raster outputs.
unc_matrix_side<- c(            # need only the AHP matrix (after 01b)
  "R/08_oat_sensitivity.R",     # RQ1
  "R/09_monte_carlo.R",         # RQ2 (weights)
  "R/11_bayesian_ahp.R"         # RQ3a (posterior + tau)
)
# 11b (prior sensitivity) needs 01b's matrix + the same Stan model as 11; it is
# a robustness check on tau, so it runs right after 11 on the matrix side, gated
# by RUN_PRIOR_SENS (it refits the model ~3 times, so it is opt-in).
unc_prior_sens<- "R/11b_prior_sensitivity.R"

unc_raster_side<- c(           # need baseline hazard / stack
  "R/10_bootstrap.R",           # RQ2 (AUC/lift CIs)    -- needs 04 + inventory
  "R/12_hazard_posterior.R"     # RQ3b (confidence map) -- needs 11 + 03 + 04
)
# 13 (SMAA) consumes 11's posterior + the reclass stack + the same class breaks
# as 12, so it runs after 12 on the raster side, gated by RUN_SMAA.
unc_smaa<- "R/13_smaa.R"

## ---- execute --------------------------------------------------------------
cat("=== Cox's Bazar LSM pipeline: run_all ===\n")
cat(sprintf("RUN_GIS = %s | RUN_UNCERTAINTY = %s\n", RUN_GIS, RUN_UNCERTAINTY))

for (s in baseline_pre)  run_step(s)
check_factor_alignment()                 # <-- guard right after the stack exists
for (s in baseline_post) run_step(s)

if (RUN_UNCERTAINTY) {
  cat("\n=================== UNCERTAINTY CHAIN ===================\n")
  for (s in unc_matrix_side) run_step(s)   # 08, 09, 11
  if (RUN_PRIOR_SENS) run_step(unc_prior_sens)   # 11b (after 11; Stan-only)
  for (s in unc_raster_side) run_step(s)   # 10, 12 (12 auto-skips if RUN_GIS=FALSE)
  if (RUN_SMAA) run_step(unc_smaa)         # 13 (after 12; auto-skips if RUN_GIS=FALSE)
} else {
  cat("\nRUN_UNCERTAINTY = FALSE -> baseline only; skipping 08-13.\n")
}

## ---- closing report -------------------------------------------------------
cat("\n=== Pipeline complete ===\n")
cat("Verification report: output/tables/input_verification.csv\n")
ver<- tryCatch(read.csv("output/tables/input_verification.csv"), error = function(e) NULL)
if (!is.null(ver)) {
  cat("\nFinal input-verification status:\n")
  print(ver[, c("check", "status")], row.names = FALSE)
  if (any(ver$status == "FAIL"))
    cat("\n>> WARNING: at least one verification FAILED — see the report.\n")
}

# Quick pointer to the uncertainty-chain outputs, when they were produced.
if (RUN_UNCERTAINTY) {
  cat("\nUncertainty-chain outputs to check:\n")
  cat("  RQ1 OAT      : output/tables/oat_sensitivity_ranked.csv, figures/oat_tornado.png\n")
  cat("  RQ2 MC       : output/tables/mc_intervals.csv, output/mc_weights.rds\n")
  cat("  RQ2 bootstrap: output/tables/bootstrap_validation.csv\n")
  cat("  RQ3a Bayesian: output/tables/bayes_intervals.csv, bayes_tau.csv, output/bayes_weights.rds\n")
  if (RUN_PRIOR_SENS)
    cat("  RQ3a prior-sens: output/tables/prior_sensitivity_tau.csv, prior_sensitivity_weights.csv,\n               output/figures/prior_sensitivity_tau.png\n")
  if (RUN_GIS)
    cat("  RQ3b confid. : output/figures/confidence_map.png, output/conf_*.rds\n")
  if (RUN_GIS && RUN_SMAA)
    cat("  RQ3c SMAA    : output/figures/smaa_acceptability.png, smaa_confidence.png,\n               output/tables/smaa_summary.csv, smaa_central_weights_by_class.csv\n")
  cat("\nKey number for RQ3 narrative: compare bayes_tau.csv (tau) against the\n")
  cat("Monte Carlo sigma (0.3) -- let the run decide the direction.\n")
  if (RUN_PRIOR_SENS)
    cat("Robustness: prior_sensitivity_tau.csv shows whether that tau finding\n  survives reasonable changes to the tau prior.\n")
}

